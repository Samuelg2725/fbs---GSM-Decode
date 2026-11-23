"""Read the LTE cell identity from one capture.

The identity in the target line is the 28-bit cellIdentity from SIB1, not the
PCI. PCI 430 only says which sync sequence was heard. Cell identity 7159822
is the broadcast name of that sector.

This is a one-shot. It does not search tuner offset and it does not loop, so
it can sit beside the direction finder and run once after a PCI lock.

A 2.4 MHz capture only contains the centre of a 5 or 10 MHz cell, so SIB1
is out of reach there. A HackRF capture at 7.68 MS/s or 15.36 MS/s covers
those cells, and the wideband path reads SIB1 from subframe 5.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from lte_coding import (
    bits_of,
    bits_to_int,
    code_block_size,
    conv_encode,
    conv_rate_match,
    conv_rate_recover,
    crc16,
    crc24a,
    gold,
    qpsk_mod,
    qpsk_soft,
    turbo_decode,
    turbo_encode,
    turbo_rate_match,
    turbo_rate_recover,
    viterbi_exact,
    viterbi_tailbiting,
    VITERBI_TABLES,
    _ones,
)
from lte_sync import CP, FFT_SIZE, FS, SYMBOL, pss_time, sss_freq

NFFT = FFT_SIZE
PBCH_BITS = 1920
FRAME_PHASES = 4
# 10 MHz first. That is the office Band 20 carrier, and a miss still tries the rest.
BANDWIDTHS = (50, 25, 100, 75, 15, 6)
_GOLD = {}
BW_MHZ = {6: 1.4, 15: 3, 5: 5, 25: 5, 50: 10, 75: 15, 100: 20}
BW_INDEX = {6: 0, 15: 1, 25: 2, 50: 3, 75: 4, 100: 5}
ANTENNA_MASK = {1: 0x0000, 2: 0xFFFF, 4: 0x5555}
# CRS is generated as if every cell were 110 resource blocks wide, then the
# live bandwidth uses the centre of that sequence. TS 36.211 6.10.1.2.
N_RB_MAX_DL = 110
SI_RNTI = 65535
# Transport sizes used by a 1.4 MHz SIB1. The CRC picks the one that was sent.
_TBS_CANDIDATES = (
    16, 24, 32, 40, 56, 72, 88, 104, 120, 136, 144, 152, 176, 208, 224, 256,
    280, 328, 376, 392, 408, 440, 488, 552, 600, 616, 648, 680, 712, 744, 776,
)


def cp_len(symbol: int) -> int:
    return 10 if symbol % 7 == 0 else CP


def symbol_useful_offset(symbol: int) -> int:
    cursor = 0
    for index in range(symbol):
        cursor += cp_len(index) + NFFT
    return cursor + cp_len(symbol)


def _k_to_bin(k: int, n_rb: int) -> int:
    """FFT bin of absolute subcarrier k when the tuner is on the cell centre."""
    delta = k - n_rb * 6
    if delta < 0:
        return delta
    return delta + 1


def _rs_residues(pci: int) -> tuple[int, int]:
    shift = pci % 6
    return shift % 6, (shift + 3) % 6


def _pbch_keys(n_rb: int):
    return [n_rb * 6 - 36 + index for index in range(72)]


def pbch_data_res(n_rb: int, pci: int):
    """PBCH resource elements, k fastest then symbol, skipping 4-port CRS."""
    low, high = _rs_residues(pci)
    elements = []
    for symbol in range(4):
        for k in _pbch_keys(n_rb):
            if symbol < 2 and (k % 6) in (low, high):
                continue
            elements.append((symbol, k))
    return elements


def _crs_init(pci: int, n_s: int, symbol: int) -> int:
    return (2 ** 10) * (7 * (n_s + 1) + symbol + 1) * (2 * pci + 1) + 2 * pci + 1


def crs_value(pci: int, n_s: int, symbol: int, k: int, port: int, n_rb: int):
    """Complex CRS for one subcarrier, or None when this port has no CRS there.

    symbol is the index inside the slot. The sequence index is the centre of
    the 110-block reference, not the start of this cell. TS 36.211 6.10.1.2.
    """
    if port == 0:
        v = 0 if symbol == 0 else 3
        if symbol not in (0, 4):
            return None
    elif port == 1:
        v = 3 if symbol == 0 else 0
        if symbol not in (0, 4):
            return None
    elif port == 2 and symbol == 1:
        v = 3 * (n_s % 2)
    elif port == 3 and symbol == 1:
        v = 3 + 3 * (n_s % 2)
    else:
        return None
    residue = (v + (pci % 6)) % 6
    if k % 6 != residue:
        return None
    m = (k - residue) // 6
    if m < 0:
        return None
    m_prime = m + N_RB_MAX_DL - int(n_rb)
    c_init = _crs_init(pci, n_s, symbol)
    sequence = _GOLD.get(c_init)
    if sequence is None or len(sequence) < 2 * (m_prime + 1):
        sequence = gold(c_init, max(2 * (m_prime + 1), 2 * 2 * N_RB_MAX_DL))
        _GOLD[c_init] = sequence
    scale = 1.0 / (2 ** 0.5)
    return scale * ((1 - 2 * sequence[2 * m_prime]) + 1j * (1 - 2 * sequence[2 * m_prime + 1]))


def _pack_mib(n_rb: int, sfn: int) -> list[int]:
    bits = bits_of(BW_INDEX[n_rb], 3)
    bits += [0, 1, 0]  # PHICH normal, Ng = 1
    bits += bits_of((sfn >> 2) & 0xFF, 8)
    bits += [0] * 10
    return bits


def _parse_mib(bits) -> dict:
    index = bits_to_int(bits[:3])
    n_rb = {value: key for key, value in BW_INDEX.items()}.get(index)
    ng = {0: "1/6", 1: "1/2", 2: "1", 3: "2"}[bits_to_int(bits[4:6])]
    return {
        "bandwidth_rb": n_rb,
        "bandwidth_mhz": None if n_rb is None else BW_MHZ[n_rb],
        "phich_duration": "extended" if bits[3] else "normal",
        "phich_resource": ng,
        "sfn_msb": bits_to_int(bits[6:14]),
    }


def _nearest_gain(gains, k):
    """Channel at this subcarrier, interpolated between the neighbouring CRS tones."""
    if not gains:
        return 1 + 0j
    if k in gains:
        return gains[k]
    keys = sorted(gains)
    if k <= keys[0]:
        return gains[keys[0]]
    if k >= keys[-1]:
        return gains[keys[-1]]
    higher = keys[0]
    for key in keys:
        if key >= k:
            higher = key
            break
    lower = keys[keys.index(higher) - 1]
    weight = (k - lower) / float(higher - lower)
    return (1.0 - weight) * gains[lower] + weight * gains[higher]


def _pbch_soft(spectra, pci: int, n_rb: int, ports: int):
    """Equalised QPSK soft bits for one PBCH occasion, or None.

    One gain for the whole 1.4 MHz block fails when the symbol window is a
    sample or two off: that error is a phase ramp across the subcarriers.
    Each reference tone sees the same ramp, so the data tone uses the nearest one.
    """
    elements = pbch_data_res(n_rb, pci)
    if len(elements) != 240:
        return None
    channel = [{}, {}]
    for port in range(min(ports, 2)):
        for k in _pbch_keys(n_rb):
            ref = crs_value(pci, 1, 0, k, port, n_rb)
            if ref is None:
                continue
            measured = spectra[0][_k_to_bin(k, n_rb)] * np.conj(ref)
            channel[port][k] = measured
        if not channel[port]:
            return None
    soft = []
    if ports == 1:
        for symbol, k in elements:
            gain = _nearest_gain(channel[0], k)
            value = spectra[symbol][_k_to_bin(k, n_rb)] * np.conj(gain)
            soft.extend(qpsk_soft(value / (np.abs(gain) ** 2 + 1e-12)))
        return soft
    if ports == 2:
        for index in range(0, len(elements), 2):
            s0, k0 = elements[index]
            s1, k1 = elements[index + 1]
            h0 = 0.5 * (_nearest_gain(channel[0], k0) + _nearest_gain(channel[0], k1))
            h1 = 0.5 * (_nearest_gain(channel[1], k0) + _nearest_gain(channel[1], k1))
            y0 = spectra[s0][_k_to_bin(k0, n_rb)]
            y1 = spectra[s1][_k_to_bin(k1, n_rb)]
            scale = np.abs(h0) ** 2 + np.abs(h1) ** 2 + 1e-12
            x0 = (np.conj(h0) * y0 + h1 * np.conj(y1)) / scale
            x1 = (np.conj(h0) * y1 - h1 * np.conj(y0)) / scale
            soft.extend(qpsk_soft(x0))
            soft.extend(qpsk_soft(x1))
        return soft
    if ports == 4:
        other = [{}, {}]
        for port, store in ((2, other[0]), (3, other[1])):
            for k in _pbch_keys(n_rb):
                ref = crs_value(pci, 1, 1, k, port, n_rb)
                if ref is None:
                    continue
                store[k] = spectra[1][_k_to_bin(k, n_rb)] * np.conj(ref)
            if not store:
                return None
        for index in range(0, len(elements), 4):
            group = elements[index:index + 4]
            gains = (
                _nearest_gain(channel[0], group[0][1]),
                _nearest_gain(channel[1], group[0][1]),
                _nearest_gain(other[0], group[2][1]),
                _nearest_gain(other[1], group[2][1]),
            )
            for pair in (0, 2):
                symbol0, k0 = group[pair]
                symbol1, k1 = group[pair + 1]
                y0 = spectra[symbol0][_k_to_bin(k0, n_rb)]
                y1 = spectra[symbol1][_k_to_bin(k1, n_rb)]
                h0, h1 = gains[pair], gains[pair + 1]
                scale = np.abs(h0) ** 2 + np.abs(h1) ** 2 + 1e-12
                x0 = (np.conj(h0) * y0 + h1 * np.conj(y1)) / scale
                x1 = (np.conj(h0) * y1 - h1 * np.conj(y0)) / scale
                soft.extend(qpsk_soft(x0))
                soft.extend(qpsk_soft(x1))
        return soft
    return None


def _score_message(message, phase: int):
    """Return the legal MIB whose checksum is closest, and that distance.

    The distance and the body have to be the same candidate. A one-bit checksum
    miss on a real bandwidth is the cell. An illegal bandwidth is ignored.
    """
    payload, parity = message[:24], message[24:]
    found = bits_to_int(parity)
    wanted = bits_to_int(crc16(payload))
    closest = 16
    best = None
    for ports, mask in ANTENNA_MASK.items():
        distance = _ones(found ^ mask ^ wanted)
        body = _parse_mib(payload)
        if body["bandwidth_rb"] is None or distance > closest:
            continue
        body["antennas"] = ports
        body["sfn"] = (body["sfn_msb"] << 2) | phase
        closest = distance
        best = body
    return best, closest


def _mib_from_soft(soft, pci: int, phase: int):
    """Return (mib or None, best CRC Hamming distance)."""
    scrambled = gold(pci, PBCH_BITS)
    start = phase * (PBCH_BITS // FRAME_PHASES)
    chunk = scrambled[start:start + len(soft)]
    llr = [metric if bit == 0 else -metric for metric, bit in zip(soft, chunk)]
    streams = conv_rate_recover(llr, 40)
    closest = 16
    best = None
    best_table = None
    for table in VITERBI_TABLES:
        body, distance = _score_message(viterbi_tailbiting(streams, table), phase)
        if distance < closest:
            closest = distance
            best = body
            best_table = table
        if distance == 0 and body is not None:
            return body, 0
    if closest <= 2 and best_table is not None:
        body, distance = _score_message(viterbi_exact(streams, best_table), phase)
        if distance == 0 and body is not None:
            return body, 0
        if distance < closest:
            closest = distance
    return None, closest


def _pbch_gold(pci: int):
    cached = _GOLD.get(("pbch", pci))
    if cached is None:
        cached = gold(pci, PBCH_BITS)
        _GOLD[("pbch", pci)] = cached
    return cached


def _descrambled_streams(soft, pci: int, phase: int):
    start = phase * (PBCH_BITS // FRAME_PHASES)
    scrambled = _pbch_gold(pci)[start:start + len(soft)]
    llr = [metric if bit == 0 else -metric for metric, bit in zip(soft, scrambled)]
    return conv_rate_recover(llr, 40)


def decode_mib_frames(frames, pci: int):
    """Soft-combine every PBCH occasion, then require an exact CRC-16.

    The same 1920 coded bits repeat every 40 ms. Frame i uses scramble phase
    (phase0 + i) mod 4, so frames 0, 4 and 8 are the same copy and their
    soft bits add. The capture may start on any phase, so all four alignments
    are tried. A one-bit checksum miss is not a decode.
    """
    closest = 16
    for n_rb in BANDWIDTHS:
        for ports in (2, 1, 4):
            accumulated = [None] * FRAME_PHASES
            found_any = False
            for index, spectra in enumerate(frames):
                if spectra is None:
                    continue
                soft = _pbch_soft(spectra, pci, n_rb, ports)
                if soft is None or len(soft) != 480:
                    continue
                slot = index % FRAME_PHASES
                vector = np.asarray(soft, dtype=np.float64)
                if accumulated[slot] is None:
                    accumulated[slot] = vector
                else:
                    accumulated[slot] = accumulated[slot] + vector
                found_any = True
            if not found_any:
                continue
            for phase0 in range(FRAME_PHASES):
                combined = [[0.0] * 40 for _ in range(3)]
                for slot in range(FRAME_PHASES):
                    if accumulated[slot] is None:
                        continue
                    phase = (phase0 + slot) % FRAME_PHASES
                    part = _descrambled_streams(accumulated[slot], pci, phase)
                    for stream in range(3):
                        for bit in range(40):
                            combined[stream][bit] += part[stream][bit]
                for table in VITERBI_TABLES:
                    body, distance = _score_message(viterbi_tailbiting(combined, table), phase0)
                    if distance < closest:
                        closest = distance
                    if _exact_mib(body, distance, ports, n_rb):
                        return body, 0
                    if distance <= 2:
                        body, distance = _score_message(viterbi_exact(combined, table), phase0)
                        if distance < closest:
                            closest = distance
                        if _exact_mib(body, distance, ports, n_rb):
                            return body, 0
    return None, closest


def _exact_mib(body, distance: int, ports: int, n_rb: int) -> bool:
    return (
        distance == 0
        and body is not None
        and body["antennas"] == ports
        and body["bandwidth_rb"] == n_rb
    )


def decode_mib_spectra(spectra, pci: int):
    """One PBCH occasion. Exact CRC only."""
    return decode_mib_frames([spectra], pci)


def _write_pbch(spectra, pci: int, n_rb: int, ports: int, payload, phase: int, h0=1 + 0j, h1=1 + 0j):
    parity = crc16(payload)
    masked = bits_to_int(parity) ^ ANTENNA_MASK[ports]
    streams = conv_encode(payload + bits_of(masked, 16))
    coded = conv_rate_match(streams, PBCH_BITS)
    scrambled = [bit ^ mask for bit, mask in zip(coded, gold(pci, PBCH_BITS))]
    start = phase * 480
    symbols = qpsk_mod(scrambled[start:start + 480])
    elements = pbch_data_res(n_rb, pci)
    if ports == 1:
        for (symbol, k), value in zip(elements, symbols):
            spectra[symbol][_k_to_bin(k, n_rb)] = h0 * value
    elif ports == 2:
        for index in range(0, len(symbols), 2):
            x0, x1 = symbols[index], symbols[index + 1]
            s0, k0 = elements[index]
            s1, k1 = elements[index + 1]
            spectra[s0][_k_to_bin(k0, n_rb)] = h0 * x0 + h1 * (-np.conj(x1))
            spectra[s1][_k_to_bin(k1, n_rb)] = h0 * x1 + h1 * np.conj(x0)
    else:
        h2, h3 = 0.2 - 0.5j, 0.6 + 0.1j
        for index in range(0, len(symbols), 4):
            block = symbols[index:index + 4]
            places = elements[index:index + 4]
            pairs = ((block[0], block[1], h0, h1), (block[2], block[3], h2, h3))
            for pair_index, (x0, x1, ha, hb) in enumerate(pairs):
                s0, k0 = places[pair_index * 2]
                s1, k1 = places[pair_index * 2 + 1]
                spectra[s0][_k_to_bin(k0, n_rb)] = ha * x0 + hb * (-np.conj(x1))
                spectra[s1][_k_to_bin(k1, n_rb)] = ha * x1 + hb * np.conj(x0)
        for k in _pbch_keys(n_rb):
            for port, gain in ((2, h2), (3, h3)):
                ref = crs_value(pci, 1, 1, k, port, n_rb)
                if ref is None:
                    continue
                spectra[1][_k_to_bin(k, n_rb)] += gain * ref
    for k in _pbch_keys(n_rb):
        for port, gain in ((0, h0), (1, h1)):
            if port >= ports:
                continue
            ref = crs_value(pci, 1, 0, k, port, n_rb)
            if ref is None:
                continue
            spectra[0][_k_to_bin(k, n_rb)] += gain * ref


def _ofdm(spec, symbol: int) -> np.ndarray:
    wave = np.fft.ifft(spec)
    length = cp_len(symbol)
    return np.concatenate((wave[-length:], wave))


def _symbol_fft(samples, useful_at: int):
    spec = np.fft.fft(samples[useful_at:useful_at + NFFT])
    # Bin 0 is the unused DC subcarrier. A HackRF spike sits here and leaks.
    spec[0] = 0.5 * (spec[1] + spec[-1])
    return spec


def _pss_match(samples, at: int, nid2: int) -> float:
    """How well a 128-sample window matches the known PSS, after the FFT."""
    if at < 0 or at + NFFT > len(samples):
        return -1.0
    from lte_sync import _pss_freq
    spec = np.fft.fft(samples[at:at + NFFT])
    observed = np.empty(62, dtype=np.complex128)
    observed[:31] = spec[-31:]
    observed[31:] = spec[1:32]
    reference = _pss_freq(nid2)
    return abs(np.vdot(reference, observed)) / (np.linalg.norm(observed) * np.sqrt(62) + 1e-12)


def _refine_pss(samples, guess: int, nid2: int) -> int:
    """Slide the useful-symbol start until the known PSS sequence matches.

    The raw correlation peak can sit several samples off the true symbol
    indoors. The cell search already does this slide. The broadcast reader
    has to use the same window or the pilots are in the wrong bins.
    """
    best_at = guess
    best = -1.0
    for delta in range(-12, 13):
        score = _pss_match(samples, guess + delta, nid2)
        if score > best:
            best = score
            best_at = guess + delta
    return best_at


def locate_subframe0(samples, pci: int, cfo_hz: float = 0.0):
    """Sample index where a subframe 0 starts, using one PSS/SSS correlation."""
    if cfo_hz:
        ramp = np.exp(-1j * 2 * np.pi * cfo_hz * np.arange(len(samples)) / FS)
        samples = samples * ramp
    nid2 = pci % 3
    template = pss_time(nid2)
    template = template / np.linalg.norm(template)
    size = len(samples) + len(template) - 1
    nfft = 1 << (size - 1).bit_length()
    corr = np.fft.ifft(np.fft.fft(samples, nfft) * np.conj(np.fft.fft(template, nfft)))
    mag = np.abs(corr[:len(samples)])
    order = np.argsort(mag)[::-1]
    nid1 = pci // 3
    seen = []
    for peak in order[:12]:
        peak = _refine_pss(samples, int(peak), nid2)
        if peak in seen:
            continue
        seen.append(peak)
        sss_at = peak - SYMBOL
        if sss_at < 0 or peak + NFFT > len(samples):
            continue
        spec = np.fft.fft(samples[sss_at:sss_at + NFFT])
        observed = np.empty(62, dtype=np.complex128)
        observed[:31] = spec[-31:]
        observed[31:] = spec[1:32]
        best = None
        best_score = -1.0
        for subframe5 in (False, True):
            ref = sss_freq(nid1, nid2, subframe5)
            score = abs(np.vdot(ref, observed)) / (np.linalg.norm(observed) * np.sqrt(62) + 1e-12)
            if score > best_score:
                best_score = score
                best = subframe5
        if best_score < 0.4:
            continue
        pss_in_subframe = symbol_useful_offset(6)
        if best:
            start = peak - pss_in_subframe - 5 * 1920
        else:
            start = peak - pss_in_subframe
        # Step back to the first subframe 0 still inside the buffer.
        while start >= 19200:
            start -= 19200
        if start < 0:
            start += 19200
        if start + 1920 <= len(samples):
            return int(start), samples
    return None, samples


def _spectra_at(samples, subframe_start: int, symbols):
    return [_symbol_fft(samples, subframe_start + symbol_useful_offset(symbol)) for symbol in symbols]


def _pdsch_elements(n_rb: int, pci: int, cfi: int):
    """Data resource elements of subframe 5, k fastest, sync and CRS removed."""
    low, high = _rs_residues(pci)
    dc = n_rb * 6
    elements = []
    for symbol in range(cfi, 14):
        for k in range(n_rb * 12):
            if symbol in (0, 4, 7, 11) and (k % 6) in (low, high):
                continue
            if symbol in (5, 6) and 1 <= abs(k - dc + 0.5) <= 31.5:
                # Centre 62 subcarriers. abs(k - (dc - 0.5)) does that.
                if dc - 31 <= k <= dc + 30:
                    continue
            elements.append((symbol, k))
    return elements


def _build_sib1_bits(cell_id: int, mcc=(2, 3, 4), mnc=(1, 5), tac: int = 1) -> list[int]:
    bits = bits_of(0, 3)  # one PLMN in the list
    bits.append(1)  # MCC is present
    for digit in mcc:
        bits += bits_of(digit, 4)
    bits.append(0 if len(mnc) == 2 else 1)
    for digit in mnc:
        bits += bits_of(digit, 4)
    bits.append(0)  # not reserved
    bits += bits_of(tac, 16)
    bits += bits_of(cell_id, 28)
    return bits


def _parse_sib1(bits):
    """Return the 28-bit cell identity, trying a leading extension bit and not."""
    for start in (0, 1):
        found = _parse_sib1_at(bits, start)
        if found is not None:
            return found
    return None


def _parse_sib1_at(bits, start: int):
    cursor = start

    def take(count):
        nonlocal cursor
        if cursor + count > len(bits):
            raise ValueError("short")
        chunk = bits[cursor:cursor + count]
        cursor += count
        return chunk

    try:
        count = bits_to_int(take(3)) + 1
        if not 1 <= count <= 6:
            return None
        for _ in range(count):
            if take(1)[0]:
                for _digit in range(3):
                    if bits_to_int(take(4)) > 9:
                        return None
            wide = take(1)[0]
            for _digit in range(3 if wide else 2):
                if bits_to_int(take(4)) > 9:
                    return None
            take(1)
        take(16)
        return bits_to_int(take(28))
    except ValueError:
        return None


def _encode_transport(info_bits, tbs: int):
    padded = list(info_bits) + [0] * (tbs - len(info_bits))
    crc = crc24a(padded)
    block = padded + crc
    size = code_block_size(len(block))
    filler = size - len(block)
    streams = turbo_encode([0] * filler + block)
    streams[0] = [None] * filler + streams[0][filler:]
    return streams, size, filler


def _decode_transport(soft, tbs: int):
    size = code_block_size(tbs + 24)
    filler = size - (tbs + 24)
    streams = turbo_rate_recover(soft, size, filler=filler)
    for index in range(filler):
        streams[0][index] = 25.0
    decoded = turbo_decode(streams)
    block = decoded[filler:]
    payload, parity = block[:tbs], block[tbs:tbs + 24]
    if parity != crc24a(payload):
        return None
    return payload


def _pdsch_soft(samples, subframe_start: int, pci: int, n_rb: int, cfi: int):
    spectra = _spectra_at(samples, subframe_start, range(14))
    # One gain from the CRS on symbol 0, enough for a flat capture.
    gains = []
    for k in range(n_rb * 12):
        ref = crs_value(pci, 10, 0, k, 0, n_rb)
        if ref is None:
            continue
        gains.append(spectra[0][_k_to_bin(k, n_rb)] * np.conj(ref))
    gain = np.mean(gains) if gains else 1 + 0j
    soft = []
    for symbol, k in _pdsch_elements(n_rb, pci, cfi):
        value = spectra[symbol][_k_to_bin(k, n_rb)] * np.conj(gain)
        soft.extend(qpsk_soft(value / (np.abs(gain) ** 2 + 1e-12)))
    n_s = 10  # first slot of subframe 5
    c_init = SI_RNTI * 2 ** 14 + (n_s // 2) * 2 ** 9 + pci
    mask = gold(c_init, len(soft))
    return [metric if bit == 0 else -metric for metric, bit in zip(soft, mask)]


def _write_pdsch(samples, subframe_start: int, pci: int, n_rb: int, cfi: int, info_bits, tbs: int):
    streams, _size, _filler = _encode_transport(info_bits, tbs)
    elements = _pdsch_elements(n_rb, pci, cfi)
    coded = turbo_rate_match(streams, len(elements) * 2)
    n_s = 10
    c_init = SI_RNTI * 2 ** 14 + (n_s // 2) * 2 ** 9 + pci
    scrambled = [bit ^ mask for bit, mask in zip(coded, gold(c_init, len(coded)))]
    symbols = qpsk_mod(scrambled)
    spectra = [np.zeros(NFFT, dtype=np.complex128) for _ in range(14)]
    for (symbol, k), value in zip(elements, symbols):
        spectra[symbol][_k_to_bin(k, n_rb)] = value
    for k in range(n_rb * 12):
        ref = crs_value(pci, 10, 0, k, 0, n_rb)
        if ref is not None:
            spectra[0][_k_to_bin(k, n_rb)] += ref
        ref4 = crs_value(pci, 10, 4, k, 0, n_rb)
        if ref4 is not None:
            spectra[4][_k_to_bin(k, n_rb)] += ref4
    for symbol, spec in enumerate(spectra):
        wave = _ofdm(spec, symbol)
        at = subframe_start + symbol_useful_offset(symbol) - cp_len(symbol)
        samples[at:at + len(wave)] += wave


def _capture_covers(sample_rate: float, n_rb: int) -> bool:
    """True when this sample rate can see every subcarrier of the cell."""
    occupied = n_rb * 12 * 15_000
    return float(sample_rate) > occupied * 1.05


def _cp_cfo(samples, sample_rate: float) -> float:
    """Hertz offset from the cyclic-prefix correlation, or 0 when it is not clear."""
    nfft = 128 if abs(sample_rate - FS) < 1 else int(round(sample_rate / 15_000))
    if nfft < 128:
        return 0.0
    cp = 9 * nfft // 128
    step = nfft + cp
    acc = 0j
    used = 0
    limit = min(len(samples) - step, step * 14)
    for origin in range(0, limit, step):
        acc += np.vdot(samples[origin:origin + cp], samples[origin + nfft:origin + nfft + cp])
        used += 1
    if used == 0 or abs(acc) <= 0:
        return 0.0
    power = float(np.mean(np.abs(samples[:limit]) ** 2)) + 1e-12
    if abs(acc) / used < power * cp * 0.2:
        return 0.0
    hertz = float(np.angle(acc) * sample_rate / (2 * np.pi * nfft))
    if abs(hertz) < 300 or abs(hertz) > 7500:
        return 0.0
    return hertz


def _rate_nfft(sample_rate: float):
    """FFT size when this rate lands on the 15 kHz LTE grid. None for 1.92 MS/s."""
    bins = float(sample_rate) / 15_000.0
    nfft = int(round(bins))
    if nfft not in (256, 512, 1024, 1536, 2048):
        return None
    if abs(bins - nfft) > 0.02:
        return None
    return nfft


def _cp_wide(symbol: int, nfft: int) -> int:
    return (10 if symbol % 7 == 0 else 9) * (nfft // 128)


def _useful_wide(symbol: int, nfft: int) -> int:
    cursor = 0
    for index in range(symbol):
        cursor += _cp_wide(index, nfft) + nfft
    return cursor + _cp_wide(symbol, nfft)


def _subframe_wide(nfft: int) -> int:
    return sum(_cp_wide(symbol, nfft) + nfft for symbol in range(14))


def _center_seq(samples, at: int, nfft: int):
    if at < 0 or at + nfft > len(samples):
        return None
    spec = np.fft.fft(samples[at:at + nfft])
    observed = np.empty(62, dtype=np.complex128)
    observed[:31] = spec[-31:]
    observed[31:] = spec[1:32]
    return observed


def _seq_score(observed, reference) -> float:
    if observed is None:
        return -1.0
    return abs(np.vdot(reference, observed)) / (np.linalg.norm(observed) * np.sqrt(62) + 1e-12)


def _refine_wide(samples, guess: int, reference, nfft: int) -> int:
    """Slide across a few samples at this rate until the PSS sequence matches."""
    step = max(1, nfft // 128)
    best_at = guess
    best = -1.0
    for delta in range(-12 * step, 12 * step + 1):
        score = _seq_score(_center_seq(samples, guess + delta, nfft), reference)
        if score > best:
            best = score
            best_at = guess + delta
    return best_at


def _locate_wide(samples, pci: int, nfft: int):
    """Subframe-0 start at the capture's own sample rate."""
    from lte_sync import _pss_freq
    nid2 = int(pci) % 3
    nid1 = int(pci) // 3
    reference = _pss_freq(nid2)
    spec = np.zeros(nfft, dtype=np.complex128)
    spec[-31:] = reference[:31]
    spec[1:32] = reference[31:]
    template = np.fft.ifft(spec)
    template = template / np.linalg.norm(template)
    size = len(samples) + nfft - 1
    fft_len = 1 << (size - 1).bit_length()
    corr = np.fft.ifft(np.fft.fft(samples, fft_len) * np.conj(np.fft.fft(template, fft_len)))
    mag = np.abs(corr[:len(samples)])
    pss_at = _useful_wide(6, nfft)
    gap = nfft + _cp_wide(6, nfft)
    sub = _subframe_wide(nfft)
    frame = 10 * sub
    seen = []
    for peak in np.argsort(mag)[::-1][:8]:
        peak = _refine_wide(samples, int(peak), reference, nfft)
        if peak in seen:
            continue
        seen.append(peak)
        observed = _center_seq(samples, peak - gap, nfft)
        if observed is None:
            continue
        best_5 = False
        best_score = -1.0
        for subframe5 in (False, True):
            score = _seq_score(observed, sss_freq(nid1, nid2, subframe5))
            if score > best_score:
                best_score = score
                best_5 = subframe5
        if best_score < 0.4:
            continue
        start = peak - pss_at - (5 * sub if best_5 else 0)
        while start >= frame:
            start -= frame
        if start < 0:
            start += frame
        if start + sub <= len(samples):
            return int(start)
    return None


def _spectra_wide(samples, subframe_start: int, nfft: int):
    spectra = []
    for symbol in (7, 8, 9, 10):
        at = subframe_start + _useful_wide(symbol, nfft)
        if at < 0 or at + nfft > len(samples):
            return None
        spec = np.fft.fft(samples[at:at + nfft])
        spec[0] = 0.5 * (spec[1] + spec[-1])
        spectra.append(spec)
    return spectra


def _decode_mib_native(samples, sample_rate: float, pci: int):
    """MIB from the capture's own 15 kHz grid, without shrinking to 1.92 MS/s."""
    nfft = _rate_nfft(sample_rate)
    if nfft is None:
        return None, 16, None
    start = _locate_wide(samples, pci, nfft)
    if start is None:
        return None, 16, None
    sub = _subframe_wide(nfft)
    frame = 10 * sub
    pss_at = _useful_wide(6, nfft)
    from lte_sync import _pss_freq
    reference = _pss_freq(int(pci) % 3)
    grouped = []
    cursor = start
    for _ in range(12):
        if cursor + frame > len(samples) and cursor + pss_at + nfft > len(samples):
            break
        guess = cursor + pss_at
        if guess < 0 or guess + nfft > len(samples):
            grouped.append(None)
        else:
            aligned = _refine_wide(samples, guess, reference, nfft) - pss_at
            grouped.append(_spectra_wide(samples, aligned, nfft))
        cursor += frame
        if cursor + pss_at + nfft > len(samples):
            break
    found, distance = decode_mib_frames(grouped, pci)
    return found, distance, start


def _derotate(samples, cfo_hz: float, sample_rate: float):
    """Remove a constant Hertz offset from the whole capture."""
    if abs(cfo_hz) < 1:
        return samples
    index = np.arange(len(samples))
    ramp = np.exp(-1j * 2 * np.pi * float(cfo_hz) * index / float(sample_rate))
    return samples * ramp


def decode_cell_identity(samples, pci: int, cfo_hz: float = 0.0, sample_rate: float = FS) -> dict:
    """Decode one capture. pci is already known. cfo_hz is the residual tuner offset."""
    wide = np.asarray(samples, dtype=np.complex128)
    rate = float(sample_rate)
    if abs(cfo_hz) < 1:
        cfo_hz = _cp_cfo(wide, rate)
    wide = _derotate(wide, cfo_hz, rate)
    samples = wide
    if abs(rate - FS) > 1:
        from lte_pci_search import resample_to_lte
        samples = resample_to_lte(wide, rate).astype(np.complex128)
    start, samples = locate_subframe0(samples, pci, 0.0)
    result = {
        "pci": int(pci),
        "cell_id": None,
        "enb_id": None,
        "sector_id": None,
        "bandwidth_rb": None,
        "bandwidth_mhz": None,
        "phich_duration": None,
        "phich_resource": None,
        "sfn": None,
        "antennas": None,
        "cell_id_status": "no LTE frame timing for PCI %d" % pci,
    }
    mib = None
    mib_start = None
    closest = 16
    if _rate_nfft(rate) is not None:
        found, distance, _wide_start = _decode_mib_native(wide, rate, int(pci))
        closest = distance
        if found is not None:
            mib = found
    if mib is None and start is None:
        return result
    if mib is None and start is not None:
        need = symbol_useful_offset(10) + NFFT
        frames = []
        cursor = start
        while cursor + need <= len(samples):
            if cursor >= 0:
                frames.append(cursor)
            cursor += 19200
        # 120 ms is 12 frames: three repeats of the 40 ms PBCH word.
        frames = frames[:12]
        nid2 = int(pci) % 3
        pss_in_subframe = symbol_useful_offset(6)
        grouped = []
        aligned = []
        for origin in frames:
            guess = origin + pss_in_subframe
            if guess < 0 or guess + NFFT > len(samples):
                grouped.append(None)
                aligned.append(origin)
                continue
            pss_at = _refine_pss(samples, guess, nid2)
            aligned_start = pss_at - pss_in_subframe
            aligned.append(aligned_start)
            if aligned_start < 0 or aligned_start + need > len(samples):
                grouped.append(None)
            else:
                grouped.append(_spectra_at(samples, aligned_start, (7, 8, 9, 10)))
        found, distance = decode_mib_frames(grouped, pci)
        closest = distance
        if found is not None:
            mib = found
            mib_start = aligned[0]
    if mib is None:
        result["cell_id_status"] = "MIB did not decode (timing locked, best checksum distance %d)" % closest
        return result
    result.update(mib)
    n_rb = mib["bandwidth_rb"]
    if n_rb > 6:
        if not _capture_covers(rate, n_rb):
            result["cell_id_status"] = (
                "SIB1 is spread across %s MHz. This capture only covers the centre 1.4 MHz, so the cell identity is not in it."
                % mib["bandwidth_mhz"]
            )
            return result
        from lte_sib1 import decode_wideband_sib1
        identity = decode_wideband_sib1(wide, rate, int(pci), mib, 0.0)
        result.update(identity)
        if result.get("cell_id") is None and not result.get("cell_id_status"):
            result["cell_id_status"] = "SIB1 did not decode"
        return result
    if mib["sfn"] % 2:
        result["cell_id_status"] = "this frame is odd, SIB1 is sent on even frames"
        return result
    sib_start = mib_start + 5 * 1920
    if sib_start + 1920 > len(samples):
        result["cell_id_status"] = "the capture ends before SIB1"
        return result
    for cfi in (1, 2, 3):
        soft = _pdsch_soft(samples, sib_start, pci, n_rb, cfi)
        for tbs in _TBS_CANDIDATES:
            if tbs + 24 > code_block_size(10 ** 9) if False else tbs > 700:
                continue
            try:
                if len(soft) < tbs or tbs > 776:
                    continue
                payload = _decode_transport(soft, tbs)
            except Exception:
                continue
            if payload is None:
                continue
            cell_id = _parse_sib1(payload)
            if cell_id is None:
                continue
            result["cell_id"] = cell_id
            result["enb_id"] = cell_id >> 8
            result["sector_id"] = cell_id & 0xFF
            result["cell_id_status"] = "decoded"
            return result
    result["cell_id_status"] = "SIB1 did not decode"
    return result


def synthesize_cell(pci: int, cell_id: int, n_rb: int = 6, ports: int = 1) -> np.ndarray:
    """One 10 ms frame at 1.92 MHz with a PBCH and, for a 1.4 MHz cell, a SIB1."""
    samples = np.zeros(19200, dtype=np.complex128)
    spectra = [np.zeros(NFFT, dtype=np.complex128) for _ in range(4)]
    h1 = 0.4 + 0.3j if ports == 2 else 0j
    _write_pbch(spectra, pci, n_rb, ports, _pack_mib(n_rb, sfn=8), phase=0, h0=1 + 0j, h1=h1)
    for index, spec in enumerate(spectra):
        symbol = 7 + index
        wave = _ofdm(spec, symbol)
        at = symbol_useful_offset(symbol) - cp_len(symbol)
        samples[at:at + len(wave)] += wave
    # PSS and SSS so frame timing can be found.
    nid2 = pci % 3
    nid1 = pci // 3
    for subframe, flag in ((0, False), (5, True)):
        base = subframe * 1920
        for symbol, freq in ((5, sss_freq(nid1, nid2, flag)), (6, None)):
            spec = np.zeros(NFFT, dtype=np.complex128)
            if freq is None:
                from lte_sync import _zc63, PSS_ROOTS
                zc = _zc63(PSS_ROOTS[nid2])
                freq = np.concatenate((zc[:31], zc[32:]))
            spec[-31:] = freq[:31]
            spec[1:32] = freq[31:]
            wave = _ofdm(np.fft.ifft(np.fft.fft(spec)), symbol) if False else _ofdm(spec, symbol)
            at = base + symbol_useful_offset(symbol) - cp_len(symbol)
            samples[at:at + len(wave)] += wave
    if n_rb == 6:
        info = _build_sib1_bits(cell_id)
        _write_pdsch(samples, 5 * 1920, pci, n_rb, 2, info, 72)
    return samples


def format_result(result: dict) -> str:
    lines = ["PCI %d" % result["pci"]]
    if result.get("bandwidth_rb"):
        lines.append(
            "MIB: %d RB (%s MHz), PHICH %s Ng=%s, SFN %s, antennas %s"
            % (
                result["bandwidth_rb"],
                result["bandwidth_mhz"],
                result["phich_duration"],
                result["phich_resource"],
                result["sfn"],
                result["antennas"],
            )
        )
    if result.get("cell_id") is not None:
        lines.append("CellID: %d" % result["cell_id"])
        lines.append("eNB ID: %d" % result["enb_id"])
        lines.append("local cell: %d" % result["sector_id"])
        if result.get("mcc") is not None:
            lines.append(
                "PLMN: %s-%s %s  TAC: %s"
                % (result["mcc"], result["mnc"], result.get("operator") or "", result.get("tac"))
            )
    else:
        lines.append("CellID: " + result["cell_id_status"])
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Decode the LTE cell identity from one Kraken capture.")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--pci", type=int, help="physical cell ID already found by lte_pci_search")
    parser.add_argument("--cfo-khz", type=float, default=0.0)
    parser.add_argument("--file", help="complex64 capture at --rate, instead of the live Kraken")
    parser.add_argument("--rate", type=float, default=2_400_000.0)
    args = parser.parse_args(argv)
    if args.self_test:
        wide = synthesize_cell(430, 7159822, n_rb=50, ports=2)
        narrow = synthesize_cell(430, 7159822, n_rb=6, ports=1)
        mib = decode_cell_identity(wide, 430)
        full = decode_cell_identity(narrow, 430)
        ok = (
            mib["bandwidth_rb"] == 50
            and mib["cell_id"] is None
            and full["cell_id"] == 7159822
            and full["enb_id"] == (7159822 >> 8)
        )
        print("self-test %s" % ("OK" if ok else "FAILED"))
        print(format_result(full))
        return 0 if ok else 1
    if args.pci is None:
        print("Pass --pci from the direction finder. This decoder does not search for it.", file=sys.stderr)
        return 1
    print("Reading one capture for PCI %d" % args.pci, flush=True)
    if args.file:
        samples = np.fromfile(args.file, dtype=np.complex64)
        rate = args.rate
    else:
        from lte_pci_search import capture_channel0
        samples, rate, _centre = capture_channel0()
    print("Decoding the MIB", flush=True)
    result = decode_cell_identity(samples, args.pci, cfo_hz=args.cfo_khz * 1000.0, sample_rate=rate)
    print(format_result(result))
    return 0 if result.get("bandwidth_rb") else 1


if __name__ == "__main__":
    sys.exit(main())
