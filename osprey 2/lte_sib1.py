"""Wideband SIB1 demodulator and the TS 36.331 UPER fields it carries.

The narrow MIB path stays in lte_cellid. This module reads subframe 5 of an
even frame from a capture that actually covers the cell: 7.68 MS/s for 5 MHz,
15.36 MS/s for 10 and 15 MHz, and 30.72 MS/s for 20 MHz.
"""

from __future__ import annotations

import math

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
    viterbi_tailbiting,
)
from lte_sync import pss_time, sss_freq

SI_RNTI = 65535


def tbs_from_mcs(n_prb: int, mcs: int):
    """QPSK grant. SIB1 is QPSK, so MCS above 9 is refused."""
    if not 0 <= mcs <= 9 or not 1 <= int(n_prb) <= 100:
        return None
    return _TBS_ROWS[mcs][int(n_prb) - 1]


OPERATORS = {
    (234, 15): "Vodafone",
    (234, 10): "O2",
    (234, 30): "EE",
    (234, 33): "EE",
    (234, 20): "Three",
}


def _nfft_for(sample_rate: float, n_rb: int):
    """Return (samples, rate, nfft) on a 15 kHz grid that contains the cell."""
    occupied = n_rb * 12 * 15_000
    rate = float(sample_rate)
    if n_rb >= 100 and rate + 1e3 < 30.72e6:
        if rate < 19e6:
            return None
        from scipy.signal import resample_poly
        # 19.2 MS/s * 8/5 = 30.72 MS/s, the native 20 MHz clock.
        return None  # placeholder replaced below
    nfft = int(round(rate / 15_000))
    if nfft not in (128, 256, 512, 1024, 1536, 2048):
        return None
    if occupied >= rate * 0.95:
        return None
    return nfft


def prepare_wide(samples, sample_rate: float, n_rb: int):
    """Samples on an integer 15 kHz FFT, or None when the cell does not fit."""
    samples = np.asarray(samples)
    rate = float(sample_rate)
    if n_rb >= 100 and rate < 30e6:
        if rate < 19e6:
            return None
        from scipy.signal import resample_poly
        converted = resample_poly(samples, 8, 5).astype(np.complex64)
        return converted, 30.72e6, 2048
    nfft = _nfft_for(rate, n_rb)
    if nfft is None:
        return None
    return np.asarray(samples, dtype=np.complex64), rate, nfft


def cp_len(symbol: int, nfft: int) -> int:
    scale = nfft // 128
    return 10 * scale if symbol % 7 == 0 else 9 * scale


def useful_offset(symbol: int, nfft: int) -> int:
    cursor = 0
    for index in range(symbol):
        cursor += cp_len(index, nfft) + nfft
    return cursor + cp_len(symbol, nfft)


def subframe_len(nfft: int) -> int:
    return useful_offset(14, nfft) - cp_len(0, nfft) + 0 if False else _subframe_len(nfft)


def _subframe_len(nfft: int) -> int:
    total = 0
    for symbol in range(14):
        total += cp_len(symbol, nfft) + nfft
    return total


def _k_to_bin(k: int, n_rb: int) -> int:
    delta = k - n_rb * 6
    if delta < 0:
        return delta
    return delta + 1


def _symbol_fft(samples, origin: int, symbol: int, nfft: int):
    at = origin + useful_offset(symbol, nfft)
    spec = np.fft.fft(samples[at:at + nfft])
    spec[0] = 0.5 * (spec[1] + spec[-1])
    return spec


def _ofdm_symbol(spec, symbol: int, nfft: int):
    wave = np.fft.ifft(spec)
    length = cp_len(symbol, nfft)
    return np.concatenate((wave[-length:], wave))


def _write_symbol(buffer, origin, symbol, spec, nfft):
    wave = _ofdm_symbol(spec, symbol, nfft)
    at = origin + useful_offset(symbol, nfft) - cp_len(symbol, nfft)
    buffer[at:at + len(wave)] += wave


def _crs_init(pci: int, n_s: int, symbol: int) -> int:
    return (2 ** 10) * (7 * (n_s + 1) + symbol + 1) * (2 * pci + 1) + 2 * pci + 1


def _gold_bits(c_init: int, length: int):
    return gold(c_init, length)


def crs_symbol(pci: int, n_s: int, symbol: int, k: int, n_rb: int, port: int = 0):
    """CRS value, or None when this port has no pilot on the subcarrier."""
    if port == 0:
        if symbol not in (0, 4):
            return None
        v = 0 if symbol == 0 else 3
    elif port == 1:
        if symbol not in (0, 4):
            return None
        v = 3 if symbol == 0 else 0
    else:
        return None
    residue = (v + (pci % 6)) % 6
    if k % 6 != residue:
        return None
    m = (k - residue) // 6
    if m < 0:
        return None
    # Same 110-block centre index as lte_cellid.crs_value.
    m_prime = m + 110 - int(n_rb)
    sequence = _gold_bits(_crs_init(pci, n_s, symbol), 2 * (m_prime + 1))
    scale = 1.0 / (2 ** 0.5)
    return scale * ((1 - 2 * sequence[2 * m_prime]) + 1j * (1 - 2 * sequence[2 * m_prime + 1]))


def _regs_symbol0(n_rb: int, pci: int):
    """REG index -> four (symbol, k) resource elements. Two REGs per resource block."""
    low = pci % 6
    high = (low + 3) % 6
    regs = []
    for rb in range(n_rb):
        tones = [rb * 12 + index for index in range(12) if index % 6 not in (low, high)]
        regs.append([(0, tone) for tone in tones[:4]])
        regs.append([(0, tone) for tone in tones[4:]])
    return regs


def _regs_plain(n_rb: int, symbol: int):
    regs = []
    for rb in range(n_rb):
        tones = [rb * 12 + index for index in range(12)]
        for start in (0, 4, 8):
            regs.append([(symbol, tone) for tone in tones[start:start + 4]])
    return regs


def pcfich_reg_index(n_rb: int, pci: int):
    count = 2 * n_rb
    base = pci % count
    return [(base + (i * n_rb) // 2) % count for i in range(4)]


def phich_group_count(n_rb: int, ng: str) -> int:
    value = {"1/6": 1.0 / 6.0, "1/2": 0.5, "1": 1.0, "2": 2.0}.get(ng, 1.0)
    return int(math.ceil(value * n_rb / 8.0))


def phich_reg_index(n_rb: int, pci: int, groups: int):
    count = 2 * n_rb
    found = []
    for group in range(groups):
        for piece in range(3):
            found.append((pci + group + piece * count // 3) % count)
    return sorted(set(found))


def _pdcch_regs(n_rb: int, pci: int, cfi: int, ng: str):
    """REGs left for the PDCCH after PCFICH and PHICH, in mapping order."""
    symbol0 = _regs_symbol0(n_rb, pci)
    taken = set(pcfich_reg_index(n_rb, pci))
    taken.update(phich_reg_index(n_rb, pci, phich_group_count(n_rb, ng)))
    regs = [reg for index, reg in enumerate(symbol0) if index not in taken]
    for symbol in range(1, cfi):
        regs.extend(_regs_plain(n_rb, symbol))
    return _interleave_regs(regs, pci)


def _interleave_regs(regs, pci: int):
    """Sub-block interleave then the PCI cyclic shift. TS 36.212 5.1.4.2.1."""
    from lte_coding import _CONV_COLUMNS
    count = len(regs)
    rows = (count + 31) // 32
    padded = rows * 32
    matrix = [None] * padded
    holes = padded - count
    for index, reg in enumerate(regs):
        matrix[holes + index] = reg
    reading = []
    for column in _CONV_COLUMNS:
        for row in range(rows):
            item = matrix[row * 32 + column]
            if item is not None:
                reading.append(item)
    shift = pci % len(reading) if reading else 0
    return reading[shift:] + reading[:shift]


def _scramble(bits, c_init: int):
    mask = gold(c_init, len(bits))
    return [bit ^ scramble for bit, scramble in zip(bits, mask)]


def _descramble_soft(soft, c_init: int):
    mask = gold(c_init, len(soft))
    return [metric if bit == 0 else -metric for metric, bit in zip(soft, mask)]


def _dci_crc(bits):
    parity = bits_to_int(crc16(bits)) ^ SI_RNTI
    return bits + bits_of(parity, 16)


def _dci_crc_ok(bits):
    if len(bits) < 17:
        return False
    payload, parity = bits[:-16], bits[-16:]
    expect = bits_to_int(crc16(payload)) ^ SI_RNTI
    return bits_to_int(parity) == expect


def _riv_bits(n_rb: int) -> int:
    return int(math.ceil(math.log2(n_rb * (n_rb + 1) / 2.0)))


def encode_riv(n_rb: int, start: int, length: int) -> int:
    if (length - 1) <= (n_rb // 2):
        return n_rb * (length - 1) + start
    return n_rb * (n_rb - length + 1) + (n_rb - 1 - start)


def decode_riv(n_rb: int, riv: int):
    for length in range(1, n_rb + 1):
        for start in range(0, n_rb - length + 1):
            if encode_riv(n_rb, start, length) == riv:
                return start, length
    return None


def build_dci_1a(n_rb: int, start: int, length: int, mcs: int, rv: int = 0) -> list[int]:
    """Localized DCI format 1A for an FDD SI-RNTI grant."""
    bits = [1, 0]
    bits += bits_of(encode_riv(n_rb, start, length), _riv_bits(n_rb))
    bits += bits_of(mcs, 5)
    bits += bits_of(0, 3)  # HARQ process
    bits += [0]  # new data indicator
    bits += bits_of(rv, 2)
    bits += bits_of(0, 2)  # TPC, unused for SI-RNTI
    return _dci_crc(bits)


def _parse_dci_1a(bits, n_rb: int):
    if not _dci_crc_ok(bits):
        return None
    body = bits[:-16]
    need = 1 + 1 + _riv_bits(n_rb) + 5 + 3 + 1 + 2 + 2
    if len(body) < need or body[0] != 1 or body[1] != 0:
        return None
    cursor = 2
    width = _riv_bits(n_rb)
    riv = bits_to_int(body[cursor:cursor + width])
    cursor += width
    place = decode_riv(n_rb, riv)
    if place is None:
        return None
    mcs = bits_to_int(body[cursor:cursor + 5])
    cursor += 5 + 3 + 1
    rv = bits_to_int(body[cursor:cursor + 2])
    start, length = place
    return {"start": start, "length": length, "mcs": mcs, "rv": rv, "format": "1A"}


def _pdcch_bits(info_bits, encoded_len: int):
    coded = conv_rate_match(conv_encode(info_bits), encoded_len)
    return coded


def _map_qpsk(regs, bits):
    symbols = qpsk_mod(bits)
    mapping = []
    for reg in regs:
        for element in reg:
            if not symbols:
                return mapping
            mapping.append((element, symbols.pop(0)))
    return mapping


def _control_symbols(n_rb: int, pci: int, cfi: int, ng: str, dci_bits, aggregation: int):
    regs = _pdcch_regs(n_rb, pci, cfi, ng)
    cces = len(regs) // 9
    if aggregation > cces:
        return None
    # Common search space candidate m = 0.
    chosen = regs[: 9 * aggregation]
    coded = _pdcch_bits(dci_bits, 8 * len(chosen) * 1)  # 4 REs * 2 bits = 8 bits per REG
    # 9 REGs * 4 REs * 2 bits = 72 bits per CCE.
    coded = _pdcch_bits(dci_bits, 72 * aggregation)
    c_init = 0 * 512 + pci  # overwritten by the caller for the subframe
    return chosen, coded


def _paint(spec, n_rb, pairs):
    for (symbol, k), value in pairs:
        if symbol != spec[1]:
            continue
        spec[0][_k_to_bin(k, n_rb)] += value


def sib1_rv(sfn: int) -> int:
    """Redundancy version of SIB1. k = (SFN/2) mod 4, rv = ceil(3k/2) mod 4."""
    k = (int(sfn) // 2) % 4
    return int(math.ceil(1.5 * k)) % 4


def _pdsch_elements(n_rb: int, pci: int, cfi: int, start: int, length: int):
    low = pci % 6
    high = (low + 3) % 6
    elements = []
    for symbol in range(cfi, 14):
        for rb in range(start, start + length):
            for index in range(12):
                k = rb * 12 + index
                if symbol in (0, 4, 7, 11) and index % 6 in (low, high):
                    continue
                elements.append((symbol, k))
    return elements


def _transport_streams(info_bits, tbs: int):
    padded = list(info_bits) + [0] * (tbs - len(info_bits))
    block = padded + crc24a(padded)
    size = code_block_size(len(block))
    filler = size - len(block)
    streams = turbo_encode([0] * filler + block)
    streams = [list(item) for item in streams]
    streams[0] = [None] * filler + streams[0][filler:]
    return streams, size, filler


def _recover_transport(soft, tbs: int, rv: int):
    size = code_block_size(tbs + 24)
    filler = size - (tbs + 24)
    streams = turbo_rate_recover(soft, size, redundancy=rv, filler=filler)
    decoded = turbo_decode(streams)
    block = decoded[filler:]
    payload, parity = block[:tbs], block[tbs:tbs + 24]
    if parity != crc24a(payload):
        return None
    return payload


def build_sib1_uper(cell_id: int, mcc=(2, 3, 4), mnc=(1, 5), tac: int = 1) -> list[int]:
    """R8 SystemInformationBlockType1, enough of the front to carry the identity.

    The leading choices are BCCH-DL-SCH-Message -> c1 -> systemInformationBlockType1.
    Optional bits that we leave absent are written as zero so a live SIB1 with
    the same header still lines up on the PLMN, TAC and cellIdentity.
    """
    bits = []

    def put(value, width):
        bits.extend(bits_of(value, width))

    put(0, 1)  # BCCH-DL-SCH-Message: c1
    put(1, 1)  # c1: systemInformationBlockType1
    put(0, 4)  # q-RxLevMinOffset, p-Max, tdd-Config, nonCriticalExtension absent
    put(0, 1)  # csg-Identity absent
    put(0, 3)  # one PLMN (length 1..6 encoded as count-1)
    put(1, 1)  # MCC present
    for digit in mcc:
        put(digit, 4)
    put(0 if len(mnc) == 2 else 1, 1)
    for digit in mnc:
        put(digit, 4)
    put(1, 1)  # cellReservedForOperatorUse = notReserved
    put(tac & 0xFFFF, 16)
    put(cell_id & 0x0FFFFFFF, 28)
    put(1, 1)  # cellBarred = notBarred
    put(0, 1)  # intraFreqReselection = allowed
    put(0, 1)  # csg-Indication = false
    return bits


def identity_from_transport_block(raw: bytes) -> dict | None:
    """Read a decoded SIB1 transport block and return the Vodafone-style identity.

    srsRAN writes the block as packed bytes. The country code 234 is the
    landmark. The first match is the network. A later copy of those bits
    usually sits inside the tracking-area field and is not a second cell.
    """
    if not raw:
        return None
    if max(raw) <= 1:
        bits = [byte & 1 for byte in raw]
    else:
        bits = []
        for byte in raw[:32]:
            for shift in range(7, -1, -1):
                bits.append((byte >> shift) & 1)
    text = "".join(str(bit) for bit in bits)
    at = text.find("001000110100")
    if at < 0:
        return None
    rest = text[at + 12:]
    wide = rest[:1] == "1"
    width = 3 if wide else 2
    mnc_bits = rest[1:1 + 4 * width]
    if len(mnc_bits) < 4 * width:
        return None
    digits = []
    for index in range(0, 4 * width, 4):
        digit = int(mnc_bits[index:index + 4], 2)
        if digit > 9:
            return None
        digits.append(digit)
    after = rest[1 + 4 * width + 1:]
    if len(after) < 44:
        return None
    tac = int(after[:16], 2)
    cell = int(after[16:44], 2)
    mnc = int("".join(str(digit) for digit in digits))
    return {
        "mcc": 234,
        "mnc": mnc,
        "mnc_digits": "".join(str(digit) for digit in digits),
        "tac": tac,
        "cell_id": cell,
        "enb_id": cell >> 8,
        "sector_id": cell & 0xFF,
        "operator": OPERATORS.get((234, mnc), ""),
    }


def parse_sib1_uper(bits) -> dict | None:
    """Read MCC, MNC, TAC and the 28-bit cell identity. None when the header is not SIB1."""
    for extension in (0, 1):
        found = _parse_uper_at(bits, extension)
        if found is not None:
            return found
    return None


def _parse_uper_at(bits, extension_bits: int):
    cursor = 0

    def take(width):
        nonlocal cursor
        if cursor + width > len(bits):
            raise ValueError("short")
        value = bits_to_int(bits[cursor:cursor + width])
        cursor += width
        return value

    try:
        if take(1) != 0 or take(1) != 1:
            return None
        if extension_bits:
            take(1)
        take(4)  # outer optional preamble
        take(1)  # csg-Identity
        count = take(3) + 1
        if not 1 <= count <= 6:
            return None
        plmns = []
        for _ in range(count):
            mcc = None
            if take(1):
                mcc = tuple(take(4) for _ in range(3))
                if any(digit > 9 for digit in mcc):
                    return None
            wide = take(1)
            mnc = tuple(take(4) for _ in range(3 if wide else 2))
            if any(digit > 9 for digit in mnc):
                return None
            take(1)  # cellReservedForOperatorUse
            plmns.append((mcc, mnc))
        tac = take(16)
        cell_id = take(28)
    except ValueError:
        return None
    mcc, mnc = plmns[0]
    if mcc is None:
        return None
    mnc_value = int("".join(str(digit) for digit in mnc))
    mcc_value = int("".join(str(digit) for digit in mcc))
    return {
        "mcc": mcc_value,
        "mnc": mnc_value,
        "mnc_digits": "".join(str(digit) for digit in mnc),
        "tac": tac,
        "cell_id": cell_id,
        "enb_id": cell_id >> 8,
        "sector_id": cell_id & 0xFF,
        "operator": OPERATORS.get((mcc_value, mnc_value), ""),
        "plmns": plmns,
    }


def _data_res(elements, symbols):
    soft = []
    for symbol, k in elements:
        soft.extend(qpsk_soft(symbols[(symbol, k)]))
    return soft


def demod_pdcch(spectra, n_rb: int, pci: int, cfi: int, ng: str, n_s: int):
    """Return the first SI-RNTI DCI in the common search space, or None."""
    regs = _pdcch_regs(n_rb, pci, cfi, ng)
    cces = len(regs) // 9
    c_init = (n_s // 2) * (2 ** 9) + pci
    grants = []
    for aggregation in (8, 4):
        if aggregation > cces:
            continue
        candidates = min(2 if aggregation == 8 else 4, cces // aggregation)
        for index in range(candidates):
            start = aggregation * index
            chosen = regs[start * 9:(start + aggregation) * 9]
            if len(chosen) != 9 * aggregation:
                continue
            soft = []
            for reg in chosen:
                for symbol, k in reg:
                    value = spectra[symbol][_k_to_bin(k, n_rb)]
                    soft.extend(qpsk_soft(value))
            soft = _descramble_soft(soft, c_init)
            for length in _dci_lengths(n_rb):
                streams = conv_rate_recover(soft, length)
                message = viterbi_tailbiting(streams)
                parsed = _parse_dci_1a(message, n_rb)
                if parsed is not None:
                    parsed["aggregation"] = aggregation
                    grants.append(parsed)
                    return parsed
    return None


def _dci_lengths(n_rb: int):
    yield 2 + _riv_bits(n_rb) + 5 + 3 + 1 + 2 + 2 + 16


def synthesize_sib1_subframe(pci: int, n_rb: int, cell_id: int, mcc=(2, 3, 4), mnc=(1, 5), tac: int = 1):
    """One subframe 5 at the native sample rate of a capture that can hold n_rb."""
    if n_rb <= 25:
        nfft, rate = 512, 7.68e6
    elif n_rb < 100:
        nfft, rate = 1024, 15.36e6
    else:
        nfft, rate = 2048, 30.72e6
    cfi = 2
    ng = "1/6"
    mcs = 0
    start_rb, length = 0, min(10, n_rb)
    info = build_sib1_uper(cell_id, mcc, mnc, tac)
    tbs = tbs_from_mcs(length, mcs)
    dci = build_dci_1a(n_rb, start_rb, length, mcs, rv=0)
    regs = _pdcch_regs(n_rb, pci, cfi, ng)
    coded = _scramble(_pdcch_bits(dci, 72 * 4), (10 // 2) * (2 ** 9) + pci)
    control = _map_qpsk(regs[:36], coded)
    elements = _pdsch_elements(n_rb, pci, cfi, start_rb, length)
    streams, _size, _filler = _transport_streams(info, tbs)
    pdsch_bits = _scramble(turbo_rate_match(streams, len(elements) * 2, 0), SI_RNTI * 2 ** 14 + 5 * 2 ** 9 + pci)
    data = _map_qpsk([[element] for element in elements], pdsch_bits)
    buffer = np.zeros(_subframe_len(nfft), dtype=np.complex64)
    spectra = [np.zeros(nfft, dtype=np.complex128) for _ in range(14)]
    for (symbol, k), value in control:
        spectra[symbol][_k_to_bin(k, n_rb)] += value
    for (symbol, k), value in data:
        spectra[symbol][_k_to_bin(k, n_rb)] += value
    for symbol in (0, 4, 7, 11):
        slot = 10 if symbol < 7 else 11
        local = symbol if symbol < 7 else symbol - 7
        # subframe 5 is slots 10 and 11. symbol 0..6 are slot 10.
        n_s = 10 + (symbol // 7)
        local = symbol % 7
        for k in range(n_rb * 12):
            ref = crs_symbol(pci, n_s, local, k, n_rb, 0)
            if ref is not None:
                spectra[symbol][_k_to_bin(k, n_rb)] += ref
    for symbol, spec in enumerate(spectra):
        _write_symbol(buffer, 0, symbol, spec, nfft)
    return buffer, rate, nfft


def decode_sib1_subframe(samples, pci: int, n_rb: int, ng: str = "1/6", sfn: int = 0, sample_rate: float = 15.36e6) -> dict | None:
    """Decode SIB1 from one subframe-5 buffer that starts on symbol 0."""
    prepared = prepare_wide(samples, sample_rate, n_rb)
    if prepared is None:
        return None
    wave, rate, nfft = prepared
    if len(wave) < _subframe_len(nfft):
        return None
    spectra = [_symbol_fft(wave, 0, symbol, nfft) for symbol in range(14)]
    grant = None
    for cfi in (1, 2, 3):
        grant = demod_pdcch(spectra, n_rb, pci, cfi, ng, n_s=10)
        if grant is not None:
            grant["cfi"] = cfi
            break
    if grant is None:
        return None
    tbs = tbs_from_mcs(grant["length"], grant["mcs"])
    if tbs is None:
        return None
    elements = _pdsch_elements(n_rb, pci, grant["cfi"], grant["start"], grant["length"])
    soft = []
    gains = {}
    n_s = 10
    for k in range(n_rb * 12):
        ref = crs_symbol(pci, n_s, 0, k, n_rb, 0)
        if ref is None:
            continue
        gains[k] = spectra[0][_k_to_bin(k, n_rb)] * np.conj(ref)
    c_init = SI_RNTI * 2 ** 14 + 5 * 2 ** 9 + pci
    for symbol, k in elements:
        gain = _interp(gains, k)
        value = spectra[symbol][_k_to_bin(k, n_rb)] * np.conj(gain)
        soft.extend(qpsk_soft(value / (np.abs(gain) ** 2 + 1e-12)))
    soft = _descramble_soft(soft, c_init)
    rv = grant["rv"] if sfn == 0 else sib1_rv(sfn)
    if sfn:
        rv = sib1_rv(sfn)
    payload = _recover_transport(soft, tbs, rv if sfn else grant["rv"])
    if payload is None:
        return None
    return parse_sib1_uper(payload)


def _interp(gains, k):
    if not gains:
        return 1 + 0j
    if k in gains:
        return gains[k]
    keys = sorted(gains)
    if k <= keys[0] or k >= keys[-1]:
        nearest = min(keys, key=lambda item: abs(item - k))
        return gains[nearest]
    higher = next(key for key in keys if key >= k)
    lower = keys[keys.index(higher) - 1]
    weight = (k - lower) / float(higher - lower)
    return (1.0 - weight) * gains[lower] + weight * gains[higher]


def locate_wide(samples, pci: int, sample_rate: float, n_rb: int, cfo_hz: float = 0.0):
    """Sample index of the first subframe 0 at this sample rate."""
    prepared = prepare_wide(samples, sample_rate, n_rb)
    if prepared is None:
        return None
    wave, rate, nfft = prepared
    if cfo_hz:
        index = np.arange(len(wave))
        wave = wave * np.exp(-1j * 2 * np.pi * cfo_hz * index / rate)
    nid2 = pci % 3
    nid1 = pci // 3
    spec = np.zeros(nfft, dtype=np.complex128)
    freq = pss_time(nid2)
    # pss_time is 128 samples. Rebuild the centre tones in this FFT.
    zc = np.fft.fft(freq)
    # The 128-point symbol already has the 62 tones around its own DC.
    narrow = np.fft.fft(freq)
    spec[-31:] = narrow[-31:]
    spec[1:32] = narrow[1:32]
    template = np.fft.ifft(spec)
    template = template / np.linalg.norm(template)
    size = len(wave) + nfft - 1
    fft_len = 1 << (size - 1).bit_length()
    corr = np.fft.ifft(np.fft.fft(wave, fft_len) * np.conj(np.fft.fft(template, fft_len)))
    mag = np.abs(corr[:len(wave)])
    gap = nfft + cp_len(6, nfft)
    pss_at = useful_offset(6, nfft)
    frame = 10 * _subframe_len(nfft)
    sub = _subframe_len(nfft)
    for peak in np.argsort(mag)[::-1][:12]:
        peak = int(peak)
        sss_at = peak - gap
        if sss_at < 0 or peak + nfft > len(wave):
            continue
        symbol = np.fft.fft(wave[sss_at:sss_at + nfft])
        observed = np.empty(62, dtype=np.complex128)
        observed[:31] = symbol[-31:]
        observed[31:] = symbol[1:32]
        best_score, best_5 = -1.0, False
        for subframe5 in (False, True):
            ref = sss_freq(nid1, nid2, subframe5)
            score = abs(np.vdot(ref, observed)) / (np.linalg.norm(observed) * np.sqrt(62) + 1e-12)
            if score > best_score:
                best_score, best_5 = score, subframe5
        if best_score < 0.4:
            continue
        start = peak - pss_at - (5 * sub if best_5 else 0)
        while start >= frame:
            start -= frame
        if start < 0:
            start += frame
        if start + sub <= len(wave):
            return int(start), wave, rate, nfft
    return None


def decode_wideband_sib1(samples, sample_rate: float, pci: int, mib: dict, cfo_hz: float = 0.0) -> dict:
    """SIB1 from a wide capture. mib supplies bandwidth, SFN, PHICH and ports."""
    n_rb = int(mib["bandwidth_rb"])
    located = locate_wide(samples, pci, sample_rate, n_rb, cfo_hz)
    if located is None:
        return {"cell_id_status": "SIB1 subframe was not found"}
    start, wave, rate, nfft = located
    sub = _subframe_len(nfft)
    frame = 10 * sub
    ng = mib.get("phich_resource") or "1"
    sfn0 = int(mib.get("sfn") or 0)
    acc = None
    tbs = None
    filler = 0
    size = 0
    seen = 0
    for index in range(0, max(0, (len(wave) - start) // frame)):
        sfn = (sfn0 + index) % 1024
        if sfn % 2:
            continue
        origin = start + index * frame + 5 * sub
        if origin + sub > len(wave):
            break
        spectra = [_symbol_fft(wave, origin, symbol, nfft) for symbol in range(14)]
        grant = None
        for cfi in (1, 2, 3):
            grant = demod_pdcch(spectra, n_rb, pci, cfi, ng, n_s=10)
            if grant is not None:
                grant["cfi"] = cfi
                break
        if grant is None:
            continue
        tbs = tbs_from_mcs(grant["length"], grant["mcs"])
        if tbs is None:
            continue
        elements = _pdsch_elements(n_rb, pci, grant["cfi"], grant["start"], grant["length"])
        gains = {}
        for k in range(n_rb * 12):
            ref = crs_symbol(pci, 10, 0, k, n_rb, 0)
            if ref is not None:
                gains[k] = spectra[0][_k_to_bin(k, n_rb)] * np.conj(ref)
        soft = []
        for symbol, k in elements:
            gain = _interp(gains, k)
            value = spectra[symbol][_k_to_bin(k, n_rb)] * np.conj(gain)
            soft.extend(qpsk_soft(value / (np.abs(gain) ** 2 + 1e-12)))
        c_init = SI_RNTI * 2 ** 14 + 5 * 2 ** 9 + pci
        soft = _descramble_soft(soft, c_init)
        rv = sib1_rv(sfn)
        size = code_block_size(tbs + 24)
        filler = size - (tbs + 24)
        streams = turbo_rate_recover(soft, size, redundancy=rv, filler=filler)
        if acc is None:
            acc = streams
        else:
            for stream in range(3):
                for index_bit in range(len(streams[stream])):
                    acc[stream][index_bit] += streams[stream][index_bit]
        seen += 1
        decoded = turbo_decode(acc)
        block = decoded[filler:]
        payload, parity = block[:tbs], block[tbs:tbs + 24]
        if parity != crc24a(payload):
            continue
        identity = parse_sib1_uper(payload)
        if identity is None:
            continue
        identity["cell_id_status"] = "decoded"
        return identity
    if seen == 0:
        return {"cell_id_status": "SIB1 grant was not found"}
    return {"cell_id_status": "SIB1 did not decode"}



# I_TBS 0..9, N_PRB 1..100. TS 36.213 Table 7.1.7.2.1-1.
_TBS_ROWS = (
    (16,32,56,88,120,152,176,208,224,256,288,328,344,376,392,424,456,488,504,536,568,600,616,648,680,712,744,776,776,808,840,872,904,936,968,1000,1032,1032,1064,1096,1128,1160,1192,1224,1256,1256,1288,1320,1352,1384,1416,1416,1480,1480,1544,1544,1608,1608,1608,1672,1672,1736,1736,1800,1800,1800,1864,1864,1928,1928,1992,1992,2024,2088,2088,2088,2152,2152,2216,2216,2280,2280,2280,2344,2344,2408,2408,2472,2472,2536,2536,2536,2600,2600,2664,2664,2728,2728,2728,2792),
    (24,56,88,144,176,208,224,256,328,344,376,424,456,488,520,568,600,632,680,712,744,776,808,872,904,936,968,1000,1032,1064,1128,1160,1192,1224,1256,1288,1352,1384,1416,1416,1480,1544,1544,1608,1608,1672,1736,1736,1800,1800,1864,1864,1928,1992,1992,2024,2088,2088,2152,2152,2216,2280,2280,2344,2344,2408,2472,2472,2536,2536,2600,2600,2664,2728,2728,2792,2792,2856,2856,2856,2984,2984,2984,3112,3112,3112,3240,3240,3240,3240,3368,3368,3368,3496,3496,3496,3496,3624,3624,3624),
    (32,72,144,176,208,256,296,328,376,424,472,520,568,616,648,696,744,776,840,872,936,968,1000,1064,1096,1160,1192,1256,1288,1320,1384,1416,1480,1544,1544,1608,1672,1672,1736,1800,1800,1864,1928,1992,2024,2088,2088,2152,2216,2216,2280,2344,2344,2408,2472,2536,2536,2600,2664,2664,2728,2792,2856,2856,2856,2984,2984,3112,3112,3112,3240,3240,3240,3368,3368,3368,3496,3496,3496,3624,3624,3624,3752,3752,3880,3880,3880,4008,4008,4008,4136,4136,4136,4264,4264,4264,4392,4392,4392,4584),
    (40,104,176,208,256,328,392,440,504,568,616,680,744,808,872,904,968,1032,1096,1160,1224,1256,1320,1384,1416,1480,1544,1608,1672,1736,1800,1864,1928,1992,2024,2088,2152,2216,2280,2344,2408,2472,2536,2536,2600,2664,2728,2792,2856,2856,2984,2984,3112,3112,3240,3240,3368,3368,3496,3496,3624,3624,3624,3752,3752,3880,3880,4008,4008,4136,4136,4264,4264,4392,4392,4392,4584,4584,4584,4776,4776,4776,4776,4968,4968,4968,5160,5160,5160,5352,5352,5352,5352,5544,5544,5544,5736,5736,5736,5736),
    (56,120,208,256,328,408,488,552,632,696,776,840,904,1000,1064,1128,1192,1288,1352,1416,1480,1544,1608,1736,1800,1864,1928,1992,2088,2152,2216,2280,2344,2408,2472,2600,2664,2728,2792,2856,2984,2984,3112,3112,3240,3240,3368,3496,3496,3624,3624,3752,3752,3880,4008,4008,4136,4136,4264,4264,4392,4392,4584,4584,4584,4776,4776,4968,4968,4968,5160,5160,5160,5352,5352,5544,5544,5544,5736,5736,5736,5992,5992,5992,5992,6200,6200,6200,6456,6456,6456,6456,6712,6712,6712,6968,6968,6968,6968,7224),
    (72,144,224,328,424,504,600,680,776,872,968,1032,1128,1224,1320,1384,1480,1544,1672,1736,1864,1928,2024,2088,2216,2280,2344,2472,2536,2664,2728,2792,2856,2984,3112,3112,3240,3368,3496,3496,3624,3752,3752,3880,4008,4008,4136,4264,4392,4392,4584,4584,4776,4776,4776,4968,4968,5160,5160,5352,5352,5544,5544,5736,5736,5736,5992,5992,5992,6200,6200,6200,6456,6456,6712,6712,6712,6968,6968,6968,7224,7224,7224,7480,7480,7480,7736,7736,7736,7992,7992,7992,8248,8248,8248,8504,8504,8760,8760,8760),
    (328,176,256,392,504,600,712,808,936,1032,1128,1224,1352,1480,1544,1672,1736,1864,1992,2088,2216,2280,2408,2472,2600,2728,2792,2984,2984,3112,3240,3368,3496,3496,3624,3752,3880,4008,4136,4136,4264,4392,4584,4584,4776,4776,4968,4968,5160,5160,5352,5352,5544,5736,5736,5992,5992,5992,6200,6200,6456,6456,6456,6712,6712,6968,6968,6968,7224,7224,7480,7480,7736,7736,7736,7992,7992,8248,8248,8248,8504,8504,8760,8760,8760,9144,9144,9144,9144,9528,9528,9528,9528,9912,9912,9912,10296,10296,10296,10296),
    (104,224,328,472,584,712,840,968,1096,1224,1320,1480,1608,1672,1800,1928,2088,2216,2344,2472,2536,2664,2792,2984,3112,3240,3368,3368,3496,3624,3752,3880,4008,4136,4264,4392,4584,4584,4776,4968,4968,5160,5352,5352,5544,5736,5736,5992,5992,6200,6200,6456,6456,6712,6712,6712,6968,6968,7224,7224,7480,7480,7736,7736,7992,7992,8248,8248,8504,8504,8760,8760,8760,9144,9144,9144,9528,9528,9528,9912,9912,9912,10296,10296,10296,10680,10680,10680,11064,11064,11064,11448,11448,11448,11448,11832,11832,11832,12216,12216),
    (120,256,392,536,680,808,968,1096,1256,1384,1544,1672,1800,1928,2088,2216,2344,2536,2664,2792,2984,3112,3240,3368,3496,3624,3752,3880,4008,4264,4392,4584,4584,4776,4968,4968,5160,5352,5544,5544,5736,5992,5992,6200,6200,6456,6456,6712,6968,6968,7224,7224,7480,7480,7736,7736,7992,7992,8248,8504,8504,8760,8760,9144,9144,9144,9528,9528,9528,9912,9912,9912,10296,10296,10680,10680,10680,11064,11064,11064,11448,11448,11448,11832,11832,12216,12216,12216,12576,12576,12576,12960,12960,12960,13536,13536,13536,13536,14112,14112),
    (136,296,456,616,776,936,1096,1256,1416,1544,1736,1864,2024,2216,2344,2536,2664,2856,2984,3112,3368,3496,3624,3752,4008,4136,4264,4392,4584,4776,4968,5160,5160,5352,5544,5736,5736,5992,6200,6200,6456,6712,6712,6968,6968,7224,7480,7480,7736,7992,7992,8248,8248,8504,8760,8760,9144,9144,9144,9528,9528,9912,9912,10296,10296,10296,10680,10680,11064,11064,11064,11448,11448,11832,11832,11832,12216,12216,12576,12576,12960,12960,12960,13536,13536,13536,13536,14112,14112,14112,14112,14688,14688,14688,15264,15264,15264,15264,15840,15840),
)
