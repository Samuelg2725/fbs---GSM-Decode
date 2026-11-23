"""LTE physical-cell-ID search from complex baseband.

This is the same PSS/SSS detection LTE Cell Scanner performs. It reads
samples that have already been resampled to 1.92 MHz. It does not open a
radio and it does not decode subscriber data.
"""

from __future__ import annotations

import numpy as np

FS = 1_920_000
FFT_SIZE = 128
CP = 9
SYMBOL = FFT_SIZE + CP  # 137 samples at 1.92 MHz
PSS_ROOTS = (25, 29, 34)  # N_ID_2 = 0, 1, 2
# A Kraken tuner can sit tens of kHz off the LTE centre. The PSS is one
# 15 kHz bin wide in time, so the search has to try the offset itself.
PSS_MIN_RATIO = 6.0
SSS_MIN_SCORE = 0.55
_PERIOD = 9600  # one PSS every 5 ms at 1.92 MHz
_EXT_GAP = FFT_SIZE + 32  # extended cyclic prefix
_COARSE_HZ = tuple(range(-50000, 50001, 5000))
_FINE_HZ = tuple(range(-4000, 4001, 1000))


def _zc63(root: int) -> np.ndarray:
    n = np.arange(63)
    return np.exp(-1j * np.pi * root * n * (n + 1) / 63.0)


def _map_62(seq62: np.ndarray) -> np.ndarray:
    """Place 62 subcarriers around DC and return the 128-sample time symbol."""
    spec = np.zeros(FFT_SIZE, dtype=np.complex128)
    spec[-31:] = seq62[:31]
    spec[1:32] = seq62[31:]
    return np.fft.ifft(spec)


def pss_time(nid2: int) -> np.ndarray:
    zc = _zc63(PSS_ROOTS[nid2])
    seq = np.concatenate((zc[:31], zc[32:]))
    return _map_62(seq)


def _mseq(taps) -> np.ndarray:
    x = [0, 0, 0, 0, 1]
    for i in range(26):
        nxt = 0
        for tap in taps:
            nxt ^= x[i + tap]
        x.append(nxt)
    return np.array([1 - 2 * bit for bit in x[:31]], dtype=np.float64)


_S = _mseq((2, 0))
_C = _mseq((3, 0))
_Z = _mseq((4, 2, 1, 0))


def _m0_m1(nid1: int):
    # 3GPP TS 36.211 6.11.2. m' is not reduced by 30*q; m0 is m' modulo 31.
    q_prime = nid1 // 30
    q = (nid1 + q_prime * (q_prime + 1) // 2) // 30
    m_prime = nid1 + q * (q + 1) // 2
    m0 = m_prime % 31
    m1 = (m0 + m_prime // 31 + 1) % 31
    return m0, m1


def sss_freq(nid1: int, nid2: int, subframe5: bool) -> np.ndarray:
    """62 SSS subcarriers, lowest frequency first. 3GPP TS 36.211 6.11.2."""
    m0, m1 = _m0_m1(nid1)
    s0 = np.take(_S, (np.arange(31) + m0) % 31)
    s1 = np.take(_S, (np.arange(31) + m1) % 31)
    c0 = np.take(_C, (np.arange(31) + nid2) % 31)
    c1 = np.take(_C, (np.arange(31) + nid2 + 3) % 31)
    z0 = np.take(_Z, (np.arange(31) + (m0 % 8)) % 31)
    z1 = np.take(_Z, (np.arange(31) + (m1 % 8)) % 31)
    d = np.empty(62, dtype=np.complex128)
    if not subframe5:
        d[0::2] = s0 * c0
        d[1::2] = s1 * c1 * z0
    else:
        d[0::2] = s1 * c0
        d[1::2] = s0 * c1 * z1
    return d


def sss_time(nid1: int, nid2: int, subframe5: bool) -> np.ndarray:
    return _map_62(sss_freq(nid1, nid2, subframe5))


def synthesize_pair(pci: int, start: int, length: int) -> np.ndarray:
    """A clean SSS symbol followed by a PSS symbol, for tests."""
    nid2 = pci % 3
    nid1 = pci // 3
    buf = np.zeros(length, dtype=np.complex128)
    sss = sss_time(nid1, nid2, False)
    pss = pss_time(nid2)
    buf[start:start + CP] = sss[-CP:]
    buf[start + CP:start + SYMBOL] = sss
    buf[start + SYMBOL:start + SYMBOL + CP] = pss[-CP:]
    buf[start + SYMBOL + CP:start + 2 * SYMBOL] = pss
    return buf


def _sss_from_symbol(symbol: np.ndarray) -> np.ndarray:
    spec = np.fft.fft(symbol)
    seq = np.empty(62, dtype=np.complex128)
    seq[:31] = spec[-31:]
    seq[31:] = spec[1:32]
    return seq


def _shifted(template: np.ndarray, offset_hz: float) -> np.ndarray:
    if offset_hz == 0:
        return template
    ramp = np.exp(1j * 2 * np.pi * offset_hz * np.arange(len(template)) / FS)
    return template * ramp


def _correlate_mag(spectrum: np.ndarray, template: np.ndarray, length: int) -> np.ndarray:
    corr = np.fft.ifft(spectrum * np.conj(np.fft.fft(template, spectrum.size)))
    return np.abs(corr[:length])


def _strongest(mag: np.ndarray):
    at = int(np.argmax(mag))
    return at, float(mag[at])


def _usable_peak(mag: np.ndarray):
    """Strongest PSS peak that still has the preceding SSS symbol in the buffer."""
    work = np.array(mag, copy=True)
    for _ in range(8):
        at, peak = _strongest(work)
        if at >= SYMBOL and at + FFT_SIZE <= len(work):
            return at, peak
        lo = max(0, at - FFT_SIZE)
        hi = min(len(work), at + FFT_SIZE + 1)
        work[lo:hi] = 0.0
        if peak <= 0.0:
            break
    return None


def _derotate(chunk: np.ndarray, start: int, offset_hz: float) -> np.ndarray:
    if offset_hz == 0:
        return chunk
    index = np.arange(start, start + len(chunk))
    return chunk * np.exp(-1j * 2 * np.pi * offset_hz * index / FS)


def _pss_freq(nid2: int) -> np.ndarray:
    zc = _zc63(PSS_ROOTS[nid2])
    return np.concatenate((zc[:31], zc[32:]))


def _extract_seq(samples: np.ndarray, start: int, offset_hz: float):
    if start < 0 or start + FFT_SIZE > len(samples):
        return None
    return _sss_from_symbol(_derotate(samples[start:start + FFT_SIZE], start, offset_hz))


def _mean_match(observed_list, reference: np.ndarray) -> float:
    if not observed_list:
        return 0.0
    scale = np.sqrt(62.0)
    total = 0.0
    for observed in observed_list:
        denom = float(np.linalg.norm(observed) * scale + 1e-12)
        total += abs(np.vdot(reference, observed)) / denom
    return total / len(observed_list)


def _blank_comb(mag: np.ndarray, anchor: int, radius: int = 250) -> None:
    """Remove one cell's PSS repeats so the next peak can belong to another cell.

    Two cells on the same EARFCN share a root only when their PCI differs by a
    multiple of 3. Their 5 ms combs sit at different sample phases. Blanking
    a short window on one comb leaves the other comb in place. The window is
    wider than the PSS main lobe and narrower than the gap between unsynchronised
    cells, so a second lock is not erased with the first.
    """
    span = len(mag)
    first = -int(anchor // _PERIOD) - 2
    last = int((span - anchor) // _PERIOD) + 2
    for step in range(first, last + 1):
        center = anchor + step * _PERIOD
        lo = max(0, center - radius)
        hi = min(span, center + radius + 1)
        mag[lo:hi] = 0.0


def _pss_times(mag: np.ndarray, anchor: int):
    """PSS locations on the 5 ms grid, tracked from the strongest peak."""
    times = []
    span = int(len(mag))
    first = -int(anchor // _PERIOD) - 1
    last = int((span - anchor) // _PERIOD) + 1
    floor = 0.5 * float(mag[anchor])
    for step in range(first, last + 1):
        center = anchor + step * _PERIOD
        lo = center - 30
        hi = center + 31
        if lo < 0 or hi > span:
            continue
        local = int(np.argmax(mag[lo:hi])) + lo
        if float(mag[local]) >= floor and local >= SYMBOL:
            times.append(local)
    if not times and SYMBOL <= anchor < span:
        times.append(anchor)
    return times


def _confirm_sss(samples: np.ndarray, mag: np.ndarray, anchor: int, offset_hz: float, nid2: int):
    """Line the FFT window up on the PSS, then average every SSS in the capture.

    A one-symbol SSS score sits near the noise floor when the cell is weak or
    the window is a few samples off. The PSS sequence is known once NID2 is
    known, so it picks the timing. Subframe 0 and subframe 5 alternate every
    5 ms and are averaged separately.
    """
    times = _pss_times(mag, anchor)
    reference = _pss_freq(nid2)
    best_delta = 0
    best_pss = -1.0
    for delta in range(-10, 11):
        observed = []
        for at in times:
            seq = _extract_seq(samples, at + delta, offset_hz)
            if seq is not None:
                observed.append(seq)
        score = _mean_match(observed, reference)
        if score > best_pss:
            best_pss = score
            best_delta = delta

    best_nid1 = 0
    best_sss = -1.0
    origin = times[0]
    for gap in (SYMBOL, _EXT_GAP):
        for delta in (best_delta - 1, best_delta, best_delta + 1):
            groups = ([], [])
            for at in times:
                seq = _extract_seq(samples, at - gap + delta, offset_hz)
                if seq is None:
                    continue
                slot = int(round((at - origin) / float(_PERIOD)))
                groups[slot % 2].append(seq)
            if not groups[0] and not groups[1]:
                continue
            for nid1 in range(168):
                for swap in (False, True):
                    scores = []
                    for index, subframe5 in ((0, swap), (1, not swap)):
                        if groups[index]:
                            scores.append(_mean_match(groups[index], sss_freq(nid1, nid2, subframe5)))
                    score = sum(scores) / len(scores)
                    if score > best_sss:
                        best_sss = score
                        best_nid1 = nid1
    return best_nid1, float(best_sss), float(best_pss)


def search_cell(samples: np.ndarray, max_cells: int = 6) -> dict:
    """PSS/SSS search, including a tuner-offset sweep.

    best_ratio is the strongest PSS peak divided by the RMS of the capture.
    A real cell is reported only when that peak clears PSS_MIN_RATIO and the
    SSS symbols on the same 5 ms grid agree.
    """
    samples = np.asarray(samples, dtype=np.complex128)
    length = int(samples.size)
    scale = float(np.sqrt(np.mean(np.abs(samples) ** 2) + 1e-12))
    nfft = 1 << (length + FFT_SIZE - 2).bit_length()
    spectrum = np.fft.fft(samples, nfft)
    bases = [pss_time(nid2) for nid2 in range(3)]
    bases = [template / np.linalg.norm(template) for template in bases]

    best_ratio = 0.0
    best_offset = 0.0
    coarse = [None] * 3

    def remember(nid2, offset, mag):
        nonlocal best_ratio, best_offset
        _at, peak = _strongest(mag)
        ratio = peak / scale
        if ratio > best_ratio:
            best_ratio = ratio
            best_offset = float(offset)
        current = coarse[nid2]
        if current is None or peak > current[0]:
            coarse[nid2] = (peak, float(offset), mag)

    for offset in _COARSE_HZ:
        for nid2, base in enumerate(bases):
            remember(nid2, offset, _correlate_mag(spectrum, _shifted(base, offset), length))

    refined = [None] * 3
    for nid2, base in enumerate(bases):
        if coarse[nid2] is None:
            continue
        winner = coarse[nid2]
        center = winner[1]
        for delta in _FINE_HZ:
            offset = center + delta
            if delta == 0:
                continue
            mag = _correlate_mag(spectrum, _shifted(base, offset), length)
            _at, peak = _strongest(mag)
            ratio = peak / scale
            if ratio > best_ratio:
                best_ratio = ratio
                best_offset = float(offset)
            if peak > winner[0]:
                winner = (peak, float(offset), mag)
        refined[nid2] = winner

    hits = []
    best_sss = 0.0
    best_pss_seq = 0.0
    for nid2, winner in enumerate(refined):
        if winner is None:
            continue
        _peak, offset, mag = winner
        work = np.array(mag, copy=True)
        # Three combs per root: the strongest cell, a weaker one on the same
        # root, and one spare for a sidelobe that fails SSS and is blanked.
        for _attempt in range(3):
            usable = _usable_peak(work)
            if usable is None:
                break
            at, peak = usable
            ratio = peak / scale
            if ratio < PSS_MIN_RATIO:
                break
            nid1, score, pss_seq = _confirm_sss(samples, work, at, offset, nid2)
            _blank_comb(work, at)
            if pss_seq > best_pss_seq:
                best_pss_seq = pss_seq
            if score > best_sss:
                best_sss = score
            if score < SSS_MIN_SCORE or pss_seq < SSS_MIN_SCORE:
                continue
            hits.append({
                "pci": 3 * nid1 + nid2,
                "nid1": nid1,
                "nid2": nid2,
                "pss_peak": peak,
                "pss_ratio": ratio,
                "pss_seq": pss_seq,
                "sss_score": score,
                "offset_hz": offset,
            })
    hits.sort(key=lambda item: item["pss_ratio"], reverse=True)
    return {
        "hits": hits[:max_cells],
        "best_ratio": best_ratio,
        "best_offset_hz": best_offset,
        "min_ratio": PSS_MIN_RATIO,
        "best_sss": best_sss,
        "best_pss_seq": best_pss_seq,
        "sss_min": SSS_MIN_SCORE,
    }


def find_pci(samples: np.ndarray, max_cells: int = 6) -> list[dict]:
    """Return every physical cell ID that clears the sync thresholds, strongest first."""
    return search_cell(samples, max_cells)["hits"]
