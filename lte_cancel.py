#!/usr/bin/env python3
"""Interference cancellation for LTE: remove a decoded cell's reference signals (CRS) from a recording, so
weaker cells on the same carrier can be found and decoded underneath it.

Every LTE cell sends cell-specific reference signals on a fixed pattern of subcarriers in every slot, whether or
not it carries traffic. Their values follow from the cell's PCI (36.211 6.10.1), so once a cell is decoded they
are fully known. On a lightly loaded carrier they are most of what a strong cell puts on top of a weak one.

For each slot of the decoded cell (its own frame timing):
  1. take the OFDM symbols that carry CRS (0 and 4 of each slot, and 1 with four antenna ports);
  2. measure the cell's channel on its CRS (received / known value), smoothed over neighbouring pilots in
     frequency (after removing the phase slope of any timing error) and in time (the same symbol of the
     neighbouring slots), so the estimate follows the strong cell and not the weak cells or the noise;
  3. rebuild the strong cell's CRS as received (known value x channel) and subtract it, cyclic prefix included.
The weak cells' own signals are untouched apart from those resource elements' share of the noise.
The cell's sync signals (PSS and SSS, also fully known from the PCI) are removed the same way: they repeat every
5 ms like a weak cell's and would otherwise pull the weak cell's timing search towards the strong cell.

  python3 lte_cancel.py      # self-check with synthetic cells (prints the suppression achieved)
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if os.path.join(_HERE, "osprey 2") not in sys.path:
    sys.path.insert(0, os.path.join(_HERE, "osprey 2"))         # lte_sync (PSS / SSS sequences)

N_RB_MAX = 110
_SEQ = {}


def gold_bits(c_init, length):
    """Pseudo-random sequence c(n) of 36.211 7.2 (length bits)."""
    nc = 1600
    total = nc + length
    x1 = np.zeros(total + 31, np.uint8)
    x2 = np.zeros(total + 31, np.uint8)
    x1[0] = 1
    for i in range(31):
        x2[i] = (c_init >> i) & 1
    for n in range(total):
        x1[n + 31] = x1[n + 3] ^ x1[n]
        x2[n + 31] = x2[n + 3] ^ x2[n + 2] ^ x2[n + 1] ^ x2[n]
    return (x1[nc:nc + length] ^ x2[nc:nc + length]).astype(np.int8)


def crs_sequence(pci, ns, l):
    """r(m') for m' = 0 .. 2*110-1 (normal cyclic prefix), cached."""
    key = (pci, ns, l)
    seq = _SEQ.get(key)
    if seq is None:
        c_init = (2 ** 10) * (7 * (ns + 1) + l + 1) * (2 * pci + 1) + 2 * pci + 1
        c = gold_bits(c_init, 4 * N_RB_MAX).astype(np.float32)
        seq = ((1 - 2 * c[0::2]) + 1j * (1 - 2 * c[1::2])).astype(np.complex64) / np.sqrt(2)
        _SEQ[key] = seq
    return seq


def crs_ports_at(l, ports):
    """Antenna ports that carry CRS on slot symbol l."""
    if l in (0, 4):
        return [0] if ports == 1 else [0, 1]
    if l == 1 and ports == 4:
        return [2, 3]
    return []


def crs_shift(port, l, ns):
    """Frequency offset v of 36.211 6.10.1.2."""
    if port == 0:
        return 0 if l == 0 else 3
    if port == 1:
        return 3 if l == 0 else 0
    if port == 2:
        return 3 * (ns % 2)
    return 3 + 3 * (ns % 2)


def grid(rate):
    """FFT size, cyclic prefix lengths (first, other) and slot length in samples for a sample rate."""
    nfft = int(round(rate / 15000.0))
    cp0 = 160 * nfft // 2048
    cp = 144 * nfft // 2048
    return nfft, cp0, cp, 7 * nfft + cp0 + 6 * cp


def symbol_offset(l, nfft, cp0, cp):
    """Start of the useful part of slot symbol l, from the start of the slot."""
    if l == 0:
        return cp0
    return cp0 + nfft + (l - 1) * (cp + nfft) + cp


def bins_for(kk, n_rb, nfft):
    """FFT bins of LTE subcarriers kk (0 .. 12*n_rb-1), DC unused."""
    d = kk - 6 * n_rb
    return np.where(d < 0, nfft + d, d + 1)


def _smooth_freq(h, width=5):
    """Moving average across pilots, after removing each row's phase slope (timing error)."""
    from scipy.ndimage import uniform_filter1d
    slope = np.angle(np.sum(h[:, 1:] * np.conj(h[:, :-1]), axis=1))
    ramp = np.exp(1j * slope[:, None] * np.arange(h.shape[1])[None, :]).astype(np.complex64)
    flat = h * np.conj(ramp)
    sm = uniform_filter1d(flat.real, width, axis=1, mode="nearest") + \
        1j * uniform_filter1d(flat.imag, width, axis=1, mode="nearest")
    return (sm * ramp).astype(np.complex64)


def _smooth_time(h):
    """0.25 / 0.5 / 0.25 across consecutive slots (rows)."""
    if h.shape[0] < 3:
        return h
    out = h.copy()
    out[1:-1] = 0.25 * h[:-2] + 0.5 * h[1:-1] + 0.25 * h[2:]
    return out


_PBCH = {}
ANTENNA_MASK = {1: 0x0000, 2: 0xFFFF, 4: 0x5555}
_CONV_PERM = [1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31,
              0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30]


def conv_rate_match_index(d, e):
    """Positions in the three concatenated encoder streams (3*d bits) that make up the e rate-matched bits:
    36.212 5.1.4.2 - sub-block interleaver with 32 columns (inter-column permutation, dummy bits in front), the
    three interleaved streams one after another, read circularly skipping the dummy bits. (The library's own
    conv_rate_match does not follow this; checked bit for bit against srsRAN's transmitter.)"""
    rows = -(-d // 32)
    nd = rows * 32 - d
    w = []
    for stream in range(3):
        y = [None] * nd + [stream * d + i for i in range(d)]
        w += [y[r * 32 + _CONV_PERM[c]] for c in range(32) for r in range(rows)]
    w = [x for x in w if x is not None]
    return np.array([w[i % len(w)] for i in range(e)])


def pbch_symbols(pci, ports, mib_hex):
    """QPSK symbols of the cell's PBCH for every SFN: array (256, 4, 240) indexed [SFN // 4, SFN % 4, element].
    Only the SFN field of the MIB changes from one 40 ms period to the next; the rest is taken from the decoded
    MIB (36.212 5.3.1: CRC16 with the antenna-port mask, tail-biting convolutional code, rate matching to 1920
    bits; 36.211 6.6: scrambling with c_init = PCI, 480 bits per frame)."""
    key = (pci, ports, mib_hex)
    if key in _PBCH:
        return _PBCH[key]
    from lte_coding import bits_of, bits_to_int, conv_encode, crc16, gold
    template = bits_of(int(mib_hex, 16), 24)
    scr = np.array(gold(pci, 1920), np.uint8)
    idx = conv_rate_match_index(40, 1920)
    out = np.zeros((256, 4, 240), np.complex64)
    for b in range(256):
        bits = template[:6] + bits_of(b, 8) + template[14:]
        masked = bits_to_int(crc16(bits)) ^ ANTENNA_MASK[ports]
        streams = np.concatenate([np.array(st, np.uint8) for st in conv_encode(bits + bits_of(masked, 16))])
        coded = streams[idx] ^ scr
        sym = ((1 - 2 * coded[0::2].astype(np.float32)) + 1j * (1 - 2 * coded[1::2].astype(np.float32))) / np.sqrt(2)
        out[b] = sym.reshape(4, 240)
    _PBCH[key] = out
    return out


def _precode(x, h0, h1, ports):
    """Received PBCH for symbols x (..., 240) through per-element channels h0/h1 (..., 240): one port, or
    space-frequency block code over pairs of elements for two ports (36.211 6.3.4.3)."""
    if ports == 1:
        return h0 * x
    r = np.empty(np.broadcast(x, h0).shape, np.complex64)
    x0, x1 = x[..., 0::2], x[..., 1::2]
    r[..., 0::2] = h0[..., 0::2] * x0 - h1[..., 0::2] * np.conj(x1)
    r[..., 1::2] = h0[..., 1::2] * x1 + h1[..., 1::2] * np.conj(x0)
    return r


def _interp_channel(h_pilots, k_pilots, k_targets):
    """Linear interpolation of per-row pilot channels (rows, P) at pilot subcarriers k_pilots onto k_targets."""
    i1 = np.clip(np.searchsorted(k_pilots, k_targets), 1, len(k_pilots) - 1)
    i0 = i1 - 1
    w1 = ((k_targets - k_pilots[i0]) / (k_pilots[i1] - k_pilots[i0])).astype(np.float32)
    return h_pilots[:, i0] * (1 - w1) + h_pilots[:, i1] * w1


def _cancel_pbch(y, js, slot_starts, pci, n_rb, nfft, cp0, cp, ports, mib_hex, slot1_h):
    """Remove the cell's PBCH (slot 1 of subframe 0, symbols 0-3, centre 72 subcarriers), in place. The SFN of each
    frame is found by matching the received PBCH against all 1024 possibilities in the first frames (it then counts
    up by one per frame); the channel comes from the cell's reference signals in the same slot; amplitude and
    phase are fitted per frame (least squares), so the PBCH power offset need not be known."""
    import scipy.fft as sf
    from lte_cellid import pbch_data_res
    rows = np.nonzero((js % 20) == 1)[0]
    if len(rows) == 0 or ports not in (1, 2) or not slot1_h:
        return None
    elements = pbch_data_res(n_rb, pci)
    sym_l = np.array([e[0] for e in elements])
    kk = np.array([e[1] for e in elements])
    b = bins_for(kk, n_rb, nfft)
    spec = []
    for l in range(4):
        useful = slot_starts[rows] + symbol_offset(l, nfft, cp0, cp)
        spec.append(sf.fft(y[useful[:, None] + np.arange(nfft)[None, :]], axis=1, workers=-1))
    yv = np.stack([spec[l][:, b[i]] for i, l in enumerate(sym_l)], axis=1)          # (frames, 240)
    hs = []
    for port in range(ports):
        k_p, h_p = slot1_h[port]
        hs.append(_interp_channel(h_p, k_p, kk))
    h0 = hs[0]
    h1 = hs[1] if ports == 2 else None
    table = pbch_symbols(pci, ports, mib_hex)                                       # (256, 4, 240)
    frame_no = (js[rows] - js[rows][0]) // 20
    # SFN of the first frame: best match over all 1024 candidates in up to 6 frames, then a vote
    votes = {}
    cand = table.reshape(1024, 240)                                                   # index = SFN
    for i in range(min(6, len(rows))):
        r = _precode(cand, h0[i][None, :], None if h1 is None else h1[i][None, :], ports)
        metric = np.abs(np.sum(np.conj(r) * yv[i][None, :], axis=1)) ** 2 / (np.sum(np.abs(r) ** 2, axis=1) + 1e-12)
        sfn0 = (int(np.argmax(metric)) - int(frame_no[i])) % 1024
        votes[sfn0] = votes.get(sfn0, 0) + 1
    sfn0 = max(votes, key=votes.get)
    sfns = (sfn0 + frame_no) % 1024
    x = cand[sfns]                                                                    # (frames, 240)
    r = _precode(x, h0, h1, ports)
    alpha = np.sum(np.conj(r) * yv, axis=1) / (np.sum(np.abs(r) ** 2, axis=1) + 1e-12)
    fit = alpha[:, None] * r
    for l in range(4):
        sel = sym_l == l
        delta = np.zeros((len(rows), nfft), np.complex64)
        delta[:, b[sel]] = fit[:, sel]
        d = sf.ifft(delta, axis=1, workers=-1).astype(np.complex64)
        useful = slot_starts[rows] + symbol_offset(l, nfft, cp0, cp)
        cpl = cp0 if l == 0 else cp
        y[useful[:, None] + np.arange(nfft)[None, :]] -= d
        y[(useful - cpl)[:, None] + np.arange(cpl)[None, :]] -= d[:, nfft - cpl:]
    return sfn0, int(js[rows][0]) // 20                 # SFN of the frame that holds the first slot-1 row


def regen_sib1(binary, pci, n_rb, ports, mib_hex, dci, payload_hex):
    """The cell's SIB1 subframe rebuilt by srsRAN's transmitter (sib1_regen): array (4, ports, 14, 12*n_rb), one
    subframe per redundancy-version slot k = (SFN/2) mod 4. dci: the decoder's SIB1dci fields. None on failure."""
    import subprocess
    import tempfile
    bits = int(mib_hex, 16)
    phich_ext = (bits >> 20) & 1
    phich_res = (bits >> 18) & 3
    with tempfile.NamedTemporaryFile(suffix=".cf32", delete=False, dir="/dev/shm") as fh:
        out = fh.name
    try:
        cmd = [str(binary), "-c", str(pci), "-p", str(n_rb), "-P", str(ports), "-g", str(phich_res),
               "-d", str(phich_ext), "-f", str(dci["cfi"]), "-L", str(dci["L"]), "-N", str(dci["ncce"]),
               "-r", str(dci["riv"]), "-M", str(dci["mode"]), "-a", str(dci["nprb1a"]), "-G", str(dci["ngap"]),
               "-m", str(dci["mcs"]), "-x", payload_hex, "-o", out]
        if subprocess.run(cmd, capture_output=True, timeout=30).returncode != 0:
            return None
        g = np.fromfile(out, np.complex64)
        if g.size != 4 * ports * 14 * 12 * n_rb:
            return None
        return g.reshape(4, ports, 14, 12 * n_rb)
    except (OSError, subprocess.SubprocessError):
        return None
    finally:
        try:
            os.remove(out)
        except OSError:
            pass


def _cancel_sib1(y, js, slot_starts, n_rb, nfft, cp0, cp, ports, grids, sf5_h, sfn_info):
    """Remove the cell's SIB1 subframes (PCFICH, the SIB1 control message and SIB1 itself, rebuilt exactly by
    srsRAN), in place. For every frame with an even SFN: the rebuilt grid of its redundancy version, through the
    channel measured on the cell's reference signals in that subframe (interpolated in frequency, and in time
    between the four reference-signal symbols), amplitude and phase fitted by least squares."""
    import scipy.fft as sf
    sfn0, frame0 = sfn_info
    sc = np.arange(12 * n_rb)
    b = bins_for(sc, n_rb, nfft)
    t_pilot = np.array([0, 4, 7, 11])
    frames = sorted(set(int(j) // 20 for j in js if (j % 20) in (10, 11)))
    for fr in frames:
        sfn = (sfn0 + fr - frame0) % 1024
        if sfn % 2:
            continue
        k = (sfn // 2) % 4
        rows = {ns: np.nonzero(js == fr * 20 + ns)[0] for ns in (10, 11)}
        if any(len(r) == 0 for r in rows.values()):
            continue
        r10, r11 = int(rows[10][0]), int(rows[11][0])
        expected = np.zeros((14, len(sc)), np.complex64)
        ok = True
        for port in range(ports):
            pil = []
            for ns, l in ((10, 0), (10, 4), (11, 0), (11, 4)):
                entry = sf5_h.get((ns, l, port))
                if entry is None:
                    ok = False
                    break
                kk, h_rows, row_index = entry
                pos = np.nonzero(row_index == (r10 if ns == 10 else r11))[0]
                if len(pos) == 0:
                    ok = False
                    break
                pil.append(_interp_channel(h_rows[pos], kk, sc)[0])
            if not ok:
                break
            pil = np.array(pil)                                                # (4, subcarriers)
            h = np.empty((14, len(sc)), np.complex64)
            for t in range(14):
                i1 = int(np.clip(np.searchsorted(t_pilot, t), 1, 3))
                w = np.clip((t - t_pilot[i1 - 1]) / (t_pilot[i1] - t_pilot[i1 - 1]), 0, 1)
                h[t] = pil[i1 - 1] * (1 - w) + pil[i1] * w
            expected += h * grids[k, port]
        if not ok:
            continue
        starts = []
        for t in range(14):
            ns, l = (10, t) if t < 7 else (11, t - 7)
            starts.append(int(slot_starts[r10 if ns == 10 else r11]) + symbol_offset(l, nfft, cp0, cp))
        starts = np.array(starts)
        rx = sf.fft(y[starts[:, None] + np.arange(nfft)[None, :]], axis=1, workers=-1)[:, b]
        m = np.abs(expected) > 0
        alpha = np.vdot(expected[m], rx[m]) / (np.vdot(expected[m], expected[m]).real + 1e-12)
        delta = np.zeros((14, nfft), np.complex64)
        delta[:, b] = alpha * expected
        d = sf.ifft(delta, axis=1, workers=-1).astype(np.complex64)
        for t in range(14):
            cpl = cp0 if t % 7 == 0 else cp
            y[starts[t]:starts[t] + nfft] -= d[t]
            y[starts[t] - cpl:starts[t]] -= d[t][nfft - cpl:]


def _cancel_sync(y, js, slot_starts, pci, n_rb, nfft, cp0, cp):
    """Remove the cell's PSS (symbol 6) and SSS (symbol 5) of slots 0 and 10, in place. Both are known from the
    PCI; the channel is measured on the sequence itself and smoothed across frequency. A strong cell's sync signals
    repeat every 5 ms like a weak cell's, so they otherwise pull the weak cell's timing search towards them."""
    import scipy.fft as sf
    import lte_sync as S
    nid1, nid2 = pci // 3, pci % 3
    kk = 6 * n_rb - 31 + np.arange(62)
    b = bins_for(kk, n_rb, nfft)
    pss = S._pss_freq(nid2).astype(np.complex64)
    for l, seqs in ((6, {0: pss, 10: pss}),
                    (5, {0: S.sss_freq(nid1, nid2, False).astype(np.complex64),
                         10: S.sss_freq(nid1, nid2, True).astype(np.complex64)})):
        for ns, ref in seqs.items():
            rows = np.nonzero((js % 20) == ns)[0]
            if len(rows) == 0:
                continue
            useful = slot_starts[rows] + symbol_offset(l, nfft, cp0, cp)
            idx = useful[:, None] + np.arange(nfft)[None, :]
            spec = sf.fft(y[idx], axis=1, workers=-1)
            h = _smooth_time(_smooth_freq(spec[:, b] / ref[None, :], width=7))
            delta = np.zeros_like(spec)
            delta[:, b] = ref[None, :] * h
            d = sf.ifft(delta, axis=1, workers=-1).astype(np.complex64)
            y[idx] -= d
            y[(useful - cp)[:, None] + np.arange(cp)[None, :]] -= d[:, nfft - cp:]


def cancel_crs(samples, rate, start, pci, n_rb, ports=2, cfo_hz=0.0, sync=True, mib_hex=None, sib1_grids=None):
    """Return a copy of samples with cell pci's CRS removed.

    start:  sample index of the start of a subframe 0 of that cell (its frame timing)
    n_rb:   cell width in resource blocks (limited to what the recording's rate can hold)
    ports:  number of transmit antenna ports (1, 2 or 4, from the MIB decode)
    cfo_hz: the cell's frequency offset in the recording
    sync:   also remove the cell's PSS and SSS
    mib_hex: the cell's decoded MIB (6 hex digits): also remove its PBCH (1- and 2-port cells)
    sib1_grids: the cell's SIB1 subframes from regen_sib1 (needs mib_hex): also remove SIB1 with its control message
    """
    import scipy.fft as sf
    nfft, cp0, cp, slot_len = grid(rate)
    n_rb = int(min(n_rb, {128: 6, 256: 15, 512: 25, 1024: 50, 1536: 75, 2048: 100}.get(nfft, n_rb)))
    ports = ports if ports in (1, 2, 4) else 2
    x = np.asarray(samples, np.complex64)
    if abs(cfo_hz) >= 1:
        ph = np.exp(-2j * np.pi * float(cfo_hz) * np.arange(len(x)) / float(rate)).astype(np.complex64)
        y = x * ph
    else:
        ph = None
        y = x.copy()
    first = -(start // slot_len)
    last = (len(y) - start) // slot_len - 1
    js = np.arange(first, last + 1)
    slot_starts = start + js * slot_len
    ok = (slot_starts >= 0) & (slot_starts + slot_len <= len(y))
    js, slot_starts = js[ok], slot_starts[ok]
    if len(js) == 0:
        return x.copy()
    ns_all = (js % 20).astype(int)
    vshift = pci % 6
    m = np.arange(2 * n_rb)
    mprime = m + N_RB_MAX - n_rb
    slot1_h = {}                                     # port -> (pilot subcarriers, channel in slot 1 of subframe 0)
    sf5_h = {}                                       # (slot, symbol, port) -> (pilot subcarriers, channel, row index)
    for l in (0, 1, 4):
        port_list = crs_ports_at(l, ports)
        if not port_list:
            continue
        cpl = cp0 if l == 0 else cp
        useful = slot_starts + symbol_offset(l, nfft, cp0, cp)
        idx = useful[:, None] + np.arange(nfft)[None, :]
        spec = sf.fft(y[idx], axis=1, workers=-1)
        delta = np.zeros_like(spec)
        xtab = np.stack([crs_sequence(pci, ns, l)[mprime] for ns in range(20)])     # (20, 2*n_rb)
        for port in port_list:
            for parity in ((0, 1) if port >= 2 else (None,)):
                rows = np.arange(len(js)) if parity is None else np.nonzero(ns_all % 2 == parity)[0]
                if len(rows) == 0:
                    continue
                v = crs_shift(port, l, parity if parity is not None else 0)
                kk = 6 * m + (v + vshift) % 6
                b = bins_for(kk, n_rb, nfft)
                ref = xtab[ns_all[rows]]
                h = spec[rows][:, b] / ref
                h = _smooth_time(_smooth_freq(h))
                delta[np.ix_(rows, b)] += ref * h
                if l == 0 and parity is None and port < 2:
                    sel = ns_all[rows] == 1
                    slot1_h[port] = (kk, h[sel])
                if sib1_grids is not None and parity is None and port < 2:
                    for ns in (10, 11):
                        sel = np.nonzero(ns_all[rows] == ns)[0]
                        sf5_h[(ns, l, port)] = (kk, h[sel], rows[sel])
        d = sf.ifft(delta, axis=1, workers=-1).astype(np.complex64)
        y[idx] -= d
        cp_idx = (useful - cpl)[:, None] + np.arange(cpl)[None, :]
        y[cp_idx] -= d[:, nfft - cpl:]
    if sync:
        _cancel_sync(y, js, slot_starts, pci, n_rb, nfft, cp0, cp)
    sfn_info = None
    if mib_hex:
        sfn_info = _cancel_pbch(y, js, slot_starts, pci, n_rb, nfft, cp0, cp, ports, mib_hex, slot1_h)
    if sib1_grids is not None and sfn_info is not None and ports in (1, 2):
        _cancel_sib1(y, js, slot_starts, n_rb, nfft, cp0, cp, ports, sib1_grids, sf5_h, sfn_info)
    if ph is not None:
        y = y * np.conj(ph)
    return y.astype(np.complex64)


# ----------------------------------------------------------------------------------------------- self-check
def _synth(rate, pci, n_rb, ports, seconds, start, gain, cfo_hz=0.0, data_load=0.0, seed=0):
    """CRS (+ optional random QPSK data on other resource elements) of one cell, through a 2-tap channel."""
    rng = np.random.default_rng(seed)
    nfft, cp0, cp, slot_len = grid(rate)
    n = int(rate * seconds)
    out = np.zeros(n + slot_len, np.complex64)
    vshift, m = pci % 6, np.arange(2 * n_rb)
    mprime = m + N_RB_MAX - n_rb
    taps = [(0, 1.0), (3, 0.4j)]
    j = 0
    while start + (j + 1) * slot_len <= n:
        ns, s0 = j % 20, start + j * slot_len
        for l in range(7):
            spec = np.zeros(nfft, np.complex64)
            if data_load:
                kk = np.arange(12 * n_rb)
                q = (rng.integers(0, 2, len(kk)) * 2 - 1 + 1j * (rng.integers(0, 2, len(kk)) * 2 - 1)) / np.sqrt(2)
                q[rng.random(len(kk)) > data_load] = 0
                spec[bins_for(kk, n_rb, nfft)] = q
            for port in crs_ports_at(l, ports):
                kk = 6 * m + (crs_shift(port, l, ns) + vshift) % 6
                spec[bins_for(kk, n_rb, nfft)] = crs_sequence(pci, ns, l)[mprime] * (1.0 if port % 2 == 0 else 0.7j)
            t = np.fft.ifft(spec)
            for delay, a in taps:
                t2 = np.roll(t, delay) * a
                cpl = cp0 if l == 0 else cp
                at = s0 + symbol_offset(l, nfft, cp0, cp) - cpl
                out[at:at + cpl + nfft] += gain * np.concatenate((t2[-cpl:], t2))
        j += 1
    out = out[:n]
    if cfo_hz:
        out *= np.exp(2j * np.pi * cfo_hz * np.arange(n) / rate).astype(np.complex64)
    return out


def _crs_power(x, rate, start, pci, n_rb, ports, cfo_hz=0.0):
    """Mean power on one cell's CRS resource elements (to measure suppression)."""
    nfft, cp0, cp, slot_len = grid(rate)
    if cfo_hz:
        x = x * np.exp(-2j * np.pi * cfo_hz * np.arange(len(x)) / rate).astype(np.complex64)
    p, cnt = 0.0, 0
    for j in range(0, (len(x) - start) // slot_len - 1, 3):
        s0 = start + j * slot_len
        Y = np.fft.fft(x[s0 + symbol_offset(0, nfft, cp0, cp): s0 + symbol_offset(0, nfft, cp0, cp) + nfft])
        kk = 6 * np.arange(2 * n_rb) + (0 + pci % 6) % 6
        p += float(np.mean(np.abs(Y[bins_for(kk, n_rb, nfft)]) ** 2))
        cnt += 1
    return p / max(1, cnt)


if __name__ == "__main__":
    # a strong cell (PCI 311) and a cell 20 dB weaker (PCI 199) on one carrier, different frame timing, both with
    # 10% data load; the weak cell's signal-to-interference on its own reference signals, before and after
    rate, n_rb = 15.36e6, 50
    strong = _synth(rate, 311, n_rb, 2, 0.1, 1234, 1.0, cfo_hz=150.0, data_load=0.1, seed=1)
    weak = _synth(rate, 199, n_rb, 2, 0.1, 5000, 0.1, cfo_hz=150.0, data_load=0.1, seed=2)
    noise = (np.random.default_rng(3).normal(size=len(strong)) +
             1j * np.random.default_rng(4).normal(size=len(strong))).astype(np.complex64) * 0.002
    mix = strong + weak + noise
    clean = cancel_crs(mix, rate, 1234, 311, n_rb, 2, cfo_hz=150.0)
    removed = 10 * np.log10(_crs_power(mix, rate, 1234, 311, n_rb, 2, 150.0) /
                            _crs_power(clean, rate, 1234, 311, n_rb, 2, 150.0))
    sig = _crs_power(weak, rate, 5000, 199, n_rb, 2, 150.0)
    sir_before = 10 * np.log10(sig / _crs_power(mix - weak, rate, 5000, 199, n_rb, 2, 150.0))
    sir_after = 10 * np.log10(sig / _crs_power(clean - weak, rate, 5000, 199, n_rb, 2, 150.0))
    # damage to the weak cell: cancel PCI 311 from a recording that holds only the weak cell and noise
    only_weak = weak + noise
    damage = 10 * np.log10(_crs_power(cancel_crs(only_weak, rate, 1234, 311, n_rb, 2, 150.0) - only_weak,
                                      rate, 5000, 199, n_rb, 2, 150.0) / sig)
    print("strong cell's reference signals removed: %.1f dB" % removed)
    print("weak cell, signal-to-interference on its reference signals: %.1f dB -> %.1f dB" % (sir_before, sir_after))
    print("what the cancellation takes from the weak cell: %.1f dB below its signal" % -damage)
    for load in (0.0, 0.3):
        s2 = _synth(rate, 311, n_rb, 2, 0.1, 1234, 1.0, cfo_hz=150.0, data_load=load, seed=1)
        w2 = _synth(rate, 199, n_rb, 2, 0.1, 5000, 0.1, cfo_hz=150.0, data_load=load, seed=2)
        m2 = s2 + w2 + noise
        c2 = cancel_crs(m2, rate, 1234, 311, n_rb, 2, cfo_hz=150.0)
        sg = _crs_power(w2, rate, 5000, 199, n_rb, 2, 150.0)
        print("  with %2.0f%% network load: signal-to-interference %.1f dB -> %.1f dB" % (
            100 * load, 10 * np.log10(sg / _crs_power(m2 - w2, rate, 5000, 199, n_rb, 2, 150.0)),
            10 * np.log10(sg / _crs_power(c2 - w2, rate, 5000, 199, n_rb, 2, 150.0))))
    print("OK" if removed > 15 and sir_after > sir_before + 1 and damage < -15 else "CHECK")
