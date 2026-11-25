#!/usr/bin/env python3
"""Faster drop-in versions of two functions in "osprey 2"/lte_sync.py (the LTE PSS/SSS cell search).
Same maths, same answers (checked against the originals: python3 lte_fastsync.py), less time:

  _confirm_sss   the original rebuilds the 336 SSS reference sequences ~20,000 times per search and scores
                 them one by one in Python; here they are built once per PSS group and all 168 x 2
                 hypotheses are scored with one matrix product per (gap, delta).
  _correlate_mag the original uses numpy's single-threaded FFT; here scipy.fft on all CPU cores.

Owen's file is not edited: install() swaps the functions in at runtime (lte_decode.py --slow-search skips it).
"""
import numpy as np

_REF = {}


def _refs(nid2):
    """(R0, R1): 168 x 62 SSS references for subframe 0 / subframe 5, built once per PSS group."""
    import lte_sync as S
    if nid2 not in _REF:
        _REF[nid2] = (np.array([S.sss_freq(n1, nid2, False) for n1 in range(168)]),
                      np.array([S.sss_freq(n1, nid2, True) for n1 in range(168)]))
    return _REF[nid2]


def _seqs(samples, starts, offset_hz):
    """SSS/PSS-style 62-subcarrier sequences for many window starts at once (None-equivalent windows dropped).
    Mirrors S._extract_seq: derotation uses the absolute sample index."""
    import lte_sync as S
    starts = np.asarray(starts, dtype=np.int64)
    ok = (starts >= 0) & (starts + S.FFT_SIZE <= len(samples))
    starts = starts[ok]
    if starts.size == 0:
        return np.zeros((0, 62), dtype=np.complex128), ok
    idx = starts[:, None] + np.arange(S.FFT_SIZE)[None, :]
    chunk = np.asarray(samples)[idx]
    if offset_hz != 0:
        chunk = chunk * np.exp(-1j * 2 * np.pi * offset_hz * idx / S.FS)
    spec = np.fft.fft(chunk, axis=1)
    seq = np.empty((spec.shape[0], 62), dtype=np.complex128)
    seq[:, :31] = spec[:, -31:]
    seq[:, 31:] = spec[:, 1:32]
    return seq, ok


def _match_rows(seqs, refs):
    """mean over rows of |vdot(ref, seq)| / (|seq| sqrt 62) for every ref row -> (n_refs,). = S._mean_match."""
    if seqs.shape[0] == 0:
        return np.zeros(refs.shape[0])
    norm = np.linalg.norm(seqs, axis=1) * np.sqrt(62.0) + 1e-12
    m = np.abs(seqs @ np.conj(refs).T) / norm[:, None]          # (n_obs, n_refs); vdot(ref, obs)
    return m.mean(axis=0)


def confirm_sss_fast(samples, mag, anchor, offset_hz, nid2):
    import lte_sync as S
    times = S._pss_times(mag, anchor)
    reference = S._pss_freq(nid2)[None, :]
    best_delta, best_pss = 0, -1.0
    for delta in range(-10, 11):
        seqs, _ok = _seqs(samples, [t + delta for t in times], offset_hz)
        score = float(_match_rows(seqs, reference)[0]) if seqs.shape[0] else 0.0
        if score > best_pss:
            best_pss, best_delta = score, delta
    R0, R1 = _refs(nid2)
    best_nid1, best_sss = 0, -1.0
    origin = times[0]
    slots = [int(round((t - origin) / float(S._PERIOD))) % 2 for t in times]
    for gap in (S.SYMBOL, S._EXT_GAP):
        for delta in (best_delta - 1, best_delta, best_delta + 1):
            seqs, ok = _seqs(samples, [t - gap + delta for t in times], offset_hz)
            par = np.array(slots)[ok]
            g = [seqs[par == 0], seqs[par == 1]]
            if g[0].shape[0] == 0 and g[1].shape[0] == 0:
                continue
            # score[nid1, swap]: swap False -> group0 vs subframe0 refs, group1 vs subframe5 refs
            parts_f, parts_t = [], []
            if g[0].shape[0]:
                parts_f.append(_match_rows(g[0], R0))
                parts_t.append(_match_rows(g[0], R1))
            if g[1].shape[0]:
                parts_f.append(_match_rows(g[1], R1))
                parts_t.append(_match_rows(g[1], R0))
            sc = np.stack([np.mean(parts_f, axis=0), np.mean(parts_t, axis=0)], axis=1)   # (168, 2)
            flat = sc.reshape(-1)                        # nid1-major, swap-minor = the original loop order
            k = int(np.argmax(flat))
            if flat[k] > best_sss:
                best_sss, best_nid1 = float(flat[k]), k // 2
    return best_nid1, float(best_sss), float(best_pss)


_TPL = {}


def correlate_mag_fast(spectrum, template, length):
    """The template FFTs are the same for every search (3 PSS x the same offset grid at the same FFT size),
    so they are cached; the per-search work is one inverse FFT per (root, offset)."""
    import scipy.fft as sf
    key = (spectrum.size, template.tobytes())
    t = _TPL.get(key)
    if t is None:
        t = np.conj(sf.fft(template, spectrum.size, workers=-1))
        if len(_TPL) > 400:
            _TPL.clear()
        _TPL[key] = t
    corr = sf.ifft(spectrum * t, workers=-1)
    return np.abs(corr[:length])


ORIG = {}


def install():
    import lte_sync as S
    if ORIG:
        return
    ORIG["confirm"] = S._confirm_sss
    ORIG["corr"] = S._correlate_mag
    S._confirm_sss = confirm_sss_fast
    S._correlate_mag = correlate_mag_fast


def uninstall():
    import lte_sync as S
    if ORIG:
        S._confirm_sss = ORIG.pop("confirm")
        S._correlate_mag = ORIG.pop("corr")


if __name__ == "__main__":
    # self-check: fast == original on synthetic cells + noise, and the speed-up
    import os
    import sys
    import time
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(here, "osprey 2"))
    from lte_cellid import synthesize_cell
    from lte_pci_search import search_samples
    rng = np.random.default_rng(1)

    def up(f):
        F = np.fft.fft(f)
        n = len(f)
        u = np.zeros(n * 8, complex)
        u[:n // 2] = F[:n // 2]
        u[-n // 2:] = F[n // 2:]
        return np.fft.ifft(u) * 8
    cases = []
    for pcis, amps, noise in (((1,), (1.0,), 0.003), ((1, 0), (1.0, 0.5), 0.003), ((21, 216), (1.0, 0.3), 0.01),
                              ((272,), (0.05,), 0.015), ((), (), 0.02)):
        x = (rng.normal(size=153600 * 20) + 1j * rng.normal(size=153600 * 20)) * noise
        for k, (p, a) in enumerate(zip(pcis, amps)):
            x = x + a * np.roll(up(np.tile(synthesize_cell(p, 1, n_rb=6), 20)), 7000 * k)
        cases.append(x.astype(np.complex64)[:int(15.36e6 * 0.08)])
    ok = True
    t_slow = t_fast = 0.0
    for x in cases:
        uninstall()
        t0 = time.time()
        a = search_samples(x, 15.36e6)
        t_slow += time.time() - t0
        install()
        t0 = time.time()
        b = search_samples(x, 15.36e6)
        t_fast += time.time() - t0
        ka = [(h["pci"], round(h["pss_ratio"], 6), round(h["sss_score"], 6)) for h in a["hits"]]
        kb = [(h["pci"], round(h["pss_ratio"], 6), round(h["sss_score"], 6)) for h in b["hits"]]
        same = ka == kb and abs(a["best_ratio"] - b["best_ratio"]) < 1e-6
        ok &= same
        print("%s  original %s  fast %s" % ("SAME" if same else "DIFFERENT", ka, kb))
    print("time: original %.2fs  fast %.2fs  (%.1fx)  -> %s" % (t_slow, t_fast, t_slow / max(t_fast, 1e-9),
                                                              "OK" if ok else "MISMATCH"))
    sys.exit(0 if ok else 1)
