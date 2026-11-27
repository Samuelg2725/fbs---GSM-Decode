#!/usr/bin/env python3
"""5G NR cell identity decoder for the BladeRF.

  python3 nr_decode.py                         # UK 5G bands n78 (3.3-3.8 GHz) and n28 (700 MHz)
  python3 nr_decode.py --bands n78             # one band (n78 n77 n28 n1 n3 n7 n8 n20 n40)
  python3 nr_decode.py --freq-mhz 3549.12      # one SSB frequency (or a few)
  python3 nr_decode.py --dual                  # both antennas (search and decode combined per cell)
  python3 nr_decode.py --loop                  # again and again; cells already identified get a quick check
  python3 nr_decode.py --selftest              # synthetic 5G cells through the whole decoder (no radio)

What it does, per band:
  1. records the band in ~22 MHz steps at 30.72 Msps (gain set automatically, remembered per step);
  2. searches every synchronisation raster position (GSCN) in each recording for the 5G sync signal (PSS),
     with both antennas added together with --dual; every peak is checked with the SSS, which gives the PCI
     (0-1007) - several cells on one SSB frequency and every beam (SSB index) of a cell are kept;
  3. decodes each cell's MIB (srsRAN's PBCH decoder, via the nr_decode helper): frame number, beam, whether the
     cell broadcasts a SIB1 and where its control channel (CORESET 0) is;
  4. for cells with a SIB1 (standalone 5G): finds the SI-RNTI control message in CORESET 0 in every slot of the
     recording and decodes SIB1 -> NR Cell ID (NCI, 36 bits), TAC and operator (MCC-MNC).
Cells without a SIB1 are non-standalone 5G (added on top of 4G): their PCI, beam and signal are shown, but they
broadcast no Cell ID - no receiver (phones included) can read one over the air.

Results: a table, ~/field/nr_cells.csv, and the MIB / SIB1 messages in ~/field/pcap/nr_<time>.pcap (Wireshark
decodes them as NR RRC). Run from the repository folder with LD_LIBRARY_PATH=/usr/local/lib.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", str(max(1, (os.cpu_count() or 4) - 1)))
os.environ.setdefault("OMP_NUM_THREADS", os.environ["OPENBLAS_NUM_THREADS"])
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                                   # noqa: E402
import lte_decode as L                               # noqa: E402  (radio, gain, srsRAN location)
import pcap_out                                      # noqa: E402

FS = 30.72e6                                         # recording rate
# usable part of a recording: the radio's filter is 0.9 x the rate wide and the radio is tuned --low-if-mhz (0.6) above
# the centre; an SS/PBCH block (7.2 MHz at 30 kHz, 3.6 MHz at 15 kHz) must lie fully inside, with 0.2 MHz to spare
USABLE_HALF = 0.9 * 30.72e6 / 2 - 0.6e6 - 0.2e6       # 13.0 MHz either side of the centre, before the SSB's width
TUNE_STEP = 2 * (USABLE_HALF - 3.6e6)                # 18.8 MHz: every 30 kHz SSB is fully inside one recording
# band: SSB subcarrier spacing (kHz), SSB pattern case, max beams, GSCN range, duplex
BANDS = {
    "n78": dict(scs=30, case="C", lmax=8, gscn=(7711, 8051), tdd=True),
    "n77": dict(scs=30, case="C", lmax=8, gscn=(7711, 8329), tdd=True),
    "n40": dict(scs=30, case="C", lmax=4, gscn=(5762, 5989), tdd=True),
    "n28": dict(scs=15, case="A", lmax=4, gscn=(1901, 2002), tdd=False),
    "n20": dict(scs=15, case="A", lmax=4, gscn=(1982, 2047), tdd=False),
    "n8": dict(scs=15, case="A", lmax=4, gscn=(2318, 2395), tdd=False),
    "n3": dict(scs=15, case="A", lmax=4, gscn=(4517, 4693), tdd=False),
    "n1": dict(scs=15, case="A", lmax=4, gscn=(5279, 5419), tdd=False),
    "n7": dict(scs=15, case="A", lmax=4, gscn=(6554, 6718), tdd=False),
}
UK_DEFAULT = ["n78", "n28"]
# TS 38.213 tables 13-1..13-4: CORESET 0 (resource blocks, symbols, offset in RBs) per {SSB, PDCCH} kHz and index
CORESET0 = {
    (15, 15): [(24, 2, 0), (24, 2, 2), (24, 2, 4), (24, 3, 0), (24, 3, 2), (24, 3, 4), (48, 1, 12), (48, 1, 16),
               (48, 2, 12), (48, 2, 16), (48, 3, 12), (48, 3, 16), (96, 1, 38), (96, 2, 38), (96, 3, 38), None],
    (15, 30): [(24, 2, 5), (24, 2, 6), (24, 2, 7), (24, 2, 8), (24, 3, 5), (24, 3, 6), (24, 3, 7), (24, 3, 8),
               (48, 1, 18), (48, 1, 20), (48, 2, 18), (48, 2, 20), (48, 3, 18), (48, 3, 20), None, None],
    (30, 15): [(48, 1, 2), (48, 1, 6), (48, 2, 2), (48, 2, 6), (48, 3, 2), (48, 3, 6), (96, 1, 28), (96, 2, 28),
               (96, 3, 28), None, None, None, None, None, None, None],
    (30, 30): [(24, 2, 0), (24, 2, 1), (24, 2, 2), (24, 2, 3), (24, 2, 4), (24, 3, 0), (24, 3, 1), (24, 3, 2),
               (24, 3, 3), (24, 3, 4), (48, 1, 12), (48, 1, 14), (48, 1, 16), (48, 2, 12), (48, 2, 14), (48, 2, 16)],
}
PCAP = [None]
MIB_RE = re.compile(r"MIB grid=(\d+) pci=(\d+) ssb_idx=(\d+) hrf=(\d+) sfn=(\d+) scs_common=(\d+) kssb=(\d+) "
                    r"coreset0=(\d+) ss0=(\d+) typeA=(\d+) barred=(\d+) intra=(\d+) corr=([-\d.na]+) "
                    r"cfo=([-\d.na]+) bch=([0-9a-f]+)")


# ------------------------------------------------------------------------------------------ raster / frequencies
def gscn_to_mhz(g):
    """Synchronisation raster (TS 38.104 5.4.3.1): GSCN -> SSB centre frequency (MHz)."""
    if g < 7499:
        n, m = (g + 1) // 3, (g + 1) % 3 * 2 + 1          # GSCN = 3N + (M-3)/2, M in {1,3,5}
        return (n * 1200e3 + m * 50e3) / 1e6
    return 3000.0 + (g - 7499) * 1.44


def nr_arfcn(mhz):
    """Global frequency raster (TS 38.104 5.4.2.1): frequency -> NR-ARFCN."""
    if mhz < 3000:
        return int(round(mhz * 1e3 / 5))
    return int(round(600000 + (mhz - 3000) * 1e3 / 15))


def band_of(mhz):
    for name, b in BANDS.items():
        if gscn_to_mhz(b["gscn"][0]) - 2 <= mhz <= gscn_to_mhz(b["gscn"][1]) + 2:
            return name
    return None


# ------------------------------------------------------------------------------------------ PSS / SSS (TS 38.211)
def _mseq(taps, init):
    x = list(init)
    for i in range(127 - 7):
        x.append((x[i + taps[0]] + x[i]) % 2)
    return np.array(x)


_XPSS = _mseq((4,), [0, 1, 1, 0, 1, 1, 1])
_X0 = _mseq((4,), [1, 0, 0, 0, 0, 0, 0])
_X1 = _mseq((1,), [1, 0, 0, 0, 0, 0, 0])


def pss_seq(nid2):
    n = np.arange(127)
    return 1.0 - 2.0 * _XPSS[(n + 43 * nid2) % 127]


def sss_seq(nid1, nid2):
    n = np.arange(127)
    m0 = 15 * (nid1 // 112) + 5 * nid2
    m1 = nid1 % 112
    return (1.0 - 2.0 * _X0[(n + m0) % 127]) * (1.0 - 2.0 * _X1[(n + m1) % 127])


SSS_ALL = {n2: np.array([sss_seq(n1, n2) for n1 in range(336)]) for n2 in range(3)}
PSS_BINS = np.arange(127) - 64                       # PSS/SSS on SSB subcarriers 56..182, SSB centre = 120


def pss_template(nid2, nfft):
    """PSS symbol (useful part, no CP) in the time domain at nfft samples per symbol."""
    spec = np.zeros(nfft, complex)
    spec[PSS_BINS % nfft] = pss_seq(nid2)
    return np.fft.ifft(spec) * np.sqrt(nfft)


# ------------------------------------------------------------------------------------------ OFDM timing
def cp_len(n, j, mu, rate=FS):
    """Cyclic prefix (samples) of symbol j of a slot: TS 38.211 5.3.1 (long CP every 0.5 ms)."""
    base = int(round(n * 144 / 2048))
    return base + int(round(16 * rate / 30.72e6)) if j % (7 * 2 ** mu) == 0 else base


def symbol_starts(n, mu, rate=FS):
    """Useful-part start (samples from the slot start) of the 14 symbols of a slot, and the slot length."""
    out, t = [], 0
    for j in range(14):
        c = cp_len(n, j, mu, rate)
        out.append(t + c)
        t += c + n
    return out, t


def ssb_first_symbols(case, lmax):
    """First symbol of each SSB candidate in the half frame (TS 38.213 4.1)."""
    if case == "C":
        base = [2, 8]
        return [b + 14 * k for k in range(lmax // 2) for b in base][:lmax]
    return [b + 14 * k for k in range(lmax // 2) for b in [2, 8]][:lmax]   # case A


def halfframe_offset(l, n, mu, rate=FS):
    """Samples from the half-frame start to the useful start of symbol l of the half frame."""
    starts, slot_len = symbol_starts(n, mu, rate)
    return (l // 14) * slot_len + starts[l % 14]


# ------------------------------------------------------------------------------------------ helper program
def helper():
    path = L.pdsch_ue().parent / "nr_decode"
    if not path.is_file():
        sys.exit("The nr_decode helper is not built (%s).\nRe-run the srsRAN patch and build it:\n"
                 "  python3 ~/fbs-upload/\"osprey 2\"/srsran_watch/print_sib1.py && cd ~/srsRAN_4G/build && "
                 "cmake .. && make -j4 pdsch_ue sib1_regen nr_decode" % path)
    return str(path)


def run_helper(args, timeout=60):
    try:
        out = subprocess.run([helper()] + [str(a) for a in args], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             timeout=timeout, preexec_fn=L._srsran_child)
        return out.stdout.decode(errors="replace")
    except subprocess.TimeoutExpired:
        return ""


def write_cf32(arrays):
    fh = tempfile.NamedTemporaryFile(prefix="nr_", suffix=".cf32", delete=False, dir="/dev/shm")
    for a in arrays:
        fh.write(np.asarray(a, np.complex64).tobytes())
    fh.close()
    return fh.name


# ------------------------------------------------------------------------------------------ search
class Band(object):
    """One SSB frequency's slice of a recording: spectrum bins around it, decimated signal on demand."""

    def __init__(self, spectra, nx, offset_hz, scs_khz):
        self.scs = scs_khz * 1e3
        self.nfft = 256                              # 256 subcarriers: 7.68 Msps (30 kHz) / 3.84 Msps (15 kHz)
        self.rate = self.nfft * self.scs
        m = int(round(nx * self.rate / FS))
        b0 = int(round(offset_hz / FS * nx))
        idx = (b0 + np.round(np.fft.fftfreq(m) * m).astype(int)) % nx
        self.spec = [S[idx] for S in spectra]        # per antenna
        self.m, self.nx = m, nx
        self._y = None

    def y(self):
        if self._y is None:
            self._y = [np.fft.ifft(s) * (self.m / self.nx) for s in self.spec]
        return self._y


_TPL = {}


def correlate(band):
    """PSS correlation of each N_ID2 (antennas added) -> list of (nid2, power array)."""
    out = []
    for nid2 in range(3):
        key = (band.m, nid2)
        if key not in _TPL:
            t = np.zeros(band.m, complex)
            t[:band.nfft] = pss_template(nid2, band.nfft)
            _TPL[key] = np.conj(np.fft.fft(t))
        p = 0
        for s in band.spec:
            p = p + np.abs(np.fft.ifft(s * _TPL[key])) ** 2
        out.append((nid2, p))
    return out


def find_peaks(p, thr, sym, max_peaks=24):
    """Local maxima of a correlation power above thr x mean, at least one symbol apart -> [(index, ratio)]."""
    mean = float(np.mean(p)) + 1e-30
    cand = np.nonzero(p > thr * mean)[0]
    if not len(cand):
        return []
    order = cand[np.argsort(p[cand])[::-1]]
    peaks = []
    for i in order:
        if all(abs(int(i) - j) > sym for j, _ in peaks):
            peaks.append((int(i), float(p[i] / mean)))
            if len(peaks) >= max_peaks:
                break
    return peaks


def ssb_grids(band, t, f_abs):
    """SS/PBCH block grids (4 x 240) of each antenna, PSS useful part at decimated sample t, with the per-symbol
    upconversion phase removed (TS 38.211 5.4: each symbol's phase restarts; compensated with the absolute
    frequency of this slice's centre)."""
    n = band.nfft
    cp = int(round(n * 144 / 2048))
    w = cp // 2                                       # FFT window half a CP early: no inter-symbol interference
    grids = []
    kk = np.arange(240) - 120
    early = np.exp(2j * np.pi * kk * w / n)            # undo the window starting w samples early
    for y in band.y():
        g = np.zeros((4, 240), complex)
        for l in range(4):
            s = t + l * (n + cp) - w
            if s < 0 or s + n > len(y):
                return None
            spec = np.fft.fft(y[s:s + n]) / np.sqrt(n)
            ph = np.exp(2j * np.pi * f_abs * ((t + l * (n + cp)) / band.rate))
            g[l] = spec[kk % n] * early * ph
        grids.append(g)
    return grids


def sss_detect(grids, nid2, keep=2):
    """N_ID1 from the SSS (equalised with the PSS) of the antennas added -> [(nid1, metric)] best first."""
    acc = np.zeros(336)
    k = np.arange(127) + 56
    d = pss_seq(nid2)
    for g in grids:
        h = g[0, k] * d                                # channel from the PSS (real +-1 sequence)
        z = g[2, k] * np.conj(h)
        acc += np.abs(SSS_ALL[nid2] @ z) ** 2
    order = np.argsort(acc)[::-1]
    best = acc[order[0]] + 1e-30
    out = [(int(order[0]), 1.0)]
    for i in order[1:keep]:
        if acc[i] > 0.35 * best:
            out.append((int(i), float(acc[i] / best)))
    return out


def combine(grids, nid1, nid2):
    """Two-antenna combining weights from PSS + SSS (wideband) -> (combined grid, RX1 share)."""
    if len(grids) == 1:
        return grids[0], None
    k = np.arange(127) + 56
    hs = []
    for g in grids:
        h = np.vdot(pss_seq(nid2), g[0, k]) + np.vdot(sss_seq(nid1, nid2), g[2, k])
        hs.append(h / 254.0)
    norm = np.sqrt(sum(abs(h) ** 2 for h in hs)) + 1e-30
    out = sum(np.conj(h) / norm * g for h, g in zip(hs, grids))
    return out, float(abs(hs[0]) ** 2 / norm ** 2)


def search_recording(x_list, fc, band_name, gscns, args, show):
    """All SSB frequencies of one recording -> list of candidate detections (not yet MIB-checked)."""
    b = BANDS[band_name]
    nx = len(x_list[0]) // 8 * 8
    spectra = [np.fft.fft(x[:nx]) for x in x_list]
    found = []
    for g in gscns:
        f = gscn_to_mhz(g) * 1e6
        off = f - fc
        if abs(off) > USABLE_HALF - (3.6e6 if b["scs"] == 30 else 1.8e6):
            continue
        band = Band(spectra, nx, off, b["scs"])
        sym = band.nfft + int(round(band.nfft * 144 / 2048))
        for nid2, p in correlate(band):
            for t, ratio in find_peaks(p, args.pss_min, sym):
                found.append(dict(gscn=g, mhz=f / 1e6, band=band, nid2=nid2, t=t, ratio=ratio))
    return found


def check_candidates(cands, args, b):
    """SSS + MIB of each candidate peak -> decoded detections (one per PCI and SSB frequency, beams merged)."""
    items, meta = [], []
    for c in cands:
        grids = ssb_grids(c["band"], c["t"], c["mhz"] * 1e6)
        if grids is None:
            continue
        for nid1, rel in sss_detect(grids, c["nid2"]):
            pci = 3 * nid1 + c["nid2"]
            g, share = combine(grids, nid1, c["nid2"])
            items.append(g)
            meta.append(dict(c, pci=pci, share=share, sss_rel=rel))
    if not items:
        return []
    path = write_cf32(items)
    try:
        text = run_helper(["mib", path, b["lmax"], b["scs"]] + [m["pci"] for m in meta])
    finally:
        os.unlink(path)
    cells = {}
    for m in MIB_RE.finditer(text):
        i = int(m.group(1))
        d = meta[i]
        mib = dict(ssb_idx=int(m.group(3)), hrf=int(m.group(4)), sfn=int(m.group(5)), scs_common=int(m.group(6)),
                   kssb=int(m.group(7)), coreset0=int(m.group(8)), ss0=int(m.group(9)), typeA=int(m.group(10)),
                   barred=int(m.group(11)), intra=int(m.group(12)), cfo=float(m.group(14)), bch=m.group(15))
        key = (d["gscn"], d["pci"])
        beam = dict(t=d["t"], ratio=d["ratio"], share=d["share"], **mib)
        if key not in cells:
            cells[key] = dict(gscn=d["gscn"], mhz=d["mhz"], pci=d["pci"], band_obj=d["band"], beams=[])
        cells[key]["beams"].append(beam)
    out = []
    for cell in cells.values():
        cell["beams"].sort(key=lambda bm: -bm["ratio"])
        best = cell["beams"][0]
        cell.update(mib=best, strength=best["ratio"], share=best["share"])
        out.append(cell)
    return out


# ------------------------------------------------------------------------------------------ SIB1
def coreset0_layouts(cell, b):
    """Where CORESET 0 lies (TS 38.213 13): [(lowest subcarrier Hz, resource blocks, symbols, scs Hz)] - with a
    30 kHz common grid the MIB leaves a 180 kHz ambiguity, so both positions are returned."""
    mib = cell["mib"]
    scs_ssb, scs_c = b["scs"], mib["scs_common"]
    table = CORESET0.get((scs_ssb, scs_c))
    if mib["kssb"] > 23 or table is None or table[mib["coreset0"]] is None:
        return []
    nrb, nsym, offset = table[mib["coreset0"]]
    f_ssb0 = cell["mhz"] * 1e6 - 120 * scs_ssb * 1e3           # SSB subcarrier 0
    rbw = 12 * scs_c * 1e3
    out = []
    for j in ((0, 1) if scs_c == 30 else (0,)):
        grid0 = f_ssb0 - mib["kssb"] * 15e3 - j * 180e3        # a common resource block boundary
        crb_lo = grid0 + np.floor((f_ssb0 - grid0) / rbw + 1e-9) * rbw
        out.append((crb_lo - offset * rbw, nrb, nsym, scs_c * 1e3))
    return out


def slot_grids(x_list, fc, frame0, f_low, nrb, scs, cfo, weights, nslots_max=200):
    """Resource grids (14 x nrb*12) of every complete slot in the recording, lowest subcarrier = f_low.
    frame0: sample index of a frame start; returns [(slot_index_in_frame, sample_start, grid)]."""
    mu = 1 if scs == 30e3 else 0
    n = int(round(FS / scs))
    starts, slot_len = symbol_starts(n, mu)
    nsc = nrb * 12
    f_dc = f_low + nsc / 2 * scs                                  # grid centre on DC
    nslot_frame = 10 * 2 ** mu
    t_all = np.arange(len(x_list[0]))
    shift = np.exp(-2j * np.pi * ((f_dc - fc) + cfo) * t_all / FS).astype(np.complex64)
    xs = [x * shift for x in x_list]
    first = frame0 - int(np.floor(frame0 / slot_len)) * slot_len  # earliest slot start >= 0 on the frame's grid
    out = []
    kk = np.arange(nsc) - nsc // 2
    k = kk % n
    w = int(round(n * 144 / 2048)) // 2
    early = np.exp(2j * np.pi * kk * w / n)            # undo the window starting w samples early
    s0 = first
    while s0 + slot_len <= len(xs[0]) and len(out) < nslots_max:
        idx = int(round((s0 - frame0) / slot_len)) % nslot_frame
        grid = np.zeros((14, nsc), complex)
        for x, wgt in zip(xs, weights):
            for l in range(14):
                st = s0 + starts[l] - w
                spec = np.fft.fft(x[st:st + n]) / np.sqrt(n)
                ph = np.exp(2j * np.pi * f_dc * ((s0 + starts[l]) / FS))
                grid[l] += wgt * spec[k] * early * ph
        out.append((idx, s0, grid))
        s0 += slot_len
    return out


def decode_sib1(cell, x_list, fc, b, args):
    """SIB1 of a cell whose MIB says it has one -> fills cell["ident"] (or a reason in cell["sib1_status"])."""
    mib = cell["mib"]
    layouts = coreset0_layouts(cell, b)
    if not layouts:
        cell["sib1_status"] = "no SIB1 (MIB: %s)" % ("non-standalone cell" if mib["kssb"] > 23 else
                                                     "CORESET 0 table missing")
        return
    band = cell["band_obj"]
    mu_ssb = 1 if b["scs"] == 30 else 0
    n_ssb = int(round(FS / (b["scs"] * 1e3)))
    lfirst = ssb_first_symbols(b["case"], b["lmax"])
    # frame start (full-rate samples) from the strongest beam: its PSS, its SSB index and half frame
    t_full = mib["t"] * FS / band.rate
    hf_start = t_full - halfframe_offset(lfirst[mib["ssb_idx"]], n_ssb, mu_ssb)
    frame0 = int(round(hf_start - mib["hrf"] * FS * 0.005))
    share = cell.get("share")
    weights = [1.0] if len(x_list) == 1 or share is None else [np.sqrt(share), np.sqrt(1 - share)]
    if len(x_list) > 1 and share is not None:
        # phase of each antenna from the SSB (combine() weights): recompute from the strongest beam's grids
        grids = ssb_grids(band, mib["t"], cell["mhz"] * 1e6)
        k = np.arange(127) + 56
        hs = [np.vdot(pss_seq(cell["pci"] % 3), g[0, k]) for g in grids]
        norm = np.sqrt(sum(abs(h) ** 2 for h in hs)) + 1e-30
        weights = [np.conj(h) / norm for h in hs]
    tried = []
    for f_low, nrb, nsym, scs in layouts:
        # srsRAN reports the frequency offset with the opposite sign (measured: cfo field = -offset)
        slots = slot_grids(x_list, fc, frame0, f_low, nrb, scs, -mib["cfo"], weights)
        if not slots:
            continue
        # SS/PBCH block REs inside CORESET 0's band, in the slots where this cell's beams were seen
        k_ssb0 = (cell["mhz"] * 1e6 - 120 * b["scs"] * 1e3 - f_low) / scs
        rb0 = int(np.floor(k_ssb0 / 12))
        rbn = int(np.ceil((k_ssb0 + 240 * b["scs"] * 1e3 / scs) / 12)) - rb0
        ssb_syms = {}
        for bm in cell["beams"]:
            t0 = bm["t"] * FS / band.rate
            for l in range(4):
                ts = t0 + l * (n_ssb + int(round(n_ssb * 144 / 2048)))
                for (idx, s0, _g) in slots:
                    starts, slot_len = symbol_starts(int(round(FS / scs)), 1 if scs == 30e3 else 0)
                    if s0 <= ts < s0 + slot_len:
                        j = int(np.argmin([abs(s0 + st - ts) for st in starts]))
                        ssb_syms.setdefault(s0, 0)
                        ssb_syms[s0] |= 1 << j
        spec = ",".join("%d:%d:%d:%x" % (idx, max(rb0, 0), rbn if ssb_syms.get(s0) else 0, ssb_syms.get(s0, 0))
                        for idx, s0, _g in slots)
        path = write_cf32([g for _i, _s, g in slots])
        try:
            text = run_helper(["sib1", path, cell["pci"], int(scs / 1e3), nrb, nsym, mib["typeA"], spec], timeout=120)
        finally:
            os.unlink(path)
        tried.append(text)
        m = re.search(r"CELL slot=(\d+) nci=(\d+) tac=(\S+) plmn=(\S+) reserved=(\d)", text)
        if m:
            mcc, mnc = (m.group(4).split(",")[0].split("-") + ["?"])[:2]
            cell["ident"] = dict(nci=int(m.group(2)), tac=int(m.group(3)) if m.group(3).isdigit() else None,
                                 mcc=mcc, mnc=mnc, plmns=m.group(4), operator=C.operator(mcc, mnc) or "")
            hx = re.search(r"SIB1 slot=\d+ hex=([0-9a-f]+)", text)
            cell["sib1_hex"] = hx.group(1) if hx else None
            return
    dci = sum(t.count("DCI slot=") for t in tried)
    ko = sum(t.count("crc=KO") for t in tried)
    cell["sib1_status"] = ("SIB1 control message found, data failed (%d tries)" % ko if dci else
                           "SIB1 control message not found")


# ------------------------------------------------------------------------------------------ bands
def tune_plan(band_name, freqs=None):
    """(tune centre Hz, [GSCNs inside]) covering the band (or the given SSB frequencies)."""
    b = BANDS[band_name]
    if freqs:
        return [(f * 1e6, [g for g in range(b["gscn"][0], b["gscn"][1] + 1) if abs(gscn_to_mhz(g) - f) < 0.0011]
                 or [None]) for f in freqs]
    lo, hi = gscn_to_mhz(b["gscn"][0]) * 1e6, gscn_to_mhz(b["gscn"][1]) * 1e6
    plan, fc = [], lo + TUNE_STEP / 2 - 1e6
    while fc - TUNE_STEP / 2 < hi:
        gs = [g for g in range(b["gscn"][0], b["gscn"][1] + 1)
              if abs(gscn_to_mhz(g) * 1e6 - fc) <= USABLE_HALF - (3.6e6 if b["scs"] == 30 else 1.8e6)]
        if gs:
            plan.append((fc, gs))
        fc += TUNE_STEP
    return plan


def record(args, fc):
    """Recording at fc (gain from memory, adjusted automatically) -> list of antenna signals or None."""
    mem = L.carrier_memory().get("nr_%d" % int(round(fc / 1e5)), {})
    args.gain_now = float(mem.get("gain", args.gain))
    args.gain_now2 = float(mem.get("gain2", args.gain_now))
    x1, x2 = L.record(args, fc / 1e6, FS)
    x1, x2 = L.auto_gain(args, fc / 1e6, FS, x1, x2)
    if x1 is None:
        return None
    L.remember_carrier("nr_%d" % int(round(fc / 1e5)), 0, args.gain_now, args.gain_now2 if args.dual else None)
    return [x1] if x2 is None else [x1, x2]


def scan_band(args, band_name, freqs=None, signal=None):
    """Search + MIB + SIB1 for one band -> list of cells."""
    b = BANDS[band_name]
    plan = tune_plan(band_name, freqs)
    C.heading("5G %s   %s   %d recording%s" % (band_name, "SSB %s MHz" % ", ".join("%.2f" % f for f in freqs) if freqs
                                               else "%.0f-%.0f MHz" % (gscn_to_mhz(b["gscn"][0]),
                                                                        gscn_to_mhz(b["gscn"][1])),
                                               len(plan), "" if len(plan) == 1 else "s"))
    cells = []
    for fc, gscns in plan:
        t0 = time.time()
        C.status("  %.1f MHz: recording..." % (fc / 1e6))
        x_list = signal(fc) if signal else record(args, fc)
        if x_list is None:
            C.warn("  %.1f MHz: recording failed" % (fc / 1e6))
            continue
        C.status("  %.1f MHz: searching %d SSB positions..." % (fc / 1e6, len(gscns)))
        if gscns == [None]:                                   # --freq-mhz off the raster: search that frequency
            gscns = []
            cands = search_offraster(x_list, fc, b, args)
        else:
            cands = search_recording(x_list, fc, band_name, gscns, args, None)
        C.status("  %.1f MHz: %d sync peak%s, checking..." % (fc / 1e6, len(cands), "" if len(cands) == 1 else "s"))
        found = check_candidates(cands, args, b) if cands else []
        C.clear_status()
        seen = {(c["gscn"], c["pci"]) for c in cells}
        new = [c for c in found if (c["gscn"], c["pci"]) not in seen]
        for cell in new:
            cell["fc"] = fc
            cell["band"] = band_name
            if not args.no_sib1:
                C.status("  PCI %d: decoding SIB1..." % cell["pci"])
                decode_sib1(cell, x_list, fc, b, args)
                C.clear_status()
            show_cell(cell)
            write_pcap(cell)
        cells += new
        if not new:
            print(C.c("  %.1f MHz: %s  (%.1f s)" % (fc / 1e6, "no 5G cell" if not cands else
                                                    "%d sync peak%s, no MIB" % (len(cands), "" if len(cands) == 1
                                                                                 else "s"), time.time() - t0),
                      "grey"), flush=True)
    return cells


def search_offraster(x_list, fc, b, args):
    """--freq-mhz not on the raster: the given frequency is the tune centre; search it directly."""
    nx = len(x_list[0]) // 8 * 8
    spectra = [np.fft.fft(x[:nx]) for x in x_list]
    band = Band(spectra, nx, 0.0, b["scs"])
    sym = band.nfft + int(round(band.nfft * 144 / 2048))
    out = []
    for nid2, p in correlate(band):
        for t, ratio in find_peaks(p, args.pss_min, sym):
            out.append(dict(gscn=0, mhz=fc / 1e6, band=band, nid2=nid2, t=t, ratio=ratio))
    return out


# ------------------------------------------------------------------------------------------ output
def show_cell(c):
    mib = c["mib"]
    beams = sorted({bm["ssb_idx"] for bm in c["beams"]})
    rx = ("  RX1 %.0f%% / RX2 %.0f%%" % (100 * c["share"], 100 * (1 - c["share"]))) if c.get("share") is not None else ""
    head = "  PCI %-4d %8.2f MHz  beam%s %s  strength %.0f" % (c["pci"], c["mhz"], "" if len(beams) == 1 else "s",
                                                              ",".join(map(str, beams)), c["strength"])
    if c.get("ident"):
        x = c["ident"]
        print("%s  %s  NCI %d  TAC %s  %s %s-%s%s" % (head, C.c("decoded", "green"), x["nci"], x["tac"],
                                                     x["operator"] or "?", x["mcc"], x["mnc"], C.c(rx, "grey")),
              flush=True)
    else:
        why = c.get("sib1_status") or ("MIB only" if mib else "")
        print("%s  %s  SFN %d%s%s" % (head, C.c(why, "yellow"), mib["sfn"], "  barred" if mib["barred"] else "",
                                      C.c(rx, "grey")), flush=True)


def write_pcap(c):
    if PCAP[0] is None:
        return
    PCAP[0].nr_rrc("nr-rrc.bcch.bch", bytes.fromhex(c["mib"]["bch"]))
    if c.get("sib1_hex"):
        PCAP[0].nr_rrc("nr-rrc.bcch.dl.sch", bytes.fromhex(c["sib1_hex"]))


def summary(cells, t_all):
    t = C.Table(["SSB MHz", "NR-ARFCN", "Band", "PCI", "Beams", "Operator", "MCC-MNC", "TAC", "NCI", "Type",
                 "Strength", "Antennas"],
                align={k: "right" for k in ("SSB MHz", "NR-ARFCN", "PCI", "TAC", "NCI", "Strength")})
    for c in sorted(cells, key=lambda c: (c["mhz"], c["pci"])):
        x = c.get("ident") or {}
        beams = ",".join(map(str, sorted({bm["ssb_idx"] for bm in c["beams"]})))
        kind = ("standalone (SIB1)", "green") if x else (
            ("non-standalone (no SIB1)", "cyan") if c["mib"]["kssb"] > 23 else (c.get("sib1_status", "MIB only"),
                                                                              "yellow"))
        t.add(["%.2f" % c["mhz"], nr_arfcn(c["mhz"]), c["band"], c["pci"], beams,
               C.op_cell(x.get("operator")) if x else "-", "%s-%s" % (x["mcc"], x["mnc"]) if x else "-",
               x.get("tac", "-") if x else "-", (x["nci"], "green+bold") if x else "-", kind,
               "%.0f" % c["strength"],
               "RX1 %.0f%% / RX2 %.0f%%" % (100 * c["share"], 100 * (1 - c["share"])) if c.get("share") is not None
               else "RX1"])
    t.show("5G cells found: %d, Cell ID decoded: %d   (%.0f s)" % (len(cells), sum(1 for c in cells if c.get("ident")),
                                                                    time.time() - t_all))


def write_csv(path, cells):
    path = os.path.expanduser(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new:
            fh.write("utc,ssb_mhz,nr_arfcn,band,pci,beams,sfn,cell_type,operator,mcc,mnc,tac,nci,strength,status\n")
        for c in cells:
            x = c.get("ident") or {}
            fh.write("%s,%.3f,%d,%s,%d,%s,%d,%s,%s,%s,%s,%s,%s,%.1f,%s\n" % (
                C.utc_iso(), c["mhz"], nr_arfcn(c["mhz"]), c["band"], c["pci"],
                " ".join(str(bm["ssb_idx"]) for bm in c["beams"]), c["mib"]["sfn"],
                "SA" if x else ("NSA" if c["mib"]["kssb"] > 23 else "?"), x.get("operator", ""), x.get("mcc", ""),
                x.get("mnc", ""), x.get("tac", ""), x.get("nci", ""), c["strength"],
                (c.get("sib1_status") or "decoded").replace(",", ";")))
    C.info("results added to %s" % path)


# ------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description="5G NR cell identity decoder (BladeRF + srsRAN).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--bands", nargs="+", choices=sorted(BANDS), help="5G bands (default: %s)" % " ".join(UK_DEFAULT))
    ap.add_argument("--freq-mhz", type=float, nargs="+", help="SSB centre frequencies to decode (instead of a scan)")
    ap.add_argument("--dual", action="store_true", help="use both antennas (RX1 + RX2), combined per cell")
    ap.add_argument("--gain", type=float, default=40.0, help="starting gain in dB (default 40; adjusted automatically)")
    ap.add_argument("--no-auto-gain", action="store_true", help="keep --gain fixed")
    ap.add_argument("--secs", type=float, default=0.055,
                    help="seconds recorded per step (default 0.055: the radio's first 10 ms are dropped, leaving "
                         "45 ms = two 20 ms sync bursts)")
    ap.add_argument("--pss-min", type=float, default=18.0,
                    help="sync peak threshold, x the average (default 18: noise alone peaks at ~14, a cell whose MIB "
                         "still decodes is at 40+; lower = more false peaks to check, slower)")
    ap.add_argument("--no-sib1", action="store_true", help="stop at the MIB (no Cell ID)")
    ap.add_argument("--loop", action="store_true", help="scan again and again (Ctrl-C to stop)")
    ap.add_argument("--low-if-mhz", type=float, default=0.6, help="radio tuned this far above (default 0.6)")
    ap.add_argument("--csv", default="~/field/nr_cells.csv", help="CSV file results are added to ('' = none)")
    ap.add_argument("--selftest", action="store_true", help="decode synthetic 5G cells (no radio needed)")
    pcap_out.add_args(ap)
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    if args.no_color:
        C.set_colour(False)
    if args.selftest:
        return selftest(args)
    helper()
    L.clean_leftovers()
    C.be_nice()
    args.file = None
    try:
        L.RADIO[0] = L.LibRadio(args.dual)
    except Exception as e:
        C.warn("could not keep the radio open (%s) - using bladeRF-cli" % e)
    PCAP[0] = pcap_out.open_from_args(args, "nr", nr=True)
    bands = args.bands or ([band_of(f) for f in args.freq_mhz] if args.freq_mhz else UK_DEFAULT)
    C.stop_on_first_ctrl_c()
    t_all, cells, passes = time.time(), {}, 0
    try:
        while True:
            passes += 1
            if args.loop:
                C.heading("Pass %d   (Ctrl-C to stop)" % passes)
            for band_name in sorted(set(bands), key=bands.index):
                if band_name is None:
                    C.warn("frequency not in a known 5G band")
                    continue
                freqs = [f for f in (args.freq_mhz or []) if band_of(f) == band_name] or None
                for c in scan_band(args, band_name, freqs):
                    key = (round(c["mhz"], 3), c["pci"])
                    if key not in cells or (c.get("ident") and not cells[key].get("ident")):
                        cells[key] = c
                if args.csv:
                    pass
            if not args.loop:
                break
    except KeyboardInterrupt:
        print("")
        C.warn("stopped - showing the results so far")
    finally:
        C.ignore_ctrl_c()
    summary(list(cells.values()), t_all)
    if args.csv and cells:
        write_csv(args.csv, list(cells.values()))
    pcap_out.report(PCAP[0])
    return 0


# ------------------------------------------------------------------------------------------ self-test
def selftest(args):
    """Synthetic recordings (made with srsRAN's own 5G transmitter code, placed in time and frequency the way a
    gNB transmits) through the whole decoder."""
    gen = L.pdsch_ue().parent / "nr_testgen"
    if not gen.is_file():
        sys.exit("nr_testgen is not built (re-run the srsRAN patch and build nr_testgen)")
    rng = np.random.default_rng(1)
    tmp = tempfile.mkdtemp(prefix="nrtest_", dir="/dev/shm")
    results = []
    tests = [
        # name, band, cells: (pci, gscn, kssb, coreset0, sfn, nci, tac, mnc, snr_db, cfo_hz, beams, sib1 slot)
        ("n78: two cells on one SSB frequency, 8 beams", "n78",
         [(123, 7880, 0, 12, 512, 6843456789, 57379, "15", 10, 1500.0, 8, 6),
          (124, 7880, 2, 13, 700, 6843456790, 57379, "15", 4, 1500.0, 4, 8)]),
        ("n78: weak cell, other SSB frequency, odd kssb", "n78",
         [(400, 7890, 7, 10, 33, 12345678901, 3372, "20", -2, -2500.0, 2, 6)]),
        ("n78: non-standalone cell (no SIB1)", "n78",
         [(77, 7880, 30, 0, 100, 0, 0, "10", 8, 300.0, 8, 6)]),
        ("n28: 15 kHz FDD cell", "n28",
         [(301, 1950, 4, 8, 222, 987654321, 50496, "10", 6, 800.0, 4, 3)]),
        ("n78: SIB1 in a slot with two sync blocks (PDSCH around them)", "n78",
         [(222, 7895, 0, 12, 300, 1234567, 4242, "30", 8, 700.0, 8, 1)]),
        ("n78: both antennas, weak cell, RX1 much weaker", "n78",
         [(509, 7870, 4, 13, 900, 34359738367, 1, "15", -5, -900.0, 8, 6)], (0.35, 1.0 * np.exp(1.2j))),
    ]
    for test in tests:
        name, band_name, cells = test[:3]
        gains = test[3] if len(test) > 3 else (1.0,)
        b = BANDS[band_name]
        fc = gscn_to_mhz(cells[0][1]) * 1e6 + 3.3e6
        x_list = synth(gen, tmp, b, fc, cells, rng, gains=gains)
        args.dual, args.pss_min, args.no_sib1 = len(x_list) > 1, 18.0, False
        C.heading("Self-test: %s" % name)
        found = selftest_scan(args, band_name, fc, x_list)
        for spec in cells:
            pci, nci = spec[0], spec[5]
            got = next((c for c in found if c["pci"] == pci), None)
            ok = got is not None and (spec[2] > 23 and not got.get("ident") or
                                      got.get("ident", {}).get("nci") == nci and spec[2] <= 23)
            results.append((name, pci, ok, got))
    t = C.Table(["Test", "PCI", "Result"])
    for name, pci, ok, got in results:
        t.add([name, pci, ("OK", "green") if ok else ("FAILED: %s" % ((got or {}).get("sib1_status") or
                                                                        ("not found" if got is None else "wrong ID")),
                                                        "red")])
    t.show("5G self-test")
    return 0 if all(r[2] for r in results) else 1


def selftest_scan(args, band_name, fc, x_list):
    b = BANDS[band_name]
    gscns = [g for g in range(b["gscn"][0], b["gscn"][1] + 1) if abs(gscn_to_mhz(g) * 1e6 - fc) <= TUNE_STEP / 2]
    cands = search_recording(x_list, fc, band_name, gscns, args, None)
    found = check_candidates(cands, args, b)
    for cell in found:
        cell["fc"], cell["band"] = fc, band_name
        decode_sib1(cell, x_list, fc, b, args)
        show_cell(cell)
    return found


def synth(gen, tmp, b, fc, cells, rng, secs=0.045, gains=(1.0,)):
    """A recording centred at fc holding the given cells: SS/PBCH bursts every 20 ms and SIB1 every 20 ms.
    One signal per antenna gain given (each antenna: the cells times its gain, plus its own noise)."""
    nsamp = int(FS * secs)
    x = np.zeros(nsamp, complex)
    scs_ssb = b["scs"] * 1e3
    mu = 1 if b["scs"] == 30 else 0
    n = int(round(FS / scs_ssb))
    lfirst = ssb_first_symbols(b["case"], b["lmax"])
    starts, slot_len = symbol_starts(n, mu)
    frame0 = int(0.0031 * FS)                                     # a frame starts 3.1 ms into the recording
    noise_p = 0.0
    for (pci, gscn, kssb, c0, sfn, nci, tac, mnc, snr_db, cfo, beams, sib1_slot) in cells:
        out = os.path.join(tmp, "c%d" % pci)
        table = CORESET0[(b["scs"], b["scs"])]
        entry = table[c0] if kssb <= 23 else (48, 1, 0)
        nrb, nsym = entry[0], entry[1]
        f_ssb = gscn_to_mhz(gscn) * 1e6
        extra = []
        if kssb <= 23:                              # sync blocks in the SIB1 slot: the PDSCH goes around them
            cell0 = dict(mhz=f_ssb / 1e6, mib=dict(kssb=kssb, coreset0=c0, scs_common=b["scs"]))
            k0 = (f_ssb - 120 * scs_ssb - coreset0_layouts(cell0, b)[0][0]) / scs_ssb
            mask = 0
            for i in range(beams):
                if lfirst[i] // 14 == sib1_slot:
                    mask |= 0xF << (lfirst[i] % 14)
            if mask:
                rb0 = int(np.floor(k0 / 12))
                extra = [rb0, int(np.ceil((k0 + 240) / 12)) - rb0, "%x" % mask]
        r = subprocess.run([str(a) for a in (gen, out, pci, b["scs"], b["lmax"], sfn, kssb, c0, 0, 2, nrb, nsym,
                                             sib1_slot, nci, tac, "234", mnc)] + [str(e) for e in extra],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if r.returncode:
            sys.exit("nr_testgen failed: %s" % r.stdout.decode())
        ssb = np.fromfile(out + "_ssb.cf32", np.complex64).reshape(b["lmax"], 4, 240)
        slot = np.fromfile(out + "_slot.cf32", np.complex64).reshape(14, nrb * 12)
        amp = 1.0
        sig = np.zeros(nsamp, complex)

        def put(grid_sym, f_low, nsc, scs, t_useful, ncp):
            """One OFDM symbol: REs from f_low upward, useful part starting at sample t_useful (TS 38.211 5.4)."""
            nn = int(round(FS / scs))
            spec = np.zeros(nn, complex)
            f_c = f_low + nsc / 2 * scs                           # reference: block centre
            spec[(np.arange(nsc) - nsc // 2) % nn] = grid_sym
            sym = np.fft.ifft(spec) * np.sqrt(nn)
            sym = np.concatenate([sym[-ncp:], sym])
            t0 = t_useful - ncp
            if t0 < 0 or t0 + len(sym) > nsamp:
                return
            tt = np.arange(t0, t0 + len(sym))
            sig[t0:t0 + len(sym)] += sym * np.exp(-2j * np.pi * f_c * (t_useful / FS)) * \
                np.exp(2j * np.pi * (f_c - fc) * tt / FS)

        for frame in range(-1, 3):
            f0 = frame0 + frame * int(FS * 0.01)
            if frame % 2:                                         # SSB and SIB1 every 20 ms
                continue
            for i in range(beams):
                for l in range(4):
                    lsym = lfirst[i] + l
                    t_u = f0 + halfframe_offset(lsym, n, mu)
                    put(ssb[i, l], f_ssb - 120 * scs_ssb, 240, scs_ssb, t_u,
                        cp_len(n, lsym % 14, mu))
            if kssb <= 23:
                cell = dict(mhz=f_ssb / 1e6, mib=dict(kssb=kssb, coreset0=c0, scs_common=b["scs"]))
                f_low = coreset0_layouts(cell, b)[0][0]
                for l in range(14):
                    t_u = f0 + sib1_slot * slot_len + starts[l]
                    put(slot[l], f_low, nrb * 12, scs_ssb, t_u, cp_len(n, l, mu))
        p = np.mean(np.abs(sig[sig != 0]) ** 2)
        sig *= amp / np.sqrt(p)
        sig *= np.exp(2j * np.pi * cfo * np.arange(nsamp) / FS)
        x += sig
        noise_p = max(noise_p, 10 ** (-snr_db / 10))
    out = []
    for g in gains:
        noise = np.sqrt(max(noise_p, 1e-3) / 2) * (rng.standard_normal(nsamp) + 1j * rng.standard_normal(nsamp))
        out.append(((g * x + noise) / 4).astype(np.complex64))
    return out


if __name__ == "__main__":
    sys.exit(main())
