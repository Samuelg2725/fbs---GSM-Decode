#!/usr/bin/env python3
"""LTE (4G) cell identity decoder for the BladeRF.

  Decode one or more carriers (by frequency or by EARFCN):
      python3 lte_decode.py --freq-mhz 806 816 952.5
      python3 lte_decode.py --earfcn 6300 6400 3725 223 347
  Use both antennas (RX1 + RX2 combined per cell):
      python3 lte_decode.py --freq-mhz 816 --dual

For every carrier it records the signal, finds every cell it can hear (PCI), and decodes each cell's
identity with srsRAN: Cell ID, eNB, sector, TAC and network (MCC-MNC). It is not limited to the strongest
cell: weaker cells on the same carrier are searched for at the carrier's own frequency, and the other
sectors of every decoded mast are tried using that mast's frame timing. Results are printed as a table.

Run from the repository folder (it uses the "osprey 2" folder next to it). Needs srsRAN 4G built with the
SIB1 patch (see README). Only one program can use the BladeRF at a time.
"""
import argparse
import os
import re
import select
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", str(max(1, (os.cpu_count() or 4) - 1)))
os.environ.setdefault("OMP_NUM_THREADS", os.environ["OPENBLAS_NUM_THREADS"])
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
for _cand in (os.path.join(HERE, "osprey 2"), os.path.expanduser("~/fbs-upload/osprey 2")):
    if os.path.isfile(os.path.join(_cand, "lte_cellid.py")):
        sys.path.insert(0, _cand)
        break
else:
    sys.exit("The 'osprey 2' folder (LTE decoder library) was not found next to this script.")

import common as C                                   # noqa: E402
import lte_fastsync                                  # noqa: E402
import pcap_out                                      # noqa: E402
from lte_pci_search import search_samples            # noqa: E402

BLADE_LIBS = "/usr/local/lib"
SC16 = "/dev/shm/_lte.sc16"
# capture rate needed for a cell width (resource blocks): 10 MHz -> 15.36, 15 MHz -> 23.04, 20 MHz -> 30.72 Msps
RATES = [(15.36e6, 50), (23.04e6, 75), (30.72e6, 100)]
BW_MHZ = {6: 1.4, 15: 3, 25: 5, 50: 10, 75: 15, 100: 20}
NATIVE_RATE = {6: 1.92e6, 15: 3.84e6, 25: 7.68e6, 50: 15.36e6, 75: 23.04e6, 100: 30.72e6}   # srsRAN file rates
SIB1_RE = re.compile(
    r"SIB1 cell identity (\d+)\s+eNB (\d+)\s+sector (\d+)\s+TAC (\d+)((?:[ \t]+PLMN \d+-\d+)*)[ \t]*\r?\n")
MIB_PRB_RE = re.compile(r"PRB:\s+(\d+)")
# measurements srsRAN prints every frame once it has the MIB
DIAG_RE = {"snr": re.compile(r"SNR:\s*([+-]?\d+(?:\.\d+)?) dB"),
           "cfo": re.compile(r"CFO:\s*([+-]?\d+(?:\.\d+)?) Hz"),
           "pdcch_miss": re.compile(r"PDCCH-Miss:\s*(\d+(?:\.\d+)?)%"),
           "bler": re.compile(r"PDSCH-BLER:\s*(\d+(?:\.\d+)?)%")}


def diagnosis(text):
    """srsRAN's last reported SNR (dB), frequency error (Hz), control-channel miss rate and SIB1 error rate
    (%) from one run, or {} when it printed none (no MIB)."""
    out = {}
    for key, rx in DIAG_RE.items():
        vals = rx.findall(text)
        if vals:
            out[key] = float(vals[-1])
    return out


def why_not(r):
    """Plain-language reason a cell was not decoded."""
    text = _why_not(r)
    if r.get("sync_with") is not None and text:
        text += "; frame-synchronised with PCI %d, so its MIB and control channel sit on top of PCI %d's" % (
            r["sync_with"], r["sync_with"])
    return text


def _why_not(r):
    st, d = r.get("status", ""), r.get("diag") or {}
    if st in ("not found", "frame timing not found"):
        return "its sync signal was too weak to find in this recording"
    if st == "MIB not decoded":
        return "timing found but even the MIB failed - very weak, or covered by a stronger cell"
    if "MIB only" not in st:
        return ""
    parts = []
    if "snr" in d:
        parts.append("SNR %+.1f dB" % d["snr"])
    if "pdcch_miss" in d:
        parts.append("control channel missed %.0f%%" % d["pdcch_miss"])
    if "bler" in d:
        parts.append("SIB1 data errors %.0f%%" % d["bler"])
    if "cfo" in d and abs(d["cfo"]) > 1000:
        parts.append("frequency error %.0f Hz" % d["cfo"])
    if not d:
        return "MIB read, srsRAN gave no measurements"
    snr = d.get("snr")
    if snr is not None and snr < 3:
        head = "too weak for SIB1"
    elif d.get("pdcch_miss", 0) >= 90:
        head = "SIB1's control message not found"
    elif d.get("bler", 0) >= 90:
        head = "SIB1 found but its data failed"
    else:
        head = "SIB1 not decoded in this recording"
    return head + " (" + ", ".join(parts) + ")"
LAST_PEAK = [0.0]
PCAP = [None]            # pcap_out.PcapWriter for the decoded MIB / SIB messages, or None
BCH_HEX_RE = re.compile(r"BCCH-BCH hex ([0-9a-f]{6})")
DLSCH_HEX_RE = re.compile(r"BCCH-DL-SCH hex ([0-9a-f]+)")
SIB1DCI_RE = re.compile(r"SIB1dci cfi=(\d+) L=(\d+) ncce=(\d+) riv=(\d+) mode=(\d+) nprb1a=(\d+) ngap=(\d+) mcs=(\d+)")
PPM = []                 # frequency error of the radio (offset / carrier) measured on decoded cells
UK_BANDS = (20, 8, 3, 1, 7, 28, 32, 38, 40)


def carrier_rate(n_rb):
    """Best recording rate once a carrier's width is known: just wide enough for it, so a stronger neighbouring
    carrier is left outside the radio's filter and the gain can rise for weak cells. 5 MHz and narrower ->
    7.68 Msps (with the 0.6 MHz low IF the cell still fits), 10 -> 15.36, 15 -> 23.04, 20 -> 30.72 Msps."""
    return 7.68e6 if n_rb <= 25 else needed_rate(n_rb)


def needed_rate(n_rb):
    for rate, max_rb in RATES:
        if n_rb <= max_rb:
            return rate
    return RATES[-1][0]


# ------------------------------------------------------------------------------------------ channels
def earfcn_to_mhz(earfcn):
    from cellular_resolver import lte_carrier
    c = lte_carrier(earfcn)
    return None if c is None else c.downlink_mhz


def mhz_to_earfcn(mhz):
    """(earfcn, band) of a downlink centre frequency, preferring bands used in the UK."""
    from cellular_resolver import _LTE_ROWS
    rows = sorted(_LTE_ROWS, key=lambda r: UK_BANDS.index(r[0]) if r[0] in UK_BANDS else 99)
    for r in rows:
        band, _duplex, dl_low, dl_off, dl_min, dl_max = r[:6]
        e = dl_off + int(round(mhz * 10 - dl_low))
        if dl_min <= e <= dl_max:
            return e, band
    return None, None


# ------------------------------------------------------------------------------------------- capture
def _bladerf_cli(cmd):
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = BLADE_LIBS + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    return subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def _cli_raw(freq_hz, rate, secs, gain, dual):
    """bladeRF-cli capture (opens the radio each time) -> int16 array, one row per sample:
    [I1, Q1] or [I1, Q1, I2, Q2]. Includes the tuning transient. None on failure."""
    n = int(rate * (secs + 0.01))
    chans = 2 if dual else 1
    if os.path.exists(SC16):
        os.remove(SC16)
    cmd = ["bladeRF-cli", "-e", "set agc off" if dual else "set agc rx off",
           "-e", "set frequency rx %d" % int(freq_hz),
           "-e", "set samplerate rx %d" % int(rate),
           "-e", "set bandwidth rx %d" % int(min(rate * 0.9, 56e6))]
    if dual:
        # with two channels, n= counts the samples of both; gains must be set after 'rx config' to reach both
        cmd += ["-e", "rx config file=%s format=bin n=%d channel=1,2 samples=16384 buffers=64 xfers=16 "
                "timeout=20000" % (SC16, 2 * n),
                "-e", "set gain rx1 %d" % int(gain), "-e", "set gain rx2 %d" % int(gain)]
    else:
        cmd += ["-e", "set gain rx %d" % int(gain),
                "-e", "rx config file=%s format=bin n=%d samples=16384 buffers=64 xfers=16 timeout=20000"
                % (SC16, n)]
    p = _bladerf_cli(cmd + ["-e", "rx start", "-e", "rx wait %d" % int(secs * 3000 * chans + 15000)])
    size = os.path.getsize(SC16) if os.path.exists(SC16) else 0
    if size < n * 4 * chans * 0.9:
        C.warn("  capture failed (%d of %d bytes); bladeRF-cli said: %s" % (
            size, n * 4 * chans, p.stdout.decode(errors="replace")[-300:]))
        return None
    raw = np.fromfile(SC16, np.int16)
    os.remove(SC16)
    return raw[:raw.size // (2 * chans) * 2 * chans].reshape(-1, 2 * chans)


class LibRadio(object):
    """Keeps the BladeRF open between captures, through libbladeRF itself (the same library bladeRF-cli uses),
    saving the time bladeRF-cli needs to open the radio on every capture. Used by lte_track.py --radio lib."""

    RX_X1, RX_X2, SC16_Q11, GAIN_MGC = 0, 2, 0, 1

    def __init__(self, dual):
        import ctypes
        import ctypes.util
        self.ct = ctypes
        path = os.path.join(BLADE_LIBS, "libbladeRF.so")
        if not os.path.exists(path):
            path = ctypes.util.find_library("bladeRF") or "libbladeRF.so.2"
        lib = self.lib = ctypes.CDLL(path)
        lib.bladerf_open.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p]
        lib.bladerf_close.argtypes = [ctypes.c_void_p]
        lib.bladerf_close.restype = None
        lib.bladerf_strerror.restype = ctypes.c_char_p
        lib.bladerf_set_frequency.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_uint64]
        lib.bladerf_set_sample_rate.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_uint32,
                                                ctypes.POINTER(ctypes.c_uint32)]
        lib.bladerf_set_bandwidth.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_uint32,
                                              ctypes.POINTER(ctypes.c_uint32)]
        lib.bladerf_set_gain_mode.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
        lib.bladerf_set_gain.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
        lib.bladerf_sync_config.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_uint,
                                            ctypes.c_uint, ctypes.c_uint, ctypes.c_uint]
        lib.bladerf_enable_module.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_bool]
        lib.bladerf_sync_rx.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_uint, ctypes.c_void_p,
                                        ctypes.c_uint]
        self.dev = ctypes.c_void_p()
        self._ok(lib.bladerf_open(ctypes.byref(self.dev), None), "open")
        self.chans = [0, 2] if dual else [0]                  # BLADERF_CHANNEL_RX(0), BLADERF_CHANNEL_RX(1)
        self.rate = self.freq = None

    def _ok(self, ret, what):
        if ret < 0:
            raise RuntimeError("libbladeRF %s failed: %s" % (what, self.lib.bladerf_strerror(ret).decode()))

    def raw(self, freq_hz, rate, secs, gain):
        ct, lib, dev = self.ct, self.lib, self.dev
        try:
            if rate != self.rate:
                for ch in self.chans:
                    self._ok(lib.bladerf_set_sample_rate(dev, ch, int(rate), None), "set sample rate")
                    self._ok(lib.bladerf_set_bandwidth(dev, ch, int(min(rate * 0.9, 56e6)), None), "set bandwidth")
                self.rate = rate
            if freq_hz != self.freq:
                for ch in self.chans:
                    self._ok(lib.bladerf_set_frequency(dev, ch, int(freq_hz)), "set frequency")
                self.freq = freq_hz
            # turning the receiver off at the end of a capture also releases the stream, so it is set up again
            # for every capture (cheap compared with re-opening the radio)
            self._ok(lib.bladerf_sync_config(dev, self.RX_X2 if len(self.chans) == 2 else self.RX_X1,
                                             self.SC16_Q11, 64, 16384, 16, 5000), "configure stream")
            for ch in self.chans:
                self._ok(lib.bladerf_set_gain_mode(dev, ch, self.GAIN_MGC), "set manual gain")
                self._ok(lib.bladerf_set_gain(dev, ch, int(gain)), "set gain")
            n = int(rate * (secs + 0.01))
            nch = len(self.chans)
            buf = np.empty((n, 2 * nch), np.int16)             # with two channels: I1 Q1 I2 Q2 per sample
            for ch in self.chans:
                self._ok(lib.bladerf_enable_module(dev, ch, True), "start receiving")
            try:
                self._ok(lib.bladerf_sync_rx(dev, buf.ctypes.data_as(ct.c_void_p), n * nch, None,
                                             int(secs * 1000 + 5000)), "receive")
            finally:
                for ch in self.chans:
                    lib.bladerf_enable_module(dev, ch, False)
            return buf
        except RuntimeError as e:
            C.warn("  capture failed: %s" % e)
            return None

    def close(self):
        if self.dev:
            self.lib.bladerf_close(self.dev)
            self.dev = None


RADIO = [None]           # a LibRadio when the radio is kept open, else bladeRF-cli is used
LAST_RAW = [None]        # int16 samples of the last capture (tuning transient removed), for saving
KEEP_RAW = [False]       # lte_track.py sets this; lte_decode.py does not need the raw copy


def capture_raw(freq_hz, rate, secs, gain, dual):
    """One capture -> int16 array [I1, Q1(, I2, Q2)] per row, tuning transient removed, or None."""
    if RADIO[0] is not None:
        raw = RADIO[0].raw(freq_hz, rate, secs, gain)
    else:
        raw = _cli_raw(freq_hz, rate, secs, gain, dual)
    if raw is None:
        return None
    raw = raw[int(rate * 0.01):]
    LAST_RAW[0] = raw if KEEP_RAW[0] else None
    return raw


def _to_iq(raw, col):
    return ((raw[:, col].astype(np.float32) + 1j * raw[:, col + 1].astype(np.float32)) / 2048.0).astype(
        np.complex64)


def capture(freq_hz, rate, secs, gain, low_if_hz):
    """One-antenna capture -> conditioned complex baseband (or None)."""
    raw = capture_raw(freq_hz + low_if_hz, rate, secs, gain, False)
    if raw is None:
        return None
    iq = _to_iq(raw, 0)
    LAST_PEAK[0] = float(np.max(np.abs(raw[:, :2]))) / 2048.0 if raw.size else 0.0
    _level_line("RX1", iq, LAST_PEAK[0], gain)
    return condition(iq, rate, low_if_hz)


def capture_dual(freq_hz, rate, secs, gain, low_if_hz):
    """RX1 + RX2 together (one tuner, same gain) -> (rx1, rx2) or None."""
    raw = capture_raw(freq_hz + low_if_hz, rate, secs, gain, True)
    if raw is None:
        return None
    out, peak = [], 0.0
    for col, name in ((0, "RX1"), (2, "RX2")):
        iq = _to_iq(raw, col)
        pk = float(np.max(np.abs(raw[:, col:col + 2]))) / 2048.0
        peak = max(peak, pk)
        _level_line(name, iq, pk, gain)
        out.append(condition(iq, rate, low_if_hz))
    LAST_PEAK[0] = peak
    return out[0], out[1]


LAST_LEVEL = [0.0]
SHOW_LEVELS = [True]


def _level_line(name, iq, peak, gain):
    lvl = 10 * np.log10(float(np.mean(np.abs(iq) ** 2)) + 1e-20)
    LAST_LEVEL[0] = lvl
    if not SHOW_LEVELS[0]:
        return
    tag = C.c("  CLIPPING", "red") if peak >= 0.99 else ""
    C.info("  %s level %.1f dBFS, peak %.2f of full scale, gain %g%s" % (name, lvl, peak, gain, tag))


def _tone_period(rate, freq_hz):
    """Samples after which exp(j*2*pi*freq*n/rate) repeats exactly (e.g. 0.6 MHz at 15.36 Msps -> 128), or None."""
    from math import gcd
    r, f = int(round(rate)), int(round(freq_hz))
    if r != rate or f != freq_hz or f == 0:
        return None
    p = r // gcd(r, abs(f))
    return p if p <= 8192 else None


def condition(samples, rate, low_if_hz):
    """Remove DC. The radio is tuned low_if_hz above the carrier to keep its DC spike off the cell's centre;
    shift the cell back to baseband and notch the spike."""
    samples = np.asarray(samples, dtype=np.complex64)
    samples = samples - np.mean(samples)
    if not low_if_hz:
        return samples.astype(np.complex64)
    p = _tone_period(rate, low_if_hz)
    if p is None:                                             # general (slower) path
        idx = np.arange(len(samples))
        samples = samples * np.exp(1j * 2 * np.pi * float(low_if_hz) * idx / float(rate))
        from hackrf_capture import notch_tone
        return notch_tone(samples, float(low_if_hz), rate)
    # fast path: the shifting tone repeats every p samples, so one period is tiled over the recording
    n = len(samples)
    full = n // p * p
    tone = np.exp(1j * 2 * np.pi * float(low_if_hz) * np.arange(p) / float(rate)).astype(np.complex64)
    out = np.empty(n, np.complex64)
    out[:full] = (samples[:full].reshape(-1, p) * tone).ravel()
    out[full:] = samples[full:] * tone[:n - full]
    # notch: after the shift the DC spike is a tone at +low_if_hz, i.e. the same tone; remove its projection
    blocks = out[:full].reshape(-1, p)
    coef = (np.vdot(tone, blocks.sum(axis=0)) + np.vdot(tone[:n - full], out[full:])) / n
    blocks -= (coef * tone).astype(np.complex64)
    out[full:] -= (coef * tone[:n - full]).astype(np.complex64)
    return out


# -------------------------------------------------------------------------------------------- search
def search_cells(samples, rate, args):
    """PSS/SSS cell search. With weak search: a second pass for weaker cells on the same carrier. The stock
    search picks one frequency offset per PSS group (the highest correlation anywhere); a strong cell's false
    peak at a wrong offset can win, so a weaker cell is never tested at the right frequency. All cells on a
    carrier share its frequency, so the second pass searches all 3 PSS groups at the strongest cell's offset."""
    import lte_sync as S
    from lte_pci_search import resample_to_lte
    import scipy.fft as sf
    x_in = samples[:int(rate * args.search_secs)] if args.search_secs else samples
    report = search_samples(x_in, rate)
    hits = list(report["hits"])
    if args.no_weak_search:
        return report
    pmin = S.PSS_MIN_RATIO if args.pss_min is None else float(args.pss_min)
    off0 = float(hits[0]["offset_hz"]) if hits else float(report.get("best_offset_hz", 0.0))
    x = np.asarray(resample_to_lte(np.asarray(x_in), rate), dtype=np.complex128)
    length = int(x.size)
    scale = float(np.sqrt(np.mean(np.abs(x) ** 2) + 1e-12))
    spectrum = sf.fft(x, sf.next_fast_len(length + S.FFT_SIZE - 1), workers=C.WORKERS)
    have = {h["pci"] for h in hits}
    for nid2 in range(3):
        base = S.pss_time(nid2)
        base = base / np.linalg.norm(base)
        best = None
        for d in range(-2000, 2001, 500):
            mag = S._correlate_mag(spectrum, S._shifted(base, off0 + d), length)
            _at, peak = S._strongest(mag)
            if best is None or peak > best[0]:
                best = (peak, off0 + d, mag)
        _peak, offset, mag = best
        work = np.array(mag, copy=True)
        for _attempt in range(4):
            usable = S._usable_peak(work)
            if usable is None:
                break
            at, peak = usable
            ratio = peak / scale
            if ratio < pmin:
                break
            nid1, score, pss_seq = S._confirm_sss(x, work, at, offset, nid2)
            S._blank_comb(work, at)
            if args.debug_search:
                ok = score >= S.SSS_MIN_SCORE and pss_seq >= S.SSS_MIN_SCORE
                C.info("    search: PSS group %d peak %.1fx -> PCI %d, SSS %.2f (%s)" % (
                    nid2, ratio, 3 * nid1 + nid2, score, "accepted" if ok else "rejected"))
            if score < S.SSS_MIN_SCORE or pss_seq < S.SSS_MIN_SCORE:
                continue
            pci = 3 * nid1 + nid2
            if pci in have:
                continue
            have.add(pci)
            hits.append({"pci": pci, "pss_ratio": ratio, "sss_score": score, "offset_hz": offset, "weak": True})
    hits.sort(key=lambda h: h["pss_ratio"], reverse=True)
    report = dict(report)
    report["hits"] = hits
    return report


# ---------------------------------------------------------------------------------------- alignment
def derotate(samples, cfo_hz, rate):
    """Remove a constant frequency offset, in pieces and in single precision (the whole-recording version
    needed several GB for a 2 s, two-antenna, 20 MHz recording)."""
    samples = np.asarray(samples)
    if abs(cfo_hz) < 1:
        return samples.astype(np.complex64, copy=False)
    out = np.empty(len(samples), np.complex64)
    step = 1 << 20
    for i in range(0, len(samples), step):
        idx = np.arange(i, min(len(samples), i + step), dtype=np.float64)
        out[i:i + len(idx)] = samples[i:i + len(idx)] * np.exp(-2j * np.pi * float(cfo_hz) * idx / float(rate))
    return out


def align_wide(samples, rate, pci, cfo_hz):
    """Start of subframe 0, found with the PSS/SSS at the capture's full rate (FFT size from the rate).
    The first 40 ms are searched first (fast); the whole capture only if that fails."""
    from lte_cellid import _locate_wide
    corrected = derotate(samples, cfo_hz, rate)
    nfft = int(round(rate / 15000.0))
    start = _locate_wide(corrected[:int(rate * 0.04)], pci, nfft)
    if start is None and len(corrected) > int(rate * 0.04):
        start = _locate_wide(corrected, pci, nfft)
    if start is None:
        return None
    return np.asarray(corrected[start:], dtype=np.complex64)


def align_narrow(samples, rate, pci, cfo_hz, both=False):
    """Subframe-0 alignment for a weak cell. Full-band alignment can be drowned by a stronger cell's whole
    transmission; here the capture is reduced to the centre 1.08 MHz (where PSS/SSS live), this PCI's PSS is
    found on its 5 ms comb, and the SSS is averaged over every repeat to tell subframe 0 from subframe 5."""
    import lte_sync as S
    from lte_pci_search import resample_to_lte
    from lte_sync import sss_freq
    wide = derotate(samples, cfo_hz, rate)
    x = np.asarray(resample_to_lte(wide, rate), dtype=np.complex128)
    nid2, nid1 = pci % 3, pci // 3
    tpl = S.pss_time(nid2)
    tpl = tpl / np.linalg.norm(tpl)
    nf = 1 << (len(x) + len(tpl) - 2).bit_length()
    mag = np.abs(np.fft.ifft(np.fft.fft(x, nf) * np.conj(np.fft.fft(tpl, nf)))[:len(x)])
    period = S._PERIOD
    folded = np.zeros(period)
    for k in range(len(mag) // period):
        folded += mag[k * period:(k + 1) * period]
    phase = int(np.argmax(folded))
    times = [phase + k * period for k in range(len(mag) // period + 1)
             if S.SYMBOL <= phase + k * period and phase + k * period + S.FFT_SIZE <= len(x)]
    if len(times) < 2:
        return None
    groups = ([], [])
    for i, t in enumerate(times):
        seq = S._extract_seq(x, t - S.SYMBOL, 0.0)
        if seq is not None:
            groups[i % 2].append(seq)
    scored = []
    for swap in (False, True):
        sc = [S._mean_match(groups[g], sss_freq(nid1, nid2, (g == 1) != swap)) for g in (0, 1) if groups[g]]
        scored.append((sum(sc) / max(1, len(sc)), swap))
    scored.sort(reverse=True)
    if both and scored[0][0] < 0.3:
        # too weak to tell subframe 0 from subframe 5 with the SSS: return both positions, best first - the
        # MIB (only in subframe 0, checked by its CRC) decides
        outs = [_narrow_start(wide, rate, times, sw, nid2) for _sc, sw in scored]
        return [o for o in outs if o is not None] or None
    score, swap = scored[0]
    if score < 0.3:
        return None
    out = _narrow_start(wide, rate, times, swap, nid2)
    return [out] if (both and out is not None) else out


def _narrow_start(wide, rate, times, swap, nid2):
    """Subframe-0 aligned recording from the PSS times found by align_narrow (swap: subframe 5 came first)."""
    from lte_cellid import _useful_wide, _subframe_wide, _refine_wide, symbol_useful_offset
    from lte_sync import _pss_freq
    start_n = (times[1] if swap else times[0]) - symbol_useful_offset(6)
    nfft = int(round(rate / 15000.0))
    pss_at, sub = _useful_wide(6, nfft), _subframe_wide(nfft)
    frame = 10 * sub
    guess = int(round(start_n * rate / 1.92e6)) + pss_at
    try:
        peak = _refine_wide(wide, guess, _pss_freq(nid2), nfft)
    except Exception:
        peak = guess
    if abs(peak - guess) > nfft // 4:
        peak = guess
    start = peak - pss_at
    while start >= frame:
        start -= frame
    while start < 0:
        start += frame
    if start + sub > len(wide):
        return None
    return np.asarray(wide[start:], dtype=np.complex64)


# -------------------------------------------------------------------------------------------- srsRAN
def pdsch_ue():
    import hackrf_lte_decode as H
    return H._pdsch()


def _srsran_child():
    """In the srsRAN process before it starts: its own process group (so it can be stopped as a whole) and
    stopped automatically if this tool ends for any reason, even when killed for lack of memory."""
    os.setsid()
    try:
        import ctypes
        ctypes.CDLL("libc.so.6", use_errno=True).prctl(1, 9)     # PR_SET_PDEATHSIG = SIGKILL
    except Exception:
        pass


def clean_leftovers():
    """Stop srsRAN / LTE / GSM tools left running by an earlier run and delete their temporary files."""
    import glob
    n = C.stop_leftovers(("lte_decode", "lte_track", "gsm_arfcn_fast", "gsm_dual", "gsm_neigh"), ("pdsch_ue",))
    freed = 0
    for f in glob.glob("/dev/shm/lte_dec_*") + glob.glob("/dev/shm/tmp*.cf32") + [SC16]:
        try:
            freed += os.path.getsize(f)
            os.remove(f)
        except OSError:
            pass
    if n or freed > 50e6:
        C.warn("Cleaned up after an earlier run: stopped %d leftover process(es), freed %.0f MB of memory."
               % (n, freed / 1e6))


def run_pdsch(path, pci, n_prb, n_sf, timeout, ports=2, extra=()):
    """srsRAN pdsch_ue on an aligned capture file; stopped as soon as it prints the SIB1 identity line.
    ports: the cell's number of transmit antennas (1, 2 or 4) - in file mode srsRAN only tries this one."""
    cmd = [str(pdsch_ue()), "-d", "-i", str(path), "-p", str(n_prb), "-c", str(pci), "-P", str(ports), "-Q",
           "-n", str(n_sf)] + list(extra)
    if shutil.which("stdbuf"):
        cmd = ["stdbuf", "-oL", "-eL"] + cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, preexec_fn=_srsran_child)
    buf = b""
    deadline = time.time() + timeout
    try:
        while time.time() < deadline:
            ready, _, _ = select.select([p.stdout], [], [], 0.25)
            if not ready:
                if p.poll() is not None:
                    break
                continue
            chunk = os.read(p.stdout.fileno(), 65536)
            if not chunk:
                break
            buf += chunk
            if SIB1_RE.search(buf.decode(errors="replace")):
                break
    finally:
        if p.poll() is None:
            try:
                os.killpg(p.pid, 9)
            except OSError:
                p.kill()
        try:
            buf += p.stdout.read() or b""
        except Exception:
            pass
        p.stdout.close()
        try:
            p.wait(timeout=2)
        except Exception:
            pass
    return buf.decode(errors="replace")


def decode_cell(samples, rate, hit, args, freq_mhz, aligned=None):
    """MIB + SIB1 of one PCI with srsRAN. Returns a result dict; 'ident' is set when the identity decoded."""
    pci = int(hit["pci"])
    out = {"pci": pci, "strength": hit.get("pss_ratio", 0.0), "offset": float(hit.get("offset_hz") or 0.0),
           "mib": None, "mib_odd": False, "ident": None,
           "aligned": None}
    if aligned is None:
        aligned = align_wide(samples, rate, pci, float(hit.get("offset_hz", 0.0)))
    if aligned is None:
        cands = align_narrow(samples, rate, pci, float(hit.get("offset_hz", 0.0)), both=True)
        if cands and len(cands) > 1:
            # subframe 0 / 5 not distinguishable: try both positions, the MIB's CRC picks the right one
            r = None
            for cand in cands:
                r = decode_cell(samples, rate, hit, args, freq_mhz, aligned=cand)
                if r["ident"] or r["mib"] or r["mib_odd"]:
                    return r
            return r
        aligned = cands[0] if cands else None
    if aligned is None:
        out["status"] = "frame timing not found"
        return out
    out["aligned"] = aligned
    # frame start in the recording (to see whether cells are frame-synchronised)
    out["start"] = (len(samples) - len(aligned)) if samples is not None else hit.get("start")
    with tempfile.NamedTemporaryFile(prefix="lte_dec_", suffix=".cf32", delete=False, dir="/dev/shm") as fh:
        path = Path(fh.name)
    try:
        aligned.tofile(path)
        n_sf = max(40, len(aligned) // max(1, int(round(rate / 1000.0))))
        own = {7680000: 25, 15360000: 50, 23040000: 75, 30720000: 100}.get(int(round(rate)), 50)
        order = [own] + [w for w in (50, 25, 75, 100) if w != own]
        if hit.get("width") in order:                        # width already known from an earlier MIB
            order.remove(hit["width"])
            order.insert(0, hit["width"])
        ports = hit.get("ports") or 2
        pdus = {}                                # decoded MIB / SIB bytes, in order, for the pcap file
        diag = {}                                # srsRAN's measurements of the cell

        # channel estimate: averaged over the subframe by default (less noisy: weak cells decode that fail with
        # srsRAN's own default); srsRAN's interpolated estimate (follows fast changes, e.g. when moving) is the
        # second try
        def srsran(n_prb, ports, extra=("-R", "average")):
            C.status("  PCI %d: srsRAN decoding (%s MHz, %d transmit antenna%s, %s channel estimate)..." % (
                pci, BW_MHZ.get(n_prb, "?"), ports, "" if ports == 1 else "s", extra[-1]))
            # srsRAN reads a file at the standard rate for the width it is told (5 MHz = 7.68 Msps, ...): give
            # a narrower cell a recording reduced to its own rate
            k = rate / NATIVE_RATE[n_prb]
            if k > 1.01 and abs(k - round(k)) < 1e-6:
                from scipy.signal import resample_poly
                small = path.with_name(path.stem + "_%d.cf32" % n_prb)
                resample_poly(aligned, 1, int(round(k))).astype(np.complex64).tofile(small)
                try:
                    text = run_pdsch(small, pci, n_prb, n_sf, args.srsran_timeout, ports, extra)
                finally:
                    small.unlink(missing_ok=True)
            else:
                text = run_pdsch(path, pci, n_prb, n_sf, args.srsran_timeout, ports, extra)
            C.clear_status()
            diag.update(diagnosis(text))
            for h in BCH_HEX_RE.findall(text):
                pdus.setdefault(("bch", h), None)
            for h in DLSCH_HEX_RE.findall(text):
                pdus.setdefault(("dlsch", h), None)
            m = SIB1DCI_RE.search(text)
            if m and "sib1_dci" not in out:
                out["sib1_dci"] = dict(zip(("cfi", "L", "ncce", "riv", "mode", "nprb1a", "ngap", "mcs"),
                                           map(int, m.groups())))
            if args.srsran_log:
                d = os.path.expanduser(args.srsran_log)
                os.makedirs(d, exist_ok=True)
                with open(os.path.join(d, "pdsch_%s_%gMHz_%.0fMsps_pci%d_prb%d_ports%d.txt" % (
                        time.strftime("%H%M%S"), freq_mhz, rate / 1e6, pci, n_prb, ports)), "w") as f:
                    f.write(text)
            return text

        def got_mib(text):
            # srsRAN is stopped as soon as it prints the Cell ID, often before its "Decoded MIB" line comes out,
            # so a Cell ID counts as a decoded MIB too
            return "Decoded MIB" in text or bool(SIB1_RE.search(text))

        for i, n_prb in enumerate(order):
            text = srsran(n_prb, ports)
            # only for a clearly strong cell: a weak cell fails the MIB for lack of signal, and each extra
            # srsRAN run on a wide two-antenna recording costs up to 30 s
            if i == 0 and not got_mib(text) and not hit.get("ports") and (hit.get("pss_ratio") or 0) >= 12:
                # no MIB: the cell may send from 4 (or 1) antennas - srsRAN must be told the right number
                for p in (4, 1, 2):
                    if p != ports:
                        text = srsran(n_prb, p)
                        if got_mib(text):
                            ports = p
                            break
            if i == 0 and not got_mib(text):
                break                                # the MIB does not depend on the width: stop here
            if "Decoded MIB" in text and out["mib"] is None:
                vals = [int(v) for v in MIB_PRB_RE.findall(text) if int(v) in BW_MHZ]
                if vals:
                    out["mib"] = vals[-1]
                else:
                    out["mib_odd"] = True            # MIB decoded but its width was not understood
            if parse_ident(text, out):
                if out["mib"] is None:
                    out["mib"] = own
                break
            if out["mib"] is not None and (out["mib"] == n_prb or needed_rate(out["mib"]) > rate):
                break                                # right width tried, or the cell needs a wider capture
        if not out["ident"] and out["mib"] and NATIVE_RATE[out["mib"]] <= rate * 1.01:
            # MIB read but not SIB1: one more go with srsRAN's interpolated channel estimate, which can do better
            # when the radio or the surroundings are moving
            text = srsran(out["mib"], ports, ("-R", "interpolate"))
            if parse_ident(text, out):
                out["note"] = "interpolated channel estimate"
    finally:
        path.unlink(missing_ok=True)
    out["ports"] = ports if (out["mib"] or out["mib_odd"]) else None
    out["diag"] = diag
    out["mib_hex"] = next((h for kind, h in pdus if kind == "bch"), None)
    out["sib1_hex"] = next((h for kind, h in pdus if kind == "dlsch"), None) if out["ident"] else None
    if PCAP[0] is not None:
        for kind, h in pdus:
            PCAP[0].lte_rrc(pcap_out.LTE_RRC_BCCH_BCH if kind == "bch" else pcap_out.LTE_RRC_BCCH_DL_SCH,
                            bytes.fromhex(h), earfcn=getattr(args, "_earfcn_now", 0), pci=pci)
    if out["ident"]:
        out["status"] = "decoded"
    elif out["mib"]:
        out["status"] = "MIB only (%s MHz), no SIB1" % BW_MHZ.get(out["mib"], "?")
        if diag.get("snr") is not None:
            out["status"] += ", SNR %+.1f dB" % diag["snr"]
    elif out["mib_odd"]:
        out["status"] = "MIB only, width unclear"
    else:
        out["status"] = "MIB not decoded"
    return out


def parse_ident(text, out):
    """Fill out["ident"] from srsRAN's SIB1 line, if there is one. Returns True when found."""
    s1 = SIB1_RE.search(text)
    if not s1:
        return False
    plmns = re.findall(r"PLMN (\d+)-(\d+)", s1.group(5) or "")
    mcc, mnc = (int(plmns[0][0]), plmns[0][1]) if plmns else (None, "?")
    out["ident"] = {"cell_id": int(s1.group(1)), "enb": int(s1.group(2)), "sector": int(s1.group(3)),
                    "tac": int(s1.group(4)), "mcc": mcc, "mnc": mnc, "operator": C.operator(mcc, mnc) if mcc else ""}
    return True


def cell_weights(x1, x2, rate, pci, cfo_hz):
    """Two-antenna maximum-ratio combining weights for one cell, from that cell's PSS on its 5 ms comb:
    each antenna's complex PSS response is averaged and its conjugate used as the weight.
    Returns (w1, w2, share_rx1)."""
    import lte_sync as S
    from lte_pci_search import resample_to_lte
    tpl = S.pss_time(pci % 3)
    tpl = tpl / np.linalg.norm(tpl)
    corr, folds = [], []
    for x in (x1, x2):
        y = np.asarray(resample_to_lte(np.asarray(x[:int(rate * 0.2)]), rate), dtype=np.complex128)
        if cfo_hz:
            y = y * np.exp(-2j * np.pi * cfo_hz * np.arange(len(y)) / S.FS)
        nf = 1 << (len(y) + len(tpl) - 2).bit_length()
        cr = np.fft.ifft(np.fft.fft(y, nf) * np.conj(np.fft.fft(tpl, nf)))[:len(y)]
        corr.append(cr)
        per = len(cr) // S._PERIOD
        folds.append((np.abs(cr[:per * S._PERIOD]) ** 2).reshape(per, S._PERIOD).sum(axis=0))
    phase = int(np.argmax(folds[0] + folds[1]))
    h = [complex(np.mean(cr[np.arange(phase, len(cr), S._PERIOD)])) for cr in corr]
    norm = np.sqrt(abs(h[0]) ** 2 + abs(h[1]) ** 2) + 1e-20
    return np.conj(h[0]) / norm, np.conj(h[1]) / norm, abs(h[0]) ** 2 / (norm ** 2)


# ------------------------------------------------------------------------------------------ carrier
def record(args, freq_mhz, rate):
    """Capture at the current gain (one or two antennas) -> (rx1, rx2 or None) or (None, None)."""
    lif = args.low_if_mhz * 1e6
    if args.dual:
        pair = capture_dual(freq_mhz * 1e6, rate, args.secs, args.gain_now, lif)
        return pair if pair else (None, None)
    return capture(freq_mhz * 1e6, rate, args.secs, args.gain_now, lif), None


def auto_gain(args, freq_mhz, rate, x1, x2, rounds=3):
    """Re-record at a better gain (aiming for a peak of about half full scale) up to `rounds` times."""
    for _round in range(0 if args.no_auto_gain else rounds):
        if x1 is None:
            break
        pk, new = LAST_PEAK[0], None
        if pk >= 0.99:
            new = args.gain_now - 10
        elif 0 < pk < 0.25 and args.gain_now < 60:
            new = min(60.0, args.gain_now + max(3.0, min(20.0, round(20 * np.log10(0.5 / pk)))))
        if new is None or new == args.gain_now:
            break
        C.info("  auto gain: %g -> %g dB" % (args.gain_now, new))
        args.gain_now = new
        x1, x2 = record(args, freq_mhz, rate)
    return x1, x2


MEMORY = os.path.expanduser("~/field/lte_carriers.json")


def carrier_memory():
    """Width and gain of each carrier from earlier runs (EARFCN -> {"width": resource blocks, "gain": dB})."""
    try:
        import json
        with open(MEMORY) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def remember_carrier(earfcn, width, gain):
    import json
    mem = carrier_memory()
    mem[str(earfcn)] = {"width": int(width), "gain": float(gain)}
    try:
        os.makedirs(os.path.dirname(MEMORY), exist_ok=True)
        with open(MEMORY, "w") as fh:
            json.dump(mem, fh, indent=1, sort_keys=True)
    except OSError:
        pass


def decode_carrier(args, freq_mhz):
    """Everything for one carrier. Returns a list of result rows (decoded and not decoded)."""
    earfcn, band = mhz_to_earfcn(freq_mhz)
    args._targets_now = list(getattr(args, "targets", {}).get(earfcn, []))
    args._earfcn_now = earfcn or 0
    C.heading("%.1f MHz   EARFCN %s   band %s" % (freq_mhz, earfcn if earfcn is not None else "?", band or "?"))
    t0 = time.time()
    rate = 15.36e6
    args.gain_now = args.gain
    known = carrier_memory().get(str(earfcn), {}) if not args.no_auto_gain else {}
    if known:
        # width and gain learned on an earlier run: record at the right rate and gain straight away
        if known.get("width"):
            rate = carrier_rate(known["width"])
        if known.get("gain") is not None:
            args.gain_now = float(known["gain"])
        C.info("  from an earlier run: %s MHz wide, gain %g - recording at %.2f Msps" % (
            BW_MHZ.get(known.get("width"), "?"), args.gain_now, rate / 1e6))
    x1, x2 = auto_gain(args, freq_mhz, rate, *record(args, freq_mhz, rate))
    rows = []
    if x1 is None:
        return [{"mhz": freq_mhz, "earfcn": earfcn, "band": band, "pci": None, "status": "capture failed"}]
    rates_done = {rate}
    best = {}                                    # pci -> best result over all recordings of this carrier
    args._done_pcis = set()
    retries = args.retries
    while True:
        for r in _decode_capture(args, freq_mhz, rate, x1, x2):
            old = best.get(r["pci"])
            if old is None or (r["ident"] and not old["ident"]) or (not old["ident"] and r["mib"] and not old["mib"]):
                best[r["pci"]] = r
        args._done_pcis = {p for p, r in best.items() if r["ident"]}
        results = list(best.values())
        need = [carrier_rate(r["mib"]) for r in results if not r["ident"] and r["mib"]]
        if not need:                             # MIB read but its width not understood: try the widest
            need = [30.72e6 for r in results if not r["ident"] and r["mib_odd"]]
        need = [n for n in need if n not in rates_done]
        if need:
            rate = max(need)
            rates_done.add(rate)
            x1 = x2 = None                       # free the old recording first
            import gc
            gc.collect()
            if rate > 15.36e6:
                C.info("  a cell is wider than this capture - recording again at %.2f Msps" % (rate / 1e6))
            else:
                C.info("  the carrier is %s MHz wide - recording just that (%.2f Msps) so the gain can rise for it"
                       % (BW_MHZ.get(min(r["mib"] for r in results if r["mib"]), "?"), rate / 1e6))
            x1, x2 = auto_gain(args, freq_mhz, rate, *record(args, freq_mhz, rate))
            if x1 is None:
                break
            continue
        # a weak cell whose MIB decoded is certainly there: SIB1 often decodes on another recording
        weak = sorted(r["pci"] for r in results if not r["ident"] and (r["mib"] or r["mib_odd"]))
        if weak and retries > 0:
            retries -= 1
            C.info("  PCI %s: MIB read but not the Cell ID - recording again for another try (%d more after this)"
                   % (", ".join(map(str, weak)), retries))
            args._targets_now = sorted(set(args._targets_now) | set(weak))
            x1 = x2 = None
            x1, x2 = record(args, freq_mhz, rate)
            if x1 is None:
                break
            continue
        break
    widths = [r["mib"] for r in results if r.get("mib")]
    if widths and earfcn is not None:
        remember_carrier(earfcn, max(widths), args.gain_now)
    for r in results:
        r.update({"mhz": freq_mhz, "earfcn": earfcn, "band": band})
        rows.append(r)
    if not results:
        rows.append({"mhz": freq_mhz, "earfcn": earfcn, "band": band, "pci": None, "status": "no LTE cell found"})
    C.info("  carrier done in %.1f s" % (time.time() - t0))
    return rows


def _decode_capture(args, freq_mhz, rate, x1, x2):
    t0 = time.time()
    C.status("  searching for cells...")
    report = search_cells(x1, rate, args)
    hits = list(report["hits"])
    if x2 is not None:                       # a cell can be heard on one antenna only
        have = {h["pci"] for h in hits}
        for h in search_cells(x2, rate, args)["hits"]:
            if h["pci"] not in have:
                hits.append(dict(h, rx2_only=True))
    C.info("  cell search %.1f s: %s" % (time.time() - t0, ", ".join(
        "PCI %d (strength %.0f%s)" % (h["pci"], h["pss_ratio"], ", weak search" if h.get("weak") else "")
        for h in hits) or "no LTE cell heard (best sync %.1fx, needs %.1fx)" % (
        report["best_ratio"], report["min_ratio"])))
    results, tried = [], set()

    def show(r):
        if r["ident"]:
            x = r["ident"]
            line = "  PCI %-3d %s  Cell ID %d (eNB %d, sector %d)  TAC %d  %s %s-%s" % (
                r["pci"], C.c("decoded", "green"), x["cell_id"], x["enb"], x["sector"], x["tac"],
                x["operator"] or "?", x["mcc"], x["mnc"])
        else:
            line = "  PCI %-3d %s" % (r["pci"], C.c(r["status"], "yellow" if "MIB only" in r["status"] else "red"))
        note = ", " + r["note"] if r.get("note") else ""
        print(line + C.c("   [%s%s%s]" % (r["how"], note, ", " + r["rx"] if r.get("rx") else ""), "grey"), flush=True)
        if not r["ident"] and why_not(r):
            print(C.c("          why: " + why_not(r), "grey"), flush=True)

    tried |= set(getattr(args, "_done_pcis", ()))     # decoded from an earlier recording of this carrier
    for h in hits:
        if int(h["pci"]) in tried:
            continue
        if getattr(args, "only", False) and int(h["pci"]) not in getattr(args, "_targets_now", []):
            continue                             # --only: decode just the PCIs given with --cells
        if x2 is not None:
            w1, w2, share = cell_weights(x1, x2, rate, int(h["pci"]), float(h.get("offset_hz", 0.0)))
            r = decode_cell((w1 * x1 + w2 * x2).astype(np.complex64), rate, h, args, freq_mhz)
            r["rx"] = "RX1 %.0f%% / RX2 %.0f%%" % (100 * share, 100 * (1 - share))
        else:
            r = decode_cell(x1, rate, h, args, freq_mhz)
        r["how"] = "weak-cell search" if h.get("weak") else ("RX2 only" if h.get("rx2_only") else "cell search")
        if r["ident"]:
            PPM.append(float(h.get("offset_hz", 0.0)) / (freq_mhz * 1e6))
        else:
            r.pop("aligned", None)               # only decoded cells' timing is used again (memory)
        tried.add(r["pci"])
        results.append(r)
        show(r)
    # PCIs given with --cells (e.g. from a phone's cell list) that the search did not find: look for each one
    # directly. Knowing the PCI, its sync signal can be found by adding up every 5 ms repeat over the whole
    # recording, which works for cells too weak for the blind search.
    for pci in getattr(args, "_targets_now", []):
        if pci in tried:
            continue
        tried.add(pci)
        parent = next((p for p in results if p["ident"] and p.get("aligned") is not None
                       and p["pci"] // 3 == pci // 3), None)
        if parent is not None:                   # same mast as a decoded cell: use its frame timing first
            r = decode_cell(None, rate, {"pci": pci, "width": parent.get("mib"), "ports": parent.get("ports"),
                             "offset_hz": parent.get("offset", 0.0), "start": parent.get("start")}, args, freq_mhz, aligned=parent["aligned"])
            if r["ident"]:
                r["how"] = "looked for (--cells), timing of PCI %d" % parent["pci"]
                results.append(r)
                show(r)
                continue
        if hits:
            off = float(hits[0].get("offset_hz", 0.0))
        else:                                    # no cell here to measure it: use the radio's error seen so far
            off = float(np.median(PPM)) * freq_mhz * 1e6 if PPM else 0.0
        h = {"pci": pci, "offset_hz": off, "pss_ratio": 0.0}
        if x2 is not None:
            w1, w2, share = cell_weights(x1, x2, rate, pci, off)
            r = decode_cell((w1 * x1 + w2 * x2).astype(np.complex64), rate, h, args, freq_mhz)
            r["rx"] = "RX1 %.0f%% / RX2 %.0f%%" % (100 * share, 100 * (1 - share))
        else:
            r = decode_cell(x1, rate, h, args, freq_mhz)
        r["how"] = "looked for (--cells)"
        if not r["ident"]:
            r.pop("aligned", None)
        if r["status"] == "frame timing not found":
            r["status"] = "not found"
        results.append(r)
        show(r)
    if not args.no_siblings:
        # the sectors of one mast share its clock and are usually numbered 3*G+0, +1, +2: try the other two
        # PCIs of every decoded cell's group with that cell's frame timing (no sync signal needed; srsRAN's
        # checksum rejects a PCI that is not there)
        for parent in [r for r in list(results) if r["ident"] and r.get("aligned") is not None]:
            base = 3 * (parent["pci"] // 3)
            for sib in (base, base + 1, base + 2):
                if sib in tried or (getattr(args, "only", False) and sib not in args._targets_now):
                    continue
                tried.add(sib)
                r = decode_cell(None, rate, {"pci": sib, "width": parent.get("mib"), "ports": parent.get("ports"),
                                                "offset_hz": parent.get("offset", 0.0), "start": parent.get("start")}, args,
                                freq_mhz, aligned=parent["aligned"])
                r["how"] = "sector of PCI %d's mast" % parent["pci"]
                if r["ident"] or args.show_failed_siblings:
                    results.append(r)
                    show(r)
    if not args.no_cancel:
        results = _after_cancelling(args, freq_mhz, rate, x1, x2, results, tried, show)
    mark_synchronised(results, rate)
    for r in results:
        r.pop("aligned", None)
    return results


def mark_synchronised(results, rate):
    """Note cells that were not decoded but whose frames start (within 10 us) where a decoded cell's do. In a
    frame-synchronised network such a cell's MIB, control channel and SIB1 slot sit exactly on top of the decoded
    cell's, so they are covered by it every time - more recordings do not help."""
    frame = int(round(rate * 0.01))
    tol = rate * 10e-6
    decoded = [r for r in results if r.get("ident") and r.get("start") is not None]
    for r in results:
        if r.get("ident") or r.get("start") is None:
            continue
        for d in decoded:
            diff = (r["start"] - d["start"]) % frame
            if min(diff, frame - diff) <= tol:
                r["sync_with"] = d["pci"]
                break
    timed = [r for r in results if r.get("start") is not None]
    if len(timed) > 1:
        C.info("  frame timing (ms into the 10 ms frame): " + ",  ".join(
            "PCI %d %.3f%s" % (r["pci"], (r["start"] % frame) / rate * 1e3,
                               " (decoded)" if r.get("ident") else
                               (" = PCI %d's" % r["sync_with"] if r.get("sync_with") is not None else ""))
            for r in sorted(timed, key=lambda r: r["start"] % frame)))


def _after_cancelling(args, freq_mhz, rate, x1, x2, results, tried, show, rounds=3):
    """Interference cancellation, in rounds. Every cell whose signals can be rebuilt - decoded, or at least its MIB
    decoded (PCI, frame timing, width and antenna ports known) - has its reference signals, sync signals and MIB
    channel removed from the recording (strongest first, each estimated from what the stronger ones left). The
    cleaned recording is searched again and every cell not yet decoded is retried, with all the OTHER cells removed.
    Cells that get further (decoded, or MIB newly decoded) are removed too in the next round (lte_cancel.py)."""
    from lte_cancel import cancel_crs, regen_sib1
    only = getattr(args, "only", False)
    regen = pdsch_ue().parent / "sib1_regen"
    grids_cache = {}

    def sib1_grids(r, n_rb):
        """The decoded cell's SIB1 subframe rebuilt by srsRAN (needs the sib1_regen helper from the patch)."""
        if not (r["ident"] and r.get("sib1_dci") and r.get("sib1_hex") and r.get("mib_hex") and regen.is_file()):
            return None
        if r.get("mib") != n_rb:
            return None                          # the cell is wider than this recording
        key = r["pci"]
        if key not in grids_cache:
            grids_cache[key] = regen_sib1(regen, r["pci"], n_rb, r.get("ports") or 2, r["mib_hex"], r["sib1_dci"],
                                          r["sib1_hex"])
        return grids_cache[key]
    own_rb = {7680000: 25, 15360000: 50, 23040000: 75, 30720000: 100}.get(int(round(rate)), 50)
    t0 = time.time()
    for rnd in range(rounds):
        known = sorted([r for r in results if r.get("start") is not None and (r["ident"] or r["mib"])],
                       key=lambda r: -(r.get("strength") or 0))
        if not known:
            return results
        C.status("  interference cancellation, round %d: removing PCI %s..." % (
            rnd + 1, ", ".join(str(r["pci"]) for r in known)))
        res1, res2 = x1, x2
        own_part = {}                            # a not-yet-decoded cell's own rebuilt part, to add back for it
        for r in known:
            n_rb = min(r.get("mib") or own_rb, own_rb)
            kw = dict(rate=rate, start=int(r["start"]), pci=r["pci"], n_rb=n_rb,
                      ports=r.get("ports") or 2, cfo_hz=r.get("offset", 0.0), mib_hex=r.get("mib_hex"),
                      sib1_grids=sib1_grids(r, n_rb))
            c1 = cancel_crs(res1, **kw)
            c2 = cancel_crs(res2, **kw) if res2 is not None else None
            if not r["ident"]:
                own_part[r["pci"]] = (res1 - c1, None if c2 is None else res2 - c2)
            res1, res2 = c1, c2
        hits = list(search_cells(res1, rate, args)["hits"])
        if res2 is not None:
            have = {h["pci"] for h in hits}
            hits += [h for h in search_cells(res2, rate, args)["hits"] if h["pci"] not in have]
        decoded = {r["pci"] for r in results if r["ident"]} | set(getattr(args, "_done_pcis", ()))
        found = {int(h["pci"]): h for h in hits}
        todo = [r["pci"] for r in results if not r["ident"]]
        todo += [p for p in found if p not in decoded and p not in todo and (not only or p in args._targets_now)]
        new_found = [p for p in todo if p not in tried]
        C.clear_status()
        if not todo:
            C.info("  interference cancellation (%.1f s): no further cells under PCI %s" % (
                time.time() - t0, ", ".join(str(r["pci"]) for r in known)))
            return results
        C.info("  interference cancellation round %d (%.1f s): removed PCI %s, retrying %s%s" % (
            rnd + 1, time.time() - t0, ", ".join(str(r["pci"]) for r in known), ", ".join(map(str, todo)),
            ("; new in the cleaned search: " + ", ".join(map(str, new_found))) if new_found else ""))
        improved = False
        # cells on one carrier normally share width and antenna-port count (same operator, often the same mast):
        # a cell not decoded yet is tried with what the known cells use
        common_ports = max(set(k.get("ports") or 2 for k in known), key=[k.get("ports") or 2 for k in known].count)
        common_width = known[0].get("mib")
        for pci in todo[:8]:
            old = next((x for x in results if x["pci"] == pci), None)
            h = dict(found.get(pci) or {"pci": pci, "offset_hz": (old or known[0]).get("offset", 0.0),
                                        "pss_ratio": 0.0})
            h["ports"] = (old or {}).get("ports") or common_ports
            h["width"] = (old or {}).get("mib") or common_width
            s1, s2 = res1, res2
            if pci in own_part:                  # keep this cell's own signal in its recording
                s1 = res1 + own_part[pci][0]
                s2 = None if res2 is None else res2 + own_part[pci][1]
            if s2 is not None:
                w1, w2, share = cell_weights(s1, s2, rate, pci, float(h.get("offset_hz", 0.0)))
                r = decode_cell((w1 * s1 + w2 * s2).astype(np.complex64), rate, h, args, freq_mhz)
            else:
                r = decode_cell(s1, rate, h, args, freq_mhz)
            r["how"] = "after removing PCI %s" % ", ".join(str(k["pci"]) for k in known if k["pci"] != pci)
            if r["status"] == "frame timing not found":
                r["status"] = "not found"
            better = old is None or (r["ident"] and not old["ident"]) or \
                (not old["ident"] and r["mib"] and not old["mib"])
            if better:
                if old is not None:
                    results.remove(old)
                results.append(r)
                tried.add(pci)
                improved = improved or bool(r["ident"] or r["mib"])
                if r["ident"] or old is None or r["mib"]:
                    show(r)
            r.pop("aligned", None)
        own_part.clear()
        if not improved:
            break
    return results


# --------------------------------------------------------------------------------------------- main
def summary(rows, t_all):
    ok = [r for r in rows if r.get("ident")]
    t = C.Table(["MHz", "EARFCN", "Band", "PCI", "Operator", "MCC-MNC", "TAC", "eNB", "Sector", "Cell ID",
                 "Width", "Found by"],
                align={k: "right" for k in ("MHz", "EARFCN", "Band", "PCI", "TAC", "eNB", "Sector", "Cell ID")})
    for r in ok:
        x = r["ident"]
        t.add(["%.1f" % r["mhz"], r["earfcn"] or "?", r["band"] or "?", r["pci"], C.op_cell(x["operator"]),
               "%s-%s" % (x["mcc"], x["mnc"]), x["tac"], x["enb"], x["sector"], (x["cell_id"], "green+bold"),
               "%s MHz" % BW_MHZ.get(r.get("mib"), "?"), r["how"] + (", " + r["rx"] if r.get("rx") else "")])
    t.show("LTE cells decoded: %d   (%.0f s)" % (len(ok), time.time() - t_all))
    bad = [r for r in rows if not r.get("ident")]
    if bad:
        f = C.Table(["MHz", "EARFCN", "PCI", "Result", "Why"],
                    align={"MHz": "right", "EARFCN": "right", "PCI": "right"})
        for r in bad:
            f.add(["%.1f" % r["mhz"], r["earfcn"] or "?", r["pci"] if r.get("pci") is not None else "-",
                   (r["status"], "yellow" if "MIB only" in r["status"] else "red"), why_not(r) or "-"])
        f.show("Not decoded")


def phone_check(args, rows, passes=None):
    """Table of the cells given with --cells: decoded, found but not decoded, or not found."""
    cols = ["EARFCN", "PCI", "Result", "Operator", "Cell ID", "How"] + (["Decoded in"] if passes else [])
    t = C.Table(cols, align={"EARFCN": "right", "PCI": "right", "Cell ID": "right", "Decoded in": "right"})
    ok = total = 0
    for earfcn in args.targets:
        for pci in args.targets[earfcn]:
            total += 1
            r = next((x for x in rows if x.get("earfcn") == earfcn and x.get("pci") == pci and x.get("ident")),
                     None) or next((x for x in rows if x.get("earfcn") == earfcn and x.get("pci") == pci), None)
            if r and r.get("ident"):
                ok += 1
                x = r["ident"]
                row = [earfcn, pci, ("decoded", "green"), C.op_cell(x["operator"]),
                       (x["cell_id"], "green+bold"), r["how"]]
            elif r:
                row = [earfcn, pci, (r["status"], "yellow" if "MIB" in r["status"] else "red"), "-", "-", r["how"]]
            else:
                row = [earfcn, pci, ("not tried", "red"), "-", "-", "-"]
            if passes:
                row.append("%d of %d passes" % ((r or {}).get("decoded_passes", 0), passes))
            t.add(row)
    t.show("Cell list check: %d of %d decoded" % (ok, total))


def write_csv(path, rows):
    path = os.path.expanduser(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new:
            fh.write("utc,mhz,earfcn,band,pci,status,operator,mcc,mnc,tac,enb,sector,cell_id,width_mhz,found_by\n")
        for r in rows:
            x = r.get("ident") or {}
            fh.write("%s,%.1f,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                C.utc_iso(), r["mhz"], r.get("earfcn") or "", r.get("band") or "",
                "" if r.get("pci") is None else r["pci"], r.get("status", "").replace(",", ";"),
                x.get("operator", ""), x.get("mcc", ""), x.get("mnc", ""), x.get("tac", ""), x.get("enb", ""),
                x.get("sector", ""), x.get("cell_id", ""), BW_MHZ.get(r.get("mib"), "") if x else "",
                r.get("how", "")))
    C.info("results added to %s" % path)


def main():
    ap = argparse.ArgumentParser(description="LTE cell identity decoder (BladeRF + srsRAN).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--freq-mhz", type=float, nargs="+", help="carrier centre frequency/frequencies in MHz")
    what.add_argument("--earfcn", type=int, nargs="+", help="carrier channel number(s) (EARFCN)")
    what.add_argument("--cells", nargs="+", metavar="EARFCN:PCI,PCI",
                      help="look for these cells, e.g. from a phone's list: 6175:216,272 1617:389,461 - every "
                           "PCI is looked for directly even when the search does not find it, and a check table "
                           "is printed at the end")
    ap.add_argument("--gain", type=float, default=20.0,
                    help="starting receive gain in dB (default 20; adjusted automatically unless --no-auto-gain)")
    ap.add_argument("--no-auto-gain", action="store_true", help="keep --gain fixed")
    ap.add_argument("--secs", type=float, default=1.0,
                    help="recording length in s (default 1.0; longer gives overlapped sectors more chances)")
    ap.add_argument("--dual", action="store_true", help="use both antennas (RX1 + RX2), combined per cell")
    ap.add_argument("--low-if-mhz", type=float, default=0.6,
                    help="tune this far above the carrier to keep the radio's DC spike off it (default 0.6)")
    ap.add_argument("--search-secs", type=float, default=0.08,
                    help="length of recording used for the cell search, s (default 0.08)")
    ap.add_argument("--no-weak-search", action="store_true", help="skip the second search for weaker cells")
    ap.add_argument("--pss-min", type=float, default=None, help="weak search threshold (default 6.0)")
    ap.add_argument("--no-siblings", action="store_true",
                    help="do not try the other sectors of each decoded mast")
    ap.add_argument("--show-failed-siblings", action="store_true",
                    help="also list sector tries that did not decode")
    ap.add_argument("--loop", action="store_true",
                    help="repeat until Ctrl-C (each pass records and decodes again); the table counts the passes "
                         "each cell decoded in")
    ap.add_argument("--only", action="store_true",
                    help="with --cells: decode only the listed PCIs (other cells found are not decoded)")
    ap.add_argument("--no-cancel", action="store_true",
                    help="skip interference cancellation (removing decoded cells' reference signals to find weaker "
                         "cells underneath)")
    ap.add_argument("--retries", type=int, default=2,
                    help="extra recordings for cells whose MIB decoded but not their Cell ID (default 2)")
    ap.add_argument("--srsran-timeout", type=int, default=30, help="seconds per srsRAN attempt (default 30)")
    ap.add_argument("--srsran-log", default=None, help="save every srsRAN output to this folder (diagnosis)")
    ap.add_argument("--debug-search", action="store_true", help="print every weak-search peak and its result")
    ap.add_argument("--csv", default="~/field/lte_cells.csv",
                    help="CSV file results are added to (default ~/field/lte_cells.csv, '' = none)")
    pcap_out.add_args(ap)
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    if args.no_color:
        C.set_colour(False)
    if not pdsch_ue().is_file():
        sys.exit("srsRAN pdsch_ue was not found at %s - build it first (see README)." % pdsch_ue())
    lte_fastsync.install()
    clean_leftovers()
    C.be_nice()
    freqs = list(args.freq_mhz or [])
    args.targets = {}
    for item in args.cells or []:
        try:
            e, pcis = item.split(":")
            args.targets.setdefault(int(e), []).extend(int(p) for p in pcis.split(",") if p.strip())
        except ValueError:
            sys.exit("--cells takes EARFCN:PCI,PCI (e.g. 6175:216,272), not '%s'" % item)
    args.earfcn = list(args.earfcn or []) + [e for e in args.targets if e not in (args.earfcn or [])]
    for e in args.earfcn or []:
        f = earfcn_to_mhz(e)
        if f is None:
            sys.exit("EARFCN %d is not an LTE downlink channel" % e)
        freqs.append(f)
    PCAP[0] = pcap_out.open_from_args(args, "lte")
    C.stop_on_first_ctrl_c()
    t_all, rows, passes, finished = time.time(), [], 0, 0
    try:
        while True:
            passes += 1
            if args.loop:
                C.heading("Pass %d   (Ctrl-C to stop)" % passes)
            for f in freqs:
                new = decode_carrier(args, f)
                for r in new:
                    r["pass"] = passes
                rows += new
                if args.csv and new:
                    write_csv(args.csv, new)
            finished = passes
            if not args.loop:
                break
            C.info("pass %d done - cells decoded so far: %s" % (passes, ", ".join(
                "%s:%d" % (e, p) for e, p in sorted({(r["earfcn"], r["pci"]) for r in rows if r.get("ident")}))
                or "none"))
    except KeyboardInterrupt:
        print("")
        C.warn("stopped - showing the results so far")
    finally:
        C.ignore_ctrl_c()
    # one row per cell: its best result over all passes, and how many passes decoded it
    best = {}
    for r in rows:
        key = (r.get("earfcn"), r.get("pci"))
        old = best.get(key)
        r["decoded_passes"] = (old or {}).get("decoded_passes", 0) + (1 if r.get("ident") else 0)
        if old is None or (r.get("ident") and not old.get("ident")) or (not old.get("ident") and r.get("mib")):
            best[key] = r
        else:
            old["decoded_passes"] = r["decoded_passes"]
    rows_best = list(best.values())
    summary(rows_best, t_all)
    if args.targets:
        phone_check(args, rows_best, max(1, finished) if args.loop else None)
    pcap_out.report(PCAP[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
