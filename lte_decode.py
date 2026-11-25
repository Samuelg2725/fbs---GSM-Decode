#!/usr/bin/env python3
"""LTE (4G) cell identity with the BladeRF.

Captures one LTE carrier with bladeRF-cli and decodes it with the pure-Python decoder in
"osprey 2" (lte_pci_search -> PCI, lte_cellid/lte_sib1 -> MIB + SIB1). No srsRAN needed.

  python3 lte_decode.py --self-test                 # check the decoder works on this machine
  python3 lte_decode.py --freq-mhz 806              # one capture + decode at 806 MHz
  python3 lte_decode.py --earfcn 6300               # same, by EARFCN
  python3 lte_decode.py --freq-mhz 806 --loop       # keep decoding the same carrier (tracking)
  python3 lte_decode.py --survey 791-821            # power map, to find where the LTE blocks are
  python3 lte_decode.py --file cap.sc16 --rate 15.36e6 --freq-mhz 806   # decode a saved capture

What you get: PCI (physical cell ID), bandwidth, then from SIB1 the 28-bit Cell ID, eNB ID,
sector, PLMN (MCC-MNC = operator) and TAC. The centre frequency must be the MIDDLE of the LTE
carrier (the sync signals sit in the centre 1.08 MHz).
"""
import argparse
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (os.path.join(HERE, "osprey 2"), os.path.expanduser("~/fbs-upload/osprey 2"),
             os.path.expanduser("~/fbs---GSM-Decode/osprey 2")):
    if os.path.isfile(os.path.join(cand, "lte_cellid.py")):
        sys.path.insert(0, cand)
        break
else:
    sys.exit("can't find the 'osprey 2' folder with lte_cellid.py - run this from the repo folder")

from lte_cellid import decode_cell_identity, format_result   # noqa: E402
from lte_pci_search import search_samples                     # noqa: E402

BLADE_LIBS = "/usr/local/lib"
SC16 = "/dev/shm/_lte.sc16"
# LTE-friendly rates (15 kHz grid) and the widest cell each one fully covers
RATES = [(15.36e6, 50), (23.04e6, 75), (30.72e6, 100)]
BW_MHZ = {6: 1.4, 15: 3, 25: 5, 50: 10, 75: 15, 100: 20}


def earfcn_to_mhz(earfcn):
    from cellular_resolver import lte_carrier
    c = lte_carrier(earfcn)
    return None if c is None else c.downlink_mhz


LAST_PEAK = [0.0]
CUR_FREQ = [0.0]


def capture(freq_hz, rate, secs, gain, low_if_hz=0.0):
    """bladeRF-cli capture -> complex64 baseband. Returns (samples, capture_seconds)."""
    n = int(rate * (secs + 0.01))
    if os.path.exists(SC16):
        os.remove(SC16)
    tune = freq_hz + low_if_hz
    cmd = ["bladeRF-cli", "-e", "set agc rx off",
           "-e", "set frequency rx %d" % int(tune),
           "-e", "set samplerate rx %d" % int(rate),
           "-e", "set bandwidth rx %d" % int(min(rate * 0.9, 56e6)),
           "-e", "set gain rx %d" % int(gain),
           "-e", "rx config file=%s format=bin n=%d samples=16384 buffers=64 xfers=16 timeout=20000" % (SC16, n),
           "-e", "rx start", "-e", "rx wait %d" % int(secs * 3000 + 15000)]
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = BLADE_LIBS + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    t0 = time.time()
    p = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    took = time.time() - t0
    size = os.path.getsize(SC16) if os.path.exists(SC16) else 0
    if size < n * 4 * 0.9:
        print("capture is short (%d of %d bytes). bladeRF-cli said:" % (size, n * 4), flush=True)
        print(p.stdout.decode(errors="replace")[-600:], flush=True)
        return None, took
    raw = np.fromfile(SC16, np.int16).astype(np.float32) / 2048.0
    iq = (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)
    iq = iq[int(rate * 0.01):]                     # drop the first 10 ms (tuning transient)
    pk = float(np.max(np.abs(raw))) if raw.size else 0.0
    LAST_PEAK[0] = pk
    lvl = 10 * np.log10(float(np.mean(np.abs(iq) ** 2)) + 1e-20)
    print("  capture level %.1f dBFS, peak %.2f of full scale%s (gain %g)" % (
        lvl, pk, "  CLIPPING - use a lower --gain" if pk >= 0.99 else "", gain), flush=True)
    return condition(iq, rate, low_if_hz), took


def condition(samples, rate, low_if_hz):
    """Remove DC; if tuned off-centre (low IF), shift the cell back to baseband."""
    samples = samples - np.mean(samples)
    if low_if_hz:
        idx = np.arange(len(samples))
        samples = samples * np.exp(1j * 2 * np.pi * float(low_if_hz) * idx / float(rate))
        # the radio's DC spike is now at +low_if: notch it out (same as Owen's HackRF path)
        from hackrf_capture import notch_tone
        samples = notch_tone(samples, float(low_if_hz), rate)
    return samples.astype(np.complex64)


def load_file(path, rate):
    if path.endswith(".sc16"):
        raw = np.fromfile(path, np.int16).astype(np.float32) / 2048.0
        iq = (raw[0::2] + 1j * raw[1::2]).astype(np.complex64)
    else:
        iq = np.fromfile(path, np.complex64)
    return iq - np.mean(iq)


def decode(samples, rate, want_pci=None, all_pci=False):
    """PCI search then MIB+SIB1. Returns (list of results, timings dict)."""
    t = {}
    t0 = time.time()
    report = search_samples(samples, rate)
    t["pci_search"] = time.time() - t0
    hits = list(report["hits"])
    if want_pci is not None:
        hits = [h for h in hits if h["pci"] == want_pci] or [{"pci": want_pci, "offset_hz": 0.0, "sss_score": 0.0}]
    if not hits:
        return [], t, report
    if not all_pci:
        hits = hits[:1]
    out = []
    t0 = time.time()
    for h in hits:
        res = decode_cell_identity(samples, int(h["pci"]), cfo_hz=float(h.get("offset_hz", 0.0)), sample_rate=rate)
        res["_hit"] = h
        out.append(res)
    t["mib_sib1"] = time.time() - t0
    return out, t, report


def srsran_decode(samples, rate, freq_mhz, args):
    """Owen's proven on-air path: PCI in Python, MIB + SIB1 by srsRAN's pdsch_ue (file mode)."""
    import hackrf_lte_decode as H
    exe = H._pdsch()
    if not exe.is_file():
        print("  srsRAN pdsch_ue not found at %s - see HANDOFF 'LTE srsRAN' to build it" % exe, flush=True)
        return None
    cell_txt = os.path.expanduser("~/sib1_cell.txt")
    before = os.path.getsize(cell_txt) if os.path.exists(cell_txt) else 0
    t0 = time.time()
    # pdsch_ue can wait forever on a file it never locks to: cap each attempt (Owen tries 4 widths)
    import types
    real_run = H.subprocess.run

    def limited_run(cmd, **kw):
        kw.setdefault("timeout", args.srsran_timeout)
        print("  srsRAN attempt: %s" % " ".join(str(c) for c in cmd[1:]), flush=True)
        try:
            return real_run(cmd, **kw)
        except H.subprocess.TimeoutExpired as exc:
            print("    (gave up after %ds)" % args.srsran_timeout, flush=True)
            out = exc.stdout or ""
            err = exc.stderr or ""
            if isinstance(out, bytes):
                out = out.decode(errors="replace")
            if isinstance(err, bytes):
                err = err.decode(errors="replace")
            return types.SimpleNamespace(stdout=out, stderr=err, returncode=-1)

    H.subprocess.run = limited_run
    try:
        text = H.decode_capture(samples, rate, freq_mhz, args.earfcn, args.pci, None, None)
    finally:
        H.subprocess.run = real_run
    print(text, flush=True)
    if os.path.exists(cell_txt) and os.path.getsize(cell_txt) > before:
        with open(cell_txt) as fh:
            fh.seek(before)
            new = fh.read().strip()
        if new:
            print("  srsRAN SIB1 log (~/sib1_cell.txt):", flush=True)
            for line in new.splitlines()[-5:]:
                print("    " + line, flush=True)
    print("  srsRAN path took %.1fs" % (time.time() - t0), flush=True)
    return None


SIB1_RE = __import__("re").compile(
    r"SIB1 cell identity (\d+)\s+eNB (\d+)\s+sector (\d+)\s+TAC (\d+)((?:[ \t]+PLMN \d+-\d+)*)[ \t]*\r?\n")
MIB_PRB_RE = __import__("re").compile(r"PRB:\s+(\d+)")


SRS_EXTRA = []
# sample rate this srsRAN build reads a FILE at, per cell width (compact symbol sizes 128/256/384/768/1024/1536)
SRS_FILE_RATE = {6: 1.92e6, 15: 3.84e6, 25: 5.76e6, 50: 11.52e6, 75: 15.36e6, 100: 23.04e6}


def run_pdsch(exe, path, pci, n_prb, n_sf, timeout, early, ports=2):
    """pdsch_ue on an aligned file. early=True: stop it as soon as it prints a SIB1 identity line
    (it would otherwise read the whole file and then the next bandwidth attempt would start)."""
    import select
    import shutil
    cmd = [str(exe), "-d", "-i", str(path), "-p", str(n_prb), "-c", str(pci), "-P", str(ports), "-Q", "-n", str(n_sf)]
    if SRS_EXTRA:
        cmd = [c for c in cmd if not (c == "-Q" and "-v" in SRS_EXTRA)] + SRS_EXTRA
    if shutil.which("stdbuf"):
        cmd = ["stdbuf", "-oL", "-eL"] + cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, start_new_session=True)
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
            if early and SIB1_RE.search(buf.decode(errors="replace")):
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


def align_wide(samples, rate, pci, cfo_hz):
    """Owen's alignment (hackrf_lte_decode._align) but with the FFT size for THIS rate (his is fixed at
    1024 = 15.36 Msps; 23.04/30.72 Msps captures for 15/20 MHz cells need 1536/2048)."""
    from lte_cellid import _derotate, _locate_wide
    corrected = _derotate(np.asarray(samples), cfo_hz, rate)
    start = _locate_wide(corrected, pci, int(round(rate / 15000.0)))
    if start is None:
        return None
    return np.asarray(corrected[start:], dtype=np.complex64)


def align_narrow(samples, rate, pci, cfo_hz):
    """Subframe-0 alignment for a WEAK cell. The stock alignment correlates the PSS across the full
    bandwidth, where a stronger cell's whole transmission drowns it. Here: resample to 1.92 Msps (only the
    centre 1.08 MHz, where PSS/SSS live), find this PCI's PSS comb (every 5 ms), decide which PSS is in
    subframe 0 by AVERAGING the SSS over every repeat, then map the time back to the full-rate capture."""
    import lte_sync as S
    from lte_pci_search import resample_to_lte
    from lte_cellid import _derotate, _useful_wide, _subframe_wide, _refine_wide, symbol_useful_offset
    from lte_sync import _pss_freq, sss_freq
    wide = _derotate(np.asarray(samples), cfo_hz, rate)
    x = np.asarray(resample_to_lte(wide, rate), dtype=np.complex128)
    nid2, nid1 = pci % 3, pci // 3
    tpl = S.pss_time(nid2)
    tpl = tpl / np.linalg.norm(tpl)
    n = len(x) + len(tpl) - 1
    nf = 1 << (n - 1).bit_length()
    mag = np.abs(np.fft.ifft(np.fft.fft(x, nf) * np.conj(np.fft.fft(tpl, nf)))[:len(x)])
    # strongest 5 ms comb phase: fold the correlation over 9600 samples
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
    best = None
    for swap in (False, True):                          # which group is subframe 0
        sc = []
        for g in (0, 1):
            if groups[g]:
                sub5 = (g == 1) != swap
                sc.append(S._mean_match(groups[g], sss_freq(nid1, nid2, sub5)))
        score = sum(sc) / max(1, len(sc))
        if best is None or score > best[0]:
            best = (score, swap)
    score, swap = best
    if score < 0.3:
        return None
    t0 = times[1] if swap else times[0]                 # a PSS that sits in subframe 0
    start_n = t0 - symbol_useful_offset(6)
    nfft = int(round(rate / 15000.0))
    pss_at = _useful_wide(6, nfft)
    sub = _subframe_wide(nfft)
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


def srsran_cell(samples, rate, hit, args, aligned=None):
    """srsRAN MIB + SIB1 for ONE PCI (any PCI the search found, not only the strongest)."""
    import tempfile
    from pathlib import Path
    import hackrf_lte_decode as H
    pci = int(hit["pci"])
    out = {"pci": pci, "strength": hit.get("pss_ratio", 0.0), "sss": hit.get("sss_score", 0.0),
           "status": "", "mib": None, "ident": None}
    if aligned is None:
        aligned = align_wide(samples, rate, pci, float(hit.get("offset_hz", 0.0)))
    if aligned is None:
        # weak cell: the full-band alignment is drowned by the strong cell; align on the narrow centre instead
        aligned = align_narrow(samples, rate, pci, float(hit.get("offset_hz", 0.0)))
        out["narrow_align"] = aligned is not None
    if aligned is None:
        out["status"] = "frame start not found"
        return out
    with tempfile.NamedTemporaryFile(suffix=".cf32", delete=False, dir="/dev/shm") as fh:
        path = Path(fh.name)
    try:
        aligned.tofile(path)
        n_sf = max(40, len(aligned) // max(1, int(round(rate / 1000.0))))
        own = {15360000: 50, 23040000: 75, 30720000: 100}.get(int(round(rate)), 50)
        widths = [own] + [w for w in (50, 25, 75, 100) if w != own]
        if args.match_rate:
            # this srsRAN build reads files at its compact rates (log: "Reading 11520 samples" per subframe
            # for 50 PRB = 11.52 Msps). Write each attempt's file at the rate srsRAN expects for that width.
            widths = [w for w in (50, 75, 100, 25, 15, 6) if SRS_FILE_RATE[w] <= rate + 1]
            mib_w = hit.get("mib_width")
            if mib_w in widths:
                widths.remove(mib_w)
                widths.insert(0, mib_w)
            own = widths[0]
        written = {}

        def file_for(n_prb):
            """path + subframe count of the aligned capture at srsRAN's rate for n_prb (match-rate mode)."""
            if not args.match_rate:
                return path, n_sf
            if n_prb not in written:
                from fractions import Fraction
                from scipy.signal import resample_poly
                fr = Fraction(int(SRS_FILE_RATE[n_prb]), int(round(rate))).limit_denominator(64)
                y = aligned if fr == 1 else resample_poly(aligned, fr.numerator, fr.denominator)
                pth = Path(str(path) + ".p%d" % n_prb)
                np.asarray(y, dtype=np.complex64).tofile(pth)
                written[n_prb] = (pth, max(40, len(y) // max(1, int(round(SRS_FILE_RATE[n_prb] / 1000.0)))))
            return written[n_prb]
        # (width, ports) attempts: Owen's fixed 2 ports first; with --try-ports, a MIB failure at the
        # capture's own width is retried with 4 and 1 ports (file mode must be told the cell's port count)
        attempts = [(w, 2) for w in widths]
        if args.try_ports:
            attempts = [(own, 2), (own, 4), (own, 1)] + [(w, 2) for w in widths if w != own]
        for n_prb, ports in attempts:
            if ports != 2 and (out["mib"] is not None or out["ident"]):
                continue                   # MIB already decoded with 2 ports: port count is not the problem
            fpath, fsf = file_for(n_prb)
            text = run_pdsch(H._pdsch(), fpath, pci, n_prb, fsf, args.srsran_timeout, args.fast_srsran, ports)
            if args.srsran_log:
                d = os.path.expanduser(args.srsran_log)
                os.makedirs(d, exist_ok=True)
                with open(os.path.join(d, "pdsch_%s_%gMHz_%.0fMsps_pci%d_p%d_P%d.txt" % (
                        time.strftime("%H%M%S"), CUR_FREQ[0], rate / 1e6, pci, n_prb, ports)), "w") as fh:
                    fh.write(text)
            vals = [int(v) for v in MIB_PRB_RE.findall(text) if int(v) in (6, 15, 25, 50, 75, 100)] \
                if "Decoded MIB" in text else []
            if vals and out["mib"] is None:
                out["mib"] = vals[-1]
            elif "Decoded MIB" in text and out["mib"] is None:
                out["mib_odd"] = MIB_PRB_RE.findall(text)        # MIB decoded, width not understood
            s1 = SIB1_RE.search(text)
            if ("Decoded MIB" in text or s1) and ports != 2:
                out["ports"] = ports
            if s1:
                plmn = " ".join(x.strip() for x in (s1.group(5) or "").split("PLMN") if x.strip())
                out["ident"] = "CellID %s (eNB %s sector %s)  TAC %s  PLMN %s" % (
                    s1.group(1), s1.group(2), s1.group(3), s1.group(4), plmn or "?")
                break
            if args.fast_srsran and out["mib"] is not None and out["mib"] == n_prb:
                break                      # right width already tried and no SIB1: other widths won't help
            if out["mib"] is not None and (SRS_FILE_RATE.get(out["mib"], 0) > rate + 1 if args.match_rate
                                           else needed_rate(out["mib"]) > rate):
                break                      # cell wider than this capture: needs a wider re-capture, not more tries
    finally:
        path.unlink(missing_ok=True)
        for pth, _n in written.values():
            pth.unlink(missing_ok=True)
    out["_aligned"] = aligned
    if out.get("ports") and out["ident"]:
        out["ident"] += "  [%d antenna ports]" % out["ports"]
    if out["ident"]:
        out["status"] = "OK"
    elif out["mib"]:
        out["status"] = "MIB ok (%d RB), no SIB1" % out["mib"]
    elif out.get("mib_odd"):
        out["status"] = "MIB decoded, width unclear (srsRAN printed PRB %s) - trying a 20 MHz capture" % (
            "/".join(out["mib_odd"]) or "?")
    else:
        out["status"] = "MIB not decoded"
    return out


def search_all_cells(samples, rate, pss_min=None, debug=False):
    """Owen's sync search, plus a second pass for WEAKER cells on the same carrier.

    The stock search picks one frequency offset per PSS group (the highest correlation anywhere). A strong
    cell's false peak at a wrong offset can beat a weak cell's real peak, so the weak cell is never tested
    at the right frequency. All cells on one carrier share the carrier's frequency, so the strongest cell's
    offset is right for every cell: search all 3 PSS groups again at that offset (+-2 kHz)."""
    import lte_sync as S
    from lte_pci_search import resample_to_lte
    report = search_samples(samples, rate)
    hits = list(report["hits"])
    pmin = S.PSS_MIN_RATIO if pss_min is None else float(pss_min)
    if not hits:
        if debug:
            print("  [debug] stock search found nothing (best PSS %.1fx at %+.1f kHz) - second pass at that offset" % (
                report["best_ratio"], report["best_offset_hz"] / 1e3), flush=True)
        # nothing confirmed: still try a second pass at the best peak's offset (the capture is tuned to the
        # carrier, so a real cell's offset is small) with the (possibly lowered) threshold
        hits_seed = [{"offset_hz": float(report.get("best_offset_hz", 0.0)), "pci": -1}]
    else:
        hits_seed = hits
    if debug:
        print("  [debug] stock search: %s; carrier offset %+.1f kHz; second pass threshold %.1f, SSS needs %.2f" % (
            ", ".join("PCI %d %.1fx sss %.2f" % (h["pci"], h["pss_ratio"], h["sss_score"]) for h in hits),
            hits[0]["offset_hz"] / 1e3, pmin, S.SSS_MIN_SCORE), flush=True)
    x = np.asarray(resample_to_lte(np.asarray(samples), rate), dtype=np.complex128)
    length = int(x.size)
    scale = float(np.sqrt(np.mean(np.abs(x) ** 2) + 1e-12))
    nfft = 1 << (length + S.FFT_SIZE - 2).bit_length()
    spectrum = np.fft.fft(x, nfft)
    off0 = float(hits_seed[0]["offset_hz"])
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
                if debug:
                    print("  [debug] group %d: next peak %.1fx < %.1f - stop" % (nid2, ratio, pmin), flush=True)
                break
            nid1, score, pss_seq = S._confirm_sss(x, work, at, offset, nid2)
            S._blank_comb(work, at)
            if debug:
                ok = score >= S.SSS_MIN_SCORE and pss_seq >= S.SSS_MIN_SCORE
                print("  [debug] group %d: peak %.1fx at %d -> PCI %d  sss %.2f pss-seq %.2f  %s" % (
                    nid2, ratio, at, 3 * nid1 + nid2, score, pss_seq, "ACCEPT" if ok else "reject"), flush=True)
            if score < S.SSS_MIN_SCORE or pss_seq < S.SSS_MIN_SCORE:
                continue
            pci = 3 * nid1 + nid2
            if pci in have:
                continue
            have.add(pci)
            hits.append({"pci": pci, "nid1": nid1, "nid2": nid2, "pss_peak": peak, "pss_ratio": ratio,
                         "pss_seq": pss_seq, "sss_score": score, "offset_hz": offset, "second_pass": True})
    hits.sort(key=lambda h: h["pss_ratio"], reverse=True)
    report = dict(report)
    report["hits"] = hits
    return report


def srsran_strongest_fast(samples, rate, freq_mhz, args):
    """--fast-srsran without --all-pci: strongest cell only (as the default), but srsRAN stops at the ID."""
    t0 = time.time()
    hits = search_samples(samples, rate)["hits"]
    if not hits:
        print("  no LTE cell found at %.3f MHz" % freq_mhz, flush=True)
        return None
    res = srsran_cell(samples, rate, hits[0], args)
    print("  PCI %-3d strength %4.0f  %s" % (res["pci"], res["strength"], res["ident"] or res["status"]), flush=True)
    print("  srsRAN path took %.1fs" % (time.time() - t0), flush=True)
    return None


def srsran_all(samples, rate, freq_mhz, args):
    """--all-pci: every cell the sync search hears on this carrier, each decoded by srsRAN."""
    import hackrf_lte_decode as H
    if not H._pdsch().is_file():
        print("  srsRAN pdsch_ue not found at %s - see HANDOFF 'LTE srsRAN' to build it" % H._pdsch(), flush=True)
        return None
    t0 = time.time()
    s_in = samples[:int(rate * args.search_secs)] if args.search_secs else samples
    if args.weak_search:
        report = search_all_cells(s_in, rate, args.pss_min, args.debug_search)
    else:
        report = search_samples(s_in, rate)
    if args.debug_search and not args.weak_search:
        search_all_cells(s_in, rate, args.pss_min, True)
    hits = report["hits"]
    print("  sync search %.1fs: %d cell(s) heard on %.1f MHz: %s" % (
        time.time() - t0, len(hits), freq_mhz,
        ", ".join("PCI %d (strength %.0f)%s" % (h["pci"], h["pss_ratio"], "*" if h.get("second_pass") else "")
                  for h in hits) or "none"), flush=True)
    if any(h.get("second_pass") for h in hits):
        print("  (* = found only by --weak-search)", flush=True)
    if not hits:
        print("  (best sync peak %.1fx, needs %.1fx)" % (report["best_ratio"], report["min_ratio"]), flush=True)
    tried = {}
    for h in hits:
        t1 = time.time()
        res = srsran_cell(samples, rate, h, args)
        tried[res["pci"]] = res
        line = res["ident"] if res["ident"] else res["status"]
        print("  PCI %-3d strength %4.0f  %s  (%.1fs)" % (res["pci"], res["strength"], line, time.time() - t1),
              flush=True)
    if args.siblings:
        # sectors of one mast share its clock (same frame timing) and are usually numbered 3*G + 0/1/2:
        # decode the other two PCIs of every decoded cell's group using THAT cell's frame timing.
        # No sync signal needed; srsRAN's CRC rejects a PCI that is not there.
        for parent in [r for r in list(tried.values()) if r["ident"] and r.get("_aligned") is not None]:
            base = 3 * (parent["pci"] // 3)
            for sib in (base, base + 1, base + 2):
                if sib in tried:
                    continue
                t1 = time.time()
                res = srsran_cell(samples, rate, {"pci": sib, "pss_ratio": 0.0, "sss_score": 0.0,
                                                  "offset_hz": 0.0}, args, aligned=parent["_aligned"])
                tried[sib] = res
                line = res["ident"] if res["ident"] else res["status"]
                print("  PCI %-3d sibling of %-3d  %s  (%.1fs)" % (sib, parent["pci"], line, time.time() - t1),
                      flush=True)
    print("  srsRAN path took %.1fs" % (time.time() - t0), flush=True)
    if args.match_rate:
        need = [needed_rate(r["mib"]) for r in tried.values()
                if not r["ident"] and r["mib"] and SRS_FILE_RATE.get(r["mib"], 0) > rate + 1]
    else:
        need = [needed_rate(r["mib"]) for r in tried.values() if not r["ident"] and r["mib"] and r["mib"] > 50]
    need += [30.72e6 for r in tried.values() if not r["ident"] and not r["mib"] and r.get("mib_odd") is not None
             and rate < 30.72e6]
    return {"need_rate": max(need) if need else None}


def needed_rate(n_rb):
    for rate, max_rb in RATES:
        if n_rb <= max_rb:
            return rate
    return RATES[-1][0]


def one_line(res, freq_mhz):
    h = res["_hit"]
    s = "%.1f MHz  PCI %d (offset %+.1f kHz)" % (freq_mhz, res["pci"], h.get("offset_hz", 0.0) / 1000.0)
    if res.get("bandwidth_rb"):
        s += "  %s MHz" % res["bandwidth_mhz"]
    if res.get("cell_id") is not None:
        s += "  CellID %d (eNB %d sector %d)" % (res["cell_id"], res["enb_id"], res["sector_id"])
        if res.get("mcc") is not None:
            s += "  PLMN %s-%s %s  TAC %s" % (res["mcc"], res["mnc"], res.get("operator") or "", res.get("tac"))
    else:
        s += "  [%s]" % res.get("cell_id_status")
    return s


def run_once(args, freq_mhz, rate):
    """Capture (or load) + decode one time; re-captures at a wider rate if the cell needs it."""
    timing = {}
    CUR_FREQ[0] = freq_mhz
    if args.file:
        samples = condition(load_file(args.file, rate), rate, args.low_if_mhz * 1e6)
    else:
        samples, timing["capture"] = capture(freq_mhz * 1e6, rate, args.secs, args.gain, args.low_if_mhz * 1e6)
        if samples is None:
            return None, timing, rate
    for _round in range(3 if (args.engine == "srsran" and args.auto_gain and not args.file) else 0):
        # clipping -> 10 dB down; quiet (peak < 0.25) -> up by up to 20 dB, aiming for peak ~0.5; max 60
        pk = LAST_PEAK[0]
        new_gain = None
        if pk >= 0.99:
            new_gain = args.gain - 10
        elif 0 < pk < 0.25 and args.gain < 60:
            new_gain = min(60.0, args.gain + max(3.0, min(20.0, round(20 * np.log10(0.5 / pk)))))
        if new_gain is None or new_gain == args.gain:
            break
        if True:
            print("  auto gain: %g -> %g and capturing again" % (args.gain, new_gain), flush=True)
            args.gain = new_gain
            samples, tc = capture(freq_mhz * 1e6, rate, args.secs, args.gain, args.low_if_mhz * 1e6)
            timing["capture"] = timing.get("capture", 0) + tc
            if samples is None:
                return None, timing, rate
    if args.engine == "srsran":
        if args.all_pci or args.fast_srsran:
            if not args.all_pci:
                return srsran_strongest_fast(samples, rate, freq_mhz, args), timing, rate
            res = srsran_all(samples, rate, freq_mhz, args)
            need = res.get("need_rate") if res else None
            if need and need > rate and args.rate == "auto" and not args.file:
                print("  a cell is wider than this %.2f Msps capture - capturing again at %.2f Msps"
                      % (rate / 1e6, need / 1e6), flush=True)
                samples, tc = capture(freq_mhz * 1e6, need, args.secs, args.gain, args.low_if_mhz * 1e6)
                timing["capture"] = timing.get("capture", 0) + tc
                if samples is not None:
                    rate = need
                    srsran_all(samples, rate, freq_mhz, args)
            return None, timing, rate
        return srsran_decode(samples, rate, freq_mhz, args), timing, rate
    results, t, report = decode(samples, rate, args.pci, args.all_pci)
    timing.update(t)
    if (results and not args.file and args.rate == "auto" and results[0].get("bandwidth_rb")
            and needed_rate(results[0]["bandwidth_rb"]) > rate):
        new_rate = needed_rate(results[0]["bandwidth_rb"])
        print("  cell is %s MHz wide - re-capturing at %.2f Msps" % (results[0]["bandwidth_mhz"], new_rate / 1e6),
              flush=True)
        samples, tc = capture(freq_mhz * 1e6, new_rate, args.secs, args.gain, args.low_if_mhz * 1e6)
        timing["capture"] = timing.get("capture", 0) + tc
        if samples is not None:
            rate = new_rate
            results, t, report = decode(samples, rate, args.pci, args.all_pci)
            for k, v in t.items():
                timing[k] = timing.get(k, 0) + v
    if not results:
        print("  no LTE cell found at %.3f MHz (best PSS %.1fx at %+.1f kHz, SSS %.2f, needs %.2f)"
              % (freq_mhz, report["best_ratio"], report["best_offset_hz"] / 1000.0, report["best_sss"],
                 report["sss_min"]), flush=True)
    if args.save and not args.file and os.path.exists(SC16):
        import shutil
        dest = os.path.expanduser(args.save)
        os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
        shutil.copyfile(SC16, dest)          # /dev/shm is a different filesystem: copy, not rename
        print("  saved capture to %s (%.2f Msps sc16)" % (args.save, rate / 1e6), flush=True)
    return results, timing, rate


def survey(args):
    """Power in 0.5 MHz bins across a range, to find LTE blocks (flat raised areas 5-20 MHz wide)."""
    lo, hi = [float(x) for x in args.survey.replace(":", "-").split("-")]
    rate = 30.72e6
    usable = 24e6
    centre = lo * 1e6 + usable / 2
    rows = []
    while centre - usable / 2 < hi * 1e6:
        samples, _ = capture(centre, rate, 0.05, args.gain)
        if samples is not None:
            nfft = 2048
            nblk = len(samples) // nfft
            x = samples[:nblk * nfft].reshape(nblk, nfft) * np.hanning(nfft)
            psd = np.fft.fftshift((np.abs(np.fft.fft(x, axis=1)) ** 2).mean(axis=0))
            axis = centre + (np.arange(nfft) - nfft / 2) * rate / nfft
            for b in np.arange(centre - usable / 2, centre + usable / 2, 0.5e6):
                if lo * 1e6 <= b < hi * 1e6:
                    sel = (axis >= b) & (axis < b + 0.5e6)
                    rows.append((b, 10 * np.log10(psd[sel].mean() + 1e-20)))
        centre += usable
    if not rows:
        return 1
    floor = np.percentile([r[1] for r in rows], 10)
    print("\nPower above the quietest 10%% (each # = 1 dB), %s MHz at gain %g:" % (args.survey, args.gain))
    for b, p in rows:
        d = p - floor
        print("  %7.1f MHz %5.1f dB %s" % (b / 1e6, d, "#" * int(max(0, min(60, d)))))
    print("LTE blocks show as flat raised runs 5/10/15/20 MHz wide; use the MIDDLE of a run as --freq-mhz.")
    return 0


def self_test():
    from lte_cellid import synthesize_cell
    ok = True
    narrow = synthesize_cell(430, 7159822, n_rb=6, ports=1)
    narrow = np.tile(narrow, 3).astype(np.complex64)
    results, t, _ = decode(narrow, 1.92e6)
    got = results[0] if results else {}
    print("1.4 MHz synthetic cell: %s" % (one_line(got, 0.0) if got else "nothing found"))
    ok &= bool(got) and got.get("cell_id") == 7159822 and got.get("pci") == 430
    from hackrf_lte_decode import self_test as sib1_test
    ok2 = sib1_test()
    print("10 MHz SIB1 decode (Vodafone test vector): %s" % ("OK" if ok2 else "FAILED"))
    print("self-test %s" % ("OK" if ok and ok2 else "FAILED"))
    return 0 if ok and ok2 else 1


def main():
    ap = argparse.ArgumentParser(description="LTE cell identity with the BladeRF (pure-Python decoder).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--freq-mhz", type=float, help="centre of the LTE carrier")
    ap.add_argument("--earfcn", type=int, help="LTE channel number instead of --freq-mhz")
    ap.add_argument("--rate", default="auto",
                    help="auto (15.36e6, re-capture wider if the cell needs it) or 15.36e6/23.04e6/30.72e6")
    ap.add_argument("--secs", type=float, default=0.2, help="capture length (default 0.2 = 20 LTE frames)")
    ap.add_argument("--gain", type=float, default=40.0, help="BladeRF gain (default 40)")
    ap.add_argument("--low-if-mhz", type=float, default=0.0,
                    help="tune this far above the cell to keep the radio's DC spike off it (default 0)")
    ap.add_argument("--pci", type=int, help="decode this PCI even if the sync search is unsure")
    ap.add_argument("--all-pci", action="store_true",
                    help="decode every PCI heard, not just the strongest (python engine and srsran engine)")
    ap.add_argument("--weak-search", action="store_true",
                    help="with --all-pci: second sync pass at the strongest cell's frequency offset to find weaker "
                         "cells on the same carrier (the stock search can lock onto a false peak and miss them)")
    ap.add_argument("--pss-min", type=float, default=None,
                    help="--weak-search: sync threshold (default 6.0 = stock); lower (e.g. 4) to try weaker cells")
    ap.add_argument("--search-secs", type=float, default=None,
                    help="run the sync search on only the first N s of the capture (e.g. 0.08 = 80 ms, as Owen's "
                         "cell search; default = whole capture). srsRAN still gets the whole capture")
    ap.add_argument("--debug-search", action="store_true",
                    help="print every sync peak per PSS group and why it was accepted / rejected")
    ap.add_argument("--siblings", action="store_true",
                    help="with --all-pci: also try the other two PCIs of each decoded cell's group (3*G+0/1/2 = "
                         "usually the same mast's other sectors) using that cell's frame timing - finds sectors "
                         "whose own sync signal is too weak to detect")
    ap.add_argument("--auto-gain", action="store_true",
                    help="srsran engine: if the capture clips, -10 dB; if it is very quiet (peak < 0.15), raise the "
                         "gain (up to +20 dB) - then capture again. 800 and 2100 MHz need very different gains")
    ap.add_argument("--srsran-log", default=None,
                    help="save every srsRAN pdsch_ue output to this folder (e.g. ~/field/srslog) for diagnosis")
    ap.add_argument("--srsran-verbose", action="store_true",
                    help="run pdsch_ue with -v (and without -Q) so its sync / MIB attempts are written to the "
                         "--srsran-log files - for diagnosing cells that fail silently")
    ap.add_argument("--try-ports", action="store_true",
                    help="if the MIB fails with 2 antenna ports (Owen's fixed -P 2), retry with 4 and 1 "
                         "(srsRAN file mode must be told the port count; many 1800/2100 cells use 4)")
    ap.add_argument("--match-rate", action="store_true",
                    help="write srsRAN's input at the rate it actually reads (11.52 Msps for 10 MHz cells, 15.36 for "
                         "15 MHz, 23.04 for 20 MHz) instead of the capture rate. Its log says 'Reading 11520 samples' "
                         "per subframe = 11.52 Msps, so 15.36 Msps files were being read 25%% too slowly")
    ap.add_argument("--fast-srsran", action="store_true",
                    help="srsran engine: stop pdsch_ue as soon as the cell identity is printed instead of also "
                         "running the other 3 bandwidth attempts (default = old behaviour)")
    ap.add_argument("--loop", action="store_true", help="keep capturing + decoding until Ctrl-C")
    ap.add_argument("--file", help="decode a saved capture (.sc16 from --save, or complex64)")
    ap.add_argument("--save", help="keep the raw capture here (sc16)")
    ap.add_argument("--survey", metavar="LO-HI", help="power map of a range in MHz, e.g. 791-821")
    ap.add_argument("--engine", choices=["python", "srsran"], default="python",
                    help="python = all-Python decoder (default); srsran = Owen's HackRF path: srsRAN pdsch_ue "
                         "reads MIB + SIB1 from the capture (needs ~/srsRAN_4G built)")
    ap.add_argument("--srsran-timeout", type=int, default=30,
                    help="--engine srsran: seconds per pdsch_ue attempt before giving up (default 30)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.srsran_verbose:
        SRS_EXTRA[:] = ["-v"]
    if args.self_test:
        return self_test()
    if args.survey:
        return survey(args)
    if args.earfcn is not None:
        freq_mhz = earfcn_to_mhz(args.earfcn)
        if freq_mhz is None:
            sys.exit("EARFCN %d is not an LTE downlink channel" % args.earfcn)
    elif args.freq_mhz is not None:
        freq_mhz = args.freq_mhz
    else:
        sys.exit("give --freq-mhz, --earfcn, --survey or --self-test")
    rate = 15.36e6 if args.rate == "auto" else float(args.rate)

    print("LTE decode %.3f MHz%s  rate %s  %.2fs  gain %g%s" % (
        freq_mhz, (" (EARFCN %d)" % args.earfcn) if args.earfcn is not None else "",
        args.rate, args.secs, args.gain, "  loop" if args.loop else ""), flush=True)
    runs = ok = 0
    seen = {}
    t_start = time.time()
    try:
        while True:
            runs += 1
            t0 = time.time()
            results, timing, rate = run_once(args, freq_mhz, rate)
            total = time.time() - t0
            stamp = time.strftime("%H:%M:%S")
            if results:
                for res in results:
                    print("%s  %s" % (stamp, one_line(res, freq_mhz)), flush=True)
                    if res.get("cell_id") is not None:
                        ok += 1
                        seen[res["cell_id"]] = one_line(res, freq_mhz)
                if not args.loop:
                    print("")
                    print(format_result(results[0]))
            print("  timing: %s  total %.1fs" % ("  ".join("%s %.1fs" % kv for kv in timing.items()), total),
                  flush=True)
            if not args.loop or args.file:
                break
    except KeyboardInterrupt:
        print("")
    if args.loop:
        el = time.time() - t_start
        print("\n===== %d run(s) in %.0fs, Cell ID decoded in %d (%.0f%%) =====" % (
            runs, el, ok, 100.0 * ok / max(1, runs)))
        for line in seen.values():
            print("  " + line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
