#!/usr/bin/env python3
"""LTE multi-band scanner for the BladeRF: find every LTE carrier, decode its identity, print one table.

  python3 ~/fbs-upload/lte_scan.py                         # known UK carriers in bands 20 8 3 1 7
  python3 ~/fbs-upload/lte_scan.py --bands 20 3            # only some bands
  python3 ~/fbs-upload/lte_scan.py --earfcns 6300 6400     # just these channels (fastest)
  python3 ~/fbs-upload/lte_scan.py --discover              # also survey for carriers not in the list
  python3 ~/fbs-upload/lte_scan.py --loop                  # keep re-decoding the live carriers (Ctrl-C)

How it works, per band:
  1. survey   - (--discover only) wide 30.72 Msps snapshots, power in 100 kHz bins -> raised blocks = candidate carriers
                (plus the known UK list, which alone is the default).
  2. probe    - 0.2 s capture at the candidate (low-IF 0.6 MHz, DC spike notched), LTE sync (PSS/SSS)
                search in Python, nudging the centre +-0.3 MHz digitally. No sync -> "no LTE", next.
  3. decode   - for EVERY PCI heard (--max-pci, default 3): srsRAN pdsch_ue reads MIB + SIB1 -> Cell ID, eNB, sector, TAC, PLMN.
                If the MIB says the cell is wider than the capture, re-capture at 23.04/30.72 Msps.
Every result is appended to --csv (default ~/field/lte_scan.csv). Needs ~/srsRAN_4G built (see HANDOFF).
Only one program can use the BladeRF at a time: stop GSM tracking first.
"""
import argparse
import os
import re
import sys
import time
import types

import numpy as np

import lte_decode as L                     # capture(), condition(), osprey 2 path setup
import hackrf_lte_decode as H              # Owen's srsRAN broadcast reader
from lte_pci_search import search_samples

# band: (DL low MHz, DL high MHz, EARFCN offset)   FDD bands used in the UK
BANDS = {
    1: (2110.0, 2170.0, 0),
    3: (1805.0, 1880.0, 1200),
    7: (2620.0, 2690.0, 2750),
    8: (925.0, 960.0, 3450),
    20: (791.0, 821.0, 6150),
}
# UK carriers seen by the phone scanner at home 2026-11-25 (EARFCN -> centre). Default candidate list;
# --discover adds survey-found blocks on top (for new places).
KNOWN_EARFCN = {
    20: [6175, 6225, 6300, 6400],     # EE 793.5, EE 798.5, Vodafone 806, O2 816
    8: [3501, 3725],                  # Vodafone 930.1, O2 952.5
    3: [1363, 1617],                  # Three 1821.3, EE 1846.7
    1: [223, 347],                    # Vodafone 2132.3, O2 2144.7
    7: [],
}
OPERATOR = {"234-10": "O2", "234-15": "Vodafone", "234-20": "Three", "234-30": "EE", "234-33": "EE",
            "234-31": "EE", "234-32": "EE", "234-58": "Manx", "234-50": "JT"}
RATE_FOR_RB = [(50, 15.36e6), (75, 23.04e6), (100, 30.72e6)]
SIB1_RE = re.compile(r"SIB1 cell identity (\d+)\s+eNB (\d+)\s+sector (\d+)\s+TAC (\d+)((?:[ \t]+PLMN \d+-\d+)*)[ \t]*\r?\n")
SURVEY_RATE = 30.72e6
SURVEY_USABLE = 24e6


def band_of(mhz):
    for b, (lo, hi, _) in BANDS.items():
        if lo <= mhz <= hi:
            return b
    return None


def mhz_of(earfcn):
    for b, (lo, hi, off) in BANDS.items():
        n = (earfcn - off) / 10.0 + lo
        if 0 <= earfcn - off and lo <= n <= hi:
            return round(n, 1)
    return None


def earfcn_of(mhz):
    b = band_of(mhz)
    if b is None:
        return None
    lo, _, off = BANDS[b]
    return int(round((mhz - lo) * 10)) + off


# ---------------------------------------------------------------- 1. survey
def survey_band(lo, hi, gain, thresh_db, min_mhz):
    """Return candidate carrier centres (MHz) from raised blocks in [lo, hi]."""
    bins = []
    centre = lo * 1e6 + SURVEY_USABLE / 2
    while centre - SURVEY_USABLE / 2 < hi * 1e6:
        s, _ = L.capture(centre, SURVEY_RATE, 0.03, gain)
        if s is not None:
            nfft = 4096
            nblk = len(s) // nfft
            x = s[:nblk * nfft].reshape(nblk, nfft) * np.hanning(nfft)
            psd = np.fft.fftshift((np.abs(np.fft.fft(x, axis=1)) ** 2).mean(axis=0))
            axis = centre + (np.arange(nfft) - nfft / 2) * SURVEY_RATE / nfft
            for b in np.arange(centre - SURVEY_USABLE / 2, centre + SURVEY_USABLE / 2, 0.1e6):
                if lo * 1e6 <= b < hi * 1e6:
                    sel = (axis >= b) & (axis < b + 0.1e6)
                    bins.append((b / 1e6 + 0.05, 10 * np.log10(psd[sel].mean() + 1e-20)))
        centre += SURVEY_USABLE
    if not bins:
        return [], None
    bins.sort()
    floor = np.percentile([p for _, p in bins], 10)
    up = [(f, p - floor >= thresh_db) for f, p in bins]
    runs, start = [], None
    for i, (f, hot) in enumerate(up):
        if hot and start is None:
            start = f
        if (not hot or i == len(up) - 1) and start is not None:
            end = f if hot else up[i - 1][0]
            if end - start + 0.1 >= min_mhz:
                runs.append((start - 0.05, end + 0.05))
            start = None
    cands = []
    for a, b in runs:
        w = b - a
        if w <= 21.0:
            cands.append(round((a + b) / 2, 1))
        else:                                  # several carriers side by side: try 20 and 10 MHz tilings
            for step in (20.0, 10.0):
                c = a + step / 2
                while c + step / 2 <= b + 0.5:
                    cands.append(round(c, 1))
                    c += step
    return cands, runs


# ---------------------------------------------------------------- 2. probe
def shift(samples, rate, hz):
    if not hz:
        return samples
    n = np.arange(len(samples))
    return (samples * np.exp(-2j * np.pi * hz * n / rate)).astype(np.complex64)


def probe(samples, rate, nudges_mhz, full_retry=True):
    """LTE sync search on the first 40 ms, trying small centre corrections. -> (nudge_hz, hits).
    If nothing is found, one more try on the WHOLE capture at the nominal centre (816 was missed once
    by the 40 ms probe in a run where it decoded fine alone)."""
    t0 = time.time()
    try:
        head = samples[:int(rate * 0.04)]
        for nudge in nudges_mhz:
            rep = search_samples(shift(head, rate, nudge * 1e6), rate)
            if rep["hits"]:
                return nudge * 1e6, rep["hits"]
        if full_retry:
            rep = search_samples(samples, rate)
            if rep["hits"]:
                return 0.0, rep["hits"]
        return None, []
    finally:
        T["probe"] += time.time() - t0


EARLY = [True]
DEBUG = [None]
T = {"capture": 0.0, "probe": 0.0, "align": 0.0, "srsran": 0.0}


def run_pdsch(path, pci, n_prb, n_sf, timeout, early):
    """pdsch_ue on a file. With early=True it is stopped as soon as a SIB1 identity line appears
    (SIB1 repeats every 20 ms; the rest of the file adds nothing). Returns all output text."""
    import select
    import shutil
    import subprocess
    cmd = [str(H._pdsch()), "-d", "-i", str(path), "-p", str(n_prb), "-c", str(pci), "-P", "2", "-Q",
           "-n", str(n_sf)]
    if shutil.which("stdbuf"):
        cmd = ["stdbuf", "-oL", "-eL"] + cmd           # line-buffered so the MIB lines arrive before a stop
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, start_new_session=True)
    out = []
    deadline = time.time() + timeout
    buf = b""
    try:
        while True:
            left = deadline - time.time()
            if left <= 0:
                break
            ready, _, _ = select.select([p.stdout], [], [], min(left, 0.5))
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
        rest = b""
        if p.poll() is None:
            try:
                os.killpg(p.pid, 9)
            except OSError:
                p.kill()
        else:
            try:
                rest = p.stdout.read() or b""
            except Exception:
                pass
        p.stdout.close()
        try:
            p.wait(timeout=2)
        except Exception:
            pass
    return (buf + rest).decode(errors="replace")


# ---------------------------------------------------------------- 3. decode
def limited_subprocess(timeout):
    real = H.subprocess.run

    def run(cmd, **kw):
        kw.setdefault("timeout", timeout)
        try:
            return real(cmd, **kw)
        except H.subprocess.TimeoutExpired as exc:
            out = exc.stdout or b""
            err = exc.stderr or b""
            out = out.decode(errors="replace") if isinstance(out, bytes) else out
            err = err.decode(errors="replace") if isinstance(err, bytes) else err
            return types.SimpleNamespace(stdout=out, stderr=err, returncode=-1)
    return real, run


# srsRAN reads a file at the standard rate for the cell width (-p N_PRB)
SRS_RATE = {6: 1.92e6, 15: 3.84e6, 25: 7.68e6, 50: 15.36e6, 75: 23.04e6, 100: 30.72e6}


def decimate(samples, k):
    """Integer down-sampling with a windowed-sinc low-pass (no scipy needed)."""
    if k == 1:
        return samples
    ntap = 16 * k + 1
    n = np.arange(ntap) - (ntap - 1) / 2
    h = np.sinc(n / k * 0.92) * np.hamming(ntap)
    h /= h.sum()
    return np.convolve(samples, h, mode="same")[::k].astype(np.complex64)


def srsran_read(samples, rate, pci, cfo_hz, n_prb, timeout):
    """One pdsch_ue run on an aligned file. -> (mib dict, identity dict or None, problem text)."""
    import tempfile
    from pathlib import Path
    from lte_cellid import _derotate, _locate_wide
    t_al = time.time()
    corrected = _derotate(np.asarray(samples), cfo_hz, rate)
    start = _locate_wide(corrected, pci, int(round(rate / 15000.0)))     # FFT size for this rate
    aligned = None if start is None else np.asarray(corrected[start:], dtype=np.complex64)
    if aligned is None:
        T["align"] += time.time() - t_al
        return {}, None, "frame start not found"
    with tempfile.NamedTemporaryFile(suffix=".cf32", delete=False, dir="/dev/shm") as fh:
        path = Path(fh.name)
    try:
        aligned.tofile(path)
        T["align"] += time.time() - t_al
        n_sf = max(40, len(aligned) // max(1, int(round(rate / 1000.0))))
        t_sr = time.time()
        H._clear_sib_bits()
        text = run_pdsch(path, pci, n_prb, n_sf, timeout, EARLY[0])
        sib = H._saved_sib_bits()
        payload = sib.read_bytes() if sib is not None else None
        T["srsran"] += time.time() - t_sr
        if DEBUG[0]:
            d = os.path.expanduser(DEBUG[0])
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "pdsch_%s_pci%d_p%d.txt" % (time.strftime("%H%M%S"), pci, n_prb)), "w") as fh:
                fh.write(text)
    finally:
        path.unlink(missing_ok=True)
    mib = H._read_mib(text)
    # our pdsch_ue build (srsran_watch/print_sib1.py) prints the decoded SIB1 as a text line
    m = SIB1_RE.search(text)
    if m:
        plmns = re.findall(r"PLMN (\d+)-(\d+)", m.group(5) or "")
        mcc, mnc = plmns[0] if plmns else ("?", "?")
        ident = {"cell_id": int(m.group(1)), "enb_id": int(m.group(2)), "sector_id": int(m.group(3)),
                 "tac": int(m.group(4)), "mcc": mcc, "mnc_digits": mnc,
                 "operator": OPERATOR.get("%s-%s" % (mcc, mnc), ""),
                 "plmns": ["%s-%s" % x for x in plmns]}
        return mib, ident, ""
    if payload:                                   # Owen's bit-file variant of the patch
        ident = H.identity_from_transport_block(payload)
        if ident:
            return mib, ident, ""
    if not mib:
        return {}, None, "MIB not decoded"
    return mib, None, "MIB ok, no SIB1"


def read_cell(s, rate, mhz, pci, cfo, args):
    """MIB + SIB1 for one PCI; re-reads at the cell's own rate if it is not 10 MHz."""
    out = {"pci": pci}
    mib, ident, problem = srsran_read(s, rate, pci, cfo, 50, args.srsran_timeout)
    if mib:
        out["rb"], out["bw"] = mib.get("prb"), mib.get("mhz")
    prb = mib.get("prb") if mib else None
    if not ident and prb in SRS_RATE and prb != 50:
        want = SRS_RATE[prb]
        if want < rate:                       # 1.4/3/5 MHz: down-sample what we have
            s2, r2 = decimate(s, int(round(rate / want))), want
        else:                                 # 15/20 MHz: needs a wider capture
            print("    PCI %d is %s MHz wide - re-capturing at %.2f Msps" % (pci, mib.get("mhz"), want / 1e6),
                  flush=True)
            s2, tc = L.capture(mhz * 1e6, want, args.secs, args.gain, args.low_if_mhz * 1e6)
            T["capture"] += tc
            r2 = want
        if s2 is not None:
            _, ident, problem = srsran_read(s2, r2, pci, cfo, prb, args.srsran_timeout)
    if ident:
        out.update(ident)
        out["status"] = "OK"
    else:
        out["status"] = problem or "no SIB1"
    return out


def decode_carrier(mhz, args):
    """Capture + probe + srsRAN for every PCI heard (strongest first). Returns a list of result dicts."""
    base = {"mhz": mhz, "band": band_of(mhz), "earfcn": earfcn_of(mhz), "status": "", "t": 0.0}
    t0 = time.time()
    rate = 15.36e6
    s, tc = L.capture(mhz * 1e6, rate, args.secs, args.gain, args.low_if_mhz * 1e6)
    T["capture"] += tc
    if s is None:
        base.update(status="capture failed", t=time.time() - t0)
        return [base]
    nudge, hits = probe(s, rate, args.nudges)
    if not hits:
        base.update(status="no LTE", t=time.time() - t0)
        return [base]
    mhz_true = round(mhz + nudge / 1e6, 1)
    base.update(mhz=mhz_true, band=band_of(mhz_true), earfcn=earfcn_of(mhz_true))
    s = shift(s, rate, nudge)
    out = []
    for hit in hits[:max(1, args.max_pci)]:
        t1 = time.time()
        r = dict(base)
        r.update(read_cell(s, rate, mhz_true, int(hit["pci"]), float(hit.get("offset_hz", 0.0)), args))
        r["t"] = time.time() - (t0 if not out else t1)
        out.append(r)
    if len(hits) > args.max_pci:
        out[-1]["also_pci"] = [int(h["pci"]) for h in hits[args.max_pci:]]
    return out


# ---------------------------------------------------------------- output
HEAD = "%-5s %-8s %-6s %-4s %-4s %-14s %-7s %-6s %-8s %-6s %-10s %s" % (
    "Band", "MHz", "EARFCN", "PCI", "BW", "Operator", "PLMN", "TAC", "eNB", "Sector", "CellID", "time")


def row(r):
    if r.get("status") != "OK":
        return "%-5s %-8.1f %-6s %-4s %-4s %s  (%.1fs)" % (
            "B%s" % r["band"], r["mhz"], r["earfcn"], r.get("pci", "-"), r.get("bw") or "-", r["status"], r["t"])
    return "%-5s %-8.1f %-6s %-4s %-4s %-14s %-7s %-6s %-8s %-6s %-10s %.1fs" % (
        "B%s" % r["band"], r["mhz"], r["earfcn"], r["pci"], r.get("bw") or "-",
        (r.get("operator") or "?")[:14], "%s-%s" % (r["mcc"], r["mnc_digits"]), r["tac"], r["enb_id"],
        r["sector_id"], r["cell_id"], r["t"])


def log_csv(path, r, label):
    path = os.path.expanduser(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new:
            fh.write("utc,label,band,mhz,earfcn,pci,bw_mhz,status,operator,plmn,tac,enb,sector,cell_id,secs\n")
        ok = r.get("status") == "OK"
        fh.write("%s,%s,%s,%.1f,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%.1f\n" % (
            time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()), label or "", r["band"], r["mhz"], r["earfcn"],
            r.get("pci", ""), r.get("bw") or "", r["status"].replace(",", ";"),
            r.get("operator", "") if ok else "", ("%s-%s" % (r["mcc"], r["mnc_digits"])) if ok else "",
            r.get("tac", "") if ok else "", r.get("enb_id", "") if ok else "", r.get("sector_id", "") if ok else "",
            r.get("cell_id", "") if ok else "", r["t"]))


def discover(args):
    cands = []
    for b in args.bands:
        lo, hi, _ = BANDS[b]
        known = [mhz_of(e) for e in KNOWN_EARFCN.get(b, [])]
        if not args.discover:
            if known:
                print("Band %d: known carriers %s" % (b, ", ".join("%.1f" % k for k in known)), flush=True)
            cands += known
            continue
        t0 = time.time()
        found, runs = survey_band(lo, hi, args.gain, args.thresh_db, args.min_mhz)
        extra = [k for k in known if all(abs(k - c) > 0.35 for c in found)]
        print("Band %d (%g-%g MHz): survey %.1fs, raised blocks %s -> try %s%s" % (
            b, lo, hi, time.time() - t0,
            ", ".join("%.1f-%.1f" % r for r in runs) if runs else "none",
            ", ".join("%.1f" % c for c in found) if found else "-",
            (" + known %s" % ", ".join("%.1f" % k for k in extra)) if extra else ""), flush=True)
        cands += found + extra
    return cands


def main():
    ap = argparse.ArgumentParser(description="LTE multi-band scan (BladeRF + srsRAN).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--bands", type=int, nargs="+", default=[20, 8, 3, 1, 7], choices=sorted(BANDS))
    ap.add_argument("--freqs", type=float, nargs="+", help="decode these centres (MHz), no survey")
    ap.add_argument("--gain", type=float, default=40.0, help="BladeRF gain (default 40, as the 806/816 decodes)")
    ap.add_argument("--secs", type=float, default=0.2, help="capture per carrier (default 0.2 s = 10 SIB1s)")
    ap.add_argument("--low-if-mhz", type=float, default=0.6, help="tune off-centre to dodge DC (default 0.6)")
    ap.add_argument("--thresh-db", type=float, default=6.0, help="survey: dB above floor counts as raised")
    ap.add_argument("--min-mhz", type=float, default=4.0, help="survey: narrowest block that can be LTE")
    ap.add_argument("--discover", action="store_true",
                    help="also survey each band for carriers not in the known UK list (new places; slower)")
    ap.add_argument("--earfcns", type=int, nargs="+", help="decode these EARFCNs (e.g. 6300 1617), no survey")
    ap.add_argument("--max-pci", type=int, default=3,
                    help="decode up to this many PCIs per carrier, strongest first (default 3)")
    ap.add_argument("--nudges", type=float, nargs="+", default=[0, 0.1, -0.1],
                    help="centre corrections (MHz) tried by the sync probe")
    ap.add_argument("--no-early-stop", action="store_true",
                    help="let pdsch_ue read the whole capture even after SIB1 decoded (old behaviour, slower)")
    ap.add_argument("--timing", action="store_true", help="print where the time goes for every carrier")
    ap.add_argument("--debug-dir", help="save every pdsch_ue output here (for failed MIB/SIB1 rows)")
    ap.add_argument("--srsran-timeout", type=int, default=20, help="seconds per pdsch_ue attempt")
    ap.add_argument("--loop", action="store_true", help="repeat until Ctrl-C")
    ap.add_argument("--rediscover", type=int, default=0, help="--loop: re-survey every N passes (0 = never)")
    ap.add_argument("--label", default=None, help="tag the CSV rows (e.g. home, field_site2)")
    ap.add_argument("--csv", default="~/field/lte_scan.csv")
    args = ap.parse_args()

    EARLY[0] = not args.no_early_stop
    DEBUG[0] = args.debug_dir
    if not H._pdsch().is_file():
        sys.exit("srsRAN pdsch_ue not found at %s - build it first (HANDOFF 'LTE srsRAN')" % H._pdsch())

    cells = {}
    passes = 0
    t_all = time.time()
    try:
        if args.earfcns:
            args.freqs = (args.freqs or []) + [mhz_of(e) for e in args.earfcns if mhz_of(e)]
        cands = sorted(set(args.freqs)) if args.freqs else None
        while True:
            passes += 1
            t_pass = time.time()
            if cands is None or (args.rediscover and passes > 1 and (passes - 1) % args.rediscover == 0
                                 and not args.freqs):
                cands = sorted(set(discover(args)))
            print("\n--- pass %d: %d carrier(s) to decode ---" % (passes, len(cands)), flush=True)
            print(HEAD, flush=True)
            got = []
            for c in cands:
                before = dict(T)
                for r in decode_carrier(c, args):
                    print(row(r), flush=True)
                    log_csv(args.csv, r, args.label)
                    if r.get("also_pci"):
                        print("      also heard PCI %s (not decoded, raise --max-pci)"
                              % ", ".join(str(p) for p in r["also_pci"]), flush=True)
                    if r.get("status") == "OK":
                        cells[(r["mhz"], r["cell_id"])] = r
                        got.append(r)
                if args.timing:
                    print("      time: " + "  ".join("%s %.1fs" % (k, T[k] - before[k]) for k in T), flush=True)
            if passes == 1 and not args.freqs:
                # later passes only revisit carriers that had LTE sync (fast tracking)
                live = sorted({round(r["mhz"], 1) for r in got})
                if live:
                    cands = live
            print("pass %d done in %.1fs: %d cell(s) decoded" % (passes, time.time() - t_pass, len(got)),
                  flush=True)
            if not args.loop:
                break
    except KeyboardInterrupt:
        print("")
    el = time.time() - t_all
    print("\n===== %d unique LTE cell(s), %d pass(es), %.0fs =====" % (len(cells), passes, el))
    print("time split: " + "  ".join("%s %.1fs" % kv for kv in T.items())
          + "  other %.1fs" % max(0.0, el - sum(T.values())))
    if cells:
        print(HEAD.rsplit(" ", 1)[0])
        for k in sorted(cells):
            print(row(cells[k]).rsplit(" ", 1)[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
