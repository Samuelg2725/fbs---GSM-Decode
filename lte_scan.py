#!/usr/bin/env python3
"""LTE multi-band scanner for the BladeRF: find every LTE carrier, decode its identity, print one table.

  python3 ~/fbs-upload/lte_scan.py                         # UK bands 20 8 3 1 7, one pass
  python3 ~/fbs-upload/lte_scan.py --bands 20 3            # only some bands
  python3 ~/fbs-upload/lte_scan.py --freqs 806 816         # skip discovery, decode these centres (fastest)
  python3 ~/fbs-upload/lte_scan.py --loop                  # discover once, then keep re-decoding (Ctrl-C)
  python3 ~/fbs-upload/lte_scan.py --loop --rediscover 5   # full re-discovery every 5th pass

How it works, per band:
  1. survey   - wide 30.72 Msps snapshots, power in 100 kHz bins -> raised blocks = candidate carriers
                (plus the known UK centres for that band). Skipped with --freqs.
  2. probe    - 0.2 s capture at the candidate (low-IF 0.6 MHz, DC spike notched), LTE sync (PSS/SSS)
                search in Python, nudging the centre +-0.3 MHz digitally. No sync -> "no LTE", next.
  3. decode   - srsRAN pdsch_ue (Owen's proven path) reads MIB + SIB1 -> Cell ID, eNB, sector, TAC, PLMN.
                If the MIB says the cell is wider than the capture, re-capture at 23.04/30.72 Msps.
Every result is appended to --csv (default ~/field/lte_scan.csv). Needs ~/srsRAN_4G built (see HANDOFF).
Only one program can use the BladeRF at a time: stop GSM tracking first.
"""
import argparse
import os
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
# centres seen on air / Owen's notes (tried even if the survey misses them)
KNOWN = {20: [793.5, 798.5, 806.0, 816.0]}
RATE_FOR_RB = [(50, 15.36e6), (75, 23.04e6), (100, 30.72e6)]
SURVEY_RATE = 30.72e6
SURVEY_USABLE = 24e6


def band_of(mhz):
    for b, (lo, hi, _) in BANDS.items():
        if lo <= mhz <= hi:
            return b
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


def probe(samples, rate, nudges_mhz):
    """LTE sync search on the first 40 ms, trying small centre corrections. -> (nudge_hz, hits)."""
    head = samples[:int(rate * 0.04)]
    for nudge in nudges_mhz:
        rep = search_samples(shift(head, rate, nudge * 1e6), rate)
        if rep["hits"]:
            return nudge * 1e6, rep["hits"]
    return None, []


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
    corrected = _derotate(np.asarray(samples), cfo_hz, rate)
    start = _locate_wide(corrected, pci, int(round(rate / 15000.0)))     # FFT size for this rate
    aligned = None if start is None else np.asarray(corrected[start:], dtype=np.complex64)
    if aligned is None:
        return {}, None, "frame start not found"
    with tempfile.NamedTemporaryFile(suffix=".cf32", delete=False, dir="/dev/shm") as fh:
        path = Path(fh.name)
    real, limited = limited_subprocess(timeout)
    try:
        aligned.tofile(path)
        n_sf = max(40, len(aligned) // max(1, int(round(rate / 1000.0))))
        H.subprocess.run = limited
        text, payload = H._run_broadcast(path, pci, n_prb, n_sf)
    finally:
        H.subprocess.run = real
        path.unlink(missing_ok=True)
    mib = H._read_mib(text)
    if not mib:
        return {}, None, "MIB not decoded"
    if not payload:
        return mib, None, "MIB ok, no SIB1"
    ident = H.identity_from_transport_block(payload)
    return mib, ident, "" if ident else "SIB1 not a UK identity"


def decode_carrier(mhz, args):
    """Capture + probe + srsRAN. Returns a result dict (always, with 'status')."""
    res = {"mhz": mhz, "band": band_of(mhz), "earfcn": earfcn_of(mhz), "status": "", "t": 0.0}
    t0 = time.time()
    rate = 15.36e6
    s, _ = L.capture(mhz * 1e6, rate, args.secs, args.gain, args.low_if_mhz * 1e6)
    if s is None:
        res.update(status="capture failed", t=time.time() - t0)
        return res
    nudge, hits = probe(s, rate, args.nudges)
    if not hits:
        res.update(status="no LTE", t=time.time() - t0)
        return res
    mhz_true = round(mhz + nudge / 1e6, 1)
    res.update(mhz=mhz_true, band=band_of(mhz_true), earfcn=earfcn_of(mhz_true))
    hit = hits[0]
    pci = int(hit["pci"])
    res["pci"] = pci
    res["also_pci"] = [int(h["pci"]) for h in hits[1:]]
    s = shift(s, rate, nudge)
    cfo = float(hit.get("offset_hz", 0.0))
    # step 1: 50-PRB read at 15.36 Msps - the MIB (centre 6 RB) decodes whatever the cell width
    mib, ident, problem = srsran_read(s, rate, pci, cfo, 50, args.srsran_timeout)
    if mib:
        res["rb"], res["bw"] = mib.get("prb"), mib.get("mhz")
    prb = mib.get("prb") if mib else None
    if not ident and prb in SRS_RATE and prb != 50:
        want = SRS_RATE[prb]
        if want < rate:                       # 1.4/3/5 MHz: down-sample what we have
            s2, r2 = decimate(s, int(round(rate / want))), want
        else:                                 # 15/20 MHz: needs a wider capture
            print("    cell is %s MHz wide - re-capturing at %.2f Msps" % (mib.get("mhz"), want / 1e6), flush=True)
            s2, _ = L.capture(mhz_true * 1e6, want, args.secs, args.gain, args.low_if_mhz * 1e6)
            r2 = want
        if s2 is not None:
            mib2, ident, problem = srsran_read(s2, r2, pci, cfo, prb, args.srsran_timeout)
            mib = mib2 or mib
    if ident:
        res.update(ident)
        res["status"] = "OK"
    else:
        res["status"] = problem or "no SIB1"
    res["t"] = time.time() - t0
    return res


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
        t0 = time.time()
        found, runs = survey_band(lo, hi, args.gain, args.thresh_db, args.min_mhz)
        extra = [k for k in KNOWN.get(b, []) if all(abs(k - c) > 0.35 for c in found)]
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
    ap.add_argument("--nudges", type=float, nargs="+", default=[0, 0.1, -0.1, 0.2, -0.2, 0.3, -0.3],
                    help="centre corrections (MHz) tried by the sync probe")
    ap.add_argument("--srsran-timeout", type=int, default=20, help="seconds per pdsch_ue attempt")
    ap.add_argument("--loop", action="store_true", help="repeat until Ctrl-C")
    ap.add_argument("--rediscover", type=int, default=0, help="--loop: re-survey every N passes (0 = never)")
    ap.add_argument("--label", default=None, help="tag the CSV rows (e.g. home, field_site2)")
    ap.add_argument("--csv", default="~/field/lte_scan.csv")
    args = ap.parse_args()

    if not H._pdsch().is_file():
        sys.exit("srsRAN pdsch_ue not found at %s - build it first (HANDOFF 'LTE srsRAN')" % H._pdsch())

    cells = {}
    passes = 0
    t_all = time.time()
    try:
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
                r = decode_carrier(c, args)
                print(row(r), flush=True)
                log_csv(args.csv, r, args.label)
                if r.get("also_pci"):
                    print("      also heard PCI %s (weaker cells on the same carrier)"
                          % ", ".join(str(p) for p in r["also_pci"]), flush=True)
                if r.get("status") == "OK":
                    cells[(r["mhz"], r["cell_id"])] = r
                    got.append(r)
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
    print("\n===== %d unique LTE cell(s), %d pass(es), %.0fs =====" % (len(cells), passes, time.time() - t_all))
    if cells:
        print(HEAD.rsplit(" ", 1)[0])
        for k in sorted(cells):
            print(row(cells[k]).rsplit(" ", 1)[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
