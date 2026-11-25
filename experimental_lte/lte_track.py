#!/usr/bin/env python3
"""Live LTE tracking on one carrier: srsRAN pdsch_ue streams from the BladeRF and decodes every SIB1.

The LTE equivalent of gsm_arfcn_fast.py tracking. SIB1 (Cell ID, TAC, PLMN) is broadcast every 20 ms,
so once locked this can report the identity many times a second.

  python3 ~/fbs-upload/lte_track.py --earfcn 6300                 # Vodafone 806 at home
  python3 ~/fbs-upload/lte_track.py --freq-mhz 816 --label o2_816
  python3 ~/fbs-upload/lte_track.py --earfcn 6300 --secs 60       # stop by itself after 60 s

Prints the time to the first decode, then one line per second (decodes in that second + identity),
and a summary on exit (also a row in --summary-csv). pdsch_ue locks to the STRONGEST cell on the carrier.
Needs ~/srsRAN_4G built with srsran_watch/print_sib1.py (see HANDOFF 'LTE srsRAN').
"""
import argparse
import os
import re
import select
import signal
import subprocess
import sys
import time
from pathlib import Path

SIB1_RE = re.compile(r"SIB1 cell identity (\d+)\s+eNB (\d+)\s+sector (\d+)\s+TAC (\d+)((?:[ \t]+PLMN \d+-\d+)*)[ \t]*\r?\n")
OPERATOR = {"234-10": "O2", "234-15": "Vodafone", "234-20": "Three", "234-30": "EE", "234-33": "EE"}
BANDS = {1: (2110.0, 0, 599), 3: (1805.0, 1200, 1949), 7: (2620.0, 2750, 3449), 8: (925.0, 3450, 3799),
         20: (791.0, 6150, 6449)}


def earfcn_mhz(e):
    for lo, off, hi in BANDS.values():
        if off <= e <= hi:
            return lo + (e - off) / 10.0
    return None


def find_pdsch():
    for home in (Path.home(), Path("/home/revector")):
        p = home / "srsRAN_4G" / "build" / "lib" / "examples" / "pdsch_ue"
        if p.is_file():
            return p
    return None


def track_capture(args, mhz):
    """Lock once (PCI, frequency offset, bandwidth), then loop: short capture -> align -> pdsch_ue,
    stopped at the first SIB1. Same proven path as lte_scan.py, minus the per-capture sync search."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import lte_scan as S
    import lte_decode as L
    S.EARLY[0] = True
    lif = args.low_if_mhz * 1e6
    print("LTE track %.1f MHz%s  gain %g  capture mode (%.2fs per decode attempt)" % (
        mhz, (" EARFCN %d" % args.earfcn) if args.earfcn is not None else "", args.gain, args.cap), flush=True)
    t_start = time.time()
    lock = None
    stats = {"ok": 0, "fail": 0, "first": None, "longest": 0.0, "last": None, "relock": 0}
    idents = {}
    cycle_times = []
    fails_in_row = 0
    try:
        while not (args.secs and time.time() - t_start >= args.secs):
            if lock is None:
                s, _ = L.capture(mhz * 1e6, 15.36e6, 0.2, args.gain, lif)
                if s is None:
                    continue
                nudge, hits = S.probe(s, 15.36e6, [0, 0.1, -0.1])
                if not hits:
                    print("%s  no LTE sync at %.1f MHz - retrying" % (time.strftime("%H:%M:%S"), mhz), flush=True)
                    continue
                pci, cfo = int(hits[0]["pci"]), float(hits[0].get("offset_hz", 0.0))
                mib, _, problem = S.srsran_read(S.shift(s, 15.36e6, nudge), 15.36e6, pci, cfo, 50, args.timeout)
                if not mib:
                    print("%s  PCI %d found, MIB failed (%s) - retrying" % (time.strftime("%H:%M:%S"), pci, problem),
                          flush=True)
                    continue
                prb = mib.get("prb", 50)
                rate = S.SRS_RATE.get(prb, 15.36e6)
                cap_rate = max(rate, 15.36e6)          # 5 MHz cells: capture 15.36 and decimate
                lock = (pci, cfo, nudge, prb, rate, cap_rate)
                print("%s  locked: PCI %d  %s MHz (%d RB)  offset %+.1f kHz  -> decoding every capture" % (
                    time.strftime("%H:%M:%S"), pci, mib.get("mhz"), prb, cfo / 1000.0), flush=True)
                continue
            pci, cfo, nudge, prb, rate, cap_rate = lock
            t0 = time.time()
            s, _ = L.capture(mhz * 1e6, cap_rate, args.cap, args.gain, lif)
            if s is None:
                continue
            s = S.shift(s, cap_rate, nudge)
            if rate < cap_rate:
                s = S.decimate(s, int(round(cap_rate / rate)))
            _, ident, problem = S.srsran_read(s, rate, pci, cfo, prb, args.timeout)
            now = time.time()
            took = now - t0
            cycle_times.append(took)
            if ident:
                fails_in_row = 0
                stats["ok"] += 1
                line = "CellID %s (eNB %s sector %s)  TAC %s  PLMN %s-%s %s  PCI %d" % (
                    ident["cell_id"], ident["enb_id"], ident["sector_id"], ident["tac"], ident["mcc"],
                    ident["mnc_digits"], ident.get("operator", ""), pci)
                idents[line] = idents.get(line, 0) + 1
                if stats["first"] is None:
                    stats["first"] = now - t_start
                if stats["last"] is not None:
                    stats["longest"] = max(stats["longest"], now - stats["last"])
                gap = (now - stats["last"]) if stats["last"] is not None else (now - t_start)
                stats["last"] = now
                print("%s  %5.1fs  %s" % (time.strftime("%H:%M:%S"), gap, line), flush=True)
            else:
                stats["fail"] += 1
                fails_in_row += 1
                print("%s  miss (%s, %.1fs)" % (time.strftime("%H:%M:%S"), problem, took), flush=True)
                if fails_in_row >= 3:
                    print("  3 misses in a row - re-locking", flush=True)
                    lock = None
                    stats["relock"] += 1
                    fails_in_row = 0
    except KeyboardInterrupt:
        pass
    run = time.time() - t_start
    att = stats["ok"] + stats["fail"]
    avg = sum(cycle_times) / len(cycle_times) if cycle_times else 0.0
    print("\n===== LTE tracking summary (capture mode) =====")
    print("[%s] %.0fs  first decode after %s  decoded %d of %d attempts (%.0f%%)  avg %.1fs per attempt  "
          "longest gap %.1fs  re-locks %d" % (
              args.label, run, "%.1fs" % stats["first"] if stats["first"] is not None else "never",
              stats["ok"], att, 100.0 * stats["ok"] / max(1, att), avg, stats["longest"], stats["relock"]))
    for ident, n in idents.items():
        print("  %5d x %s" % (n, ident))
    print("time split: " + "  ".join("%s %.1fs" % kv for kv in S.T.items()))
    path = os.path.expanduser(args.summary_csv)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new:
            fh.write("start_utc,label,mhz,gain,runtime_s,first_decode_s,sib1_decoded,per_second,longest_gap_s,identities\n")
        fh.write("%s,%s,%.1f,%g,%.1f,%s,%d,%.2f,%.2f,%s\n" % (
            time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t_start)), (args.label or "") + "_capture", mhz,
            args.gain, run, "%.2f" % stats["first"] if stats["first"] is not None else "", stats["ok"],
            stats["ok"] / run if run else 0.0, stats["longest"], " | ".join(idents).replace(",", ";")))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--earfcn", type=int)
    ap.add_argument("--freq-mhz", type=float)
    ap.add_argument("--gain", type=float, default=40.0)
    ap.add_argument("--secs", type=float, default=0, help="stop after this many seconds (0 = until Ctrl-C)")
    ap.add_argument("--label", default=None)
    ap.add_argument("--summary-csv", default="~/field/lte_track_summary.csv")
    ap.add_argument("--extra", default="", help="--live: extra pdsch_ue arguments")
    ap.add_argument("--live", action="store_true",
                    help="use srsRAN's live receiver (failed at home: never locks, DC spike on centre)")
    ap.add_argument("--cap", type=float, default=0.1,
                    help="capture mode: seconds per capture once locked (default 0.1 = 5 SIB1 chances)")
    ap.add_argument("--low-if-mhz", type=float, default=0.6)
    ap.add_argument("--timeout", type=int, default=10, help="capture mode: seconds per pdsch_ue run")
    args = ap.parse_args()

    mhz = args.freq_mhz if args.freq_mhz else (earfcn_mhz(args.earfcn) if args.earfcn is not None else None)
    if mhz is None:
        sys.exit("give --earfcn or --freq-mhz")
    if not args.live:
        return track_capture(args, mhz)
    exe = find_pdsch()
    if exe is None:
        sys.exit("pdsch_ue not found under ~/srsRAN_4G/build - see HANDOFF 'LTE srsRAN'")
    env = dict(os.environ)
    rf = exe.parents[2] / "lib" / "src" / "phy" / "rf"
    env["LD_LIBRARY_PATH"] = "/usr/local/lib:%s%s" % (
        rf, (":" + env["LD_LIBRARY_PATH"]) if env.get("LD_LIBRARY_PATH") else "")
    cmd = [str(exe), "-f", "%d" % int(round(mhz * 1e6)), "-g", "%g" % args.gain] + args.extra.split()
    print("LTE track %.1f MHz%s  gain %g  (%s)" % (
        mhz, (" EARFCN %d" % args.earfcn) if args.earfcn is not None else "", args.gain, " ".join(cmd[1:])),
        flush=True)

    t_start = time.time()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env, start_new_session=True)
    buf = ""
    total = 0
    first = None
    last_t = None
    longest = 0.0
    second_count = 0
    second_mark = None
    idents = {}
    current = None
    tail = []
    try:
        while True:
            if args.secs and time.time() - t_start >= args.secs:
                break
            ready, _, _ = select.select([p.stdout], [], [], 0.2)
            now = time.time()
            if ready:
                chunk = os.read(p.stdout.fileno(), 65536)
                if not chunk:
                    break
                text = chunk.decode(errors="replace")
                keep = "".join(l for l in text.splitlines(True) if "Requested buffer size" not in l)
                if keep:
                    tail = (tail + [keep])[-40:]
                buf += text
                while True:
                    m = SIB1_RE.search(buf)
                    if not m:
                        buf = buf[-400:]
                        break
                    buf = buf[m.end():]
                    plmns = re.findall(r"PLMN (\d+)-(\d+)", m.group(5) or "")
                    plmn = "%s-%s" % plmns[0] if plmns else "?"
                    ident = "CellID %s (eNB %s sector %s)  TAC %s  PLMN %s %s" % (
                        m.group(1), m.group(2), m.group(3), m.group(4), plmn, OPERATOR.get(plmn, ""))
                    total += 1
                    second_count += 1
                    idents[ident] = idents.get(ident, 0) + 1
                    if last_t is not None:
                        longest = max(longest, now - last_t)
                    last_t = now
                    if first is None:
                        first = now - t_start
                        second_mark = now
                        print("%s  first decode after %.1fs: %s" % (time.strftime("%H:%M:%S"), first, ident), flush=True)
                        current = ident
                    elif ident != current:
                        print("%s  IDENTITY CHANGED: %s" % (time.strftime("%H:%M:%S"), ident), flush=True)
                        current = ident
            elif p.poll() is not None:
                break
            if second_mark is not None and now - second_mark >= 1.0:
                print("%s  %3d SIB1 in last %.1fs   total %d   %s" % (
                    time.strftime("%H:%M:%S"), second_count, now - second_mark, total, current), flush=True)
                second_count = 0
                second_mark = now
            if first is None and now - t_start > 30 and int(now - t_start) % 30 == 0:
                print("  no SIB1 yet after %.0fs (no lock?) last output:\n%s" % (
                    now - t_start, "".join(tail)[-500:]), flush=True)
                time.sleep(1)
    except KeyboardInterrupt:
        pass
    finally:
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        if p.poll() is None:
            try:
                os.killpg(p.pid, signal.SIGINT)
                p.wait(timeout=3)
            except Exception:
                try:
                    os.killpg(p.pid, 9)
                except OSError:
                    pass
    run = time.time() - t_start
    locked = (run - first) if first is not None else 0.0
    rate = total / locked if locked > 0 else 0.0
    print("\n===== LTE tracking summary =====")
    print("[%s] %.0fs  first decode after %s  SIB1 decoded %d (%.1f per second after lock)  longest gap %.2fs" % (
        args.label, run, "%.1fs" % first if first is not None else "never", total, rate, longest))
    for ident, n in idents.items():
        print("  %5d x %s" % (n, ident))
    if first is None:
        print("no SIB1 decoded. Last pdsch_ue output:\n" + "".join(tail)[-1500:])
    path = os.path.expanduser(args.summary_csv)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new:
            fh.write("start_utc,label,mhz,gain,runtime_s,first_decode_s,sib1_decoded,per_second,longest_gap_s,identities\n")
        fh.write("%s,%s,%.1f,%g,%.1f,%s,%d,%.2f,%.2f,%s\n" % (
            time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t_start)), args.label or "", mhz, args.gain, run,
            "%.2f" % first if first is not None else "", total, rate, longest,
            " | ".join(idents).replace(",", ";")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
