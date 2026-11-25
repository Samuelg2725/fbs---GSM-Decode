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


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--earfcn", type=int)
    ap.add_argument("--freq-mhz", type=float)
    ap.add_argument("--gain", type=float, default=40.0)
    ap.add_argument("--secs", type=float, default=0, help="stop after this many seconds (0 = until Ctrl-C)")
    ap.add_argument("--label", default=None)
    ap.add_argument("--summary-csv", default="~/field/lte_track_summary.csv")
    ap.add_argument("--extra", default="", help="extra pdsch_ue arguments, e.g. \"-a device_name=bladerf\"")
    args = ap.parse_args()

    mhz = args.freq_mhz if args.freq_mhz else (earfcn_mhz(args.earfcn) if args.earfcn is not None else None)
    if mhz is None:
        sys.exit("give --earfcn or --freq-mhz")
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
                tail = (tail + [text])[-20:]
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
