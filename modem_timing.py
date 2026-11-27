#!/usr/bin/env python3
"""Full 2G / 3G / 4G scan with the CALYX modem (ALA440), repeated for accuracy, with timings.

  sudo python3 modem_timing.py                 # 2G, 3G and 4G, each a full scan, 3 rounds
  sudo python3 modem_timing.py --rounds 5      # more rounds: weak cells that come and go are caught
  sudo python3 modem_timing.py --rat 2g 4g     # only some technologies
  sudo python3 modem_timing.py --per-band      # 4G also band by band (fast but the modem checks less)

The modem has no "accuracy" setting: each AT^NETSCAN is one snapshot, and a cell near the noise (about -100 dBm)
is in some snapshots and not in others. So every technology is scanned in full (all its bands in one command,
which found the most), and the whole set is repeated --rounds times. Every cell found in any round is kept,
with how many rounds saw it. A cell seen in every round is solid; one seen once is weak or intermittent.

Timings: each scan is timed from the moment the command is sent until the modem answers. The modem sends all
its results at the end of a scan, so the time is per scan (per technology), not per cell; "first found" is the
time into the whole run at which a cell first appeared. Results: ~/field/modem_timing.csv (one row per cell
per round). Raw answers: ~/field/modem_scan_raw.txt.
"""
import argparse
import os
import re
import select
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                                   # noqa: E402
import modem_scan as M                               # noqa: E402
from modem_probe import open_port, command          # noqa: E402

OUT = os.path.join(M.FIELD, "modem_timing.csv")
FULL = [("2G", 0, "GSM, all bands (900 + 1800)"), ("3G", 1, "UMTS, all bands"), ("4G", 2, "LTE, all bands")]
LTE_BANDS = [(20, "LTE 800 (band 20)"), (28, "LTE 700 (band 28)"), (8, "LTE 900 (band 8)"),
             (3, "LTE 1800 (band 3)"), (1, "LTE 2100 (band 1)"), (7, "LTE 2600 (band 7)"),
             (38, "LTE 2600 TDD (band 38)"), (40, "LTE 2300 TDD (band 40)")]

CELL_RE = re.compile(r'^\s*\d+\s*,\s*"?(2G|3G|4G)"?\s*,(.*)$', re.I)
PAIR_RE = re.compile(r"([A-Za-z][A-Za-z ]*?)\s*:\s*(-?[0-9A-Fa-fx]+)")
OP_RE = re.compile(r'\^NETSCAN:\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,\s*"([^"]*)"')


def num(v):
    try:
        return int(v)
    except ValueError:
        try:
            return int(v, 16)
        except ValueError:
            return v


def timed_scan(fd, mode, band, timeout):
    """Send AT^NETSCAN and read its answer line by line -> (lines [(seconds, text)], total seconds, error)."""
    cmd = "AT^NETSCAN=%d" % mode + (",%d" % band if band is not None else "")
    os.write(fd, (cmd + "\r").encode())
    t0 = time.time()
    buf, lines, err = "", [], None
    while time.time() - t0 < timeout:
        C.status("    %s ... %.1f s" % (cmd, time.time() - t0))
        r, _, _ = select.select([fd], [], [], 0.2)
        if not r:
            continue
        try:
            chunk = os.read(fd, 65536).decode(errors="replace")
        except BlockingIOError:
            continue
        now = time.time() - t0
        buf += chunk
        while "\n" in buf:
            line, buf = buf.split("\n", 1)
            line = line.strip()
            if line:
                lines.append((now, line))
        done = [l for _, l in lines if l == "OK" or "ERROR" in l]
        if done:
            if "ERROR" in done[0]:
                err = done[0]
            break
    else:
        err = "no answer within %d s" % timeout
    total = time.time() - t0
    C.clear_status()
    with open(M.RAW, "a") as fh:
        fh.write("=== %s %s (%.1f s)\n%s\n" % (C.utc_iso(), cmd, total, "\n".join("%6.1f %s" % x for x in lines)))
    return lines, total, err


def cells_from(lines):
    """-> [{tech, t, plmn, fields...}] from the timed NETSCAN lines."""
    out, plmn = [], ""
    for t, line in lines:
        m = OP_RE.search(line)
        if m:
            plmn = m.group(3) or m.group(2)
            continue
        m = CELL_RE.match(line)
        if not m:
            continue
        f = {k.strip().lower().replace(" ", ""): num(v) for k, v in PAIR_RE.findall(m.group(2))}
        out.append({"tech": m.group(1).upper(), "t": t, "plmn": plmn, "f": f})
    return out


def describe(c):
    f = c["f"]
    op = M.operator(c["plmn"]) or c["plmn"] or "?"
    chan = f.get("freq", f.get("arfcn", f.get("uarfcn", "?")))
    if c["tech"] == "4G":
        ident = "EARFCN %s  PCI %s  TAC %s  Cell ID %s" % (chan, f.get("pci", "?"), f.get("tac", "?"),
                                                           f.get("cellid", "?"))
        level = "RSRP %s" % f.get("rsrp", "?")
    elif c["tech"] == "3G":
        ident = "UARFCN %s  PSC %s  LAC %s  Cell ID %s" % (chan, f.get("psc", "?"), f.get("lac", "?"),
                                                           f.get("cellid", "?"))
        level = "RSCP %s" % f.get("rscp", f.get("rssi", "?"))
    else:
        ident = "ARFCN %s  BSIC %s  LAC %s  Cell ID %s" % (chan, f.get("bsic", "?"), f.get("lac", "?"),
                                                           f.get("cellid", "?"))
        level = "RSSI %s" % f.get("rssi", "?")
    return op, chan, ident, level


def main():
    ap = argparse.ArgumentParser(description="Full 2G/3G/4G modem scan, repeated for accuracy, with timings.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--rat", choices=["2g", "3g", "4g"], nargs="+", help="only these technologies")
    ap.add_argument("--rounds", type=int, default=3, help="full scans of each technology (default 3)")
    ap.add_argument("--per-band", action="store_true", help="4G also band by band in each round")
    ap.add_argument("--timeout", type=float, default=300, help="seconds per scan before giving up (default 300)")
    ap.add_argument("--port", help="modem AT port (default: find it)")
    args = ap.parse_args()
    if os.geteuid() != 0:
        sys.exit("Run with sudo (the modem port needs it).")
    os.makedirs(M.FIELD, exist_ok=True)
    port = M.find_port(args.port)
    fd = open_port(port)
    command(fd, "ATE0", 2)
    if "+CFUN: 1" not in command(fd, "AT+CFUN?", 3):
        command(fd, "AT+CFUN=1", 15)
        time.sleep(3)
    plan = []
    for tech, mode, label in FULL:
        if args.rat and tech.lower() not in args.rat:
            continue
        plan.append((tech, mode, None, label))
        if tech == "4G" and args.per_band:
            plan += [("4G", 2, b, l) for b, l in LTE_BANDS]
    C.info("modem on %s - %d round%s of: %s" % (port, args.rounds, "" if args.rounds == 1 else "s",
                                               ", ".join(p[3] for p in plan)))
    seen = {}                                            # (tech, identity) -> cell record
    times = {}                                           # scan label -> [seconds per round]
    rows, t_start, done_rounds = [], time.time(), 0
    try:
        for rnd in range(1, args.rounds + 1):
            C.heading("Round %d of %d  (%.0f s into the run)" % (rnd, args.rounds, time.time() - t_start))
            for tech, mode, band, label in plan:
                lines, total, err = timed_scan(fd, mode, band, args.timeout)
                times.setdefault(label, []).append(total)
                cells = [c for c in cells_from(lines) if c["tech"] == tech]
                if err:
                    why = ("band not supported by the modem" if "50" in err else
                           "nothing found (the modem gave up after %.0f s)" % total if total > 30 else err)
                    print("  %-28s %6.1f s   %s" % (label, total, C.c(why, "grey")), flush=True)
                    continue
                new, new_keys = 0, set()
                for c in cells:
                    op, chan, ident, level = describe(c)
                    key = (tech, ident)
                    rec = seen.get(key)
                    if rec is None:
                        new += 1
                        new_keys.add(key)
                        rec = seen[key] = {"tech": tech, "op": op, "ident": ident, "level": level, "rounds": set(),
                                           "first": time.time() - t_start, "scan": label}
                    rec["rounds"].add(rnd)
                    rec["level"] = level
                    rows.append([C.utc_iso(), rnd, tech, label, "%.1f" % total, op, chan, ident, level])
                print("  %-28s %6.1f s   %s%s" % (
                    label, total, C.c("%d cell%s" % (len(cells), "" if len(cells) == 1 else "s"),
                                      "green" if cells else "yellow"),
                    C.c("  (%d new)" % new, "cyan") if new and rnd > 1 else ""), flush=True)
                for c in cells:
                    op, chan, ident, level = describe(c)
                    mark = C.c("NEW ", "cyan") if rnd > 1 and (tech, ident) in new_keys else "    "
                    print("      %s%-9s %s  %s" % (mark, C.c(*C.op_cell(op)), ident, level), flush=True)
            done_rounds = rnd
    except KeyboardInterrupt:
        print("")
        C.clear_status()
    finally:
        os.close(fd)
    n = max(1, done_rounds)
    t = C.Table(["Tech", "Operator", "Cell", "Last level", "Seen in", "First found"],
                align={"First found": "right"})
    for rec in sorted(seen.values(), key=lambda r: (r["tech"], r["op"], r["ident"])):
        k = len(rec["rounds"])
        t.add([rec["tech"], C.op_cell(rec["op"]), rec["ident"], rec["level"],
               ("%d of %d rounds" % (k, n), "green" if k == n else ("yellow" if k > 1 else "red")),
               "%.0f s" % rec["first"]])
    t.show("Every cell the modem found (%d), all rounds together" % len(seen))
    s = C.Table(["Scan", "Rounds", "Average s", "Fastest s", "Slowest s", "Cells (all rounds)"],
                align={k: "right" for k in ("Rounds", "Average s", "Fastest s", "Slowest s", "Cells (all rounds)")})
    for tech, mode, band, label in plan:
        v = times.get(label, [])
        if not v:
            continue
        cells = sum(1 for r in seen.values() if r["scan"] == label or (band is None and r["tech"] == tech))
        s.add([label, len(v), "%.1f" % (sum(v) / len(v)), "%.1f" % min(v), "%.1f" % max(v), cells])
    s.show("Scan timings")
    C.info("whole run %.0f s; per round about %.0f s" % (time.time() - t_start,
                                                          (time.time() - t_start) / n))
    new_file = not os.path.exists(OUT)
    with open(OUT, "a") as fh:
        if new_file:
            fh.write("utc,round,technology,scan,scan_s,operator,channel,identity,level\n")
        for r in rows:
            fh.write(",".join('"%s"' % str(x).replace('"', "'") for x in r) + "\n")
    M.give_back(OUT)
    M.give_back(M.RAW)
    C.info("results added to %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
