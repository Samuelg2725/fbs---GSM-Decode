#!/usr/bin/env python3
"""How long the CALYX modem (ALA440) takes to scan each technology and band, and what it finds - with timings.

  sudo python3 modem_timing.py                 # 2G, then 4G band by band, then all-in-one (no 3G scan)
  sudo python3 modem_timing.py --rat 4g        # only 4G (LTE) bands
  sudo python3 modem_timing.py --no-all        # skip the final all-technologies scan

For each band it sends AT^NETSCAN=<technology>,<band> and records:
  * the time from sending the command to the first cell, to each cell, and to the end of the scan
  * every cell found: operator, channel (ARFCN / UARFCN / EARFCN), BSIC / PSC / PCI, LAC or TAC, Cell ID, level
Each line's arrival time is measured as the modem sends it (so if the modem reports operators one by one, the
timings show it; if it sends everything at the end, all cells show the same time).
Then a summary table: band by band, seconds and cells found, and the totals per technology. Finally (unless
--no-all) one AT^NETSCAN=3 scan of everything at once, for comparison.
Results: ~/field/modem_timing.csv (one row per cell, with its timings). Raw answers: ~/field/modem_scan_raw.txt.
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
# (technology, NETSCAN mode, [(band number, label)])  - UK bands; the modem says if one is not supported
PLAN = [
    # this firmware refuses a band for 2G (error 50): 2G is scanned whole. 3G is not scanned: UK 3G is
    # switched off and the modem only gives up after 60 s
    ("2G", 0, [(None, "GSM (900 + 1800)")]),
    ("4G", 2, [(20, "LTE 800 (band 20)"), (28, "LTE 700 (band 28)"), (8, "LTE 900 (band 8)"),
               (32, "LTE 1500 (band 32)"), (3, "LTE 1800 (band 3)"), (1, "LTE 2100 (band 1)"),
               (7, "LTE 2600 (band 7)"), (38, "LTE 2600 TDD (band 38)"), (40, "LTE 2300 TDD (band 40)"),
               (None, "LTE, all bands in one scan")]),
]
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
    ap = argparse.ArgumentParser(description="Time the CALYX modem's scan of each technology and band.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--rat", choices=["2g", "4g"], nargs="+", help="only these technologies")
    ap.add_argument("--no-all", action="store_true", help="skip the final scan of everything at once")
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
    C.info("modem on %s - scanning band by band, timings from the moment each scan command is sent" % port)
    summary = C.Table(["Technology", "Band", "Seconds", "First cell at", "Cells", "Operators", "Result"],
                      align={"Seconds": "right", "First cell at": "right", "Cells": "right"})
    rows, per_tech, t_start = [], {}, time.time()
    plan = [p for p in PLAN if not args.rat or p[0].lower() in args.rat]
    try:
        for tech, mode, bands in plan:
            C.heading("%s scan" % tech)
            for band, label in bands:
                lines, total, err = timed_scan(fd, mode, band, args.timeout)
                cells = [c for c in cells_from(lines) if c["tech"] == tech]
                per_tech.setdefault(tech, [0.0, 0])
                per_tech[tech][0] += total
                per_tech[tech][1] += len(cells)
                if err:
                    why = ("band not supported by the modem" if "50" in err else
                           "nothing found (the modem gave up after %.0f s)" % total if total > 30 else err)
                    print("  %-24s %6.1f s   %s" % (label, total, C.c(why, "grey")), flush=True)
                    summary.add([tech, label, "%.1f" % total, "-", 0, "-", (why, "grey")])
                    continue
                print("  %-24s %6.1f s   %s" % (label, total, C.c("%d cell%s" % (len(cells), "" if len(cells) == 1
                                                                                  else "s"),
                                                                  "green" if cells else "yellow")), flush=True)
                for c in cells:
                    op, chan, ident, level = describe(c)
                    print("      at %5.1f s  %-9s %s  %s" % (c["t"], C.c(*C.op_cell(op)), ident, level), flush=True)
                    rows.append([C.utc_iso(), tech, label, "%.1f" % c["t"], "%.1f" % total, op, chan, ident, level])
                ops = sorted({describe(c)[0] for c in cells})
                summary.add([tech, label, "%.1f" % total, "%.1f" % min(c["t"] for c in cells) if cells else "-",
                             len(cells), ", ".join(ops) or "-", ("found cells", "green") if cells else
                             ("nothing found", "yellow")])
        if not args.no_all:
            C.heading("Everything at once (AT^NETSCAN=3, all technologies, all bands)")
            lines, total, err = timed_scan(fd, 3, None, args.timeout)
            cells = cells_from(lines)
            print("  all technologies       %6.1f s   %s" % (total, err or "%d cells" % len(cells)), flush=True)
            for c in cells:
                op, chan, ident, level = describe(c)
                print("      at %5.1f s  %s  %-9s %s  %s" % (c["t"], c["tech"], C.c(*C.op_cell(op)), ident, level),
                      flush=True)
            summary.add(["all", "everything at once", "%.1f" % total,
                         "%.1f" % min(c["t"] for c in cells) if cells else "-", len(cells),
                         ", ".join(sorted({describe(c)[0] for c in cells})) or "-",
                         (err, "red") if err else ("found cells", "green") if cells else ("nothing found", "yellow")])
    except KeyboardInterrupt:
        print("")
        C.clear_status()
    finally:
        os.close(fd)
    summary.show("Modem scan timings (band by band)")
    C.info("Totals: " + ",  ".join("%s %.1f s, %d cells" % (t, v[0], v[1]) for t, v in per_tech.items())
           + "   |   whole run %.1f s" % (time.time() - t_start))
    new = not os.path.exists(OUT)
    with open(OUT, "a") as fh:
        if new:
            fh.write("utc,technology,band,cell_found_at_s,band_scan_s,operator,channel,identity,level\n")
        for r in rows:
            fh.write(",".join('"%s"' % str(x).replace('"', "'") for x in r) + "\n")
    M.give_back(OUT)
    M.give_back(M.RAW)
    C.info("results added to %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
