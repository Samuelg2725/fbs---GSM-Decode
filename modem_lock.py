#!/usr/bin/env python3
"""Point the CALYX modem (ALA440) at exact EARFCNs: lock it to one carrier at a time and list every cell (PCI)
it measures there, with Cell IDs where the modem decodes them.

  sudo python3 modem_lock.py --earfcn 6300 6400 76
  sudo python3 modem_lock.py --earfcn 6300 --secs 60 --raw

For each EARFCN:
  1. radio off (AT+CFUN=0), lock to the EARFCN (AT^CELLLOCK=1,3,0,<earfcn>), radio on (AT+CFUN=1)
  2. for --secs seconds, every 2 s: the serving cell (AT^CELLINFO=1) and the cells it measures on that
     EARFCN (AT^CELLINFO=2, "intra" = same carrier) - PCI and RSRP
  3. a band scan (AT^NETSCAN) for the Cell IDs on that EARFCN
The manual says a lock only fully clears after a modem restart, so if the modem does not move to the next
EARFCN it is restarted (AT+RESET, about 30 s) and locked again. At the end the lock is removed and the modem
restarted, so it is back to normal (--keep-lock to skip that).

Slower than modem_scan.py (about 1 minute per EARFCN), but it lists every PCI the modem measures on that carrier,
not only the best cell per operator. Results: ~/field/modem_lock.csv; raw answers in ~/field/modem_scan_raw.txt.
"""
import argparse
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                                   # noqa: E402
import modem_scan as M                               # noqa: E402
from modem_probe import open_port, command          # noqa: E402

OUT = os.path.join(M.FIELD, "modem_lock.csv")


def log(text):
    with open(M.RAW, "a") as fh:
        fh.write("=== %s\n%s\n" % (C.utc_iso(), text.strip()))


def ask(fd, cmd, timeout=5.0, raw=False):
    reply = command(fd, cmd, timeout)
    log("%s\n%s" % (cmd, reply))
    if raw:
        print("   %s -> %s" % (cmd, " | ".join(l for l in reply.splitlines() if l.strip())), flush=True)
    return reply


def serving(text):
    """AT^CELLINFO=1 LTE line -> {pci, earfcn, cell_id, tac, plmn} or None (e.g. NO SERVICE)."""
    for line in text.splitlines():
        if "^CELLINFO:" not in line or "NO SERVICE" in line.upper():
            continue
        f = [x.strip().strip('"') for x in line.split(":", 1)[1].split(",")]
        if len(f) < 12 or "LTE" not in f[0].upper():
            continue
        try:
            return {"plmn": f[2] + f[3], "cell_id": M.num(f[4]), "pci": int(f[5]), "tac": M.num(f[8]),
                    "earfcn": int(f[11])}
        except (ValueError, IndexError):
            continue
    return None


def reconnect(port, wait=120):
    """After AT+RESET the USB ports vanish and come back: wait until the modem answers again."""
    t0 = time.time()
    time.sleep(10)
    while time.time() - t0 < wait:
        C.status("  modem restarting... %.0f s" % (time.time() - t0))
        if os.path.exists(port):
            try:
                fd = open_port(port)
                if "OK" in command(fd, "AT", 1.5):
                    C.clear_status()
                    command(fd, "ATE0", 2)
                    return fd
                os.close(fd)
            except OSError:
                pass
        time.sleep(2)
    C.clear_status()
    sys.exit("The modem did not come back after the restart. Unplug and replug it, then reload the driver.")


def restart(fd, port):
    command(fd, "AT+RESET", 5)
    try:
        os.close(fd)
    except OSError:
        pass
    return reconnect(port)


def lock(fd, earfcn, raw):
    ask(fd, "AT+CFUN=0", 15, raw)
    ask(fd, "AT^CELLLOCK=0", 5, raw)
    r = ask(fd, 'AT^CELLLOCK=1,3,0,%d' % earfcn, 5, raw)
    ask(fd, "AT+CFUN=1", 15, raw)
    return "OK" in r


def watch(fd, earfcn, secs, raw):
    """Poll serving + neighbour info; -> (cells {pci: {...}}, on_earfcn)."""
    cells, on, t0 = {}, False, time.time()
    while time.time() - t0 < secs:
        C.status("  EARFCN %d: listening %.0f/%.0f s, %d PCIs so far" % (earfcn, time.time() - t0, secs, len(cells)))
        s = serving(ask(fd, "AT^CELLINFO=1", 5, raw))
        if s and s["earfcn"] == earfcn:
            on = True
            c = cells.setdefault(s["pci"], {"pci": s["pci"], "rsrp": None})
            c.update(cell_id=s["cell_id"], tac=s["tac"], plmn=s["plmn"], how="serving cell")
        text = ask(fd, "AT^CELLINFO=2", 5, raw)
        for m in M.NEIGH_RE.finditer(text):
            if int(m.group(2)) != earfcn or int(m.group(3)) > 503:
                continue
            if m.group(1) == "intra":
                on = True
            pci, rsrp = int(m.group(3)), float(m.group(5))
            c = cells.setdefault(pci, {"pci": pci, "rsrp": rsrp, "how": "measured"})
            if c.get("rsrp") is None or rsrp > c["rsrp"]:
                c["rsrp"] = rsrp
        time.sleep(2)
    C.clear_status()
    return cells, on


def main():
    ap = argparse.ArgumentParser(description="Lock the CALYX modem to each EARFCN and list the cells there.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--earfcn", type=int, nargs="+", required=True, help="LTE carriers to visit")
    ap.add_argument("--secs", type=float, default=30, help="seconds to listen on each EARFCN (default 30)")
    ap.add_argument("--port", help="modem AT port (default: find it)")
    ap.add_argument("--keep-lock", action="store_true", help="leave the modem locked at the end")
    ap.add_argument("--no-netscan", action="store_true", help="skip the band scan for Cell IDs")
    ap.add_argument("--raw", action="store_true", help="print every AT command and answer")
    args = ap.parse_args()
    if os.geteuid() != 0:
        sys.exit("Run with sudo (the modem port needs it).")
    os.makedirs(M.FIELD, exist_ok=True)
    port = M.find_port(args.port)
    fd = open_port(port)
    command(fd, "ATE0", 2)
    C.info("modem on %s" % port)
    results = {}
    try:
        for earfcn in args.earfcn:
            C.heading("EARFCN %d (band %s)" % (earfcn, M.band_of(earfcn) or "?"))
            if not lock(fd, earfcn, args.raw):
                C.warn("  the modem refused the lock (see --raw)")
            time.sleep(3)
            cells, on = watch(fd, earfcn, args.secs, args.raw)
            if not on:
                C.info("  modem not on EARFCN %d yet - restarting it and locking again" % earfcn)
                fd = restart(fd, port)
                lock(fd, earfcn, args.raw)
                time.sleep(3)
                more, on = watch(fd, earfcn, args.secs, args.raw)
                cells.update(more)
            if not args.no_netscan and M.band_of(earfcn):
                lines = M.scan(fd, M.RAT["lte"], M.band_of(earfcn), 300)
                for c in M.parse(lines):
                    if c["freq"] != earfcn:
                        continue
                    x = cells.setdefault(c["pci"], {"pci": c["pci"], "rsrp": c["level"], "how": "band scan"})
                    x.update(cell_id=c["cell_id"], tac=c["tac"], plmn=c["plmn"])
            results[earfcn] = (cells, on)
            print("  %s: %s" % ("on the carrier" if on else C.c("never reached this carrier", "red"),
                                " ".join("PCI %d" % p for p in sorted(cells)) or "no cells"), flush=True)
    except KeyboardInterrupt:
        print("")
        C.clear_status()
    finally:
        if not args.keep_lock:
            C.info("removing the lock and restarting the modem (back to normal)...")
            try:
                command(fd, "AT^CELLLOCK=0", 5)
                fd = restart(fd, port)
            except OSError:
                pass
        try:
            os.close(fd)
        except OSError:
            pass
    t = C.Table(["EARFCN", "PCI", "RSRP", "Operator", "TAC", "Cell ID", "From"],
                align={k: "right" for k in ("EARFCN", "PCI", "RSRP", "TAC", "Cell ID")})
    new = not os.path.exists(OUT)
    with open(OUT, "a") as fh:
        if new:
            fh.write("utc,earfcn,pci,rsrp,plmn,tac,cell_id,source\n")
        for earfcn, (cells, on) in results.items():
            if not cells:
                t.add([earfcn, "-", "-", "-", "-", "-", ("no cells" if on else "not reached", "red")])
            for c in sorted(cells.values(), key=lambda c: -(c.get("rsrp") or -999)):
                t.add([earfcn, c["pci"], "%g" % c["rsrp"] if c.get("rsrp") is not None else "-",
                       C.op_cell(M.operator(c.get("plmn", ""))) if c.get("plmn") else "-", c.get("tac", "-"),
                       (c["cell_id"], "green+bold") if c.get("cell_id") else "-", c.get("how", "")])
                fh.write("%s,%d,%d,%s,%s,%s,%s,%s\n" % (C.utc_iso(), earfcn, c["pci"], c.get("rsrp", ""),
                                                        c.get("plmn", ""), c.get("tac", ""), c.get("cell_id", ""),
                                                        c.get("how", "")))
    M.give_back(OUT)
    M.give_back(M.RAW)
    t.show("Cells the modem measured on each EARFCN")
    C.info("results added to %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
