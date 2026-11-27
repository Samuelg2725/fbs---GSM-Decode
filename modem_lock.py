#!/usr/bin/env python3
"""Point the CALYX modem (ALA440) at exact carriers or exact cells and read the Cell ID it finds there.

  sudo python3 modem_lock.py --earfcn 6300 76              # lock to each carrier: which cell does it camp on?
  sudo python3 modem_lock.py --cells 6300:21,94 76:229     # lock to each exact cell (EARFCN:PCI): its Cell ID
  sudo python3 modem_lock.py --cells 6300:21 --raw         # show every AT command and answer

For each target:
  1. radio off (AT+CFUN=0), remove any lock, lock (AT^CELLLOCK: to the EARFCN, or to the EARFCN + PCI), radio on
  2. wait until the modem camps on a cell. Without a SIM it has "limited service" only, but it still reports the
     cell it camps on in its registration status (+CEREG: TAC and Cell ID) - that is what is read here
  3. a band scan (AT^NETSCAN) of the target's band, to match the Cell ID to its PCI
With --cells the result is a check of each cell's identity by a second, independent receiver: compare the
Cell ID column with what lte_decode.py decoded for the same EARFCN:PCI (shown in the SDR column).

Without a SIM the modem does not measure other cells, so it reports one cell per lock: use --cells to visit each
PCI the SDR found. If the modem does not camp within --secs the target is reported as not camped (the cell is
not there or too weak; --restart-retry restarts the modem and tries once more). At the end the lock is removed and the modem restarted, so it is back to normal (--keep-lock skips
that). Results: ~/field/modem_lock.csv; raw answers in ~/field/modem_scan_raw.txt.
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
# +CEREG: [<n>,]<stat>,"<tac>","<ci>"[,<act>]  (TAC and Cell ID in hex) - as a reply to AT+CEREG? or unsolicited
CEREG_RE = re.compile(r'\+CEREG:\s*(?:\d+\s*,\s*)?(\d+)\s*,\s*"([0-9A-Fa-f]+)"\s*,\s*"([0-9A-Fa-f]+)"')


def log(text):
    with open(M.RAW, "a") as fh:
        fh.write("=== %s\n%s\n" % (C.utc_iso(), text.strip()))


def ask(fd, cmd, timeout=5.0, raw=False):
    reply = command(fd, cmd, timeout)
    log("%s\n%s" % (cmd, reply))
    if raw:
        print("   %s -> %s" % (cmd, " | ".join(l for l in reply.splitlines() if l.strip())), flush=True)
    return reply


def camped(text):
    """Last cell the modem reported in a +CEREG line -> {tac, cell_id} or None."""
    found = None
    for m in CEREG_RE.finditer(text):
        tac, ci = int(m.group(2), 16), int(m.group(3), 16)
        if ci and ci != 0xFFFFFFF:
            found = {"tac": tac, "cell_id": ci}
    return found


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
                    command(fd, "AT+CEREG=2", 2)
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


def lock(fd, earfcn, pci, raw):
    """Lock to the carrier (pci None) or to one cell. Returns (accepted, text of the answers)."""
    text = ask(fd, "AT+CFUN=0", 15, raw)
    text += ask(fd, "AT^CELLLOCK=0", 5, raw)
    cmd = "AT^CELLLOCK=1,3,0,%d" % earfcn if pci is None else "AT^CELLLOCK=1,3,1,%d,%d" % (earfcn, pci)
    r = ask(fd, cmd, 5, raw)
    ask(fd, "AT+CEREG=2", 3, raw)                # registration reports with TAC and Cell ID
    ask(fd, "AT+CFUN=1", 15, raw)
    return "OK" in r


def wait_camp(fd, label, secs, raw):
    """Poll until the modem reports the cell it camps on -> {tac, cell_id} or None."""
    t0, cell = time.time(), None
    while time.time() - t0 < secs:
        C.status("  %s: waiting for the modem to camp... %.0f/%.0f s" % (label, time.time() - t0, secs))
        text = ask(fd, "AT+CEREG?", 3, raw)
        cell = camped(text) or cell
        if cell:
            break
        time.sleep(1)
    C.clear_status()
    return cell


def parse_targets(args, ap):
    """-> [(earfcn, pci or None)]"""
    out = [(e, None) for e in (args.earfcn or [])]
    for item in args.cells or []:
        try:
            e, pcis = item.split(":")
            out += [(int(e), int(p)) for p in pcis.split(",") if p.strip()]
        except ValueError:
            ap.error("--cells takes EARFCN:PCI,PCI (e.g. 6300:21,94), not '%s'" % item)
    if not out:
        ap.error("give --earfcn and/or --cells")
    return out


def main():
    ap = argparse.ArgumentParser(description="Lock the CALYX modem to exact carriers / cells and read the Cell ID.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--earfcn", type=int, nargs="+", help="LTE carriers to lock to")
    ap.add_argument("--cells", nargs="+", metavar="EARFCN:PCI,PCI", help="exact cells to lock to")
    ap.add_argument("--secs", type=float, default=15,
                    help="seconds to wait for the modem to camp (default 15; a cell that is there camps in a few)")
    ap.add_argument("--restart-retry", action="store_true",
                    help="if the modem does not camp, restart it and try once more (slow, rarely helps)")
    ap.add_argument("--port", help="modem AT port (default: find it)")
    ap.add_argument("--keep-lock", action="store_true", help="leave the modem locked at the end")
    ap.add_argument("--no-netscan", action="store_true", help="skip the band scan that matches Cell IDs to PCIs")
    ap.add_argument("--raw", action="store_true", help="print every AT command and answer")
    args = ap.parse_args()
    targets = parse_targets(args, ap)
    if os.geteuid() != 0:
        sys.exit("Run with sudo (the modem port needs it).")
    os.makedirs(M.FIELD, exist_ok=True)
    port = M.find_port(args.port)
    fd = open_port(port)
    command(fd, "ATE0", 2)
    C.info("modem on %s" % port)
    results, scans = [], {}
    try:
        for earfcn, pci in targets:
            label = "EARFCN %d" % earfcn + ("  PCI %d" % pci if pci is not None else "")
            C.heading("%s (band %s)" % (label, M.band_of(earfcn) or "?"))
            if not lock(fd, earfcn, pci, args.raw):
                C.warn("  the modem refused the lock (see --raw)")
            cell = wait_camp(fd, label, args.secs, args.raw)
            if cell is None and args.restart_retry:
                C.info("  the modem has not camped - restarting it and locking again")
                fd = restart(fd, port)
                lock(fd, earfcn, pci, args.raw)
                cell = wait_camp(fd, label, args.secs, args.raw)
            band = M.band_of(earfcn)
            # the band scan only serves to find the PCI: not needed when the target names it
            if cell and pci is None and not args.no_netscan and band and band not in scans:
                scans[band] = M.parse(M.scan(fd, M.RAT["lte"], band, 300))
            match = None
            if cell:
                match = next((c for c in scans.get(M.band_of(earfcn), []) if c["cell_id"] == cell["cell_id"]), None)
            row = {"earfcn": earfcn, "pci": pci, "cell": cell,
                   "seen_pci": match["pci"] if match else None, "seen_earfcn": match["freq"] if match else None,
                   "plmn": match["plmn"] if match else "", "rsrp": match["level"] if match else None}
            results.append(row)
            if cell:
                print("  camped on Cell ID %d (eNB %d, sector %d), TAC %d%s" % (
                    cell["cell_id"], cell["cell_id"] >> 8, cell["cell_id"] & 255, cell["tac"],
                    "  = PCI %d on EARFCN %d" % (match["pci"], match["freq"]) if match else ""), flush=True)
            else:
                print("  " + C.c("the modem did not camp on this %s" % ("cell" if pci is not None else "carrier"),
                                 "red"), flush=True)
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
    sdr = M.sdr_results()
    t = C.Table(["Locked to", "Operator", "TAC", "eNB", "Sector", "Cell ID", "PCI (scan)", "RSRP", "SDR"],
                align={k: "right" for k in ("TAC", "eNB", "Sector", "Cell ID", "PCI (scan)", "RSRP")})
    new = not os.path.exists(OUT)
    with open(OUT, "a") as fh:
        if new:
            fh.write("utc,earfcn,locked_pci,tac,cell_id,scan_pci,scan_earfcn,plmn,rsrp\n")
        for r in results:
            c = r["cell"]
            locked = "%d" % r["earfcn"] + (":%d" % r["pci"] if r["pci"] is not None else "")
            pci = r["pci"] if r["pci"] is not None else r["seen_pci"]
            got = sdr.get((r["earfcn"], pci)) if pci is not None else None
            if not c:
                mark = ("-", "")
            elif got and c["cell_id"] in got:
                mark = ("same Cell ID", "green")
            elif got:
                mark = ("SDR has %s" % ",".join(map(str, sorted(got))), "yellow")
            else:
                mark = ("not decoded by SDR", "red")
            if c:
                t.add([locked, C.op_cell(M.operator(r["plmn"])) if r["plmn"] else "-", c["tac"], c["cell_id"] >> 8,
                       c["cell_id"] & 255, (c["cell_id"], "green+bold"),
                       r["seen_pci"] if r["seen_pci"] is not None else "-",
                       r["rsrp"] if r["rsrp"] is not None else "-", mark])
            else:
                t.add([locked, "-", "-", "-", "-", ("did not camp", "red"), "-", "-", mark])
            fh.write("%s,%d,%s,%s,%s,%s,%s,%s,%s\n" % (
                C.utc_iso(), r["earfcn"], "" if r["pci"] is None else r["pci"], c["tac"] if c else "",
                c["cell_id"] if c else "", "" if r["seen_pci"] is None else r["seen_pci"],
                "" if r["seen_earfcn"] is None else r["seen_earfcn"], r["plmn"],
                "" if r["rsrp"] is None else r["rsrp"]))
    M.give_back(OUT)
    M.give_back(M.RAW)
    t.show("What the modem camped on for each lock")
    C.info("results added to %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
