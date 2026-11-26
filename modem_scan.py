#!/usr/bin/env python3
"""Reference cell scan with the Teltonika CALYX (ALA440 module): every LTE/GSM cell the modem hears, with its
EARFCN, PCI, TAC and Cell ID - then a comparison with what the SDR decoder has decoded.

  sudo python3 modem_scan.py                    # scan LTE, all bands (can take a few minutes)
  sudo python3 modem_scan.py --band 20 3        # only band 20 (800) and band 3 (1800) - faster
  sudo python3 modem_scan.py --rat all          # LTE + GSM (+ 3G)
  sudo python3 modem_scan.py --loop             # repeat, adding new cells to the list

Uses the ALA440's AT^NETSCAN command (manual section 8.15). It works without a SIM and lists the cells of all
operators, like a phone's network-monitor scan, but each line is decoded by the modem from the cell itself.

The table shows for each cell whether the SDR has decoded it (from ~/field/lte_cells.csv, written by lte_decode.py
and lte_track.py) and whether the Cell IDs agree. At the end it prints the lte_decode.py command to look for the
cells the SDR has not decoded yet. Results are added to ~/field/modem_cells.csv; the modem's raw answers are kept
in ~/field/modem_scan_raw.txt (send that file if the table looks wrong).

Before the first run after a reboot the serial ports need the driver:
  sudo modprobe option && echo 1d12 0101 | sudo tee /sys/bus/usb-serial/drivers/option1/new_id
"""
import argparse
import csv
import glob
import os
import re
import select
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                                   # noqa: E402
from modem_probe import open_port, command          # noqa: E402

RAT = {"gsm": 0, "umts": 1, "lte": 2, "all": 3}
LTE_RE = re.compile(r'"?4G"?\s*,\s*Freq:\s*(\d+)\s*,\s*RSRP:\s*(-?\d+)\s*,\s*RSRQ:\s*(-?\d+)\s*,\s*PCI:\s*(\d+)\s*,'
                    r'\s*TAC:\s*(\w+)\s*,\s*CellID:\s*(\w+)(?:\s*,\s*cellBar:\s*(\d+))?', re.I)
GSM_RE = re.compile(r'"?2G"?\s*,\s*Freq:\s*(\d+)\s*,\s*RSSI:\s*(-?\d+)\s*,\s*bsic:\s*(\d+)\s*,\s*LAC:\s*(\w+)'
                    r'.*?CellID:\s*(\w+)(?:\s*,\s*cellBar:\s*(\d+))?', re.I)
OP_RE = re.compile(r'\^NETSCAN:\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,\s*"([^"]*)"')
SDR_CSV = os.path.expanduser("~/field/lte_cells.csv")
OUT_CSV = os.path.expanduser("~/field/modem_cells.csv")
RAW = os.path.expanduser("~/field/modem_scan_raw.txt")


def num(text):
    text = text.strip()
    try:
        return int(text)
    except ValueError:
        return int(text, 16)                         # some firmwares print TAC / Cell ID in hex


def find_port(given):
    if given:
        return given
    ports = sorted(glob.glob("/dev/ttyUSB*") + glob.glob("/dev/ttyACM*"))
    ports.sort(key=lambda p: p != "/dev/ttyUSB2")    # the CALYX AT port is usually ttyUSB2
    for p in ports:
        try:
            fd = open_port(p)
        except OSError:
            continue
        try:
            if "OK" in command(fd, "AT", 1.5):
                return p
        finally:
            os.close(fd)
    sys.exit("No modem AT port found. Is the CALYX plugged in? Load the driver (see --help) and run with sudo.")


def scan(fd, rat, band, timeout):
    """Send AT^NETSCAN, show progress while the modem searches, return its answer lines."""
    cmd = "AT^NETSCAN=%d" % rat + (",%d" % band if band else "")
    os.write(fd, (cmd + "\r").encode())
    buf, t0 = b"", time.time()
    while time.time() - t0 < timeout:
        C.status("  %s: modem scanning... %.0f s" % (cmd, time.time() - t0))
        r, _, _ = select.select([fd], [], [], 1.0)
        if not r:
            continue
        try:
            buf += os.read(fd, 65536)
        except BlockingIOError:
            continue
        text = buf.decode(errors="replace")
        if re.search(r"\r\n(OK|ERROR|\+CME ERROR:[^\r]*)\r\n", text):
            break
    C.clear_status()
    text = buf.decode(errors="replace")
    with open(RAW, "a") as fh:
        fh.write("=== %s %s (%.0f s)\n%s\n" % (C.utc_iso(), cmd, time.time() - t0, text.strip()))
    if "OK" not in text:
        err = re.search(r"(\+CME ERROR:[^\r\n]*|ERROR)", text)
        C.warn("  %s -> %s" % (cmd, err.group(1) if err else "no answer in %d s" % timeout))
    return text.splitlines()


def parse(lines):
    """-> list of cells {rat, freq, pci/bsic, tac/lac, cell_id, level, plmn}."""
    cells, plmn = [], ""
    for line in lines:
        m = OP_RE.search(line)
        if m:
            plmn = m.group(3) or m.group(2)
            continue
        m = LTE_RE.search(line)
        if m:
            cid = num(m.group(6))
            cells.append({"rat": "LTE", "freq": int(m.group(1)), "level": int(m.group(2)), "rsrq": int(m.group(3)),
                          "pci": int(m.group(4)), "tac": num(m.group(5)), "cell_id": cid, "enb": cid >> 8,
                          "sector": cid & 255, "barred": m.group(7) or "", "plmn": plmn})
            continue
        m = GSM_RE.search(line)
        if m:
            cells.append({"rat": "GSM", "freq": int(m.group(1)), "level": int(m.group(2)), "bsic": int(m.group(3)),
                          "tac": num(m.group(4)), "cell_id": num(m.group(5)), "barred": m.group(6) or "",
                          "plmn": plmn})
    return cells


def sdr_results():
    """(earfcn, pci) -> set of Cell IDs the SDR decoded, from lte_cells.csv."""
    out = {}
    if not os.path.exists(SDR_CSV):
        return out
    with open(SDR_CSV) as fh:
        for r in csv.DictReader(fh):
            try:
                key = (int(r["earfcn"]), int(r["pci"]))
            except (ValueError, KeyError, TypeError):
                continue
            out.setdefault(key, set())
            if r.get("cell_id"):
                out[key].add(int(r["cell_id"]))
    return out


def operator(plmn):
    plmn = re.sub(r"\D", "", plmn or "")
    if len(plmn) >= 5:
        return C.operator(plmn[:3], plmn[3:]) or plmn
    return plmn or "-"


def save(cells):
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    new = not os.path.exists(OUT_CSV)
    with open(OUT_CSV, "a") as fh:
        if new:
            fh.write("utc,rat,freq,pci_or_bsic,plmn,tac_or_lac,cell_id,enb,sector,level_dbm,rsrq,barred\n")
        for c in cells:
            fh.write("%s,%s,%d,%s,%s,%d,%d,%s,%s,%d,%s,%s\n" % (
                C.utc_iso(), c["rat"], c["freq"], c.get("pci", c.get("bsic")), c["plmn"], c["tac"], c["cell_id"],
                c.get("enb", ""), c.get("sector", ""), c["level"], c.get("rsrq", ""), c["barred"]))


def show(cells, sdr):
    lte = sorted((c for c in cells if c["rat"] == "LTE"), key=lambda c: (c["freq"], -c["level"]))
    todo = {}
    if lte:
        t = C.Table(["EARFCN", "PCI", "Operator", "TAC", "eNB", "Sector", "Cell ID", "RSRP", "RSRQ", "SDR"],
                    align={k: "right" for k in ("EARFCN", "PCI", "TAC", "eNB", "Sector", "Cell ID", "RSRP", "RSRQ")})
        ok = 0
        for c in lte:
            got = sdr.get((c["freq"], c["pci"]))
            if got and c["cell_id"] in got:
                ok += 1
                mark = ("decoded, same Cell ID", "green")
            elif got:
                mark = ("decoded, other ID %s" % ",".join(map(str, sorted(got))), "yellow")
            elif got is not None:
                mark = ("seen, not decoded", "yellow")
            else:
                mark = ("not decoded", "red")
            if not (got and c["cell_id"] in got):
                todo.setdefault(c["freq"], []).append(c["pci"])
            t.add([c["freq"], c["pci"], C.op_cell(operator(c["plmn"])), c["tac"], c["enb"], c["sector"],
                   (c["cell_id"], "green+bold"), c["level"], c["rsrq"], mark])
        t.show("LTE cells the modem decoded: %d   (SDR has decoded %d of them)" % (len(lte), ok))
    gsm = sorted((c for c in cells if c["rat"] == "GSM"), key=lambda c: (c["freq"], -c["level"]))
    if gsm:
        t = C.Table(["ARFCN", "BSIC", "Operator", "LAC", "Cell ID", "RSSI"],
                    align={k: "right" for k in ("ARFCN", "BSIC", "LAC", "Cell ID", "RSSI")})
        for c in gsm:
            t.add([c["freq"], c["bsic"], C.op_cell(operator(c["plmn"])), c["tac"], (c["cell_id"], "green+bold"),
                   c["level"]])
        t.show("GSM cells the modem decoded: %d" % len(gsm))
    if todo:
        C.info("Cells the SDR has not decoded yet - look for exactly these:")
        print("  python3 lte_decode.py --cells %s --only --loop" % " ".join(
            "%d:%s" % (f, ",".join(map(str, sorted(set(p))))) for f, p in sorted(todo.items())), flush=True)


def main():
    ap = argparse.ArgumentParser(description="Reference cell scan with the CALYX / ALA440 modem (AT^NETSCAN).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--port", help="modem AT port (default: find it, usually /dev/ttyUSB2)")
    ap.add_argument("--rat", choices=sorted(RAT), default="lte", help="what to scan (default lte)")
    ap.add_argument("--band", type=int, nargs="+", help="only these bands (LTE e.g. 20 3 1 8 7; GSM 1=900 2=1800)")
    ap.add_argument("--timeout", type=float, default=600, help="seconds to wait per scan (default 600)")
    ap.add_argument("--loop", action="store_true", help="scan again and again; the table keeps every cell seen")
    ap.add_argument("--raw", action="store_true", help="also print the modem's raw answer")
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    if args.no_color:
        C.set_colour(False)
    os.makedirs(os.path.dirname(RAW), exist_ok=True)
    port = find_port(args.port)
    fd = open_port(port)
    C.info("modem on %s" % port)
    command(fd, "ATE0", 2)
    fun = command(fd, "AT+CFUN?", 3)
    if "+CFUN: 1" not in fun:                         # radio off (flight mode): the scan needs it on
        C.info("radio was off (%s) - switching it on" % " ".join(fun.split()))
        command(fd, "AT+CFUN=1", 15)
        time.sleep(3)
    seen, n = {}, 0
    try:
        while True:
            n += 1
            for band in (args.band or [None]):
                lines = scan(fd, RAT[args.rat], band, args.timeout)
                if args.raw:
                    print("\n".join("   | " + l for l in lines if l.strip()), flush=True)
                cells = parse(lines)
                if not cells and any(l.strip() not in ("", "OK") for l in lines):
                    C.warn("  answer not understood - raw text saved in %s; run with --raw and send it" % RAW)
                save(cells)
                for c in cells:
                    key = (c["rat"], c["freq"], c.get("pci", c.get("bsic")), c["cell_id"])
                    old = seen.get(key)
                    seen[key] = c if old is None or c["level"] > old["level"] else old
            C.heading("Scan %d - %s" % (n, time.strftime("%H:%M:%S")))
            show(list(seen.values()), sdr_results())
            if not args.loop:
                break
    except KeyboardInterrupt:
        print("")
        C.clear_status()
        if seen:
            show(list(seen.values()), sdr_results())
    finally:
        os.close(fd)
    C.info("results added to %s" % OUT_CSV)
    return 0


if __name__ == "__main__":
    sys.exit(main())
