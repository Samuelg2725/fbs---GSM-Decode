#!/usr/bin/env python3
"""Side-by-side test: the CALYX modem and the SDR scan at the same time, then one table compares them.

  sudo python3 compare_scan.py
  sudo python3 compare_scan.py --sdr-args '--dual'          # SDR with both antennas
  sudo python3 compare_scan.py --earfcn 6300 6400 1617      # only these carriers

What it does:
  1. The modem does one full scan first (AT^NETSCAN, all LTE bands, ~20 s), so its carriers are known.
  2. The SDR (lte_decode.py) then goes once through every carrier: the UK list below plus every EARFCN the
     modem has reported (now or in earlier scans). Its output is shown as it runs.
  3. Meanwhile the modem keeps scanning over and over (plus its neighbour list), collecting every cell it hears.
  4. When the SDR has finished: one table of every cell either of them found, with both Cell IDs:
       both        - both decoded it (Cell IDs compared)
       modem only  - the modem decoded it, the SDR did not
       SDR only    - the SDR decoded it, the modem never listed it
     Only results from THIS run count (the SDR writes to its own CSV for the run).
Results: ~/field/compare_<time>.csv. Ctrl+C stops the SDR early and still prints the comparison.
"""
import argparse
import csv
import os
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                                   # noqa: E402
import modem_scan as M                               # noqa: E402
from modem_probe import open_port, command          # noqa: E402

# UK LTE carriers seen so far (Vodafone, O2, EE, Three) - the modem's EARFCNs are added to these
UK_EARFCN = [76, 223, 347, 1363, 1617, 1815, 2850, 3179, 3350, 3501, 3725, 6175, 6225, 6300, 6400]


class Modem(threading.Thread):
    """Keeps scanning in the background; .cells and .neigh collect everything heard."""

    def __init__(self, port, timeout):
        super().__init__(daemon=True)
        self.port, self.timeout = port, timeout
        self.cells, self.neigh, self.scans = {}, {}, 0
        self.stop, self.lock = threading.Event(), threading.Lock()
        self.first = threading.Event()

    def run(self):
        fd = open_port(self.port)
        try:
            command(fd, "ATE0", 2)
            if "+CFUN: 1" not in command(fd, "AT+CFUN?", 3):
                command(fd, "AT+CFUN=1", 15)
                time.sleep(3)
            while not self.stop.is_set():
                cells = M.parse(M.scan(fd, M.RAT["lte"], None, self.timeout, quiet=True))
                neigh = M.neighbours(fd)
                M.save(cells)
                with self.lock:
                    for c in cells:
                        key = (c["freq"], c["pci"])
                        old = self.cells.get(key)
                        if old is None or c["level"] > old["level"]:
                            c["times"] = (old or {}).get("times", 0)
                            self.cells[key] = c
                        self.cells[key]["times"] = self.cells[key].get("times", 0) + 1
                    for n in neigh:
                        key = (n["freq"], n["pci"])
                        if key not in self.neigh or n["rsrp"] > self.neigh[key]["rsrp"]:
                            self.neigh[key] = n
                    self.scans += 1
                self.first.set()
        finally:
            os.close(fd)
            self.first.set()


def sdr_rows(path):
    """(earfcn, pci) -> SDR row for this run (decoded wins over not decoded)."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh):
            try:
                key = (int(r["earfcn"]), int(r["pci"]))
            except (ValueError, KeyError, TypeError):
                continue
            if key not in out or (r.get("cell_id") and not out[key].get("cell_id")):
                out[key] = r
    return out


def compare(modem, sdr_csv, out_csv, secs):
    with modem.lock:
        cells, neigh, scans = dict(modem.cells), dict(modem.neigh), modem.scans
    sdr = sdr_rows(sdr_csv)
    keys = sorted(set(cells) | set(sdr) | {k for k in neigh if k not in cells})
    t = C.Table(["EARFCN", "PCI", "Operator", "Modem Cell ID", "SDR Cell ID", "Modem RSRP", "Result"],
                align={k: "right" for k in ("EARFCN", "PCI", "Modem Cell ID", "SDR Cell ID", "Modem RSRP")})
    count = {"both": 0, "modem only": 0, "SDR only": 0, "different ID": 0, "neither decoded": 0}
    rows = []
    for key in keys:
        m, n, s = cells.get(key), neigh.get(key), sdr.get(key)
        mid = m["cell_id"] if m else None
        sid = int(s["cell_id"]) if s and s.get("cell_id") else None
        op = M.operator(m["plmn"]) if m else (s.get("operator") if s and s.get("operator") else "-")
        rsrp = m["level"] if m else (n["rsrp"] if n else None)
        if mid and sid:
            res = ("both" if mid == sid else "different ID")
        elif mid:
            res = "modem only"
        elif sid:
            res = "SDR only"
        else:
            res = "neither decoded"
        count[res] += 1
        note = res
        if res == "modem only" and s:
            note += " (SDR: %s)" % (s.get("status") or "not decoded")
        elif res == "modem only":
            note += " (SDR did not find the PCI)"
        elif res == "neither decoded":
            note += " (modem neighbour list%s)" % ("; SDR: %s" % s.get("status") if s else ", SDR did not find it")
        colour = {"both": "green", "modem only": "red", "SDR only": "cyan", "different ID": "yellow"}.get(res, "grey")
        t.add([key[0], key[1], C.op_cell(op), mid if mid else "-", sid if sid else "-",
               "%g" % rsrp if rsrp is not None else "-", (note, colour)])
        rows.append([key[0], key[1], op, mid or "", sid or "", "" if rsrp is None else rsrp, note])
    C.heading("Modem vs SDR - same %.0f minutes, %d modem scans" % (secs / 60, scans))
    t.show("Every LTE cell found by either")
    C.info("  ".join("%s: %d" % (k, v) for k, v in count.items() if v))
    C.info("modem only = what the SDR still misses.  SDR only = cells the modem's scan never listed "
           "(it lists only the best cell per operator per carrier).")
    with open(out_csv, "w") as fh:
        w = csv.writer(fh)
        w.writerow(["earfcn", "pci", "operator", "modem_cell_id", "sdr_cell_id", "modem_rsrp", "result"])
        w.writerows(rows)
    M.give_back(out_csv)
    C.info("saved %s" % out_csv)


def main():
    ap = argparse.ArgumentParser(description="Modem and SDR scan at the same time, then compare.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--port", help="modem AT port (default: find it)")
    ap.add_argument("--earfcn", type=int, nargs="+", help="SDR carriers (default: UK list + the modem's)")
    ap.add_argument("--sdr-args", default="", help="extra lte_decode.py options, e.g. '--dual --retries 4'")
    ap.add_argument("--timeout", type=float, default=600, help="seconds per modem scan (default 600)")
    args = ap.parse_args()
    if os.geteuid() != 0:
        sys.exit("Run with sudo (the modem port needs it); the SDR part runs as your normal user.")
    os.makedirs(M.FIELD, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    sdr_csv = os.path.join(M.FIELD, "compare_%s_sdr.csv" % stamp)
    out_csv = os.path.join(M.FIELD, "compare_%s.csv" % stamp)
    port = M.find_port(args.port)
    t0 = time.time()
    modem = Modem(port, args.timeout)
    modem.start()
    C.info("modem on %s: first full scan (about 20 s)..." % port)
    while not modem.first.wait(1):
        C.status("  modem first scan... %.0f s" % (time.time() - t0))
    C.clear_status()
    with modem.lock:
        seen = sorted(modem.cells)
    C.info("modem heard %d cells: %s" % (len(seen), " ".join("%d:%d" % k for k in seen) or "-"))
    earfcns = args.earfcn
    if not earfcns:
        old = set()
        if os.path.exists(M.OUT_CSV):
            with open(M.OUT_CSV) as fh:
                for r in csv.DictReader(fh):
                    if r.get("rat") == "LTE":
                        try:
                            old.add(int(r["freq"]))
                        except (ValueError, TypeError):
                            pass
        earfcns = sorted(set(UK_EARFCN) | {k[0] for k in seen} | old)
    cmd = [sys.executable, os.path.join(HERE, "lte_decode.py"), "--earfcn"] + [str(e) for e in earfcns] + \
          ["--csv", sdr_csv] + args.sdr_args.split()
    user = os.environ.get("SUDO_USER")
    if user and user != "root":
        cmd = ["sudo", "-u", user, "env", "LD_LIBRARY_PATH=/usr/local/lib"] + cmd
    C.heading("SDR: one pass through %d carriers (the modem keeps scanning meanwhile)" % len(earfcns))
    print("  " + " ".join(cmd[cmd.index(sys.executable):]), flush=True)
    try:
        subprocess.run(cmd, env=dict(os.environ, LD_LIBRARY_PATH="/usr/local/lib"))
    except KeyboardInterrupt:
        print("")
        C.warn("SDR stopped early - comparing what it did so far")
    modem.stop.set()
    compare(modem, sdr_csv, out_csv, time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
