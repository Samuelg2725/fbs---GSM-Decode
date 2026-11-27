#!/usr/bin/env python3
"""Modem and SDR on the same cell at the same time: for each target the CALYX modem locks onto it while the
SDR decodes it, then one line compares the two Cell IDs.

  sudo python3 sync_check.py --cells 6300:21 6400:1 223:21 9360:236     # exact cells (EARFCN:PCI)
  sudo python3 sync_check.py --earfcn 6300 6400 223 9360                 # whole carriers
  sudo python3 sync_check.py --cells 9360:236 --dual                     # any lte_decode.py option is passed on

Per target both start together; the next target starts when both are done (the modem usually takes 5-10 s,
the SDR 10-30 s). With --cells the SDR looks for that PCI only; with --earfcn it decodes every cell on the
carrier and the modem shows which one it camps on. The SDR's own output goes to ~/field/sync_check_sdr.log
(--show-sdr prints it). Results: ~/field/sync_check.csv. At the end the modem is restarted to remove the lock.
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
import modem_lock as ML                              # noqa: E402
from modem_probe import open_port, command          # noqa: E402

OUT = os.path.join(M.FIELD, "sync_check.csv")
SDR_LOG = os.path.join(M.FIELD, "sync_check_sdr.log")


def run_sdr(earfcn, pci, extra, show):
    """lte_decode.py on one target, as the normal user -> (rows from its CSV, seconds)."""
    tmp = os.path.join(M.FIELD, "sync_check_tmp.csv")
    if os.path.exists(tmp):
        os.remove(tmp)
    what = ["--cells", "%d:%d" % (earfcn, pci), "--only"] if pci is not None else ["--earfcn", str(earfcn)]
    cmd = [sys.executable, os.path.join(HERE, "lte_decode.py")] + what + ["--csv", tmp, "--no-color"] + extra
    user = os.environ.get("SUDO_USER")
    if user and user != "root":
        cmd = ["sudo", "-u", user, "env", "LD_LIBRARY_PATH=/usr/local/lib"] + cmd
    t0 = time.time()
    with open(SDR_LOG, "a") as log:
        log.write("\n=== %s %s\n" % (C.utc_iso(), " ".join(what)))
        log.flush()
        subprocess.run(cmd, env=dict(os.environ, LD_LIBRARY_PATH="/usr/local/lib"),
                       stdout=None if show else log, stderr=subprocess.STDOUT if not show else None)
    rows = []
    if os.path.exists(tmp):
        with open(tmp) as fh:
            rows = list(csv.DictReader(fh))
        os.remove(tmp)
    M.give_back(SDR_LOG)
    return rows, time.time() - t0


class ModemJob(threading.Thread):
    def __init__(self, fd, earfcn, pci, secs, raw):
        super().__init__(daemon=True)
        self.args = (fd, earfcn, pci, secs, raw)
        self.cell, self.secs, self.at = None, 0.0, None

    def run(self):
        fd, earfcn, pci, secs, raw = self.args
        t0 = time.time()
        self.at = time.strftime("%H:%M:%S")
        ML.lock(fd, earfcn, pci, raw)
        while time.time() - t0 < secs + 10 and self.cell is None:    # +10 s: the lock itself takes a few
            self.cell = ML.camped(ML.ask(fd, "AT+CEREG?", 3, raw))
            if self.cell is None:
                time.sleep(1)
        self.secs = time.time() - t0


def main():
    ap = argparse.ArgumentParser(description="Modem and SDR on the same cell at the same time.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--earfcn", type=int, nargs="+", help="whole carriers")
    ap.add_argument("--cells", nargs="+", metavar="EARFCN:PCI,PCI", help="exact cells")
    ap.add_argument("--secs", type=float, default=15, help="seconds the modem may take to camp (default 15)")
    ap.add_argument("--port", help="modem AT port (default: find it)")
    ap.add_argument("--show-sdr", action="store_true", help="print the SDR's own output too")
    ap.add_argument("--raw", action="store_true", help="print the modem's AT commands and answers")
    args, extra = ap.parse_known_args()
    targets = ML.parse_targets(args, ap)
    if os.geteuid() != 0:
        sys.exit("Run with sudo (the modem port needs it); the SDR part runs as your normal user.")
    os.makedirs(M.FIELD, exist_ok=True)
    port = M.find_port(args.port)
    fd = open_port(port)
    command(fd, "ATE0", 2)
    C.info("modem on %s; SDR output in %s" % (port, SDR_LOG))
    t = C.Table(["Target", "Time", "Modem Cell ID", "Modem s", "SDR Cell ID", "SDR s", "SDR found", "Result"],
                align={k: "right" for k in ("Modem Cell ID", "Modem s", "SDR Cell ID", "SDR s")})
    lines = []
    try:
        for earfcn, pci in targets:
            label = "%d" % earfcn + (":%d" % pci if pci is not None else "")
            C.status("  %s: modem locking and SDR decoding at the same time..." % label)
            job = ModemJob(fd, earfcn, pci, args.secs, args.raw)
            job.start()
            rows, sdr_secs = run_sdr(earfcn, pci, extra, args.show_sdr)
            job.join()
            C.clear_status()
            decoded = {}
            for r in rows:
                if r.get("cell_id") and r.get("pci"):
                    decoded[int(r["pci"])] = (int(r["cell_id"]), r.get("operator", ""))
            mcell = job.cell["cell_id"] if job.cell else None
            if pci is not None:
                sdr_id, op = decoded.get(pci, (None, ""))
            else:                                    # whole carrier: the SDR cell that matches the modem's, if any
                hit = next(((p, v) for p, v in decoded.items() if v[0] == mcell), None)
                sdr_id, op = (hit[1] if hit else (None, ""))
            found = ", ".join("%d" % p for p in sorted(decoded)) or "-"
            if mcell and sdr_id:
                res = ("same Cell ID", "green") if mcell == sdr_id else ("DIFFERENT", "red+bold")
            elif mcell:
                res = ("only the modem", "yellow")
            elif sdr_id:
                res = ("only the SDR", "cyan")
            else:
                res = ("neither", "red")
            print("  %-9s modem %-10s SDR %-10s %s" % (label, mcell or "-", sdr_id or "-", C.c(*res)), flush=True)
            t.add([label, job.at, mcell or "-", "%.0f" % job.secs, sdr_id or "-", "%.0f" % sdr_secs, found, res])
            lines.append([C.utc_iso(), earfcn, "" if pci is None else pci, mcell or "", sdr_id or "", op,
                          "%.1f" % job.secs, "%.1f" % sdr_secs, found.replace(", ", " "), res[0]])
    except KeyboardInterrupt:
        print("")
        C.clear_status()
    finally:
        C.info("removing the lock and restarting the modem (back to normal)...")
        try:
            command(fd, "AT^CELLLOCK=0", 5)
            fd = ML.restart(fd, port)
            os.close(fd)
        except OSError:
            pass
    t.show("Modem and SDR, same cell, same time")
    new = not os.path.exists(OUT)
    with open(OUT, "a") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["utc", "earfcn", "pci", "modem_cell_id", "sdr_cell_id", "operator", "modem_secs",
                        "sdr_secs", "sdr_pcis", "result"])
        w.writerows(lines)
    M.give_back(OUT)
    C.info("results added to %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
