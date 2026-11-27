#!/usr/bin/env python3
"""Modem and SDR on the same cell at the same time, with the full output of both on screen.

  sudo python3 sync_check.py --cells 6300:21 6400:1 223:21 9360:236     # exact cells (EARFCN:PCI)
  sudo python3 sync_check.py --earfcn 6300 6400 223 9360                 # whole carriers
  sudo python3 sync_check.py --cells 9360:236 --dual                     # any lte_decode.py option is passed on

For each target both start together:
  modem | the CALYX locks onto the cell or carrier, camps, and reports the Cell ID, TAC, eNB, sector and its
          signal measurements (RSRP, RSRQ, SINR, RSSI). Usually done in about 5 s.
  SDR   | lte_decode.py's normal output, live: recording, gain, cell search, decode, the reasons when a cell
          does not decode, interference cancellation. Usually 10-30 s, so its lines keep coming after the modem's.
Then a side-by-side block for the target, and the next target starts. At the end, a table of all targets.
--quiet-sdr hides the SDR's own lines (only the results). Results: ~/field/sync_check.csv; the SDR's output is
also kept in ~/field/sync_check_sdr.log. The modem is restarted at the end to remove the lock.
"""
import argparse
import csv
import os
import re
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
HCSQ_RE = re.compile(r'\^HCSQ:\s*(?:\d+\s*,\s*\d+\s*,\s*)?"LTE"\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)')
PRINT = threading.Lock()
ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def say(who, text, colour):
    with PRINT:
        print("  %s %s" % (C.c("%-5s|" % who, colour), text), flush=True)


def signal(text):
    """AT^HCSQ? LTE values -> {rssi, rsrq, rsrp, sinr} in dBm / dB (None where not known)."""
    m = HCSQ_RE.search(text)
    if not m:
        return None
    rxlev, rsrq, rsrp, snr = (int(v) for v in m.groups())
    return {"rssi": None if rxlev in (99, 255) else rxlev - 111,
            "rsrq": None if rsrq == 255 else -20 + 0.5 * rsrq,
            "rsrp": None if rsrp == 255 else rsrp - 141,
            "sinr": None if snr == 255 else -20.2 + 0.2 * snr}


def signal_text(s):
    if not s:
        return "signal not reported"
    parts = [("RSRP %.0f dBm", s["rsrp"]), ("RSRQ %.1f dB", s["rsrq"]), ("SINR %.1f dB", s["sinr"]),
             ("RSSI %.0f dBm", s["rssi"])]
    return ", ".join(f % v for f, v in parts if v is not None) or "signal not reported"


class ModemJob(threading.Thread):
    """Lock, camp, read the cell and its signal - printing each step as it happens."""

    def __init__(self, fd, earfcn, pci, secs, raw):
        super().__init__(daemon=True)
        self.args = (fd, earfcn, pci, secs, raw)
        self.cell, self.sig, self.secs, self.at = None, None, 0.0, None

    def run(self):
        fd, earfcn, pci, secs, raw = self.args
        t0 = time.time()
        self.at = time.strftime("%H:%M:%S")
        target = "EARFCN %d" % earfcn + (" PCI %d" % pci if pci is not None else "")
        ok = ML.lock(fd, earfcn, pci, raw)
        say("modem", "locked to %s (radio off/on) in %.1f s%s" % (
            target, time.time() - t0, "" if ok else " - " + C.c("the modem refused the lock", "red")), "yellow")
        while time.time() - t0 < secs + 10 and self.cell is None:    # +10 s: the lock itself takes a few
            self.cell = ML.camped(ML.ask(fd, "AT+CEREG?", 3, raw))
            if self.cell is None:
                time.sleep(1)
        self.secs = time.time() - t0
        if not self.cell:
            say("modem", C.c("did not camp within %.0f s - no usable cell there for the modem" % self.secs, "red"),
                "yellow")
            return
        c = self.cell
        say("modem", "camped after %.1f s: %s  (eNB %d, sector %d)  TAC %d" % (
            self.secs, C.c("Cell ID %d" % c["cell_id"], "green+bold"), c["cell_id"] >> 8, c["cell_id"] & 255,
            c["tac"]), "yellow")
        time.sleep(1)                                    # let the first measurement come in
        self.sig = signal(ML.ask(fd, "AT^HCSQ?", 3, raw))
        say("modem", "signal: %s" % signal_text(self.sig), "yellow")


def run_sdr(earfcn, pci, extra, quiet):
    """lte_decode.py on one target as the normal user, its output streamed live -> (rows from its CSV, secs)."""
    tmp = os.path.join(M.FIELD, "sync_check_tmp.csv")
    if os.path.exists(tmp):
        os.remove(tmp)
    what = ["--cells", "%d:%d" % (earfcn, pci), "--only"] if pci is not None else ["--earfcn", str(earfcn)]
    cmd = [sys.executable, "-u", os.path.join(HERE, "lte_decode.py")] + what + ["--csv", tmp] + extra
    user = os.environ.get("SUDO_USER")
    if user and user != "root":
        cmd = ["sudo", "-u", user, "env", "LD_LIBRARY_PATH=/usr/local/lib"] + cmd
    t0 = time.time()
    with open(SDR_LOG, "a") as log:
        log.write("\n=== %s %s\n" % (C.utc_iso(), " ".join(what)))
        p = subprocess.Popen(cmd, env=dict(os.environ, LD_LIBRARY_PATH="/usr/local/lib"), stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, bufsize=1, text=True, errors="replace")
        skip_tables = False
        for line in p.stdout:
            log.write(ANSI.sub("", line))
            line = line.rstrip("\n")
            plain = ANSI.sub("", line)
            if plain.startswith("== LTE cells decoded") or plain.startswith("== Not decoded"):
                skip_tables = True                       # its end tables: the comparison below replaces them
            if plain.startswith("=="):
                if plain.startswith("== Pass") or "MHz   EARFCN" not in plain:
                    continue
            if quiet or skip_tables or not plain.strip() or plain.startswith(("+-", "|")):
                continue
            say("SDR", line.strip(), "cyan")
        p.wait()
    rows = []
    if os.path.exists(tmp):
        with open(tmp) as fh:
            rows = list(csv.DictReader(fh))
        os.remove(tmp)
    M.give_back(SDR_LOG)
    return rows, time.time() - t0


def pick_sdr(rows, pci, modem_cell):
    """The SDR row to compare: the PCI asked for, or (whole carrier) the cell matching the modem's."""
    rows = [r for r in rows if r.get("pci")]
    if pci is not None:
        mine = [r for r in rows if int(r["pci"]) == pci]
        return next((r for r in mine if r.get("cell_id")), mine[0] if mine else None)
    if modem_cell:
        return next((r for r in rows if r.get("cell_id") and int(r["cell_id"]) == modem_cell), None)
    return None


def compare_block(label, job, row, rows, sdr_secs):
    mc = job.cell
    sid = int(row["cell_id"]) if row and row.get("cell_id") else None
    if mc and sid:
        res = ("same Cell ID", "green+bold") if mc["cell_id"] == sid else ("DIFFERENT Cell ID", "red+bold")
    elif mc:
        res = ("only the modem decoded it", "yellow")
    elif sid:
        res = ("only the SDR decoded it", "cyan")
    else:
        res = ("neither decoded it", "red")
    t = C.Table(["", "Cell ID", "eNB", "Sector", "TAC", "Operator", "Signal", "Time"],
                align={k: "right" for k in ("Cell ID", "eNB", "Sector", "TAC", "Time")})
    if mc:
        t.add([("modem", "yellow"), (mc["cell_id"], "green+bold"), mc["cell_id"] >> 8, mc["cell_id"] & 255,
               mc["tac"], "-", signal_text(job.sig), "%.0f s" % job.secs])
    else:
        t.add([("modem", "yellow"), ("did not camp", "red"), "-", "-", "-", "-", "-", "%.0f s" % job.secs])
    if sid:
        t.add([("SDR", "cyan"), (sid, "green+bold"), row["enb"], row["sector"], row["tac"],
               C.op_cell(row.get("operator")), "%s MHz wide, %s" % (row.get("width_mhz") or "?", row.get("found_by")),
               "%.0f s" % sdr_secs])
    else:
        status = row.get("status") if row else ("not found" if rows else "no LTE cell found")
        t.add([("SDR", "cyan"), (status or "not decoded", "red"), "-", "-", "-", "-", "-", "%.0f s" % sdr_secs])
    t.show("%s  ->  %s" % (label, C.c(*res)))
    others = sorted({int(r["pci"]) for r in rows if r.get("pci") and r.get("cell_id")}
                    - ({int(row["pci"])} if row and row.get("pci") else set()))
    if others:
        C.info("  the SDR also decoded PCI %s on this carrier" % ", ".join(map(str, others)))
    return res, sid


def main():
    ap = argparse.ArgumentParser(description="Modem and SDR on the same cell at the same time.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--earfcn", type=int, nargs="+", help="whole carriers")
    ap.add_argument("--cells", nargs="+", metavar="EARFCN:PCI,PCI", help="exact cells")
    ap.add_argument("--secs", type=float, default=15, help="seconds the modem may take to camp (default 15)")
    ap.add_argument("--port", help="modem AT port (default: find it)")
    ap.add_argument("--quiet-sdr", action="store_true", help="hide the SDR's own lines, show only the results")
    ap.add_argument("--raw", action="store_true", help="print the modem's AT commands and answers")
    args, extra = ap.parse_known_args()
    targets = ML.parse_targets(args, ap)
    if os.geteuid() != 0:
        sys.exit("Run with sudo (the modem port needs it); the SDR part runs as your normal user.")
    os.makedirs(M.FIELD, exist_ok=True)
    port = M.find_port(args.port)
    fd = open_port(port)
    command(fd, "ATE0", 2)
    C.info("modem on %s; SDR output also kept in %s" % (port, SDR_LOG))
    summary = C.Table(["Target", "Time", "Modem Cell ID", "Modem signal", "Modem s", "SDR Cell ID", "SDR result",
                       "SDR s", "Result"],
                      align={k: "right" for k in ("Modem Cell ID", "Modem s", "SDR Cell ID", "SDR s")})
    lines = []
    try:
        for earfcn, pci in targets:
            label = "%d" % earfcn + (":%d" % pci if pci is not None else "")
            C.heading("%s  (band %s) - modem and SDR start now" % (
                "EARFCN %d PCI %d" % (earfcn, pci) if pci is not None else "EARFCN %d" % earfcn,
                M.band_of(earfcn) or "?"))
            job = ModemJob(fd, earfcn, pci, args.secs, args.raw)
            job.start()
            rows, sdr_secs = run_sdr(earfcn, pci, extra, args.quiet_sdr)
            job.join()
            row = pick_sdr(rows, pci, job.cell["cell_id"] if job.cell else None)
            res, sid = compare_block(label, job, row, rows, sdr_secs)
            mc = job.cell["cell_id"] if job.cell else None
            s = job.sig or {}
            summary.add([label, job.at, mc or "-", "%.0f dBm" % s["rsrp"] if s.get("rsrp") is not None else "-",
                          "%.0f" % job.secs, sid or "-", (row or {}).get("status") or "-", "%.0f" % sdr_secs, res])
            lines.append([C.utc_iso(), earfcn, "" if pci is None else pci, mc or "",
                          job.cell["tac"] if job.cell else "", s.get("rsrp", ""), s.get("rsrq", ""), s.get("sinr", ""),
                          sid or "", (row or {}).get("operator", ""), (row or {}).get("status", ""),
                          "%.1f" % job.secs, "%.1f" % sdr_secs, res[0]])
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
    summary.show("Modem and SDR, same cell, same time")
    new = not os.path.exists(OUT)
    with open(OUT, "a") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["utc", "earfcn", "pci", "modem_cell_id", "modem_tac", "modem_rsrp", "modem_rsrq",
                        "modem_sinr", "sdr_cell_id", "operator", "sdr_status", "modem_secs", "sdr_secs", "result"])
        w.writerows(lines)
    M.give_back(OUT)
    C.info("results added to %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
