#!/usr/bin/env python3
"""LTE tracking on one carrier: every cell (PCI) sharing that EARFCN, decoded and watched continuously,
with the raw IQ recordings saved for later decoding or direction finding.

  Track one carrier (Ctrl-C to stop and show the table):
      python3 lte_track.py --earfcn 6400
      python3 lte_track.py --freq-mhz 816
  Both antennas (recordings hold RX1 and RX2 together - what direction finding needs):
      python3 lte_track.py --earfcn 6400 --dual
  Save every recording instead of only those with a new cell:
      python3 lte_track.py --earfcn 6400 --save-iq all

Every cycle it:
  1. records the carrier: a short --watch-secs recording when it only needs to see which PCIs are present,
     the full --secs when something must be decoded or saved;
  2. searches the recording for every PCI, including weaker cells under a strong one (and searches a second
     window if a cell heard before is missing from the first);
  3. decodes with srsRAN any PCI it has not identified yet (Cell ID, eNB, sector, TAC, operator), and from
     time to time tries the other sectors of each decoded mast (they can be decoded with the mast's timing
     even when their own sync signal is too weak to be found);
  4. prints one line with every PCI heard, and saves the raw IQ when a new PCI appears (see --save-iq).

IQ recordings are saved in SigMF format in ~/field/iq/: a .sigmf-data file (16-bit I/Q, full scale 2048;
with --dual each sample is I1 Q1 I2 Q2) and a .sigmf-meta file (JSON) with the frequency, sample rate, gain,
time and, for every PCI in the recording, its subframe-0 sample position, frame length, frequency offset and
identity. The radio is tuned --low-if-mhz above the carrier, so the carrier sits at -0.6 MHz in the files.

Results: a table on exit (and every --table-every cycles), ~/field/lte_track_summary.csv (one row per PCI)
and ~/field/lte_track_log.csv (one row per PCI per cycle, with the IQ file it is in).
Run with LD_LIBRARY_PATH=/usr/local/lib.
"""
import argparse
import json
import os
import queue
import sys
import threading
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                                   # noqa: E402
import lte_decode as L                               # noqa: E402  (also sets up the "osprey 2" library path)
import lte_fastsync                                  # noqa: E402


# ------------------------------------------------------------------------------------------- cells
class Cell(object):
    """Everything known about one PCI on the tracked carrier."""

    def __init__(self, pci, cycle):
        self.pci = pci
        self.ident = None            # dict from lte_decode.decode_cell when identified
        self.width = None            # MIB resource blocks
        self.status = "not decoded yet"
        self.how = ""                # "cell search", "weak-cell search", "RX2 only", "sector of PCI n's mast"
        self.sector_only = False     # found only through its mast's timing (no sync signal of its own)
        self.first = time.time()
        self.last = None
        self.first_cycle = cycle
        self.detected = 0            # cycles in which the search found it
        self.strengths = []
        self.tries = 0               # srsRAN attempts
        self.last_try = -999
        self.checks = 0              # sector checks (sector_only cells)
        self.confirmed = 0
        self.missing = 0             # consecutive cycles not heard
        self.offset = 0.0            # frequency offset from the last search that found it
        self.lost_reported = False
        self.files = 0               # IQ recordings it appears in

    @property
    def label(self):
        if self.ident:
            return "%s %d" % (self.ident["operator"] or "?", self.ident["cell_id"])
        return "?"


def colour_of(cell):
    if cell.ident:
        return C.OPERATOR_COLOUR.get(cell.ident["operator"], "") or "green"
    return "yellow"


def result_line(cell, extra=""):
    if cell.ident:
        x = cell.ident
        text = "PCI %-3d %s  Cell ID %d (eNB %d, sector %d)  TAC %d  %s %s-%s" % (
            cell.pci, C.c("decoded", "green"), x["cell_id"], x["enb"], x["sector"], x["tac"],
            C.c(x["operator"] or "?", colour_of(cell)), x["mcc"], x["mnc"])
    else:
        text = "PCI %-3d %s" % (cell.pci, C.c(cell.status, "yellow" if "MIB only" in cell.status else "red"))
    return "           " + text + C.c("   [%s%s]" % (cell.how, extra), "grey")


# ------------------------------------------------------------------------------------------ IQ files
class IQWriter(object):
    """Writes recordings in a background thread so saving never slows the tracking down. Keeps the folder
    under --iq-max-gb by deleting the oldest recordings this tool made."""

    def __init__(self, folder, max_gb):
        self.folder = os.path.expanduser(folder)
        os.makedirs(self.folder, exist_ok=True)
        self.max_bytes = max_gb * 1e9
        self.q = queue.Queue()
        self.saved = 0
        self.thread = threading.Thread(target=self._run, daemon=True)
        self.thread.start()

    def save(self, name, raw, meta):
        self.q.put((name, raw, meta))
        return name + ".sigmf-data"

    def _run(self):
        while True:
            item = self.q.get()
            if item is None:
                return
            name, raw, meta = item
            try:
                base = os.path.join(self.folder, name)
                raw.astype("<i2").tofile(base + ".sigmf-data")
                with open(base + ".sigmf-meta", "w") as fh:
                    json.dump(meta, fh, indent=1)
                self.saved += 1
                self._trim()
            except Exception as e:                       # disk full etc.: tracking carries on
                C.warn("  could not save IQ recording %s: %s" % (name, e))

    def _trim(self):
        if self.max_bytes <= 0:
            return
        files = sorted((os.path.getmtime(os.path.join(self.folder, f)), f) for f in os.listdir(self.folder)
                       if f.startswith("lte_") and f.endswith(".sigmf-data"))
        total = sum(os.path.getsize(os.path.join(self.folder, f)) for _t, f in files)
        while files and total > self.max_bytes:
            _t, f = files.pop(0)
            path = os.path.join(self.folder, f)
            total -= os.path.getsize(path)
            for ext in (".sigmf-data", ".sigmf-meta"):
                try:
                    os.remove(path[:-len(".sigmf-data")] + ext)
                except OSError:
                    pass

    def close(self):
        if self.q.qsize():
            C.info("finishing %d IQ recording(s) still being saved..." % self.q.qsize())
        self.q.put(None)
        self.thread.join()


def sigmf_meta(args, rate, t_capture, cells_in_file):
    carrier_hz = args.freq_mhz_now * 1e6
    tuned_hz = carrier_hz + args.low_if_mhz * 1e6
    annotations = []
    for cell, info in cells_in_file:
        width_mhz = L.BW_MHZ.get(cell.width, 10)
        a = {"core:sample_start": int(info["start"]), "core:sample_count": int(round(rate * 0.001)),
             "core:label": "PCI %d subframe 0" % cell.pci,
             "core:freq_lower_edge": carrier_hz - width_mhz * 0.45e6,
             "core:freq_upper_edge": carrier_hz + width_mhz * 0.45e6,
             "lte:pci": cell.pci, "lte:frame_samples": int(round(rate * 0.01)),
             "lte:frequency_offset_hz": round(float(info["offset"]), 1),
             "lte:pss_strength": round(float(info["strength"]), 1) if info["strength"] else None,
             "lte:found_by": cell.how}
        if cell.ident:
            x = cell.ident
            a.update({"lte:cell_id": x["cell_id"], "lte:enb": x["enb"], "lte:sector": x["sector"],
                      "lte:tac": x["tac"], "lte:mcc": x["mcc"], "lte:mnc": x["mnc"], "lte:operator": x["operator"],
                      "lte:width_mhz": width_mhz})
        annotations.append(a)
    annotations.sort(key=lambda a: a["core:sample_start"])
    chans = "RX1 and RX2, interleaved per sample (I1 Q1 I2 Q2)" if args.dual else "RX1"
    return {
        "global": {
            "core:datatype": "ci16_le", "core:sample_rate": rate, "core:version": "1.0.0",
            "core:num_channels": 2 if args.dual else 1, "core:hw": "BladeRF 2.0 micro",
            "core:description": "LTE EARFCN %s (%.1f MHz). %s. 12-bit samples, full scale 2048. The radio was "
                                "tuned %.1f MHz above the carrier." % (args.earfcn_now, args.freq_mhz_now, chans,
                                                                         args.low_if_mhz),
            "lte:earfcn": args.earfcn_now, "lte:carrier_hz": carrier_hz, "lte:gain_db": args.gain_now,
            "lte:channels": chans,
        },
        "captures": [{"core:sample_start": 0, "core:frequency": tuned_hz,
                      "core:datetime": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t_capture)) +
                      ".%03dZ" % int((t_capture % 1) * 1000)}],
        "annotations": annotations,
    }


# --------------------------------------------------------------------------------------------- tables
def show_table(cells, cycles, t_all, final=False):
    t = C.Table(["PCI", "Operator", "MCC-MNC", "TAC", "eNB", "Sector", "Cell ID", "Width", "Found by",
                 "Heard", "Strength", "First seen", "Last seen"],
                align={k: "right" for k in ("PCI", "TAC", "eNB", "Sector", "Cell ID", "Strength")})
    order = sorted(cells.values(), key=lambda c: (c.ident is None, c.sector_only, c.pci))
    for cell in order:
        x = cell.ident or {}
        if cell.sector_only:
            heard = ("%d/%d checks" % (cell.confirmed, cell.checks), C.rate_colour(
                100.0 * cell.confirmed / max(1, cell.checks)))
        else:
            span = max(1, cycles - cell.first_cycle + 1)
            pct = 100.0 * cell.detected / span
            heard = ("%d/%d (%.0f%%)" % (cell.detected, span, pct), C.rate_colour(pct))
        strength = ("%.0f" % np.mean(cell.strengths)) if cell.strengths else "-"
        if cell.ident:
            cid = (x["cell_id"], "green+bold")
        else:
            cid = (cell.status, "yellow" if "MIB only" in cell.status else "red")
        t.add([cell.pci, C.op_cell(x.get("operator")) if x else ("?", "yellow"),
               "%s-%s" % (x["mcc"], x["mnc"]) if x else "-", x.get("tac", "-"), x.get("enb", "-"),
               x.get("sector", "-"), cid, "%s MHz" % L.BW_MHZ.get(cell.width, "?") if cell.width else "-",
               cell.how, heard, strength, time.strftime("%H:%M:%S", time.localtime(cell.first)),
               time.strftime("%H:%M:%S", time.localtime(cell.last)) if cell.last else "-"])
    decoded = sum(1 for c in cells.values() if c.ident)
    t.show("%s: %d PCI(s) on this carrier, %d identified   (%d cycles, %.0f s)" % (
        "Cells on the carrier" if final else "So far", len(cells), decoded, cycles, time.time() - t_all))


def write_summary_csv(path, args, cells, cycles):
    path = os.path.expanduser(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new:
            fh.write("utc,label,earfcn,mhz,pci,operator,mcc,mnc,tac,enb,sector,cell_id,width_mhz,found_by,"
                     "cycles_heard,cycles,sector_checks,sector_confirmed,mean_strength,first_seen,last_seen,"
                     "iq_files,status\n")
        for cell in sorted(cells.values(), key=lambda c: c.pci):
            x = cell.ident or {}
            fh.write(",".join(str(v) for v in (
                C.utc_iso(), args.label, args.earfcn_now, "%.1f" % args.freq_mhz_now, cell.pci,
                x.get("operator", ""), x.get("mcc", ""), x.get("mnc", ""), x.get("tac", ""), x.get("enb", ""),
                x.get("sector", ""), x.get("cell_id", ""), L.BW_MHZ.get(cell.width, "") if cell.width else "",
                cell.how.replace(",", ";"), cell.detected, cycles, cell.checks, cell.confirmed,
                "%.1f" % np.mean(cell.strengths) if cell.strengths else "", C.utc_iso(cell.first),
                C.utc_iso(cell.last) if cell.last else "", cell.files, cell.status.replace(",", ";"))) + "\n")
    C.info("summary added to %s" % path)


class CycleLog(object):
    def __init__(self, path):
        self.fh = None
        if not path:
            return
        path = os.path.expanduser(path)
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        new = not os.path.exists(path)
        self.fh = open(path, "a")
        if new:
            self.fh.write("utc,cycle,earfcn,mhz,gain_db,pci,heard_by,strength,cell_id,operator,iq_file\n")

    def add(self, args, cycle, cell, heard_by, strength, iq_file):
        if self.fh:
            x = cell.ident or {}
            self.fh.write("%s,%d,%s,%.1f,%g,%d,%s,%s,%s,%s,%s\n" % (
                C.utc_iso(), cycle, args.earfcn_now, args.freq_mhz_now, args.gain_now, cell.pci, heard_by,
                "%.1f" % strength if strength else "", x.get("cell_id", ""), x.get("operator", ""), iq_file or ""))

    def close(self):
        if self.fh:
            self.fh.close()


# ---------------------------------------------------------------------------------------------- track
def _search_at(args, rate, x1, x2, start):
    hits = list(L.search_cells(x1[start:], rate, args)["hits"])
    if x2 is not None:
        have = {h["pci"] for h in hits}
        for h in L.search_cells(x2[start:], rate, args)["hits"]:
            if h["pci"] not in have:
                hits.append(dict(h, rx2_only=True))
    return hits


def search(args, rate, x1, x2, expect=()):
    """Cell search on the start of the recording. If a cell heard before is missing (one short search window
    can miss a cell that is there), search a second window later in the same recording."""
    hits = _search_at(args, rate, x1, x2, 0)
    have = {int(h["pci"]) for h in hits}
    if set(expect) - have and len(x1) >= int(rate * args.search_secs * 2):
        for h in _search_at(args, rate, x1, x2, len(x1) // 2):
            if int(h["pci"]) not in have:
                have.add(int(h["pci"]))
                hits.append(h)
    return hits


def recording(args, rate, secs):
    keep, args.secs = args.secs, secs
    try:
        return L.record(args, args.freq_mhz_now, rate)
    finally:
        args.secs = keep


def retry_due(args, cell, cycles):
    """An unidentified PCI is tried every --retry-every cycles (default: every cycle) for --max-tries tries,
    then every 10 cycles. A weak cell is decoded on some recordings and not others, so tries pay off."""
    if cell.ident:
        return False
    every = args.retry_every if cell.tries < args.max_tries else max(args.retry_every, 10)
    return cell.tries == 0 or cycles - cell.last_try >= every


def work_due(args, cells, cycles):
    """True when this cycle needs a full-length recording: an unidentified PCI is due another srsRAN try,
    or the other sectors of a mast are due a check."""
    for cell in cells.values():
        if cell.sector_only:
            continue
        if not cell.ident and cell.missing <= args.lost_after and retry_due(args, cell, cycles):
            return True
        if cell.missing:
            continue
        if cell.ident and not args.no_siblings:
            g = cell.pci // 3
            for sib in (3 * g, 3 * g + 1, 3 * g + 2):
                other = cells.get(sib)
                if other is not None and not other.sector_only:
                    continue
                if cycles - args._sector_tried.get(sib, -999) >= args.sector_every:
                    return True
    return False


def signal_for(args, rate, x1, x2, pci, offset):
    """The recording to decode one PCI from: RX1, or both antennas combined for that cell."""
    if x2 is None:
        return x1
    w1, w2, _share = L.cell_weights(x1, x2, rate, pci, offset)
    return (w1 * x1 + w2 * x2).astype(np.complex64)


def apply_result(cell, r):
    cell.tries += 1
    if r.get("mib"):
        cell.width = r["mib"]
    if r["ident"]:
        cell.ident = r["ident"]
    cell.status = r["status"]


def track(args):
    rate = 15.36e6
    args.gain_now = args.gain
    cells, cycles, t_all = {}, 0, time.time()
    writer = IQWriter(args.iq_dir, args.iq_max_gb) if args.save_iq != "none" else None
    log = CycleLog(args.log)
    fails = 0
    regain = False
    C.info("Tracking EARFCN %s (%.1f MHz), %s, %.1f s recordings. IQ: %s. Ctrl-C to stop." % (
        args.earfcn_now, args.freq_mhz_now, "both antennas" if args.dual else "RX1",
        args.secs, {"new": "saved when a new PCI appears, in " + args.iq_dir, "all": "every recording saved, in " +
                    args.iq_dir, "none": "not saved"}[args.save_iq]))
    try:
        while args.cycles == 0 or cycles < args.cycles:
            cycles += 1
            t0 = time.time()
            full = cycles == 1 or args.save_iq == "all" or work_due(args, cells, cycles)
            C.status("  cycle %d: recording..." % cycles)
            t_capture = time.time()
            x1, x2 = recording(args, rate, args.secs if full else args.watch_secs)
            if cycles == 1 or regain:
                regain = False
                L.SHOW_LEVELS[0] = True
                x1, x2 = L.auto_gain(args, args.freq_mhz_now, rate, x1, x2)
                L.SHOW_LEVELS[0] = False
            if x1 is None:
                cycles -= 1
                fails += 1
                if fails >= 2 and L.RADIO[0] is not None:
                    C.warn("  the open radio keeps failing - switching to bladeRF-cli (slower) for the rest")
                    L.RADIO[0].close()
                    L.RADIO[0] = None
                    fails = 0
                elif fails >= 5:
                    C.warn("Stopping: the BladeRF could not be used 5 times in a row. Check that it is plugged in "
                           "and that no other program is using it (bladeRF-cli -p lists it; "
                           "ps aux | grep -i -e blade -e gsm -e lte shows other users).")
                    break
                time.sleep(1.0)
                continue
            fails = 0
            peak = L.LAST_PEAK[0]
            raw = L.LAST_RAW[0]
            t_rec = time.time() - t0
            C.status("  cycle %d: searching for cells..." % cycles)
            t1 = time.time()
            hits = search(args, rate, x1, x2, [c.pci for c in cells.values() if not c.sector_only and not c.missing])
            t_search = time.time() - t1
            if not full and ({int(h["pci"]) for h in hits} - set(cells)):
                # a new PCI: record the full length so it can be decoded and saved
                C.status("  cycle %d: new PCI - recording %.1f s to decode it..." % (cycles, args.secs))
                t_capture = time.time()
                y1, y2 = recording(args, rate, args.secs)
                if y1 is not None:
                    x1, x2, full = y1, y2, True
                    raw, peak = L.LAST_RAW[0], max(peak, L.LAST_PEAK[0])
            now = time.time()
            heard = {}                                          # pci -> (cell, offset, strength, how)
            new_pcis = []
            for h in hits:
                pci = int(h["pci"])
                cell = cells.get(pci)
                if cell is None:
                    cell = cells[pci] = Cell(pci, cycles)
                    new_pcis.append(pci)
                if cell.sector_only:                            # now strong enough to be found on its own
                    cell.sector_only = False
                    cell.first_cycle = cycles
                cell.how = cell.how if cell.how and not cell.how.startswith("sector") else (
                    "weak-cell search" if h.get("weak") else ("RX2 only" if h.get("rx2_only") else "cell search"))
                cell.detected += 1
                cell.strengths.append(float(h["pss_ratio"]))
                cell.last = now
                if cell.lost_reported:
                    C.clear_status()
                    print("%s  %s" % (C.c(C.now_stamp(), "grey"), C.c("PCI %d (%s) heard again" % (
                        pci, cell.label), "cyan")), flush=True)
                cell.missing, cell.lost_reported = 0, False
                cell.offset = float(h.get("offset_hz", 0.0))
                heard[pci] = (cell, cell.offset, float(h["pss_ratio"]), cell.how)
            for cell in cells.values():
                if cell.pci not in heard and not cell.sector_only:
                    cell.missing += 1
                    if cell.missing == args.lost_after and not cell.lost_reported:
                        cell.lost_reported = True
                        C.clear_status()
                        print("%s  %s" % (C.c(C.now_stamp(), "grey"), C.c("PCI %d (%s) not heard for %d cycles" % (
                            cell.pci, cell.label, args.lost_after), "yellow")), flush=True)
            for pci in new_pcis:
                C.clear_status()
                print("%s  %s" % (C.c(C.now_stamp(), "grey"), C.c("NEW PCI %d (strength %.0f, %s)" % (
                    pci, heard[pci][2], heard[pci][3]), "cyan+bold")), flush=True)

            # identify PCIs not decoded yet (first time at once, then every --retry-every cycles)
            aligned = {}
            want_rate = rate
            # a PCI the search missed this time but heard in the last few cycles is tried too (with its
            # last known frequency offset): the short search window misses weak cells that are still there
            todo = [(pci, v[0], v[1], v[2]) for pci, v in heard.items()]
            todo += [(c.pci, c, c.offset, None) for c in cells.values() if c.pci not in heard and not c.ident
                     and not c.sector_only and c.missing <= args.lost_after]
            for pci, cell, offset, strength in sorted(todo, key=lambda t: -(t[3] or 0)):
                if cell.ident or not full or not retry_due(args, cell, cycles):
                    continue
                C.status("  cycle %d: decoding PCI %d with srsRAN..." % (cycles, pci))
                cell.last_try = cycles
                r = L.decode_cell(signal_for(args, rate, x1, x2, pci, offset), rate,
                                  {"pci": pci, "pss_ratio": strength, "offset_hz": offset, "width": cell.width}, args, args.freq_mhz_now)
                apply_result(cell, r)
                if r.get("aligned") is not None:
                    aligned[pci] = r["aligned"]
                if cell.ident and pci not in heard:              # decoded although the search missed it
                    cell.last = time.time()
                    heard[pci] = (cell, offset, None, cell.how)
                C.clear_status()
                print(result_line(cell, ", try %d" % cell.tries if cell.tries > 1 else ""), flush=True)
                if cell.width:
                    want_rate = L.carrier_rate(cell.width)

            # the other sectors of each decoded mast (same frame timing) - new masts at once, then every
            # --sector-every cycles to confirm they are still there
            if full and not args.no_siblings:
                groups = {}
                for pci, (cell, offset, _s, _h) in heard.items():
                    if cell.ident:
                        groups.setdefault(pci // 3, (pci, offset))
                for g, (parent, offset) in sorted(groups.items()):
                    for sib in (3 * g, 3 * g + 1, 3 * g + 2):
                        if sib in heard:
                            continue
                        cell = cells.get(sib)
                        if cell is not None and not cell.sector_only:
                            continue                            # found by the search before: not a sector try
                        if cell is not None and cycles - cell.last_try < args.sector_every:
                            continue
                        if cell is None and sib in args._sector_tried and \
                                cycles - args._sector_tried[sib] < args.sector_every:
                            continue
                        args._sector_tried[sib] = cycles
                        if parent not in aligned:
                            C.status("  cycle %d: timing of PCI %d..." % (cycles, parent))
                            aligned[parent] = L.align_wide(signal_for(args, rate, x1, x2, parent, offset), rate,
                                                           parent, offset)
                        if aligned[parent] is None:
                            continue
                        C.status("  cycle %d: trying PCI %d (sector of PCI %d's mast)..." % (cycles, sib, parent))
                        r = L.decode_cell(None, rate, {"pci": sib, "width": cells[parent].width}, args, args.freq_mhz_now, aligned=aligned[parent])
                        if cell is None and not r["ident"]:
                            continue
                        if cell is None:
                            cell = cells[sib] = Cell(sib, cycles)
                            cell.sector_only = True
                            cell.how = "sector of PCI %d's mast" % parent
                        cell.last_try = cycles
                        cell.checks += 1
                        was = cell.ident
                        if r["ident"]:
                            apply_result(cell, r)
                            cell.confirmed += 1
                            cell.last = time.time()
                            heard[sib] = (cell, offset, None, cell.how)
                            aligned[sib] = aligned[parent]
                            if not was:
                                C.clear_status()
                                print("%s  %s" % (C.c(C.now_stamp(), "grey"), C.c(
                                    "NEW PCI %d (sector of PCI %d's mast)" % (sib, parent), "cyan+bold")),
                                      flush=True)
                                print(result_line(cell), flush=True)
                                new_pcis.append(sib)

            # save the raw IQ
            iq_file = None
            if writer and raw is not None and full and (args.save_iq == "all" or new_pcis):
                in_file = []
                for pci, (cell, offset, strength, _how) in heard.items():
                    a = aligned.get(pci)
                    if a is None:
                        a = L.align_wide(signal_for(args, rate, x1, x2, pci, offset), rate, pci, offset)
                    if a is None:
                        continue
                    in_file.append((cell, {"start": len(x1) - len(a), "offset": offset, "strength": strength}))
                    cell.files += 1
                name = "lte_%s_%s_c%04d" % (args.earfcn_now, time.strftime("%Y%m%dT%H%M%SZ",
                                                                           time.gmtime(t_capture)), cycles)
                iq_file = writer.save(name, raw, sigmf_meta(args, rate, t_capture, in_file))

            # one line per cycle
            parts = []
            for pci in sorted(heard, key=lambda p: (heard[p][2] is None, -(heard[p][2] or 0))):
                cell, _o, strength, _h = heard[pci]
                txt = "PCI %d %s" % (pci, cell.label if cell.ident else "?")
                txt += " [%.0f]" % strength if strength else " [sector]"
                parts.append(C.c(txt, colour_of(cell)))
            level = "%.0f dBFS%s" % (L.LAST_LEVEL[0], " CLIPPING" if peak >= 0.99 else "")
            C.clear_status()
            print("%s  %s  %s  %s%s" % (
                C.c(C.now_stamp(), "grey"), C.c("#%-4d %4.1fs  gain %2g %s" % (cycles, time.time() - t0,
                                                                              args.gain_now, level),
                                                "red" if peak >= 0.99 else "grey"),
                "  ".join(parts) if parts else C.c("no LTE cell heard", "red"),
                C.c("   IQ saved" if iq_file else "", "grey"),
                C.c("   (record %.1fs, search %.1fs%s)" % (t_rec, t_search, ", full" if full else ""), "grey")
                ), flush=True)
            for pci, (cell, _o, strength, how) in heard.items():
                log.add(args, cycles, cell, "search" if strength else "sector check", strength, iq_file)

            # keep the gain right as conditions change; switch to a wider recording for a wider cell
            if not args.no_auto_gain:
                if peak >= 0.99 and args.gain_now > 0:
                    args.gain_now = max(0.0, args.gain_now - 6)
                    C.info("  clipping - gain lowered to %g dB" % args.gain_now)
                elif 0 < peak < 0.12 and args.gain_now < 60:
                    args.gain_now = min(60.0, args.gain_now + 6)
                    C.info("  signal low - gain raised to %g dB" % args.gain_now)
            if want_rate != rate:
                rate = want_rate
                C.info("  the carrier is %s MHz wide - recording %.2f Msps from now on%s" % (
                    L.BW_MHZ.get(max(c.width or 0 for c in cells.values()), "?"), rate / 1e6,
                    " (only this carrier, so the gain can rise for weak cells)" if rate < 15.36e6 else ""))
                regain = True
                for cell in cells.values():
                    if not cell.ident:
                        cell.tries, cell.last_try = 0, -999
            if args.table_every and cycles % args.table_every == 0:
                show_table(cells, cycles, t_all)
    except KeyboardInterrupt:
        C.clear_status()
        print("")
    finally:
        C.ignore_ctrl_c()
        if writer:
            writer.close()
        log.close()
        if L.RADIO[0] is not None:
            L.RADIO[0].close()
    show_table(cells, cycles, t_all, final=True)
    if writer:
        C.info("%d IQ recording(s) saved in %s" % (writer.saved, writer.folder))
    if args.csv and cells:
        write_summary_csv(args.csv, args, cells, cycles)
    if args.log:
        C.info("per-cycle log added to %s" % os.path.expanduser(args.log))


def main():
    ap = argparse.ArgumentParser(description="Track every LTE cell on one carrier, saving raw IQ.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--earfcn", type=int, help="carrier channel number (EARFCN)")
    what.add_argument("--freq-mhz", type=float, help="carrier centre frequency in MHz")
    ap.add_argument("--dual", action="store_true", help="use both antennas (RX1 + RX2); recordings hold both")
    ap.add_argument("--gain", type=float, default=20.0, help="starting gain in dB (default 20, then automatic)")
    ap.add_argument("--no-auto-gain", action="store_true", help="keep --gain fixed")
    ap.add_argument("--secs", type=float, default=1.0,
                    help="recording length for decoding and for saved IQ, s (default 1.0)")
    ap.add_argument("--watch-secs", type=float, default=0.2,
                    help="recording length on cycles that only check which PCIs are present, s (default 0.2)")
    ap.add_argument("--cycles", type=int, default=0, help="stop after this many cycles (default 0 = until Ctrl-C)")
    ap.add_argument("--save-iq", choices=("new", "all", "none"), default="new",
                    help="save recordings: 'new' = when a new PCI appears (default), 'all' = every cycle, 'none'")
    ap.add_argument("--iq-dir", default="~/field/iq", help="folder for IQ recordings (default ~/field/iq)")
    ap.add_argument("--iq-max-gb", type=float, default=20.0,
                    help="delete the oldest recordings above this total size, GB (default 20, 0 = no limit)")
    ap.add_argument("--retry-every", type=int, default=1,
                    help="re-try decoding an unidentified PCI every N cycles (default 1 = every cycle)")
    ap.add_argument("--max-tries", type=int, default=10, help="tries at --retry-every before slowing to every 10 cycles (default 10)")
    ap.add_argument("--sector-every", type=int, default=10,
                    help="re-check the other sectors of each mast every N cycles (default 10)")
    ap.add_argument("--no-siblings", action="store_true", help="do not try the other sectors of each mast")
    ap.add_argument("--lost-after", type=int, default=3, help="report a PCI as not heard after N cycles")
    ap.add_argument("--table-every", type=int, default=20, help="print the table every N cycles (0 = only at the end)")
    ap.add_argument("--radio", choices=("lib", "cli"), default="lib",
                    help="lib (default): keep the radio open between recordings (libbladeRF); cli: start "
                         "bladeRF-cli for every recording (slower; used automatically if lib fails)")
    ap.add_argument("--low-if-mhz", type=float, default=0.6, help="tune this far above the carrier (default 0.6)")
    ap.add_argument("--search-secs", type=float, default=0.08, help="recording length used for the search, s")
    ap.add_argument("--no-weak-search", action="store_true", help="skip the search for weaker cells")
    ap.add_argument("--pss-min", type=float, default=None, help="weak search threshold (default 6.0)")
    ap.add_argument("--srsran-timeout", type=int, default=12, help="seconds per srsRAN attempt (default 12)")
    ap.add_argument("--srsran-log", default=None, help="save every srsRAN output to this folder (diagnosis)")
    ap.add_argument("--debug-search", action="store_true", help="print every weak-search peak")
    ap.add_argument("--label", default="-", help="name for this run in the summary CSV")
    ap.add_argument("--csv", default="~/field/lte_track_summary.csv", help="summary CSV ('' = none)")
    ap.add_argument("--log", default="~/field/lte_track_log.csv", help="per-cycle CSV ('' = none)")
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    args._sector_tried = {}
    if args.no_color:
        C.set_colour(False)
    if not L.pdsch_ue().is_file():
        sys.exit("srsRAN pdsch_ue was not found at %s - build it first (see README)." % L.pdsch_ue())
    lte_fastsync.install()
    if args.earfcn is not None:
        args.freq_mhz_now = L.earfcn_to_mhz(args.earfcn)
        if args.freq_mhz_now is None:
            sys.exit("EARFCN %d is not an LTE downlink channel" % args.earfcn)
        args.earfcn_now = args.earfcn
        band = L.mhz_to_earfcn(args.freq_mhz_now)[1]
    else:
        args.freq_mhz_now = args.freq_mhz
        args.earfcn_now, band = L.mhz_to_earfcn(args.freq_mhz)
    C.heading("LTE tracking   EARFCN %s   %.1f MHz   band %s" % (args.earfcn_now, args.freq_mhz_now, band or "?"))
    if args.radio == "lib":
        try:
            L.RADIO[0] = L.LibRadio(args.dual)
            C.info("radio kept open between recordings (libbladeRF)")
        except Exception as e:
            C.warn("could not keep the radio open (%s) - using bladeRF-cli" % e)
    L.KEEP_RAW[0] = True                         # raw samples are needed for the IQ files
    C.stop_on_first_ctrl_c()
    track(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
