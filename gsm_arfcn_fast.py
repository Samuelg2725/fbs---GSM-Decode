#!/usr/bin/env python3
"""GSM cell scanner and tracker for the BladeRF (gr-gsm decoder).

  Scan a whole band (finds every cell and its identity):
      python3 gsm_arfcn_fast.py --band 900
      python3 gsm_arfcn_fast.py --band 1800
  Scan chosen channels:
      python3 gsm_arfcn_fast.py --arfcn 1 6 13 22 114 116 119
  Track one channel continuously (prints every cell-ID broadcast, statistics on exit):
      python3 gsm_arfcn_fast.py --arfcn 13

Run with LD_LIBRARY_PATH=/usr/local/lib (libbladeRF). Stop with Ctrl-C at any time; the summary table is
always printed. The cell identity comes from System Information type 3, which a cell broadcasts every
0.94 s, so tracking can report a cell at most about once a second.
"""
import argparse
import os
import signal
import struct
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                     # noqa: E402
import pmt                             # noqa: E402
from gnuradio import blocks, gr, gsm   # noqa: E402
import osmosdr                         # noqa: E402

SAMPLE_RATE = 2e6
SI3_PERIOD = 4 * 51 * 120.0 / 26.0 / 1000.0          # one SI3 every 4 x 51-multiframes = 0.9415 s
BANDS = {
    "900": list(range(975, 1024)) + list(range(1, 125)),      # E-GSM 900: 925.2-959.8 MHz
    "1800": list(range(512, 886)),                              # DCS 1800: 1805.2-1879.8 MHz
}
BANDS["all"] = BANDS["900"] + BANDS["1800"]


# ------------------------------------------------------------------------------- other radio users
def _ancestors():
    """This process and its parents (timeout, shell, sudo...): never stopped by stop_other_scanners()."""
    keep, pid = set(), os.getpid()
    while pid > 1 and pid not in keep:
        keep.add(pid)
        try:
            with open("/proc/%d/stat" % pid) as fh:
                pid = int(fh.read().rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            break
    return keep


def stop_other_scanners():
    """Stop any other GSM tool still holding the BladeRF (e.g. left running in another terminal)."""
    keep, pids = _ancestors(), []
    for name in os.listdir("/proc"):
        if not name.isdigit() or int(name) in keep:
            continue
        try:
            cmd = open("/proc/%s/cmdline" % name, "rb").read().replace(b"\x00", b" ").decode(errors="replace")
        except OSError:
            continue
        if "python" in cmd and any(t in cmd for t in ("gsm_arfcn_fast", "gsm_dual", "gsm_neigh", "grgsm_scanner")):
            pids.append(int(name))
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    if pids:
        time.sleep(1.5)
        for pid in pids:
            try:
                os.kill(pid, signal.SIGKILL)
            except OSError:
                pass
    return pids


# ---------------------------------------------------------------------------------- message parsing
def _to_bytes(obj):
    if isinstance(obj, (bytes, bytearray)):
        return bytes(obj)
    if isinstance(obj, tuple):
        for part in obj:
            got = _to_bytes(part)
            if got and len(got) >= 12:
                return got
    if isinstance(obj, list) and obj and isinstance(obj[0], int):
        try:
            return bytes(bytearray(obj))
        except ValueError:
            return None
    return None


def blob_of(msg):
    """Raw bytes (GSMTAP header + L2 frame) of a gr-gsm decoder message."""
    candidates = []
    for getter in (lambda: pmt.to_python(msg), lambda: pmt.to_python(pmt.cdr(msg)),
                   lambda: pmt.blob_data(pmt.cdr(msg))):
        try:
            candidates.append(getter())
        except Exception:
            pass
    for obj in candidates:
        if obj is None:
            continue
        try:
            raw = bytes(bytearray(obj))
            if len(raw) >= 12:
                return raw
        except (TypeError, ValueError):
            pass
        raw = _to_bytes(obj)
        if raw and len(raw) >= 12:
            return raw
    return None


def _si3_body(raw):
    """The SI3 message (RR type 0x1b) inside a decoded frame, or None."""
    if not raw:
        return None
    starts = []
    if len(raw) > 1 and 1 <= raw[1] <= 12:
        starts.append(raw[1] * 4)                      # GSMTAP header length in 32-bit words
    starts.extend((16, 0))
    for off in starts:
        if off + 12 <= len(raw) and raw[off + 2] == 0x1b:
            return raw[off:]
    return None


def si3_cell(raw, arfcn, freq):
    """Cell identity from an SI3 frame -> dict(arfcn, mhz, cid, lac, mcc, mnc, operator, power) or None."""
    body = _si3_body(raw)
    if body is None:
        return None
    cid = (body[3] << 8) + body[4]
    if cid == 0:
        return None
    mcc = (body[5] & 0xF) * 100 + ((body[5] & 0xF0) >> 4) * 10 + (body[6] & 0xF)
    mnc = (body[7] & 0xF) * 10 + (body[7] >> 4)
    if ((body[6] & 0xF0) >> 4) < 10:                   # 3-digit MNC
        mnc = mnc * 10 + ((body[6] & 0xF0) >> 4)
    lac = (body[8] << 8) + body[9]
    try:
        power = struct.unpack_from("b", raw, 6)[0]      # GSMTAP signal level
    except struct.error:
        power = None
    return {"arfcn": arfcn, "mhz": freq / 1e6, "cid": cid, "lac": lac, "mcc": mcc, "mnc": mnc,
            "operator": C.operator(mcc, mnc), "power": power}


def cell_line(cell, gap=None):
    """One coloured console line for a decoded cell."""
    op = cell["operator"] or "?"
    return "%s  %s  ARFCN %4d  %7.1f MHz  %s  LAC %5d  %s  power %4s" % (
        C.c(C.now_stamp(), "grey"), C.c("%5.1fs" % gap, "grey") if gap is not None else "      ",
        cell["arfcn"], cell["mhz"], C.c("Cell ID %5d" % cell["cid"], "green+bold"), cell["lac"],
        C.c("%-9s %03d-%02d" % (op, cell["mcc"], cell["mnc"]), C.OPERATOR_COLOUR.get(op, "")),
        cell["power"] if cell["power"] is not None else "-")


# ----------------------------------------------------------------------------------------- receiver
class Receiver(gr.top_block):
    """BladeRF RX1 -> gr-gsm receiver -> BCCH/CCCH decoder, on one channel."""

    def __init__(self, freq, gain=20.0, agc=False, bandwidth=None):
        gr.top_block.__init__(self, "gsm_receiver")
        self.src = osmosdr.source(args="bladerf=0")
        self.src.set_sample_rate(SAMPLE_RATE)
        self.src.set_center_freq(freq, 0)
        self.src.set_freq_corr(0, 0)
        self.src.set_gain_mode(bool(agc), 0)
        if not agc:
            self.src.set_gain(gain, 0)
        self.src.set_bandwidth(bandwidth or SAMPLE_RATE, 0)
        self.gsm_input = gsm.gsm_input(ppm=0, osr=4, fc=freq, samp_rate_in=SAMPLE_RATE)
        self.receiver = gsm.receiver(4, ([0]), ([]))
        self.clock = gsm.clock_offset_control(freq, SAMPLE_RATE, osr=4)
        self.demapper = gsm.gsm_bcch_ccch_demapper(0)
        self.decoder = gsm.control_channels_decoder()
        self.dbg = blocks.message_debug()
        self.sync = blocks.message_debug()       # receiver emits C0 bursts only once locked (FCCH + SCH)
        self.connect((self.src, 0), (self.gsm_input, 0))
        self.connect((self.gsm_input, 0), (self.receiver, 0))
        self.msg_connect(self.clock, "ctrl", self.gsm_input, "ctrl_in")
        self.msg_connect(self.receiver, "C0", self.demapper, "bursts")
        self.msg_connect(self.receiver, "measurements", self.clock, "measurements")
        self.msg_connect(self.demapper, "bursts", self.decoder, "bursts")
        self.msg_connect(self.decoder, "msgs", self.dbg, "store")
        self.msg_connect(self.receiver, "C0", self.sync, "store")

    def synced(self):
        try:
            return self.sync.num_messages() > 0
        except Exception:
            return True


class Radio(object):
    """Opens a fresh Receiver per channel. The previous (stopped) flowgraph is kept until the new one is
    running, so gr-osmosdr re-uses the open BladeRF instead of closing and reopening it (~0.3 s vs ~2.8 s).
    A fresh decoder per channel matters: a decoder kept across many empty channels drifts in frequency and
    then misses even strong cells."""

    def __init__(self, args):
        self.args = args
        self.rx = None
        self._old = None

    def tune(self, freq):
        self.stop()
        with C.quiet():
            rx = Receiver(freq, self.args.gain, self.args.agc, self.args.bw)
            rx.start()
        self._old = self.rx
        self.rx = rx
        with C.quiet():
            self._old = None
        return rx

    def stop(self):
        if self.rx is not None:
            with C.quiet():
                try:
                    self.rx.stop()
                    self.rx.wait()
                except Exception:
                    pass

    def close(self):
        self.stop()
        with C.quiet():
            self.rx = None
            self._old = None


def downlink(arfcn):
    freq = gsm.arfcn.arfcn2downlink(arfcn)
    return freq if freq and freq > 100e6 else None


def new_messages(rx, seen):
    try:
        n = rx.dbg.num_messages()
    except Exception:
        return seen, []
    out = [rx.dbg.get_message(i) for i in range(seen, n)]
    return n, out


# --------------------------------------------------------------------------------------------- scan
def scan(args, channels):
    radio = Radio(args)
    found = {}                     # arfcn -> best cell dict (+ "seen" count)
    nosync = cellless = 0
    passes = 0
    t_all = time.time()
    C.info("Scanning %d channel(s), gain %s, up to %.1f s per channel (%.1f s if no GSM sync). Ctrl-C to stop."
           % (len(channels), "AGC" if args.agc else "%g" % args.gain, args.dwell, args.sync_secs))
    try:
        while args.passes == 0 or passes < args.passes:
            passes += 1
            t_pass = time.time()
            for arfcn in channels:
                freq = downlink(arfcn)
                if freq is None:
                    continue
                rx = radio.tune(freq)
                t0, seen, cell, synced = time.time(), 0, None, False
                while time.time() - t0 < args.dwell:
                    synced = synced or rx.synced()
                    if not synced and time.time() - t0 > args.sync_secs:
                        break
                    seen, msgs = new_messages(rx, seen)
                    for m in msgs:
                        cell = si3_cell(blob_of(m), arfcn, freq)
                        if cell:
                            break
                    if cell:
                        break
                    time.sleep(0.02)
                took = time.time() - t0
                if cell:
                    print("%s   (%.1fs)" % (cell_line(cell), took), flush=True)
                    prev = found.get(arfcn)
                    cell["seen"] = (prev["seen"] + 1) if prev else 1
                    found[arfcn] = cell
                elif synced:
                    cellless += 1
                    print("%s          ARFCN %4d  %7.1f MHz  %s   (%.1fs)" % (
                        C.c(C.now_stamp(), "grey"), arfcn, freq / 1e6, C.c("GSM signal, no cell ID decoded",
                                                                            "yellow"), took), flush=True)
                elif args.verbose:
                    nosync += 1
                    print("%s          ARFCN %4d  %7.1f MHz  %s" % (C.c(C.now_stamp(), "grey"), arfcn, freq / 1e6,
                                                                    C.c("no GSM signal", "grey")), flush=True)
                else:
                    nosync += 1
            C.info("pass %d finished in %.0f s - %d cell(s) so far" % (passes, time.time() - t_pass, len(found)))
    except KeyboardInterrupt:
        print("")
    finally:
        C.ignore_ctrl_c()
        radio.close()
    t = C.Table(["ARFCN", "MHz", "Operator", "MCC-MNC", "LAC", "Cell ID", "Power", "Seen"],
                align={"ARFCN": "right", "MHz": "right", "LAC": "right", "Cell ID": "right", "Power": "right",
                       "Seen": "right"})
    order = sorted(found, key=lambda a: (BANDS["1800"].count(a) > 0, a < 900, a))
    for a in order:
        x = found[a]
        t.add([a, "%.1f" % x["mhz"], C.op_cell(x["operator"]), "%03d-%02d" % (x["mcc"], x["mnc"]), x["lac"],
               (x["cid"], "green+bold"), x["power"] if x["power"] is not None else "-",
               "%d/%d" % (x["seen"], passes)])
    t.show("GSM cells found: %d   (%d pass(es), %.0f s)" % (len(found), passes, time.time() - t_all))
    if args.csv:
        write_scan_csv(args.csv, found)
    return 0


def write_scan_csv(path, found):
    path = os.path.expanduser(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new:
            fh.write("utc,arfcn,mhz,operator,mcc,mnc,lac,cell_id,power\n")
        for a, x in sorted(found.items()):
            fh.write("%s,%d,%.1f,%s,%d,%d,%d,%d,%s\n" % (C.utc_iso(), a, x["mhz"], x["operator"], x["mcc"], x["mnc"],
                                                          x["lac"], x["cid"], x["power"]))
    C.info("results added to %s" % path)


# -------------------------------------------------------------------------------------------- track
class TrackStats(object):
    """Decode-rate bookkeeping for single-channel tracking."""

    def __init__(self):
        self.t_start = time.time()
        self.t_first = self.t_last = None
        self.decoded = self.missed = self.restarts = 0
        self.longest = 0.0
        self.cells = {}            # cid -> cell dict
        self.power = []

    def add(self, cell):
        now = time.time()
        gap = None
        if self.t_last is None:
            self.t_first = now
        else:
            gap = now - self.t_last
            self.longest = max(self.longest, gap)
            self.missed += max(0, int(round(gap / SI3_PERIOD)) - 1)
        self.t_last = now
        self.decoded += 1
        self.cells[cell["cid"]] = cell
        if cell["power"] is not None:
            self.power.append(cell["power"])
        return gap

    def caught_pct(self):
        # decoded / (decoded + missed): some cells send their ID more often than every 0.94 s, so counting
        # gaps is exact where "expected = duration / 0.94 s" is not
        return 100.0 * self.decoded / (self.decoded + self.missed) if self.decoded else 0.0

    def line(self):
        return "%.0f s: %d decoded, %d missed (%.1f%% caught), longest gap %.1f s, restarts %d" % (
            time.time() - self.t_start, self.decoded, self.missed, self.caught_pct(), self.longest, self.restarts)


def track(args, arfcn):
    freq = downlink(arfcn)
    if freq is None:
        sys.exit("ARFCN %d is not a GSM downlink channel" % arfcn)
    radio = Radio(args)
    stats = TrackStats()
    C.info("Tracking ARFCN %d (%.1f MHz), gain %s. Ctrl-C to stop and show the summary."
           % (arfcn, freq / 1e6, "AGC" if args.agc else "%g" % args.gain))
    rx = radio.tune(freq)
    seen, last_msg = 0, time.time()
    next_stats = time.time() + args.stats_every if args.stats_every > 0 else None
    log = open_log(args.log)
    try:
        while True:
            seen_before = seen
            seen, msgs = new_messages(rx, seen)
            if seen > seen_before:
                last_msg = time.time()
            for m in msgs:
                cell = si3_cell(blob_of(m), arfcn, freq)
                if not cell:
                    continue
                gap = stats.add(cell)
                print(cell_line(cell, gap if gap is not None else time.time() - stats.t_start), flush=True)
                if log:
                    log.write("%s,%d,%d,%d,%d,%d,%s,%s\n" % (C.utc_iso(), arfcn, cell["cid"], cell["lac"], cell["mcc"],
                                                            cell["mnc"], cell["power"], "" if gap is None else "%.3f" % gap))
                    log.flush()
            if next_stats and time.time() >= next_stats:
                C.info(stats.line())
                next_stats = time.time() + args.stats_every
            if time.time() - last_msg > 2.0:
                stats.restarts += 1
                C.warn("%s  no decoder output for 2 s - restarting the receiver" % C.now_stamp())
                rx = radio.tune(freq)
                seen, last_msg = 0, time.time()
            time.sleep(0.02)
    except KeyboardInterrupt:
        print("")
    finally:
        C.ignore_ctrl_c()
        radio.close()
        if log:
            log.close()
    track_summary(args, arfcn, freq, stats)
    return 0


def open_log(path):
    if not path:
        return None
    path = os.path.expanduser(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    new = not os.path.exists(path)
    fh = open(path, "a")
    if new:
        fh.write("utc,arfcn,cell_id,lac,mcc,mnc,power,gap_s\n")
    return fh


def track_summary(args, arfcn, freq, s):
    run = time.time() - s.t_start
    pct = s.caught_pct()
    pw = ("%.1f (%d to %d)" % (sum(s.power) / len(s.power), min(s.power), max(s.power))) if s.power else "-"
    t = C.Table(["ARFCN", "MHz", "Operator", "MCC-MNC", "LAC", "Cell ID", "Decoded", "Missed", "Caught",
                 "Longest gap", "First after", "Restarts", "Power avg (min to max)"],
                align={k: "right" for k in ("ARFCN", "MHz", "LAC", "Cell ID", "Decoded", "Missed", "Caught",
                                            "Longest gap", "First after", "Restarts")})
    if not s.cells:
        t.add([arfcn, "%.1f" % (freq / 1e6), "-", "-", "-", ("none decoded", "red"), 0, 0, ("0%", "red"),
               "-", "-", s.restarts, "-"])
    for cid, x in s.cells.items():
        t.add([arfcn, "%.1f" % x["mhz"], C.op_cell(x["operator"]), "%03d-%02d" % (x["mcc"], x["mnc"]), x["lac"],
               (cid, "green+bold"), s.decoded, s.missed, ("%.1f%%" % pct, C.rate_colour(pct)),
               "%.1f s" % s.longest, "%.1f s" % (s.t_first - s.t_start) if s.t_first else "-", s.restarts, pw])
    t.show("GSM tracking summary   (%.0f s%s)" % (run, ", " + args.label if args.label else ""))
    if args.summary_csv:
        path = os.path.expanduser(args.summary_csv)
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        new = not os.path.exists(path)
        with open(path, "a") as fh:
            if new:
                fh.write("start_utc,label,arfcn,gain,runtime_s,decoded,missed,caught_pct,longest_gap_s,restarts,"
                         "first_decode_s,power_avg,cell_ids\n")
            fh.write("%s,%s,%d,%s,%.1f,%d,%d,%.1f,%.2f,%d,%s,%s,%s\n" % (
                C.utc_iso(s.t_start), args.label or "", arfcn, "agc" if args.agc else "%g" % args.gain, run,
                s.decoded, s.missed, pct, s.longest, s.restarts,
                "%.2f" % (s.t_first - s.t_start) if s.t_first else "",
                "%.1f" % (sum(s.power) / len(s.power)) if s.power else "", "/".join(str(k) for k in s.cells)))
        C.info("summary row added to %s" % path)


# --------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="GSM cell scanner and tracker (BladeRF + gr-gsm).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--band", choices=sorted(BANDS), help="scan every channel of a band (900, 1800 or all)")
    what.add_argument("--arfcn", type=int, nargs="+",
                      help="one ARFCN = track it; several = scan them")
    ap.add_argument("--gain", type=float, default=20.0,
                    help="receive gain in dB (default 20; 40 overloads the radio near strong cells)")
    ap.add_argument("--agc", action="store_true", help="use the radio's automatic gain control instead of --gain")
    ap.add_argument("--bw", type=float, default=None,
                    help="analog filter bandwidth in Hz (default 2e6). Narrow it (e.g. 300e3) when a strong "
                         "neighbouring channel overloads the radio")
    ap.add_argument("--dwell", type=float, default=2.5, help="scan: longest wait per channel in s (default 2.5)")
    ap.add_argument("--sync-secs", type=float, default=1.2,
                    help="scan: leave a channel after this long without GSM sync (default 1.2)")
    ap.add_argument("--passes", type=int, default=1, help="scan: number of passes (default 1, 0 = until Ctrl-C)")
    ap.add_argument("--verbose", action="store_true", help="scan: also list channels with no GSM signal")
    ap.add_argument("--csv", default=None, help="scan: append the cells found to this CSV file")
    ap.add_argument("--label", default=None, help="track: name for this run in the summary and CSV")
    ap.add_argument("--log", default=None, help="track: CSV file with one row per decoded cell ID")
    ap.add_argument("--summary-csv", default="~/field/track_summary.csv",
                    help="track: CSV file the summary row is added to (default ~/field/track_summary.csv)")
    ap.add_argument("--stats-every", type=float, default=60.0,
                    help="track: print a statistics line every N s (default 60, 0 = off)")
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    if args.no_color:
        C.set_colour(False)
    C.stop_on_first_ctrl_c()
    stopped = stop_other_scanners()
    if stopped:
        C.warn("stopped another scanner that was using the radio (pid %s)" % " ".join(map(str, stopped)))
    if args.band:
        return scan(args, BANDS[args.band])
    if len(args.arfcn) == 1:
        return track(args, args.arfcn[0])
    return scan(args, args.arfcn)


if __name__ == "__main__":
    sys.exit(main())
