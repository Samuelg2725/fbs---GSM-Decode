#!/usr/bin/env python3
"""Two-antenna GSM on the BladeRF: RX1 and RX2 on the same channel at different gains, one decoder each.
A cell that overloads the high-gain antenna is still decoded by the low-gain one, and a weak cell that is
lost on the low-gain antenna is picked up by the high-gain one, so no single gain has to be guessed.

  Track one channel (both antennas):
      python3 gsm_dual.py --arfcn 13
      python3 gsm_dual.py --arfcn 13 --gain1 20 --gain2 45
  Scan channels or a band with both antennas:
      python3 gsm_dual.py --band 900
      python3 gsm_dual.py --arfcn 1 6 13 22 114 116 119
  Show both antennas' signal levels only (checks the gains / clipping):
      python3 gsm_dual.py --arfcn 13 --levels

Each decoded cell ID is tagged with which antenna decoded it: [RX1+RX2], [RX1 only] or [RX2 only].
The radio's two receive channels share one tuner (same frequency) but have independent gains.
Run with LD_LIBRARY_PATH=/usr/local/lib.
"""
import argparse
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                     # noqa: E402
import gsm_arfcn_fast as G             # noqa: E402
import pcap_out                        # noqa: E402
from gnuradio import blocks, gr, gsm   # noqa: E402

SAMPLE_RATE = 2e6


class DualReceiver(gr.top_block):
    """SoapySDR BladeRF source with two channels -> one gr-gsm decoder per channel + level meters."""

    def __init__(self, freq, gains, agcs, bandwidth):
        gr.top_block.__init__(self, "gsm_dual")
        from gnuradio import soapy
        self.src = soapy.source("driver=bladerf", "fc32", 2, "", "", ["", ""], ["", ""])
        for ch in (0, 1):
            self.src.set_sample_rate(ch, SAMPLE_RATE)
            self.src.set_frequency(ch, freq)
            self.src.set_bandwidth(ch, bandwidth)
            self.src.set_gain_mode(ch, bool(agcs[ch]))
            if not agcs[ch]:
                self.src.set_gain(ch, float(gains[ch]))
        self.dbg, self.sync, self.meters, self._keep = [], [], [], []
        for ch in (0, 1):
            gi = gsm.gsm_input(ppm=0, osr=4, fc=freq, samp_rate_in=SAMPLE_RATE)
            rx = gsm.receiver(4, ([0]), ([]))
            clk = gsm.clock_offset_control(freq, SAMPLE_RATE, osr=4)
            dm = gsm.gsm_bcch_ccch_demapper(0)
            dec = gsm.control_channels_decoder()
            dbg = blocks.message_debug()
            syn = blocks.message_debug()
            self.connect((self.src, ch), (gi, 0))
            self.connect((gi, 0), (rx, 0))
            self.msg_connect(clk, "ctrl", gi, "ctrl_in")
            self.msg_connect(rx, "C0", dm, "bursts")
            self.msg_connect(rx, "measurements", clk, "measurements")
            self.msg_connect(dm, "bursts", dec, "bursts")
            self.msg_connect(dec, "msgs", dbg, "store")
            self.msg_connect(rx, "C0", syn, "store")
            # level meter: mean power over 0.1 s, and peak magnitude over 2 ms blocks (clipping check)
            mag2 = blocks.complex_to_mag_squared(1)
            avg = blocks.moving_average_ff(int(SAMPLE_RATE * 0.1), 1.0 / (SAMPLE_RATE * 0.1), 4000, 1)
            p_avg = blocks.probe_signal_f()
            mag = blocks.complex_to_mag(1)
            vec = blocks.stream_to_vector(gr.sizeof_float, 4096)
            peak = blocks.max_ff(4096, 1)
            p_peak = blocks.probe_signal_f()
            self.connect((self.src, ch), mag2, avg, p_avg)
            self.connect((self.src, ch), mag, vec, peak, p_peak)
            self.dbg.append(dbg)
            self.sync.append(syn)
            self.meters.append((p_avg, p_peak))
            self._keep += [gi, rx, clk, dm, dec, mag2, avg, mag, vec, peak]

    def level(self, ch):
        """(power dBFS, peak 0..1.4) for channel ch."""
        p_avg, p_peak = self.meters[ch]
        p = p_avg.level()
        return (10 * math.log10(p) if p > 0 else -999.0), p_peak.level()

    def synced(self, ch):
        try:
            return self.sync[ch].num_messages() > 0
        except Exception:
            return False


def open_dual(args, freq):
    with C.quiet():
        tb = DualReceiver(freq, (args.gain1, args.gain2), (args.agc1, args.agc2), args.bw)
        tb.start()
    return tb


def close_dual(tb):
    if tb is None:
        return
    with C.quiet():
        try:
            tb.stop()
            tb.wait()
        except Exception:
            pass


def gain_name(n, agc, gain):
    return "RX%d %s" % (n, "AGC" if agc else "gain %g" % gain)


def who(chans):
    return {frozenset({0, 1}): ("RX1+RX2", "green"), frozenset({0}): ("RX1 only", "cyan"),
            frozenset({1}): ("RX2 only", "magenta")}[frozenset(chans)]


# -------------------------------------------------------------------------------------------- track
def track(args, arfcn):
    freq = G.downlink(arfcn)
    if freq is None:
        sys.exit("ARFCN %d is not a GSM downlink channel" % arfcn)
    names = (gain_name(1, args.agc1, args.gain1), gain_name(2, args.agc2, args.gain2))
    C.info("Two-antenna tracking of ARFCN %d (%.1f MHz): %s, %s, filter %.0f kHz. Ctrl-C to stop."
           % (arfcn, freq / 1e6, names[0], names[1], args.bw / 1e3))
    tb = open_dual(args, freq)
    t0 = time.time()
    seen, dec, events = [0, 0], [[], []], []      # events: [time, set(channels), cell]
    pending, last_print, last_msg, restarts = None, None, time.time(), 0
    next_level = t0

    def flush(force=False):
        nonlocal pending, last_print
        if pending and (force or time.time() - pending[0] > 0.3):
            t, chans, cell = pending
            gap = t - (last_print or t0)
            last_print = t
            tag, col = who(chans)
            print("%s  %s" % (G.cell_line(cell, gap), C.c("[%s]" % tag, col)), flush=True)
            pending = None
    try:
        while True:
            now = time.time()
            if args.levels and now >= next_level:
                (l1, p1), (l2, p2) = tb.level(0), tb.level(1)
                clip = lambda p: C.c("  CLIPPING", "red") if p >= 0.99 else ""
                print("%s  RX1 %6.1f dBFS peak %.2f%s   RX2 %6.1f dBFS peak %.2f%s   difference %.1f dB" % (
                    C.c(C.now_stamp(), "grey"), l1, p1, clip(p1), l2, p2, clip(p2), l2 - l1), flush=True)
                next_level = now + 1.0
            for ch in (0, 1):
                try:
                    n = tb.dbg[ch].num_messages()
                except Exception:
                    continue
                while seen[ch] < n:
                    msg = tb.dbg[ch].get_message(seen[ch])
                    seen[ch] += 1
                    last_msg = time.time()
                    blob = G.blob_of(msg)
                    if G.PCAP[0] is not None:
                        G.PCAP[0].gsm(blob, antenna=ch)
                    cell = G.si3_cell(blob, arfcn, freq)
                    if not cell:
                        continue
                    t = time.time()
                    dec[ch].append(t)
                    if events and t - events[-1][0] < 0.3 and ch not in events[-1][1]:
                        events[-1][1].add(ch)             # the other antenna's copy of the same broadcast
                    else:
                        flush(force=True)
                        events.append([t, {ch}, cell])
                        pending = [t, events[-1][1], cell]
            if not args.levels:
                flush()
            if time.time() - last_msg > 3.0:
                restarts += 1
                C.warn("%s  no decoder output on either antenna for 3 s - restarting the receiver" % C.now_stamp())
                close_dual(tb)
                tb = open_dual(args, freq)
                seen, last_msg = [0, 0], time.time()
            time.sleep(0.02)
    except KeyboardInterrupt:
        print("")
    finally:
        C.ignore_ctrl_c()
        flush(force=True)
        levels = [tb.level(0), tb.level(1)]
        close_dual(tb)
    track_summary(args, arfcn, freq, names, t0, dec, events, restarts, levels)
    return 0


def _missed(times):
    m, longest = 0, 0.0
    for a, b in zip(times, times[1:]):
        m += max(0, int(round((b - a) / G.SI3_PERIOD)) - 1)
        longest = max(longest, b - a)
    return m, longest


def track_summary(args, arfcn, freq, names, t0, dec, events, restarts, levels):
    run = time.time() - t0
    cells = {}
    for e in events:
        cells.setdefault(e[2]["cid"], e[2])
    t = C.Table(["Receiver", "Decoded", "Missed", "Caught", "Longest gap", "Level (dBFS)", "Peak"],
                align={k: "right" for k in ("Decoded", "Missed", "Caught", "Longest gap", "Level (dBFS)", "Peak")})
    rows = []
    for name, times, lv in ((names[0], dec[0], levels[0]), (names[1], dec[1], levels[1]),
                            ("Both combined", [e[0] for e in events], None)):
        m, lg = _missed(times)
        pct = 100.0 * len(times) / (len(times) + m) if times else 0.0
        rows.append((name, len(times), m, pct, lg))
        t.add([(name, "bold" if lv is None else ""), len(times), m, ("%.1f%%" % pct, C.rate_colour(pct)),
               "%.1f s" % lg, "%.1f" % lv[0] if lv else "", (("%.2f" % lv[1]) + (" CLIP" if lv[1] >= 0.99 else ""),
                                                            "red" if lv[1] >= 0.99 else "") if lv else ""])
    t.show("Two-antenna tracking summary - ARFCN %d (%.0f s, restarts %d)" % (arfcn, run, restarts))
    ct = C.Table(["Cell ID", "Operator", "MCC-MNC", "LAC", "Broadcasts", "Both", "RX1 only", "RX2 only"],
                 align={k: "right" for k in ("Cell ID", "LAC", "Broadcasts", "Both", "RX1 only", "RX2 only")})
    for cid, x in cells.items():
        ev = [e for e in events if e[2]["cid"] == cid]
        ct.add([(cid, "green+bold"), C.op_cell(x["operator"]), "%03d-%02d" % (x["mcc"], x["mnc"]), x["lac"], len(ev),
                sum(1 for e in ev if e[1] == {0, 1}), sum(1 for e in ev if e[1] == {0}),
                sum(1 for e in ev if e[1] == {1})])
    ct.show("Cells")
    if args.summary_csv:
        path = os.path.expanduser(args.summary_csv)
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        new = not os.path.exists(path)
        with open(path, "a") as fh:
            if new:
                fh.write("start_utc,label,arfcn,rx1,rx2,runtime_s,rx1_decoded,rx1_caught_pct,rx2_decoded,"
                         "rx2_caught_pct,both_decoded,both_caught_pct,restarts,cell_ids\n")
            fh.write("%s,%s,%d,%s,%s,%.1f,%d,%.1f,%d,%.1f,%d,%.1f,%d,%s\n" % (
                C.utc_iso(t0), args.label or "", arfcn, names[0], names[1], run, rows[0][1], rows[0][3],
                rows[1][1], rows[1][3], rows[2][1], rows[2][3], restarts, "/".join(str(k) for k in cells)))
        C.info("summary row added to %s" % path)


# --------------------------------------------------------------------------------------------- scan
def scan(args, channels):
    names = (gain_name(1, args.agc1, args.gain1), gain_name(2, args.agc2, args.gain2))
    C.info("Two-antenna scan of %d channel(s): %s, %s. The radio is reopened for every channel (~4 s each). "
           "Ctrl-C to stop." % (len(channels), names[0], names[1]))
    found, t_all, tb = {}, time.time(), None
    try:
        for i, arfcn in enumerate(channels, 1):
            freq = G.downlink(arfcn)
            if freq is None:
                continue
            C.status("  checking ARFCN %d (%.1f MHz)   channel %d of %d   %d cell(s) found so far"
                     % (arfcn, freq / 1e6, i, len(channels), len(found)))
            close_dual(tb)
            tb = open_dual(args, freq)
            t0, seen, got, cell = time.time(), [0, 0], set(), None
            while time.time() - t0 < args.dwell:
                for ch in (0, 1):
                    try:
                        n = tb.dbg[ch].num_messages()
                    except Exception:
                        continue
                    while seen[ch] < n:
                        blob = G.blob_of(tb.dbg[ch].get_message(seen[ch]))
                        if G.PCAP[0] is not None:
                            G.PCAP[0].gsm(blob, antenna=ch)
                        x = G.si3_cell(blob, arfcn, freq)
                        seen[ch] += 1
                        if x:
                            got.add(ch)
                            cell = cell or x
                if got == {0, 1} or (got and time.time() - t0 > 1.0):
                    break                    # both antennas have it, or one has and the other had 1 s more
                if not got and time.time() - t0 > args.sync_secs and not tb.synced(0) and not tb.synced(1):
                    break                    # no GSM sync on either antenna
                time.sleep(0.02)
            if cell or args.verbose:
                C.clear_status()
            if cell:
                tag, col = who(got)
                print("%s  %s" % (G.cell_line(cell, time.time() - t0), C.c("[%s]" % tag, col)), flush=True)
                prev = found.get(arfcn)
                found[arfcn] = (cell, set(got) | (prev[1] if prev else set()))
            elif args.verbose:
                print("%s          ARFCN %4d  %7.1f MHz  %s" % (C.c(C.now_stamp(), "grey"), arfcn, freq / 1e6,
                                                                C.c("no cell", "grey")), flush=True)
    except KeyboardInterrupt:
        C.clear_status()
        print("")
    finally:
        C.ignore_ctrl_c()
        close_dual(tb)
    t = C.Table(["ARFCN", "MHz", "Operator", "MCC-MNC", "LAC", "Cell ID", "Power", "Decoded by"],
                align={k: "right" for k in ("ARFCN", "MHz", "LAC", "Cell ID", "Power")})
    for a in sorted(found, key=lambda a: (a >= 512 and a <= 885, a < 900, a)):
        x, chans = found[a]
        t.add([a, "%.1f" % x["mhz"], C.op_cell(x["operator"]), "%03d-%02d" % (x["mcc"], x["mnc"]), x["lac"],
               (x["cid"], "green+bold"), x["power"] if x["power"] is not None else "-", who(chans)])
    t.show("GSM cells found with two antennas: %d (%.0f s)" % (len(found), time.time() - t_all))
    n1 = sum(1 for v in found.values() if 0 in v[1])
    n2 = sum(1 for v in found.values() if 1 in v[1])
    C.info("%s decoded %d, %s decoded %d" % (names[0], n1, names[1], n2))
    return 0


def main():
    ap = argparse.ArgumentParser(description="Two-antenna GSM tracking and scanning (BladeRF RX1 + RX2).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--band", choices=sorted(G.BANDS), help="scan every channel of a band")
    what.add_argument("--arfcn", type=int, nargs="+", help="one ARFCN = track it; several = scan them")
    ap.add_argument("--gain1", type=float, default=20.0, help="RX1 gain in dB (default 20, for strong cells)")
    ap.add_argument("--gain2", type=float, default=40.0, help="RX2 gain in dB (default 40, for weak cells)")
    ap.add_argument("--agc1", action="store_true", help="RX1 uses automatic gain control")
    ap.add_argument("--agc2", action="store_true", help="RX2 uses automatic gain control")
    ap.add_argument("--bw", type=float, default=2e6,
                    help="analog filter bandwidth in Hz, shared by both antennas (default 2e6, minimum 200e3)")
    ap.add_argument("--levels", action="store_true", help="track: print both antennas' signal levels every second")
    ap.add_argument("--dwell", type=float, default=2.5, help="scan: longest wait per channel in s (default 2.5)")
    ap.add_argument("--sync-secs", type=float, default=1.2,
                    help="scan: leave a channel after this long without GSM sync on either antenna (default 1.2)")
    ap.add_argument("--verbose", action="store_true", help="scan: also list channels without a cell")
    ap.add_argument("--label", default=None, help="track: name for this run in the summary CSV")
    ap.add_argument("--summary-csv", default="~/field/dual_summary.csv",
                    help="track: CSV file the summary row is added to (default ~/field/dual_summary.csv)")
    pcap_out.add_args(ap)
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    args.bw = max(200e3, float(args.bw))
    if args.no_color:
        C.set_colour(False)
    C.stop_on_first_ctrl_c()
    G.stop_other_scanners()
    G.PCAP[0] = pcap_out.open_from_args(args, "gsm_dual")
    try:
        if args.band:
            return scan(args, G.BANDS[args.band])
        if len(args.arfcn) == 1:
            return track(args, args.arfcn[0])
        return scan(args, args.arfcn)
    finally:
        pcap_out.report(G.PCAP[0])


if __name__ == "__main__":
    sys.exit(main())
