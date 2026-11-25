#!/usr/bin/env python3
"""Two-antenna GSM tracking on the BladeRF: RX1 and RX2 on the SAME ARFCN at DIFFERENT gains,
one decoder each, results merged. Gain stops being a guess: a strong cell is caught by the
low-gain channel, a weak one by the high-gain channel.

  LD_LIBRARY_PATH=/usr/local/lib python3 ~/gsm_dual.py --arfcn 13
  LD_LIBRARY_PATH=/usr/local/lib python3 ~/gsm_dual.py --arfcn 9 --gain1 20 --gain2 50
  LD_LIBRARY_PATH=/usr/local/lib python3 ~/gsm_dual.py --arfcn 13 --agc1 --gain2 45
  LD_LIBRARY_PATH=/usr/local/lib python3 ~/gsm_dual.py --arfcn 13 --levels      # only show signal levels

Every broadcast prints once, tagged with who decoded it: [RX1+RX2], [RX1 only], [RX2 only].
The summary (every --stats-every s and on Ctrl-C, also a row in --summary-csv) gives decoded/missed
for RX1 alone, RX2 alone and both combined, so the second antenna's value is measured, not assumed.
level = signal level in dB below full scale per channel: near 0 dBFS = clipping (gain too high).

The two channels share one tuner (same frequency, same sample rate) but have separate gains.
--src picks how both channels are streamed (auto tries them in order until both deliver samples).
"""
import argparse
import os
import re
import signal
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.expanduser("~"))
import gsm_arfcn_fast as F            # noqa: E402  (decoder helpers + clean output module)
from gnuradio import blocks, gr, gsm  # noqa: E402

SAMP = 2e6
SOURCES = ["soapy", "osmo-soapy", "osmo-numchan", "osmo-nchan"]


BW = [SAMP]   # analog filter bandwidth (shared by RX1/RX2), set from --bw


def make_source(kind, freq, gains, agcs):
    """Return a block with 2 outputs (RX1, RX2), already configured."""
    if kind == "soapy":
        from gnuradio import soapy
        src = soapy.source("driver=bladerf", "fc32", 2, "", "", ["", ""], ["", ""])
        for ch in (0, 1):
            src.set_sample_rate(ch, SAMP)
            src.set_frequency(ch, freq)
            src.set_bandwidth(ch, BW[0])
            src.set_gain_mode(ch, bool(agcs[ch]))
            if not agcs[ch]:
                src.set_gain(ch, float(gains[ch]))
        return src
    import osmosdr
    args = {"osmo-soapy": "numchan=2 soapy=0,driver=bladerf",
            "osmo-numchan": "numchan=2 bladerf=0",
            "osmo-nchan": "bladerf=0,nchan=2"}[kind]
    src = osmosdr.source(args=args)
    if src.get_num_channels() < 2:
        raise RuntimeError("only %d channel(s)" % src.get_num_channels())
    src.set_sample_rate(SAMP)
    for ch in (0, 1):
        src.set_center_freq(freq, ch)
        src.set_freq_corr(0, ch)
        src.set_bandwidth(BW[0], ch)
        src.set_gain_mode(bool(agcs[ch]), ch)
        if not agcs[ch]:
            src.set_gain(float(gains[ch]), ch)
    return src


class dual(gr.top_block):
    def __init__(self, kind, freq, gains, agcs):
        gr.top_block.__init__(self, "dual")
        self.kind = kind
        self.src = make_source(kind, freq, gains, agcs)
        self._chains(freq)

    def _chains(self, freq):
        """Two fresh decoder chains (+ level meters) behind the source."""
        self.dbg, self.probe, self.keep, self.sync = [], [], [], []
        for ch in (0, 1):
            gi = gsm.gsm_input(ppm=0, osr=4, fc=freq, samp_rate_in=SAMP)
            rx = gsm.receiver(4, ([0]), ([]))
            clk = gsm.clock_offset_control(freq, SAMP, osr=4)
            dm = gsm.gsm_bcch_ccch_demapper(0)
            dec = gsm.control_channels_decoder()
            dbg = blocks.message_debug()
            self.connect((self.src, ch), (gi, 0))
            self.connect((gi, 0), (rx, 0))
            self.msg_connect(clk, "ctrl", gi, "ctrl_in")
            self.msg_connect(rx, "C0", dm, "bursts")
            self.msg_connect(rx, "measurements", clk, "measurements")
            self.msg_connect(dm, "bursts", dec, "bursts")
            self.msg_connect(dec, "msgs", dbg, "store")
            # the receiver only emits C0 bursts once it has locked to a cell's FCCH+SCH
            syn = blocks.message_debug()
            self.msg_connect(rx, "C0", syn, "store")
            self.sync.append(syn)
            # level meter: mean |x|^2 over ~0.1 s, and peak |x| (clipping check)
            mag = blocks.complex_to_mag_squared(1)
            avg = blocks.moving_average_ff(int(SAMP * 0.1), 1.0 / (SAMP * 0.1), 4000, 1)
            pr = blocks.probe_signal_f()
            pk = blocks.complex_to_mag(1)
            # peak |x| over blocks of 20000 samples (10 ms)
            st = blocks.stream_to_vector(gr.sizeof_float, 20000)
            mxv = blocks.max_ff(20000, 1)
            prp = blocks.probe_signal_f()
            self.connect((self.src, ch), mag, avg, pr)
            self.connect((self.src, ch), pk, st, mxv, prp)
            self.dbg.append(dbg)
            self.probe.append((pr, prp))
            self.keep += [gi, rx, clk, dm, dec, mag, avg, pk, st, mxv]

    def synced(self, ch):
        try:
            return self.sync[ch].num_messages() > 0
        except Exception:
            return False

    def retune(self, freq):
        """Keep the radio open; new frequency + brand-new decoders (no carried-over clock correction)."""
        try:
            self.stop()
            self.wait()
        except Exception:
            pass
        self.disconnect_all()
        for ch in (0, 1):
            if self.kind == "soapy":
                self.src.set_frequency(ch, freq)
            else:
                self.src.set_center_freq(freq, ch)
        self._chains(freq)
        self.start()

    def level(self, ch):
        import math
        pr, prp = self.probe[ch]
        p = pr.level()
        peak = prp.level()
        return (10 * math.log10(p) if p > 0 else -999.0), peak


def build(args, freq, check=True, quiet=False):
    kinds = SOURCES if args.src == "auto" else [args.src]
    gains, agcs = (args.gain1, args.gain2), (args.agc1, args.agc2)
    last_err = None
    for kind in kinds:
        muted = F.clean.silence() if not args.verbose else None
        tb = None
        try:
            tb = dual(kind, freq, gains, agcs)
            tb.start()
        except Exception as exc:
            last_err = "%s: %s" % (kind, exc)
            if muted:
                F.clean.restore(*muted)
            print("  source %-12s failed to open (%s)" % (kind, str(exc).splitlines()[0][:80]), flush=True)
            continue
        if muted:
            F.clean.restore(*muted)
        if not check:
            return tb, kind
        # both channels must actually deliver samples
        time.sleep(1.5)
        lv = [tb.level(0)[0], tb.level(1)[0]]
        if min(lv) <= -998:
            print("  source %-12s opened but channel(s) silent (levels %s) - trying next" % (
                kind, ", ".join("%.0f" % x for x in lv)), flush=True)
            try:
                tb.stop()
                tb.wait()
            except Exception:
                pass
            del tb
            continue
        return tb, kind
    sys.exit("no source could stream both channels. Last error: %s\n"
             "Try: --src soapy --verbose  to see the full error." % last_err)


def sweep(args):
    """Both antennas on each channel in turn; a channel counts as found if EITHER decodes its cell ID.
    The decoder is rebuilt per channel (same as gsm_arfcn_fast --fresh)."""
    BW[0] = max(200e3, float(args.bw))
    names = (tag("RX1", args.agc1, args.gain1), tag("RX2", args.agc2, args.gain2))
    chans = []
    for a in args.arfcn:
        f = gsm.arfcn.arfcn2downlink(a)
        if f is None or f < 100e6:
            sys.exit("ARFCN %d is not a downlink channel" % a)
        chans.append((a, f))
    print("dual sweep %d channels  %s  |  %s  |  dwell %.1fs  filter %.0f kHz" % (
        len(chans), names[0], names[1], args.dwell, BW[0] / 1e3), flush=True)
    for pid in F.other_scanners():
        try:
            os.kill(pid, signal.SIGINT)
        except OSError:
            pass
    signal.signal(signal.SIGINT, F._first_ctrl_c_only)
    found = {}          # arfcn -> [text, set(channels that ever decoded it)]
    tb, kind = None, None
    t_all = time.time()
    passes = 0
    try:
        while args.passes == 0 or passes < args.passes:
            passes += 1
            t_pass = time.time()
            for a, f in chans:
                if tb is not None and args.quick:
                    try:
                        tb.retune(f)
                    except Exception as exc:
                        print("  retune failed (%s) - reopening the radio" % str(exc)[:60], flush=True)
                        try:
                            tb.stop()
                            tb.wait()
                        except Exception:
                            pass
                        tb = None
                elif tb is not None:
                    try:
                        tb.stop()
                        tb.wait()
                    except Exception:
                        pass
                    del tb
                    tb = None
                if tb is None:
                    tb, kind = build(args, f, check=(kind is None))
                if kind and passes == 1 and a == chans[0][0]:
                    print("streaming via %s" % kind, flush=True)
                if kind and args.src == "auto":
                    args.src = kind              # keep the route that worked
                t0 = time.time()
                seen, got, text = [0, 0], set(), None
                nosync = False
                while time.time() - t0 < args.dwell:
                    for ch in (0, 1):
                        try:
                            n = tb.dbg[ch].num_messages()
                        except Exception:
                            continue
                        while seen[ch] < n:
                            msg = tb.dbg[ch].get_message(seen[ch])
                            seen[ch] += 1
                            t = F.si3_text(F.blob_of(msg), a, f)
                            if t:
                                got.add(ch)
                                text = text or t
                    if got == {0, 1} or (got and time.time() - t0 > 1.0):
                        break                     # both have it, or one has it and we gave the other 1 s
                    if (args.quick and not got and time.time() - t0 > args.sync_secs
                            and not tb.synced(0) and not tb.synced(1)):
                        nosync = True
                        break                     # no GSM sync on either antenna: nothing to wait for
                    time.sleep(0.02)
                took = time.time() - t0
                if got:
                    who = "RX1+RX2" if got == {0, 1} else ("RX1 only" if got == {0} else "RX2 only")
                    print("%s  %5.1fs  [%-8s] %s" % (F.clean.now_stamp(), took, who, text), flush=True)
                    plain = re.sub("\x1b\\[[0-9;]*m", "", text).split(", Pwr")[0]
                    found.setdefault(a, [plain, set()])[1].update(got)
                else:
                    print("%s  %5.1fs  ARFCN %d  %s" % (F.clean.now_stamp(), took, a,
                          "no GSM sync" if nosync else "no cell"), flush=True)
            print("pass %d done in %.1fs, %d cell(s) so far" % (passes, time.time() - t_pass, len(found)), flush=True)
    except KeyboardInterrupt:
        pass
    finally:
        try:
            signal.signal(signal.SIGINT, signal.SIG_IGN)
        except KeyboardInterrupt:
            pass
        if tb is not None:
            try:
                tb.stop()
            except Exception:
                pass
    r1 = sum(1 for v in found.values() if 0 in v[1])
    r2 = sum(1 for v in found.values() if 1 in v[1])
    only2 = sum(1 for v in found.values() if v[1] == {1})
    only1 = sum(1 for v in found.values() if v[1] == {0})
    print("\n===== %d unique cell(s) in %.0fs =====" % (len(found), time.time() - t_all))
    for a in sorted(found, key=lambda k: (k < 900 and k > 500, k < 900, k)):
        v = found[a]
        who = "RX1+RX2" if v[1] == {0, 1} else ("RX1 only" if v[1] == {0} else "RX2 only")
        print("  [%-8s] %s" % (who, v[0]))
    print("%s found %d   %s found %d   (RX1 only %d, RX2 only %d)" % (names[0], r1, names[1], r2, only1, only2))
    return 0


def tag(name, agc, gain):
    return "%s %s" % (name, "AGC" if agc else "g%g" % gain)


def main():
    ap = argparse.ArgumentParser(description="Two-antenna GSM tracking (BladeRF RX1+RX2, different gains).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--arfcn", type=int, nargs="+", required=True,
                    help="one ARFCN = track it; several = sweep them (both antennas on every channel)")
    ap.add_argument("--dwell", type=float, default=2.5, help="sweep: seconds per channel (default 2.5)")
    ap.add_argument("--passes", type=int, default=1, help="sweep: rounds (default 1, 0 = forever)")
    ap.add_argument("--quick", action="store_true",
                    help="sweep: keep the radio open between channels (fresh decoders each time) and move on "
                         "after --sync-secs if neither antenna gets GSM sync")
    ap.add_argument("--sync-secs", type=float, default=1.2,
                    help="--quick: seconds to wait for GSM sync before moving on (default 1.2)")
    ap.add_argument("--gain1", type=float, default=20.0, help="RX1 gain (default 20: strong cells)")
    ap.add_argument("--gain2", type=float, default=45.0, help="RX2 gain (default 45: weak cells)")
    ap.add_argument("--agc1", action="store_true", help="RX1 uses the radio's AGC instead of --gain1")
    ap.add_argument("--agc2", action="store_true", help="RX2 uses the radio's AGC instead of --gain2")
    ap.add_argument("--src", default="auto", choices=["auto"] + SOURCES)
    ap.add_argument("--bw", type=float, default=2e6,
                    help="analog filter bandwidth in Hz (default 2e6, min 200e3). Narrow it (e.g. 400e3) so a strong "
                         "neighbouring cell is filtered out before the converter and high gain can work")
    ap.add_argument("--levels", action="store_true", help="only print both channels' signal levels")
    ap.add_argument("--label", default=None)
    ap.add_argument("--stats-every", type=float, default=60.0)
    ap.add_argument("--summary-csv", default="~/field/dual_summary.csv")
    ap.add_argument("--verbose", action="store_true", help="show radio driver messages")
    args = ap.parse_args()

    BW[0] = max(200e3, float(args.bw))
    if len(args.arfcn) > 1:
        return sweep(args)
    args.arfcn = args.arfcn[0]
    freq = gsm.arfcn.arfcn2downlink(args.arfcn)
    if freq is None or freq < 100e6:
        sys.exit("ARFCN %d is not a downlink channel" % args.arfcn)
    names = (tag("RX1", args.agc1, args.gain1), tag("RX2", args.agc2, args.gain2))
    print("dual tracking ARFCN %d (%.1f MHz)  %s  |  %s  |  filter %.0f kHz" % (
        args.arfcn, freq / 1e6, names[0], names[1], BW[0] / 1e3), flush=True)
    for pid in F.other_scanners():
        try:
            os.kill(pid, signal.SIGINT)
        except OSError:
            pass
    signal.signal(signal.SIGINT, F._first_ctrl_c_only)
    tb, kind = build(args, freq)
    print("streaming via %s" % kind, flush=True)

    P = F.SI3_PERIOD
    seen = [0, 0]
    dec = [[], []]            # decode times per channel
    events = []               # merged broadcasts: [t_first, set(channels), text]
    pending = None
    last_msg = time.time()
    restarts = 0
    t0 = time.time()
    next_stats = t0 + args.stats_every if args.stats_every > 0 else None
    next_level = t0
    last_print = None

    def missed_of(times):
        m, longest = 0, 0.0
        for a, b in zip(times, times[1:]):
            m += max(0, int(round((b - a) / P)) - 1)
            longest = max(longest, b - a)
        return m, longest

    def summary(final=False):
        run = time.time() - t0
        u = [e[0] for e in events]
        rows = []
        for name, times in ((names[0], dec[0]), (names[1], dec[1]), ("COMBINED", u)):
            m, lg = missed_of(times)
            rate = 100.0 * len(times) / (len(times) + m) if times else 0.0
            rows.append((name, len(times), m, rate, lg))
        only1 = sum(1 for e in events if e[1] == {0})
        only2 = sum(1 for e in events if e[1] == {1})
        both = sum(1 for e in events if e[1] == {0, 1})
        lv = [tb.level(0), tb.level(1)]
        out = ["[%s] %.0fs  restarts %d  source %s" % (args.label, run, restarts, kind)]
        for name, n, m, rate, lg in rows:
            out.append("  %-12s decoded %4d  missed %4d  = %5.1f%% caught  longest gap %.1fs" % (name, n, m, rate, lg))
        out.append("  broadcasts: both %d  RX1 only %d  RX2 only %d   levels RX1 %.0f dBFS (peak %.2f)  RX2 %.0f dBFS (peak %.2f)"
                   % (both, only1, only2, lv[0][0], lv[0][1], lv[1][0], lv[1][1]))
        print(("\n===== dual tracking summary =====\n" if final else "") + "\n".join(out), flush=True)
        if final:
            path = os.path.expanduser(args.summary_csv)
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            new = not os.path.exists(path)
            with open(path, "a") as fh:
                if new:
                    fh.write("start_utc,label,arfcn,rx1,rx2,runtime_s,rx1_decoded,rx1_missed,rx1_pct,rx2_decoded,"
                             "rx2_missed,rx2_pct,comb_decoded,comb_missed,comb_pct,both,rx1_only,rx2_only,restarts,source\n")
                fh.write("%s,%s,%d,%s,%s,%.1f,%d,%d,%.1f,%d,%d,%.1f,%d,%d,%.1f,%d,%d,%d,%d,%s\n" % (
                    time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t0)), args.label or "", args.arfcn,
                    names[0], names[1], run, rows[0][1], rows[0][2], rows[0][3], rows[1][1], rows[1][2], rows[1][3],
                    rows[2][1], rows[2][2], rows[2][3], both, only1, only2, restarts, kind))
                print("  summary row added to %s" % path, flush=True)

    def flush_pending(force=False):
        nonlocal pending, last_print
        if pending and (force or time.time() - pending[0] > 0.3):
            t, chans, text = pending
            who = "RX1+RX2" if chans == {0, 1} else ("RX1 only" if chans == {0} else "RX2 only")
            gap = (t - last_print) if last_print else (t - t0)
            last_print = t
            print("%s  %5.1fs  [%-8s] %s" % (F.clean.now_stamp(), gap, who, text), flush=True)
            pending = None

    try:
        while True:
            now = time.time()
            if args.levels and now >= next_level:
                l1, l2 = tb.level(0), tb.level(1)
                clip = lambda l: "  CLIPPING" if l[1] >= 0.99 else ""
                print("%s  RX1 %6.1f dBFS peak %.2f%s   RX2 %6.1f dBFS peak %.2f%s   (difference %.1f dB)" % (
                    F.clean.now_stamp(), l1[0], l1[1], clip(l1), l2[0], l2[1], clip(l2), l2[0] - l1[0]), flush=True)
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
                    text = F.si3_text(F.blob_of(msg), args.arfcn, freq)
                    if not text:
                        continue
                    t = time.time()
                    dec[ch].append(t)
                    if events and t - events[-1][0] < 0.3 and ch not in events[-1][1]:
                        events[-1][1].add(ch)
                        if pending is not None and pending[0] == events[-1][0]:
                            pending[1].add(ch)
                    else:
                        flush_pending(force=True)
                        events.append([t, {ch}, text])
                        pending = [t, events[-1][1], text]
            if not args.levels:
                flush_pending()
            if next_stats and now >= next_stats:
                flush_pending(force=True)
                summary()
                next_stats = now + args.stats_every
            if time.time() - last_msg > 3.0:
                restarts += 1
                print("%s  no messages on either channel for 3 s - restarting radio" % F.clean.now_stamp(), flush=True)
                try:
                    tb.stop()
                    tb.wait()
                except Exception:
                    pass
                del tb
                tb, kind = build(args, freq)
                seen = [0, 0]
                last_msg = time.time()
            time.sleep(0.02)
    except KeyboardInterrupt:
        pass
    finally:
        try:
            signal.signal(signal.SIGINT, signal.SIG_IGN)
        except KeyboardInterrupt:
            pass
        flush_pending(force=True)
        try:
            summary(final=True)
        except BrokenPipeError:
            pass
        try:
            tb.stop()
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
