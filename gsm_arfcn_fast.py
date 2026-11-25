#!/usr/bin/env python3
import argparse
import importlib.util
import os
import signal
import struct
import sys
import json
import time
import urllib.request

import numpy
import pmt
from gnuradio import blocks, gr, gsm
import osmosdr

spec = importlib.util.spec_from_file_location("gsm_clean", "/home/revector/gsm_arfcn_clean.py")
clean = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clean)

def my_ancestors():
    """This process and every parent above it (timeout, tee's shell, sudo...). Never kill these:
    'timeout 120 python3 gsm_arfcn_fast.py' has gsm_arfcn in its command line too."""
    keep = set()
    pid = os.getpid()
    while pid > 1 and pid not in keep:
        keep.add(pid)
        try:
            with open("/proc/%d/stat" % pid) as fh:
                pid = int(fh.read().rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            break
    return keep

def other_scanners():
    keep = my_ancestors()
    found = []
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        pid = int(name)
        if pid in keep:
            continue
        try:
            cmd = open("/proc/%d/cmdline" % pid, "rb").read().replace(b"\x00", b" ").decode(errors="replace")
        except OSError:
            continue
        if "python" in cmd and ("gsm_arfcn" in cmd or "grgsm_scanner" in cmd or "gsm_sweep.py" in cmd):
            found.append(pid)
    return found

def kill_other_scanners():
    pids = other_scanners()
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    if not pids:
        return []
    time.sleep(1)
    for pid in pids:
        try:
            os.kill(pid, 0)
        except OSError:
            continue
        try:
            os.kill(pid, signal.SIGKILL)
        except OSError:
            pass
    time.sleep(1)
    return pids

def py_to_bytes(py):
    if isinstance(py, (bytes, bytearray)):
        return bytes(py)
    if isinstance(py, tuple):
        for part in py:
            got = py_to_bytes(part)
            if got and len(got) >= 12:
                return got
    if isinstance(py, list) and py and isinstance(py[0], int):
        try:
            return bytes(bytearray(py))
        except Exception:
            return None
    return None

def blob_of(msg):
    objs = []
    try:
        objs.append(pmt.to_python(msg))
    except Exception:
        pass
    try:
        objs.append(pmt.to_python(pmt.cdr(msg)))
    except Exception:
        pass
    try:
        objs.append(pmt.blob_data(pmt.cdr(msg)))
    except Exception:
        pass
    for obj in objs:
        if obj is None:
            continue
        try:
            raw = bytes(bytearray(obj))
            if len(raw) >= 12:
                return raw
        except Exception:
            pass
        raw = py_to_bytes(obj)
        if raw and len(raw) >= 12:
            return raw
    return None

def si3_body(raw):
    if not raw: 
        return None
    starts = []
    if len(raw) > 1 and 1 <= raw[1] <= 12:
        starts.append(raw[1] * 4)
    starts.extend((16, 0))
    for off in starts:
        if off + 12 <= len(raw) and raw[off + 2] == 0x1b:
            return raw[off:]
    return None

def si3_text(raw, arfcn, freq):
    body = si3_body(raw)
    if body is None:
        return None
    cid = (body[3] << 8) + body[4]
    if cid == 0:
        return None
    mcc = ((body[5] & 0xF) * 100) + (((body[5] & 0xF0) >> 4) * 10) + (body[6] & 0xF)
    mnc = ((body[7] & 0xF) * 10) + (body[7] >> 4)
    if ((body[6] & 0xF0) >> 4) < 10:
        mnc = mnc * 10 + ((body[6] & 0xF0) >> 4)
    lac = (body[8] << 8) + body[9]
    try:
        pwr = struct.unpack_from("b", raw, 6)[0]
    except Exception:
        pwr = 0
    return clean.cell_text(arfcn, freq, cid, lac, mcc, mnc, pwr)

def radio_args(radio):
    """osmosdr device string. BladeRF: plain bladerf=0 - numchan/nchan BREAK it."""
    if radio == "bladerf":
        return "bladerf=0"
    try:
        dev_args = str(gsm.device.get_default_args("driver=lime,soapy=0,nchan=1"))
    except Exception:
        dev_args = "driver=lime,soapy=0,nchan=1"
    return "numchan=1 " + dev_args

class top(gr.top_block):
    def __init__(self, freq, gain=50, radio="lime", agc=False):
        gr.top_block.__init__(self, "fast")
        samp_rate = 2e6
        self.src = osmosdr.source(args=radio_args(radio))
        self.src.set_sample_rate(samp_rate)
        self.src.set_center_freq(freq, 0)
        self.src.set_freq_corr(0, 0)
        # Lime: AGC on (Owen's setting). BladeRF: manual gain (40 clean, ~45+ clips).
        # --agc: BladeRF uses the AD9361 hardware AGC instead (test: gain 40 overloads strong cells).
        self.src.set_gain_mode(radio != "bladerf" or agc, 0)
        if radio != "bladerf" or not agc:
            self.src.set_gain(gain, 0)
        self.src.set_bandwidth(samp_rate, 0)
        self.gsm_input = gsm.gsm_input(ppm=0, osr=4, fc=freq, samp_rate_in=samp_rate)
        self.gsm_receiver = gsm.receiver(4, ([0]), ([]))
        self.gsm_clock = gsm.clock_offset_control(freq, samp_rate, osr=4)
        self.demapper = gsm.gsm_bcch_ccch_demapper(0)
        self.decoder = gsm.control_channels_decoder()
        self.dbg = blocks.message_debug()
        self.connect((self.src, 0), (self.gsm_input, 0))
        self.connect((self.gsm_input, 0), (self.gsm_receiver, 0))
        self.msg_connect(self.gsm_clock, "ctrl", self.gsm_input, "ctrl_in")
        self.msg_connect(self.gsm_receiver, "C0", self.demapper, "bursts")
        self.msg_connect(self.gsm_receiver, "measurements", self.gsm_clock, "measurements")
        self.msg_connect(self.demapper, "bursts", self.decoder, "bursts")
        self.msg_connect(self.decoder, "msgs", self.dbg, "store")

    def retune(self, freq):
        try:
            self.stop()
        except Exception:
            pass
        time.sleep(0.2)
        self.src.set_center_freq(freq, 0)
        try:
            self.gsm_input.set_fc(freq)
        except Exception:
            pass
        self.start()


KRAKEN = "http://192.168.2.20:8090"
kraken_offset = 0.0
kraken_mhz = None
kraken_sync_at = 0.0
df_cellid = None
df_arfcn = None

def kraken_call(path, timeout):
    with urllib.request.urlopen(KRAKEN + path, timeout=timeout) as resp:
        return json.loads(resp.read().decode())

def kraken_sync():
    global kraken_offset, kraken_sync_at
    best = None
    for _ in range(5):
        t1 = time.time()
        try:
            remote = float(kraken_call("/time", 1.0)["time"])
        except Exception:
            return False
        t2 = time.time()
        delay = t2 - t1
        offset = remote - ((t1 + t2) / 2.0)
        if best is None or delay < best[0]:
            best = (delay, offset)
    kraken_offset = best[1]
    kraken_sync_at = time.time()
    return True

def kraken_tune(freq_hz):
    global kraken_mhz
    kraken_mhz = freq_hz / 1e6
    print("kraken tune %.3f MHz" % kraken_mhz, flush=True)
    try:
        info = kraken_call("/tune?mhz=%.6f" % kraken_mhz, 8)
    except Exception as exc:
        print("kraken: %s" % exc, flush=True)
        return
    if info.get("error"):
        print("kraken: %s" % info["error"], flush=True)
    else:
        print("kraken tuned %.3f MHz" % kraken_mhz, flush=True)

def emit_df(item):
    def dye(code, text):
        text = str(text)
        try:
            import sys
            if not sys.stdout.isatty():
                return text
        except Exception:
            return text
        return "\033[%sm%s\033[0m" % (code, text)
    try:
        info = kraken_call("/at?t=%.6f&mhz=%.6f" % (item["t"] + kraken_offset, kraken_mhz), 1)
    except Exception as exc:
        info = {"bearing": None, "error": str(exc)}
    bearing = info.get("bearing")
    if bearing is None:
        if info.get("reason") == "not calibrated":
            note, note_code = "not calibrated", "91"
        elif info.get("error"):
            note, note_code = "no link", "91"
        else:
            note, note_code = "unmatched", "93"
        angle, conf, tail = "  -", "    -", "%-14s" % note
        angle_code, conf_code, tail_code = "93", "93", note_code
    else:
        angle = "%3.0f" % float(bearing)
        conf = "%5.2f" % float(info.get("confidence") or 0)
        delta = info.get("delta")
        tail = "%-14s" % ("%.2fs" % float(delta) if isinstance(delta, (int, float)) else "-")
        angle_code, conf_code, tail_code = "93", "93", "96"
    extra = "Dir: %s, Conf: %s, Dt: %s" % (
        dye(angle_code, angle), dye(conf_code, conf), dye(tail_code, tail),
    )
    gap = dye("96", "%6.1fs" % item["took"])
    print("%s  %s  %s, %s" % (item["stamp"], gap, item["text"], extra), flush=True)
    if df_cellid is not None and item["cid"] != df_cellid:
        print("ARFCN %d expected %d, got %d" % (df_arfcn, df_cellid, item["cid"]), flush=True)

def flush_df(pending):
    now = time.time()
    if now - kraken_sync_at > 30:
        kraken_sync()
    while pending and now - pending[0]["t"] >= 1.5:
        emit_df(pending.pop(0))

SI3_PERIOD = 4 * 51 * 120.0 / 26.0 / 1000.0   # SI3 every 4 x 51-multiframes = 0.9415 s

def parse_cell(text):
    """Pull CID/LAC/MCC/MNC/Pwr out of a cell_text line (colour codes stripped)."""
    import re
    plain = re.sub("\x1b\\[[0-9;]*m", "", text)
    out = {}
    for key in ("CID", "LAC", "MCC", "MNC", "Pwr"):
        m = re.search(key + r":\s*(-?\d+)", plain)
        out[key] = int(m.group(1)) if m else None
    m = re.search(r"MNC:\s*-?\d+\s+([^,]+)", plain)
    out["op"] = m.group(1).strip() if m else ""
    return out

class TrackStats(object):
    """Decode-rate bookkeeping for single-channel tracking (--label / --log)."""
    def __init__(self, args, arfcn):
        self.args = args
        self.arfcn = arfcn
        self.t_listen = time.time()
        self.t_first = None
        self.t_last = None
        self.decodes = 0
        self.missed = 0
        self.longest = 0.0
        self.restarts = 0
        self.wrong = 0
        self.cids = {}
        self.pwr = []
        self.logf = None
        if args.log:
            path = os.path.expanduser(args.log)
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            new = not os.path.exists(path)
            self.logf = open(path, "a")
            if new:
                self.logf.write("utc,label,arfcn,cid,lac,mcc,mnc,operator,pwr,gap_s,missed_before,restarts\n")

    def decoded(self, text):
        now = time.time()
        c = parse_cell(text)
        gap = None
        missed_here = 0
        if self.t_last is None:
            self.t_first = now
        else:
            gap = now - self.t_last
            self.longest = max(self.longest, gap)
            missed_here = max(0, int(round(gap / SI3_PERIOD)) - 1)
            self.missed += missed_here
        self.t_last = now
        self.decodes += 1
        self.cids[c["CID"]] = self.cids.get(c["CID"], 0) + 1
        if self.args.cellid is not None and c["CID"] != self.args.cellid:
            self.wrong += 1
        if c["Pwr"] is not None:
            self.pwr.append(c["Pwr"])
        if self.logf:
            self.logf.write("%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%d,%d\n" % (
                clean.now_stamp(), self.args.label, self.arfcn, c["CID"], c["LAC"], c["MCC"], c["MNC"],
                c["op"], c["Pwr"], "" if gap is None else "%.3f" % gap, missed_here, self.restarts))
            self.logf.flush()

    def numbers(self):
        end = time.time()
        if self.t_first is None:
            expected = 0
        else:
            expected = int(round((self.t_last - self.t_first) / SI3_PERIOD)) + 1
        rate = 100.0 * self.decodes / expected if expected else 0.0
        return dict(runtime=end - self.t_listen, expected=expected, rate=rate,
                    first=(self.t_first - self.t_listen) if self.t_first else None)

    def line(self):
        n = self.numbers()
        pw = ("avg %.1f  min %d  max %d" % (sum(self.pwr) / len(self.pwr), min(self.pwr), max(self.pwr))
              if self.pwr else "-")
        return ("[%s] %.0fs  decoded %d / expected %d = %.1f%%  missed %d  longest gap %.1fs  "
                "restarts %d  first after %s  pwr %s%s"
                % (self.args.label, n["runtime"], self.decodes, n["expected"], n["rate"], self.missed,
                   self.longest, self.restarts, "-" if n["first"] is None else "%.1fs" % n["first"], pw,
                   ("  WRONG CID %d" % self.wrong) if self.wrong else ""))

    def finish(self):
        n = self.numbers()
        path = os.path.expanduser(self.args.summary_csv)
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        new = not os.path.exists(path)
        with open(path, "a") as fh:
            if new:
                fh.write("start_utc,label,arfcn,radio,gain,runtime_s,decoded,expected,rate_pct,missed,"
                         "longest_gap_s,restarts,first_decode_s,pwr_avg,pwr_min,pwr_max,cids,wrong_cid\n")
            fh.write("%s,%s,%d,%s,%s,%.1f,%d,%d,%.1f,%d,%.2f,%d,%s,%s,%s,%s,%s,%d\n" % (
                time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(self.t_listen)), self.args.label, self.arfcn,
                self.args.radio, "agc" if self.args.agc else "%g" % self.args.gain, n["runtime"], self.decodes, n["expected"], n["rate"],
                self.missed, self.longest, self.restarts,
                "" if n["first"] is None else "%.2f" % n["first"],
                "%.1f" % (sum(self.pwr) / len(self.pwr)) if self.pwr else "",
                min(self.pwr) if self.pwr else "", max(self.pwr) if self.pwr else "",
                "/".join(str(k) for k in self.cids), self.wrong))
        if self.logf:
            self.logf.close()
        # print after the CSV is safe on disk, so a closed pipe (tee killed by Ctrl-C) loses nothing
        try:
            print("")
            print("===== tracking summary =====", flush=True)
            print(self.line(), flush=True)
            if len(self.cids) > 1:
                print("  more than one CID seen: %s" % ", ".join("%s x%d" % kv for kv in self.cids.items()),
                      flush=True)
            print("  summary row added to %s" % path, flush=True)
        except BrokenPipeError:
            pass

def main():
    global df_cellid, df_arfcn
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--arfcn", type=int, nargs="+", required=True)
    parser.add_argument("--cellid", type=int)
    parser.add_argument("--df", action="store_true")
    parser.add_argument("--radio", choices=["lime", "bladerf"], default="lime",
                        help="lime (default, Owen's rig) or bladerf")
    parser.add_argument("--gain", type=float, default=None,
                        help="default 50 for lime, 40 for bladerf (home test: 40 overloads, 20 = 100%% decode)")
    parser.add_argument("--agc", action="store_true",
                        help="bladerf: use the radio's hardware AGC instead of the fixed --gain")
    parser.add_argument("--dwell", type=float, default=1.5,
                        help="multi-channel: seconds to wait for a cell ID per channel (default 1.5)")
    parser.add_argument("--passes", type=int, default=0,
                        help="multi-channel: stop after this many rounds (default 0 = forever)")
    parser.add_argument("--fresh", action="store_true",
                        help="rebuild the decoder on every retune/restart (no carried-over frequency "
                             "correction); test against the default when long sweeps miss cells")
    parser.add_argument("--label", default=None,
                        help="single channel: name this run; prints decode-rate stats and adds a row to "
                             "--summary-csv on Ctrl-C (e.g. --label gain40_roof)")
    parser.add_argument("--log", default=None,
                        help="single channel: CSV file with one row per decoded cell ID")
    parser.add_argument("--summary-csv", default="~/field/track_summary.csv",
                        help="where --label runs add their summary row (default ~/field/track_summary.csv)")
    parser.add_argument("--stats-every", type=float, default=60.0,
                        help="single channel with --label: print a stats line every N s (0 = only at the end)")
    args = parser.parse_args()
    if args.gain is None:
        args.gain = 40.0 if args.radio == "bladerf" else 50.0
    if args.df and len(args.arfcn) != 1:
        print("df uses one channel")
        return 1
    if args.df:
        df_cellid = args.cellid
        df_arfcn = args.arfcn[0]
    if args.cellid is not None and len(args.arfcn) != 1:
        print("use one channel with --cellid")
        return 1
    arfcns = args.arfcn
    cleared = kill_other_scanners()
    if cleared:
        print("cleared leftover scan %s" % " ".join(str(pid) for pid in cleared), flush=True)
    freqs = []
    for arfcn in arfcns:
        freq = gsm.arfcn.arfcn2downlink(arfcn)
        if freq is None or freq < 100e6:
            print("ARFCN %d is not a downlink channel" % arfcn)
            return 1
        freqs.append(freq)
    print("opening ARFCN %s  radio %s  gain %s" % (" ".join(str(a) for a in arfcns), args.radio,
          "AGC" if args.agc else "%g" % args.gain), flush=True)
    muted = clean.silence()
    tb = None
    try:
        tb = top(freqs[0], args.gain, args.radio, args.agc)
        tb.start()
    finally:
        clean.restore(*muted)
    print("listening", flush=True)
    if args.df:
        kraken_tune(freqs[0])
        if not kraken_sync():
            print("kraken clock not read", flush=True)
    last_print = {a: time.time() for a in arfcns}
    found_cells = {}

    def reopen(tb, freq):
        """--fresh: throw the whole decoder away and build a new one for this channel, so no
        frequency/clock correction learned on the previous (possibly empty) channel carries over.
        Default: Owen's retune (same decoder, new frequency)."""
        if not args.fresh:
            tb.retune(freq)
            return tb
        try:
            tb.stop()
            tb.wait()
        except Exception:
            pass
        del tb
        muted = clean.silence()
        try:
            new = top(freq, args.gain, args.radio, args.agc)
            new.start()
        finally:
            clean.restore(*muted)
        return new
    try:
        if len(arfcns) == 1:
            arfcn, freq = arfcns[0], freqs[0]
            seen = 0
            last_msg = time.time()
            shown = False
            pending = []
            stats = TrackStats(args, arfcn) if (args.label or args.log) else None
            if stats and not args.label:
                args.label = "-"
            next_stats = time.time() + args.stats_every if args.stats_every > 0 else None
            while True:
                if args.df:
                    flush_df(pending)
                n = tb.dbg.num_messages()
                if n > seen:
                    last_msg = time.time()
                while seen < n:
                    text = si3_text(blob_of(tb.dbg.get_message(seen)), arfcn, freq)
                    seen += 1
                    if not text:
                        continue
                    now = time.time()
                    took = now - last_print[arfcn]
                    last_print[arfcn] = now
                    cid = int(__import__("re").sub("\x1b\\[[0-9;]*m", "", text).split("CID:")[1].split(",")[0])
                    if stats:
                        stats.decoded(text)
                    if args.df:
                        pending.append({"t": time.time(), "stamp": clean.now_stamp(), "took": took, "text": text, "cid": cid})
                    elif args.cellid is None:
                        print("%s  %6.1fs  %s" % (clean.now_stamp(), took, text), flush=True)
                    elif cid == args.cellid:
                        if not shown:
                            print("%s  %6.1fs  %s" % (clean.now_stamp(), took, text), flush=True)
                            print("watching ARFCN %d for CID %d" % (arfcn, args.cellid), flush=True)
                            shown = True
                    else:
                        print("%s  %6.1fs  %s" % (clean.now_stamp(), took, text), flush=True)
                        print("ARFCN %d expected %d, got %d" % (arfcn, args.cellid, cid), flush=True)
                if stats and next_stats and time.time() >= next_stats:
                    print(stats.line(), flush=True)
                    next_stats = time.time() + args.stats_every
                if time.time() - last_msg > 2:
                    if stats:
                        stats.restarts += 1
                    if args.cellid is None:
                        print("%s  restarting radio" % clean.now_stamp(), flush=True)
                    tb = reopen(tb, freq)
                    last_msg = time.time()
                    try:
                        seen = tb.dbg.num_messages()
                    except Exception:
                        pass
                time.sleep(0.02)
        else:
            seen = 0
            done_passes = 0
            pass_start = time.time()
            while args.passes == 0 or done_passes < args.passes:
                for arfcn, freq in zip(arfcns, freqs):
                    print("tune %d" % arfcn, flush=True)
                    tb = reopen(tb, freq)
                    time.sleep(0.1)
                    try:
                        seen = tb.dbg.num_messages()
                    except Exception:
                        pass
                    t0 = time.time()
                    text = None
                    while time.time() - t0 < args.dwell:
                        try:
                            n = tb.dbg.num_messages()
                        except Exception:
                            break
                        while seen < n:
                            got = si3_text(blob_of(tb.dbg.get_message(seen)), arfcn, freq)
                            seen += 1
                            if got:
                                text = got
                                break
                        if text:
                            break
                        time.sleep(0.02)
                    now = time.time()
                    took = now - last_print[arfcn]
                    last_print[arfcn] = now
                    if text:
                        print("%s  %6.1fs  %s" % (clean.now_stamp(), took, text), flush=True)
                        found_cells[arfcn] = __import__("re").sub("\x1b\\[[0-9;]*m", "", text).split(", Pwr")[0]
                    else:
                        print("%s  %6.1fs  ARFCN %d  no cell" % (clean.now_stamp(), took, arfcn), flush=True)
                done_passes += 1
                print("pass %d done in %.1fs, %d cell(s) so far"
                      % (done_passes, time.time() - pass_start, len(found_cells)), flush=True)
                pass_start = time.time()
    except KeyboardInterrupt:
        pass
    finally:
        # Ctrl-C under `timeout` arrives twice (terminal + timeout forwards it); ignore the 2nd so cleanup finishes
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        if tb is not None:
            try:
                tb.stop()
            except Exception:
                pass
        if len(arfcns) == 1 and locals().get("stats"):
            stats.finish()
    if len(arfcns) > 1:
        try:
            print("")
            print("===== %d unique cell(s) =====" % len(found_cells))
            for a in sorted(found_cells, key=lambda k: (k < 900 and k > 500, k < 900, k)):
                print("  " + found_cells[a])
            sys.stdout.flush()
        except BrokenPipeError:
            pass
    try:
        sys.stdout.close()
    except BrokenPipeError:
        pass
    return 0

if __name__ == "__main__":
    sys.exit(main())