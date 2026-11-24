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

def other_scanners():
    me = os.getpid()
    found = []
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        pid = int(name)
        if pid == me:
            continue
        try:
            cmd = open("/proc/%d/cmdline" % pid, "rb").read().replace(b"\x00", b" ").decode(errors="replace")
        except OSError:
            continue
        if "gsm_arfcn" in cmd or "grgsm_scanner" in cmd or "gsm_sweep.py" in cmd:
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
            os.kill(pid, signal.SIGKIL)
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

class top(gr.top_block):
    def __init__(self, freq, gain=50):
        gr.top_block.__init__(self, "fast")
        samp_rate = 2e6
        try:
            dev_args = str(gsm.device.get_default_args("driver=lime,soapy=0,nchan=1"))
        except Exception:
            dev_args = "driver=lime,soapy=0,nchan=1"
        self.src = osmosdr.source(args="numchan=1 " + dev_args)
        self.src.set_sample_rate(samp_rate)
        self.src.set_center_freq(freq, 0)
        self.src.set_freq_corr(0, 0)
        self.src.set_gain_mode(True, 0)
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

def main():
    global df_cellid, df_arfcn
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--arfcn", type=int, nargs="+", required=True)
    parser.add_argument("--cellid", type=int)
    parser.add_argument("--df", action="store_true")
    args = parser.parse_args()
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
    print("opening ARFCN %s" % " ".join(str(a) for a in arfcns), flush=True)
    muted = clean.silence()
    tb = None
    try:
        tb = top(freqs[0])
        tb.start()
    finally:
        clean.restore(*muted)
    print("listening", flush=True)
    if args.df:
        kraken_tune(freqs[0])
        if not kraken_sync():
            print("kraken clock not read", flush=True)
    last_print = {a: time.time() for a in arfcns}
    try:
        if len(arfcns) == 1:
            arfcn, freq = arfcns[0], freqs[0]
            seen = 0
            last_msg = time.time()
            shown = False
            pending = []
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
                if time.time() - last_msg > 2:
                    if args.cellid is None:
                        print("%s  restarting radio" % clean.now_stamp(), flush=True)
                    tb.retune(freq)
                    last_msg = time.time()
                    try:
                        seen = tb.dbg.num_messages()
                    except Exception:
                        pass
                time.sleep(0.02)
        else:
            seen = 0
            while True:
                for arfcn, freq in zip(arfcns, freqs):
                    print("tune %d" % arfcn, flush=True)
                    tb.retune(freq)
                    time.sleep(0.1)
                    try:
                        seen = tb.dbg.num_messages()
                    except Exception:
                        pass
                    t0 = time.time()
                    text = None
                    while time.time() - t0 < 1.5:
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
                    else:
                        print("%s  %6.1fs  ARFCN %d  no cell" % (clean.now_stamp(), took, arfcn), flush=True)
    except KeyboardInterrupt:
        print("")
    finally:
        if tb is not None:
            try:
                tb.stop()
            except Exception:
                pass
    return 0

if __name__ == "__main__":
    sys.exit(main())