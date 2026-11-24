#!/usr/bin/env python3
import gc
import importlib.util
import os
import struct
import sys
import time
import pmt
from gnuradio import blocks, filter, gr, gsm
from gnuradio.filter import firdes
import osmosdr

spec = importlib.util.spec_from_file_location("gsm_clean", "/home/revector/gsm_arfcn_clean.py")
clean = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clean)

RATE = 40e6
CENTER = 942.5e6
SECS = 2.0
GAIN = 50.0
PATH = os.path.expanduser("~/gsm_wide.cfile")

def channels():
    found = []
    for n in range(975, 1024):
        freq = 935e6 + 0.2e6 * (n - 1024)
        if 925e6 <= freq <= 960e6:
            found.append((n, freq))
    for n in range(1, 125):
        freq = 935e6 + 0.2e6 * n
        if 925e6 <= freq <= 960e6:
            found.append((n, freq))
    return found

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

def capture():
    src = osmosdr.source(args="numchan=1 driver=lime,soapy=0,nchan=1")
    src.set_sample_rate(RATE)
    got = float(src.get_sample_rate())
    if got <= 0:
        got = RATE
    if abs(got - RATE) / RATE > 0.02:
        print("lime stayed at %.1f MS/s, wanted 40" % (got / 1e6), flush=True)
        return 0
    src.set_center_freq(CENTER, 0)
    src.set_freq_corr(0, 0)
    src.set_gain(GAIN, 0)
    src.set_bandwidth(RATE, 0)
    src.set_antenna("LNAL", 0)
    n = int(got * SECS)
    tb = gr.top_block()
    head = blocks.head(gr.sizeof_gr_complex, n)
    sink = blocks.file_sink(gr.sizeof_gr_complex, PATH, False)
    tb.connect(src, head, sink)
    print("capturing %.1fs at %.1f MS/s, centre %.1f MHz" % (SECS, got / 1e6, CENTER / 1e6), flush=True)
    started = time.time()
    tb.start()
    tb.wait()
    elapsed = time.time() - started
    size = os.path.getsize(PATH)
    print("wrote %d bytes in %.1fs" % (size, elapsed), flush=True)
    del tb
    if size < n * 8 * 0.9:
        print("capture is short", flush=True)
        return 0
    if elapsed > SECS * 1.8:
        print("lime did not keep up, the recording has gaps", flush=True)
    return size // 8

class one(gr.top_block):
    def __init__(self, n_in, offset, fc):
        gr.top_block.__init__(self, "wide")
        src = blocks.file_source(gr.sizeof_gr_complex, PATH, False)
        head = blocks.head(gr.sizeof_gr_complex, int(n_in))
        taps = firdes.low_pass(1.0, RATE, 500e3, 700e3)
        xlate = filter.freq_xlating_fir_filter_ccc(10, taps, offset, RATE)
        taps2 = firdes.low_pass(1.0, RATE / 10.0, 100e3, 50e3)
        low = filter.fir_filter_ccc(2, taps2)
        self.gsm_input = gsm.gsm_input(ppm=0, osr=4, fc=fc, samp_rate_in=RATE / 20.0)
        self.gsm_receiver = gsm.receiver(4, ([0]), ([]))
        self.gsm_clock = gsm.clock_offset_control(fc, RATE / 20.0, osr=4)
        self.demapper = gsm.gsm_bcch_ccch_demapper(0)
        self.decoder = gsm.control_channels_decoder()
        self.dbg = blocks.message_debug()
        self.connect(src, head, xlate, low, self.gsm_input, self.gsm_receiver)
        self.msg_connect(self.gsm_clock, "ctrl", self.gsm_input, "ctrl_in")
        self.msg_connect(self.gsm_receiver, "C0", self.demapper, "bursts")
        self.msg_connect(self.gsm_receiver, "measurements", self.gsm_clock, "measurements")
        self.msg_connect(self.demapper, "bursts", self.decoder, "bursts")
        self.msg_connect(self.decoder, "msgs", self.dbg, "store")

def decode(n_in, arfcn, freq):
    tb = one(n_in, freq - CENTER, freq)
    tb.start()
    tb.wait()
    time.sleep(0.1)
    seen = tb.dbg.num_messages()
    lines = []
    shown = set()
    for i in range(seen):
        text = si3_text(blob_of(tb.dbg.get_message(i)), arfcn, freq)
        if not text or text in shown:
            continue
        shown.add(text)
        lines.append(text)
    del tb
    gc.collect()
    return lines

def main():
    global PATH
    if len(sys.argv) > 1:
        PATH = os.path.abspath(os.path.expanduser(sys.argv[1]))
        if not os.path.isfile(PATH):
            print("no file", PATH)
            return 1
        n_in = os.path.getsize(PATH) // 8
        print("using", PATH, "samples", n_in, flush=True)
    else:
        n_in = capture()
        if not n_in:
            return 1
    chans = channels()
    print("%d channels in the %.0f MS/s capture, one after another on this Pi" % (len(chans), RATE / 1e6), flush=True)
    cells = 0
    for index, (arfcn, freq) in enumerate(chans, 1):
        print("%d/%d  ARFCN %d  %.1f MHz" % (index, len(chans), arfcn, freq / 1e6), flush=True)
        lines = decode(n_in, arfcn, freq)
        if not lines:
            print("  nothing on this channel", flush=True)
            continue
        for text in lines:
            cells += 1
            print(text, flush=True)
    print("cells", cells, flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
