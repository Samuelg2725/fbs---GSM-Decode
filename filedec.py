#!/usr/bin/env python3
import sys, time
import numpy, pmt
from gnuradio import blocks, gr, gsm

CFILE = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "/tmp/blade13.cfile"
ARFCN = 13
CONJ = "--conj" in sys.argv
for a in sys.argv:
    if a.startswith("--arfcn="):
        ARFCN = int(a.split("=")[1])
freq = gsm.arfcn.arfcn2downlink(ARFCN)
samp_rate = 2e6

def py_to_bytes(py):
    if isinstance(py, (bytes, bytearray)):
        return bytes(py)
    if isinstance(py, tuple):
        for part in py:
            g = py_to_bytes(part)
            if g and len(g) >= 12:
                return g
    if isinstance(py, list) and py and isinstance(py[0], int):
        try:
            return bytes(bytearray(py))
        except Exception:
            return None
    return None

def blob_of(msg):
    for getter in (lambda: pmt.to_python(msg),
                   lambda: pmt.to_python(pmt.cdr(msg)),
                   lambda: pmt.blob_data(pmt.cdr(msg))):
        try:
            obj = getter()
        except Exception:
            continue
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

def si3(raw):
    if not raw:
        return None
    starts = []
    if len(raw) > 1 and 1 <= raw[1] <= 12:
        starts.append(raw[1] * 4)
    starts.extend((16, 0))
    for off in starts:
        if off + 12 <= len(raw) and raw[off + 2] == 0x1b:
            b = raw[off:]
            cid = (b[3] << 8) + b[4]
            if cid == 0:
                return None
            mcc = ((b[5] & 0xF) * 100) + (((b[5] & 0xF0) >> 4) * 10) + (b[6] & 0xF)
            mnc = ((b[7] & 0xF) * 10) + (b[7] >> 4)
            if ((b[6] & 0xF0) >> 4) < 10:
                mnc = mnc * 10 + ((b[6] & 0xF0) >> 4)
            lac = (b[8] << 8) + b[9]
            return "CID %d  LAC %d  MCC %d  MNC %d" % (cid, lac, mcc, mnc)
    return None

class dec(gr.top_block):
    def __init__(self):
        gr.top_block.__init__(self, "filedec")
        self.src = blocks.file_source(gr.sizeof_gr_complex, CFILE, True)
        self.thr = blocks.throttle(gr.sizeof_gr_complex, samp_rate * 4)
        self.gsm_input = gsm.gsm_input(ppm=0, osr=4, fc=freq, samp_rate_in=samp_rate)
        self.rcv = gsm.receiver(4, ([0]), ([]))
        self.clk = gsm.clock_offset_control(freq, samp_rate, osr=4)
        self.dem = gsm.gsm_bcch_ccch_demapper(0)
        self.dcd = gsm.control_channels_decoder()
        self.dbg = blocks.message_debug()
        tail = (self.src, 0)
        if CONJ:
            self.conj = blocks.conjugate_cc()
            self.connect(tail, (self.conj, 0)); tail = (self.conj, 0)
        self.connect(tail, (self.thr, 0))
        self.connect((self.thr, 0), (self.gsm_input, 0))
        self.connect((self.gsm_input, 0), (self.rcv, 0))
        self.msg_connect(self.clk, "ctrl", self.gsm_input, "ctrl_in")
        self.msg_connect(self.rcv, "C0", self.dem, "bursts")
        self.msg_connect(self.rcv, "measurements", self.clk, "measurements")
        self.msg_connect(self.dem, "bursts", self.dcd, "bursts")
        self.msg_connect(self.dcd, "msgs", self.dbg, "store")

print("decoding %s  arfcn %d (%.1f MHz)  conj=%s" % (CFILE, ARFCN, freq/1e6, CONJ))
tb = dec()
tb.start()
seen = 0
t0 = time.time()
hits = 0
try:
    while time.time() - t0 < 25:
        n = tb.dbg.num_messages()
        while seen < n:
            info = si3(blob_of(tb.dbg.get_message(seen)))
            seen += 1
            if info:
                hits += 1
                print("  DECODED  " + info, flush=True)
        time.sleep(0.2)
except KeyboardInterrupt:
    pass
tb.stop()
tb.wait()
print("total control-channel messages: %d   SI3 cells: %d" % (seen, hits))
print("RESULT:", "IQ DECODES - cells found" if hits else ("decoder ran, %d raw msgs" % seen if seen else "NO decode - 0 messages"))
