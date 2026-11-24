#!/usr/bin/env python3
import sys
import numpy as np
import pmt
from gnuradio import blocks, gr, gsm
from gnuradio import filter as grfilter
from gnuradio.filter import firdes

sc16 = sys.argv[1] if len(sys.argv) > 1 else "/dev/shm/wide.sc16"
fc_mhz = float(sys.argv[2]) if len(sys.argv) > 2 else 942.5
RATE = float(sys.argv[3]) if len(sys.argv) > 3 else 15.36e6
FC = fc_mhz * 1e6
DECIM = 8
CHANRATE = RATE / DECIM
CFILE = "/dev/shm/_wide.cfile"

print("loading %s ..." % sc16, flush=True)
r = np.fromfile(sc16, np.int16).astype(np.float32) / 2048.0
iq = (r[0::2] + 1j * r[1::2]).astype(np.complex64)
iq.tofile(CFILE)
N = iq.size
print("%d samples  (%.2f s at %.2f Msps, centre %.3f MHz)" % (N, N / RATE, RATE / 1e6, fc_mhz), flush=True)

Nf = 4096
nblk = min(N // Nf, 4000)
w = np.hanning(Nf)
acc = np.zeros(Nf)
for i in range(nblk):
    acc += np.abs(np.fft.fftshift(np.fft.fft(iq[i * Nf:(i + 1) * Nf] * w))) ** 2
psd = acc / max(nblk, 1)
binhz = RATE / Nf
axis = FC + (np.arange(Nf) - Nf / 2) * binhz
noise = np.median(psd)
half = max(1, int(100e3 / binhz))
lo = FC - RATE / 2 + 0.3e6
hi = FC + RATE / 2 - 0.3e6
cands = []
for n in list(range(1, 125)) + list(range(975, 1024)):
    f = gsm.arfcn.arfcn2downlink(n)
    if f is None or f < lo or f > hi:
        continue
    k = int(np.argmin(np.abs(axis - f)))
    p = psd[max(0, k - half):k + half].mean()
    snr = 10 * np.log10(p / noise + 1e-12)
    if snr > 6.0:
        cands.append((n, f, snr))
cands.sort(key=lambda x: -x[2])
print("live channels found: %d" % len(cands), flush=True)
for n, f, snr in cands:
    print("  ARFCN %d  %.1f MHz  +%.0f dB" % (n, f / 1e6, snr), flush=True)

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
            return (cid, lac, mcc, mnc)
    return None

taps = firdes.low_pass(1.0, RATE, 100e3, 30e3)
found = {}
for n, f, snr in cands:
    tb = gr.top_block()
    src = blocks.file_source(gr.sizeof_gr_complex, CFILE, False)
    xlat = grfilter.freq_xlating_fir_filter_ccf(DECIM, taps, f - FC, RATE)
    gi = gsm.gsm_input(ppm=0, osr=4, fc=f, samp_rate_in=CHANRATE)
    rcv = gsm.receiver(4, ([0]), ([]))
    clk = gsm.clock_offset_control(f, CHANRATE, osr=4)
    dem = gsm.gsm_bcch_ccch_demapper(0)
    dcd = gsm.control_channels_decoder()
    dbg = blocks.message_debug()
    tb.connect(src, xlat, gi, rcv)
    tb.msg_connect(clk, "ctrl", gi, "ctrl_in")
    tb.msg_connect(rcv, "C0", dem, "bursts")
    tb.msg_connect(rcv, "measurements", clk, "measurements")
    tb.msg_connect(dem, "bursts", dcd, "bursts")
    tb.msg_connect(dcd, "msgs", dbg, "store")
    tb.run()
    got = None
    for i in range(dbg.num_messages()):
        info = si3(blob_of(dbg.get_message(i)))
        if info:
            got = info
            break
    if got:
        cid, lac, mcc, mnc = got
        if cid not in found:
            found[cid] = (n, f, cid, lac, mcc, mnc)
            print("  DECODED  ARFCN %d  CID %d  LAC %d  MCC %d  MNC %d" % (n, cid, lac, mcc, mnc), flush=True)

print("\n===== %d cell(s) decoded =====" % len(found))
for cid in sorted(found):
    n, f, cid2, lac, mcc, mnc = found[cid]
    print("  ARFCN %d  %.1f MHz  CID %d  LAC %d  MCC %d  MNC %d" % (n, f / 1e6, cid, lac, mcc, mnc))
