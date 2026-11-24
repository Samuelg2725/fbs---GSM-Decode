#!/usr/bin/env python3
"""Full-band GSM scanner for the BladeRF. Captures 925-960 MHz in wide chunks,
finds live channels, decodes them, prints a merged cell list.
  python3 bandscan.py                 # single gain 40, whole band
  python3 bandscan.py --gain-sweep    # gains 35,40,45 merged (max accuracy)
  python3 bandscan.py --gains 38,45 --secs 10 --log ~/scan.txt
"""
import sys, os, argparse, subprocess
import numpy as np
import pmt
from gnuradio import blocks, gr, gsm
from gnuradio import filter as grfilter
from gnuradio.filter import firdes

LIBS = "/usr/local/lib"
SC16 = "/dev/shm/_bandscan.sc16"
CFILE = "/dev/shm/_bandscan.cfile"
DECIM = 8
OPERATORS = {2: "O2", 10: "O2", 15: "Vodafone", 20: "Three", 30: "EE", 33: "EE"}

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

def main():
    ap = argparse.ArgumentParser(description="Wideband GSM band scanner (BladeRF)")
    ap.add_argument("--chunks", default="931,943,955", help="MHz chunk centres")
    ap.add_argument("--gains", default="40", help="comma list, e.g. 35,40,45")
    ap.add_argument("--gain-sweep", action="store_true", help="shortcut for 35,40,45")
    ap.add_argument("--rate", type=float, default=15.36e6)
    ap.add_argument("--secs", type=float, default=8.0)
    ap.add_argument("--snr", type=float, default=6.0)
    ap.add_argument("--log", default=None)
    args = ap.parse_args()
    gains = [35.0, 40.0, 45.0] if args.gain_sweep else [float(g) for g in args.gains.split(",") if g.strip()]
    chunks = [float(c) * 1e6 for c in args.chunks.split(",") if c.strip()]
    rate = args.rate
    chanrate = rate / DECIM
    logf = open(args.log, "a") if args.log else None
    def say(s):
        print(s, flush=True)
        if logf:
            logf.write(s + "\n"); logf.flush()
    def capture(center_hz, gain):
        n = int(rate * args.secs)
        cmd = ["bladeRF-cli", "-e", "set agc rx off",
               "-e", "set frequency rx %d" % int(center_hz),
               "-e", "set samplerate rx %d" % int(rate),
               "-e", "set bandwidth rx %d" % int(rate * 0.9),
               "-e", "set gain rx %d" % int(gain),
               "-e", "rx config file=%s format=bin n=%d samples=16384 buffers=64 xfers=16 timeout=20000" % (SC16, n),
               "-e", "rx start", "-e", "rx wait 60000"]
        env = dict(os.environ, LD_LIBRARY_PATH=LIBS)
        try:
            subprocess.run(cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=180)
        except Exception as exc:
            say("  capture error: %s" % exc); return False
        return os.path.exists(SC16) and os.path.getsize(SC16) > 100000
    def live_channels(iq, center_hz):
        Nf = 4096; N = iq.size; nblk = min(N // Nf, 4000)
        if nblk < 1:
            return []
        w = np.hanning(Nf); acc = np.zeros(Nf)
        for i in range(nblk):
            acc += np.abs(np.fft.fftshift(np.fft.fft(iq[i * Nf:(i + 1) * Nf] * w))) ** 2
        psd = acc / nblk; binhz = rate / Nf
        axis = center_hz + (np.arange(Nf) - Nf / 2) * binhz
        noise = np.median(psd); half = max(1, int(100e3 / binhz))
        lo = center_hz - rate / 2 + 0.3e6; hi = center_hz + rate / 2 - 0.3e6
        cands = []
        for n in list(range(1, 125)) + list(range(975, 1024)):
            f = gsm.arfcn.arfcn2downlink(n)
            if f is None or f < lo or f > hi:
                continue
            k = int(np.argmin(np.abs(axis - f)))
            p = psd[max(0, k - half):k + half].mean()
            snr = 10 * np.log10(p / noise + 1e-12)
            if snr >= args.snr:
                cands.append((n, f, snr))
        cands.sort(key=lambda x: -x[2]); return cands
    taps = firdes.low_pass(1.0, rate, 100e3, 30e3)
    def decode_channel(center_hz, f):
        tb = gr.top_block()
        src = blocks.file_source(gr.sizeof_gr_complex, CFILE, False)
        xlat = grfilter.freq_xlating_fir_filter_ccf(DECIM, taps, f - center_hz, rate)
        gi = gsm.gsm_input(ppm=0, osr=4, fc=f, samp_rate_in=chanrate)
        rcv = gsm.receiver(4, ([0]), ([]))
        clk = gsm.clock_offset_control(f, chanrate, osr=4)
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
        for i in range(dbg.num_messages()):
            info = si3(blob_of(dbg.get_message(i)))
            if info:
                return info
        return None
    say("band scan: chunks %s MHz | gains %s | %.1f Msps | %.0fs each"
        % (",".join("%.0f" % (c / 1e6) for c in chunks), ",".join("%g" % g for g in gains), rate / 1e6, args.secs))
    found = {}
    for gain in gains:
        for center in chunks:
            say("=== chunk %.1f MHz  gain %g ===" % (center / 1e6, gain))
            if not capture(center, gain):
                say("  (no capture - check radio/libbladeRF)"); continue
            r = np.fromfile(SC16, np.int16).astype(np.float32) / 2048.0
            iq = (r[0::2] + 1j * r[1::2]).astype(np.complex64)
            iq.tofile(CFILE)
            cands = live_channels(iq, center)
            say("  %d live channels" % len(cands))
            for n, f, snr in cands:
                if n in found:
                    continue
                info = decode_channel(center, f)
                if info:
                    cid, lac, mcc, mnc = info
                    op = OPERATORS.get(mnc, "?")
                    found[n] = (f, cid, lac, mcc, mnc, op)
                    say("  DECODED  ARFCN %d  %.1f MHz  CID %d  LAC %d  PLMN %d-%02d [%s]"
                        % (n, f / 1e6, cid, lac, mcc, mnc, op))
    say("")
    say("===== FINAL: %d unique cell(s) across 925-960 MHz =====" % len(found))
    for n in sorted(found):
        f, cid, lac, mcc, mnc, op = found[n]
        say("  ARFCN %3d  %.1f MHz  CID %-6d  LAC %-5d  PLMN %d-%02d  [%s]"
            % (n, f / 1e6, cid, lac, mcc, mnc, op))
    if logf:
        logf.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
