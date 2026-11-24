#!/usr/bin/env python3
# Capture raw IQ from the bladeRF to a GNU Radio .cfile (complex float32),
# so it can be decoded offline with grgsm_decode - exactly the format gr-gsm reads.
import sys, numpy as np, SoapySDR
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CF32

freq = float(sys.argv[1]) if len(sys.argv) > 1 else 937.6e6   # ARFCN 13 downlink
rate = 2e6
secs = 8.0
gain = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
out  = sys.argv[3] if len(sys.argv) > 3 else "/tmp/blade13.cfile"

dev = SoapySDR.Device("driver=bladerf")
dev.setSampleRate(SOAPY_SDR_RX, 0, rate)
dev.setFrequency(SOAPY_SDR_RX, 0, freq)
try:
    dev.setGainMode(SOAPY_SDR_RX, 0, False)
except Exception:
    pass
dev.setGain(SOAPY_SDR_RX, 0, gain)
st = dev.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CF32)
dev.activateStream(st)

total = int(rate * secs)
chunk = np.zeros(1 << 16, np.complex64)
f = open(out, "wb")
got = 0
stalls = 0
while got < total:
    sr = dev.readStream(st, [chunk], len(chunk), timeoutUs=1000000)
    if sr.ret > 0:
        chunk[:sr.ret].astype(np.complex64).tofile(f)
        got += sr.ret
        stalls = 0
    else:
        stalls += 1
        if stalls > 20:
            break
f.close()
dev.deactivateStream(st)
dev.closeStream(st)
pwr = None
try:
    a = np.fromfile(out, np.complex64, count=1 << 20)
    if a.size:
        pwr = 10 * np.log10(np.mean(np.abs(a) ** 2) + 1e-12)
except Exception:
    pass
print("wrote %d samples (complex float32, interleaved I/Q) to %s" % (got, out))
print("file is %.1f MB" % (got * 8 / 1e6))
if pwr is not None:
    print("mean power in capture: %.1f dB  (should be well above -60 if signal present)" % pwr)
