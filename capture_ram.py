#!/usr/bin/env python3
import sys, numpy as np, SoapySDR
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CF32, SOAPY_SDR_OVERFLOW
freq = 937.6e6; rate = 2e6; secs = 8.0; gain = 40.0
out = sys.argv[1] if len(sys.argv) > 1 else "/dev/shm/blade13.cfile"
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
total = int(rate * secs); chunk = np.zeros(1 << 16, np.complex64)
f = open(out, "wb"); got = 0; over = 0; err = 0
while got < total:
    sr = dev.readStream(st, [chunk], len(chunk), timeoutUs=1000000)
    if sr.ret > 0:
        chunk[:sr.ret].tofile(f); got += sr.ret
    elif sr.ret == SOAPY_SDR_OVERFLOW:
        over += 1
    else:
        err += 1
        if err > 50:
            break
f.close(); dev.deactivateStream(st); dev.closeStream(st)
print("wrote %d samples to %s" % (got, out))
print("overflow (dropped) events: %d    other errors: %d" % (over, err))
if over == 0 and err == 0:
    print(">>> CLEAN capture - no dropped samples")
else:
    print(">>> %d drop/error events - samples were lost (this breaks GSM timing)" % (over + err))
