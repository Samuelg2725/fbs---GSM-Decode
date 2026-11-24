#!/usr/bin/env python3
import sys, numpy as np
rate = 2e6
lime = sys.argv[1] if len(sys.argv) > 1 else "/tmp/lime13.cfile"
blade = sys.argv[2] if len(sys.argv) > 2 else "/tmp/blade13.cfile"

def psd(fn, N=4096, maxs=8_000_000):
    a = np.fromfile(fn, np.complex64, count=maxs)
    if a.size < N:
        return None
    nch = a.size // N
    w = np.hanning(N)
    acc = np.zeros(N)
    for i in range(nch):
        acc += np.abs(np.fft.fftshift(np.fft.fft(a[i * N:(i + 1) * N] * w))) ** 2
    return acc / nch

L = psd(lime)
B = psd(blade)
if L is None or B is None:
    print("could not read one of the files"); sys.exit(1)
N = len(L)
binhz = rate / N

def peak_khz(P):
    return (int(np.argmax(P)) - N / 2) * binhz / 1e3

print("lime  strongest carrier at %+.1f kHz from centre" % peak_khz(L))
print("blade strongest carrier at %+.1f kHz from centre" % peak_khz(B))
Ln = (L - L.mean()) / (L.std() + 1e-12)
Bn = (B - B.mean()) / (B.std() + 1e-12)
xc = np.correlate(Bn, Ln, mode="full")
lag = int(np.argmax(xc)) - (N - 1)
xci = np.correlate(Bn, Ln[::-1], mode="full")
inverted = xci.max() > xc.max() * 1.05
print("\n--- how blade differs from lime ---")
print("frequency shift:  blade is %+.1f kHz vs lime" % (lag * binhz / 1e3))
print("spectrum match:   normal %.0f   inverted %.0f   -> %s"
      % (xc.max(), xci.max(), "INVERTED (mirror-flipped)" if inverted else "same orientation"))
