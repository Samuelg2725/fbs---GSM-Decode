#!/usr/bin/env python3
# Show the spectrum of a captured .cfile (complex float32) as an ASCII chart,
# so you can SEE the signal in the IQ. No decode, just the raw captured data.
import sys, numpy as np

cfile = sys.argv[1] if len(sys.argv) > 1 else "/tmp/blade13.cfile"
rate = 2e6
fc = float(sys.argv[2]) if len(sys.argv) > 2 else 937.6e6

a = np.fromfile(cfile, np.complex64)
print("file: %s" % cfile)
print("samples: %d   (%.2f s at %.1f Msps)" % (a.size, a.size / rate, rate / 1e6))
if a.size:
    print("overall power: %.1f dB" % (10 * np.log10(np.mean(np.abs(a) ** 2) + 1e-12)))

N = 2048
nch = a.size // N
if nch < 1:
    print("not enough samples")
    sys.exit(0)
nch = min(nch, 3000)
w = np.hanning(N)
acc = np.zeros(N)
for i in range(nch):
    seg = a[i * N:(i + 1) * N] * w
    acc += np.abs(np.fft.fftshift(np.fft.fft(seg))) ** 2
psd = 10 * np.log10(acc / nch + 1e-12)

bins = 60
step = N // bins
vals = [float(psd[i * step:(i + 1) * step].max()) for i in range(bins)]
freqs = [((i - bins / 2 + 0.5) * (rate / bins) + fc) / 1e6 for i in range(bins)]
floor = min(vals)
peak = max(vals)
print("\nspectrum across %.1f MHz (centre %.3f MHz).  floor %.0f dB  peak %.0f dB\n"
      % (rate / 1e6, fc / 1e6, floor, peak))
for f, v in zip(freqs, vals):
    bar = "#" * int(round(v - floor))
    mark = "  <-- strongest" if abs(v - peak) < 0.1 else ""
    print("%9.3f MHz | %s%s" % (f, bar, mark))
