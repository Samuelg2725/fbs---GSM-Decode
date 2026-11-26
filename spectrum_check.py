#!/usr/bin/env python3
"""Antenna check: how much signal the antennas bring in at each frequency - no decoding, just power.

  python3 spectrum_check.py                          # the default list of UK LTE carriers
  python3 spectrum_check.py --freq-mhz 806 1846.7 2662.9

For each frequency it records 0.1 s at 30.72 Msps on both antennas (RX1 + RX2) and shows the spectrum
(+-13 MHz around the carrier) and the power in the carrier's central 9 MHz.

Then it asks you to UNPLUG BOTH ANTENNAS and measures again with the same gains. With nothing connected the
radio only sees its own noise, so:
    antenna effect = power with antenna - power without antenna
  * 10 dB or more: the antenna brings in real signal at that frequency
  * under 3 dB:    nothing comes through the antenna there - only the radio's own noise is being measured
Use --no-compare to skip the second step. Run with LD_LIBRARY_PATH=/usr/local/lib.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                                   # noqa: E402
import lte_decode as L                               # noqa: E402

RATE = 30.72e6
DEFAULT = [(806.0, "800 MHz band 20 (EARFCN 6300)"), (952.5, "900 MHz band 8 (EARFCN 3725)"),
           (1821.3, "1800 MHz band 3 (EARFCN 1363)"), (1846.7, "1800 MHz band 3 (EARFCN 1617)"),
           (2132.3, "2100 MHz band 1 (EARFCN 223)"), (2630.0, "2600 MHz band 7 (EARFCN 2850)"),
           (2662.9, "2600 MHz band 7 (EARFCN 3179)")]
BARS = " ▁▂▃▄▅▆▇█"


def psd(raw, col):
    """Averaged power spectrum (dBFS per bin, 1024 bins, DC in the middle) of one antenna."""
    x = (raw[:, col].astype(np.float32) + 1j * raw[:, col + 1].astype(np.float32)) / 2048.0
    n = 1024
    seg = x[:len(x) // n * n].reshape(-1, n) * np.hanning(n)
    p = np.mean(np.abs(np.fft.fftshift(np.fft.fft(seg, axis=1), axes=1)) ** 2, axis=0) / n
    return 10 * np.log10(p + 1e-20)


def measure(freq_mhz, gain):
    """-> {"RX1": psd, "RX2": psd, "peak": peak} or None."""
    raw = L._cli_raw(freq_mhz * 1e6, RATE, 0.1, gain, True)
    if raw is None:
        return None
    raw = raw[int(RATE * 0.01):]
    return {"RX1": psd(raw, 0), "RX2": psd(raw, 2), "peak": float(np.max(np.abs(raw))) / 2048.0}


def band_power(p, half_mhz=4.5):
    """Mean power (dB) across the carrier's centre +-half_mhz, without the DC spike."""
    f = (np.arange(len(p)) - len(p) // 2) * RATE / len(p) / 1e6
    sel = (np.abs(f) <= half_mhz) & (np.abs(f) > 0.15)
    return 10 * np.log10(np.mean(10 ** (p[sel] / 10)))


def spectrum_line(p, lo, hi, cols=52):
    """One-line text spectrum from -13 to +13 MHz, scaled between lo and hi dB."""
    f = (np.arange(len(p)) - len(p) // 2) * RATE / len(p) / 1e6
    edges = np.linspace(-13, 13, cols + 1)
    out = ""
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (f >= a) & (f < b) & (np.abs(f) > 0.15)
        v = 10 * np.log10(np.mean(10 ** (p[sel] / 10))) if sel.any() else lo
        out += BARS[int(np.clip((v - lo) / max(1e-6, hi - lo) * (len(BARS) - 1), 0, len(BARS) - 1))]
    return out


def pick_gain(freq_mhz, start):
    """Highest gain (from start down) without clipping, so weak signal is not hidden below the noise."""
    g = start
    while True:
        m = measure(freq_mhz, g)
        if m is None or m["peak"] < 0.7 or g <= 0:
            return g, m
        g = max(0, g - 8)


def main():
    ap = argparse.ArgumentParser(description="Antenna check: signal power per frequency, with and without the "
                                             "antennas connected.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--freq-mhz", type=float, nargs="+", help="frequencies to check (default: UK LTE list)")
    ap.add_argument("--gain", type=float, default=60, help="highest gain to use, dB (default 60)")
    ap.add_argument("--no-compare", action="store_true", help="skip the measurement with the antennas unplugged")
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    if args.no_color:
        C.set_colour(False)
    L.clean_leftovers()
    freqs = [(f, "") for f in args.freq_mhz] if args.freq_mhz else DEFAULT
    C.heading("Step 1: antennas connected")
    rows = []
    for f, label in freqs:
        C.status("  measuring %.1f MHz..." % f)
        g, m = pick_gain(f, args.gain)
        if m is None:
            C.warn("  %.1f MHz: recording failed" % f)
            continue
        rows.append({"f": f, "label": label, "gain": g, "with": m})
    if not rows:
        return 1
    lo = min(min(np.median(r["with"]["RX1"]), np.median(r["with"]["RX2"])) for r in rows) - 3
    for r in rows:
        C.clear_status()
        print("  %7.1f MHz  gain %2g  %s" % (r["f"], r["gain"], C.c(r["label"], "grey")))
        for ant in ("RX1", "RX2"):
            p = r["with"][ant]
            print("     %s  -13 MHz %s +13 MHz   centre %.1f dB" % (
                ant, spectrum_line(p, lo, max(lo + 25, float(np.max(p)))), band_power(p)), flush=True)
    if not args.no_compare:
        C.heading("Step 2: antennas unplugged")
        print(C.c("  Unplug BOTH antennas from the BladeRF (leave the ports empty, or fit 50-ohm terminators),\n"
                  "  then press Enter.", "yellow+bold"), flush=True)
        try:
            input()
        except (EOFError, KeyboardInterrupt):
            print("")
            return 1
        for r in rows:
            C.status("  measuring %.1f MHz without antennas..." % r["f"])
            r["without"] = measure(r["f"], r["gain"])
        C.clear_status()
        print(C.c("  Plug the antennas back in.", "yellow"), flush=True)
    t = C.Table(["MHz", "Band", "Gain", "RX1 with", "RX1 without", "RX1 effect", "RX2 with", "RX2 without",
                 "RX2 effect", "Result"],
                align={k: "right" for k in ("MHz", "Gain", "RX1 with", "RX1 without", "RX1 effect", "RX2 with",
                                            "RX2 without", "RX2 effect")})
    for r in rows:
        cells, best = [], None
        for ant in ("RX1", "RX2"):
            w = band_power(r["with"][ant])
            if r.get("without"):
                wo = band_power(r["without"][ant])
                d = w - wo
                best = d if best is None else max(best, d)
                cells += ["%.1f" % w, "%.1f" % wo, ("%+.1f dB" % d, "green" if d >= 10 else
                                                    ("yellow" if d >= 3 else "red"))]
            else:
                cells += ["%.1f" % w, "-", "-"]
        if best is None:
            result = ("not compared", "grey")
        elif best >= 10:
            result = ("signal received", "green")
        elif best >= 3:
            result = ("weak signal", "yellow")
        else:
            result = ("nothing through antenna", "red+bold")
        t.add(["%.1f" % r["f"], r["label"].split(" (")[0] or "-", "%g" % r["gain"]] + cells + [result])
    t.show("Antenna check (power in the carrier's centre 9 MHz, dB relative to full scale)")
    C.info("effect = with antenna minus without: how much the antenna adds above the radio's own noise")
    return 0


if __name__ == "__main__":
    sys.exit(main())
