#!/usr/bin/env python3
"""Turn the tracker's IQ recordings into files any LTE decoder can use, one per cell.

  python3 iq_export.py ~/field/iq/earfcn76                    # every recording in the folder
  python3 iq_export.py ~/field/iq/earfcn76/lte_76_..._c0019.sigmf-data --pci 21
  python3 iq_export.py ~/field/iq/earfcn76 --newest 3         # only the last 3 recordings

The tracker stores the radio's raw samples: 16-bit integers (12-bit, full scale 2048), recorded 0.6 MHz above the
carrier (to keep the radio's DC spike off the cell), starting at an arbitrary point. For each recording this writes:

  <name>_centred.cf32        the whole recording, shifted back so the carrier is at 0 Hz, DC spike removed,
                             complex float32 (I,Q interleaved, little-endian) at the recording's sample rate
  <name>_pci<P>.cf32         the same, starting exactly at subframe 0 of cell P (its frame timing from the
                             tracker) with its frequency offset removed, at the standard LTE rate for its width
                             (5 MHz 7.68, 10 MHz 15.36, 15 MHz 23.04, 20 MHz 30.72 Msps) - what srsRAN pdsch_ue
                             in file mode and most offline decoders expect
  each with a .sigmf-meta, and README.txt with the exact decode command for every file.

With two antennas (--dual recordings) RX1 is used; --antenna 2 picks RX2. Run on the Pi (needs the repository).
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lte_decode as L                               # noqa: E402  (signal conditioning and alignment helpers)

RB_OF_MHZ = {1.4: 6, 3: 15, 5: 25, 10: 50, 15: 75, 20: 100}


def load(meta_path, antenna):
    with open(meta_path) as fh:
        meta = json.load(fh)
    g = meta["global"]
    rate = float(g["core:sample_rate"])
    nch = int(g.get("core:num_channels", 1))
    raw = np.fromfile(meta_path[:-len(".sigmf-meta")] + ".sigmf-data", "<i2")
    raw = raw[:raw.size // (2 * nch) * 2 * nch].reshape(-1, 2 * nch)
    col = 2 * (min(antenna, nch) - 1)
    iq = ((raw[:, col].astype(np.float32) + 1j * raw[:, col + 1].astype(np.float32)) / 2048.0).astype(np.complex64)
    carrier = float(g.get("lte:carrier_hz") or meta["captures"][0]["core:frequency"])
    low_if = float(meta["captures"][0]["core:frequency"]) - carrier
    return meta, rate, carrier, L.condition(iq, rate, low_if)


def write_sigmf(path, x, rate, carrier, when, description, extra=None):
    x.astype(np.complex64).tofile(path)
    meta = {"global": {"core:datatype": "cf32_le", "core:sample_rate": rate, "core:version": "1.0.0",
                       "core:description": description},
            "captures": [{"core:sample_start": 0, "core:frequency": carrier, "core:datetime": when}],
            "annotations": []}
    if extra:
        meta["global"].update(extra)
    with open(path[:-len(".cf32")] + ".sigmf-meta", "w") as fh:
        json.dump(meta, fh, indent=1)


def main():
    ap = argparse.ArgumentParser(description="Export tracker IQ recordings as decoder-ready files.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("inputs", nargs="+", help="recordings (.sigmf-data or .sigmf-meta) or folders of them")
    ap.add_argument("--out", default="~/field/iq/export", help="output folder (default ~/field/iq/export)")
    ap.add_argument("--pci", type=int, nargs="+", help="only these cells")
    ap.add_argument("--newest", type=int, default=0, help="only the newest N recordings")
    ap.add_argument("--antenna", type=int, choices=(1, 2), default=1, help="antenna for two-antenna recordings")
    args = ap.parse_args()
    metas = []
    for p in args.inputs:
        p = os.path.expanduser(p)
        if os.path.isdir(p):
            metas += sorted(glob.glob(os.path.join(p, "*.sigmf-meta")))
        else:
            metas.append(p[:-len(".sigmf-data")] + ".sigmf-meta" if p.endswith(".sigmf-data") else p)
    metas = [m for m in metas if "_joined" not in os.path.basename(m)]
    if args.newest:
        metas = metas[-args.newest:]
    if not metas:
        sys.exit("no recordings found")
    out = os.path.expanduser(args.out)
    os.makedirs(out, exist_ok=True)
    readme = ["Decoder-ready LTE recordings", "=" * 28, "",
              "Format: complex float32, I and Q interleaved, little-endian (.cf32). The carrier is at 0 Hz.",
              "Each .sigmf-meta holds the sample rate, carrier frequency and time; *_pci*.cf32 files start at",
              "subframe 0 of that cell, at the standard LTE sample rate for its width.", "",
              "srsRAN 4G (file mode): pdsch_ue -i <file> -p <PRB> -c <PCI> -P <ports> -Q -d",
              "  PRB: 5 MHz = 25, 10 MHz = 50, 15 MHz = 75, 20 MHz = 100. Ports: try 2, then 4 or 1.",
              "  It prints 'Decoded MIB' and the PDCCH/PDSCH statistics; SIB1 needs a decoder that prints it.", ""]
    for mp in metas:
        meta, rate, carrier, x = load(mp, args.antenna)
        name = os.path.basename(mp)[:-len(".sigmf-meta")]
        when = meta["captures"][0].get("core:datetime", "")
        earfcn = meta["global"].get("lte:earfcn", "?")
        f = os.path.join(out, name + "_centred.cf32")
        write_sigmf(f, x, rate, carrier, when, "LTE EARFCN %s, centred on the carrier, %.2f Msps" % (earfcn, rate / 1e6),
                    {"lte:earfcn": earfcn})
        readme += ["%s  (EARFCN %s, %.1f MHz, %.2f Msps, %s)" % (name, earfcn, carrier / 1e6, rate / 1e6, when),
                   "  %s_centred.cf32   whole recording, not aligned (for decoders with their own cell search)" % name]
        print("%s_centred.cf32  %.2f Msps, %.1f s" % (name, rate / 1e6, len(x) / rate))
        for a in meta.get("annotations", []):
            pci = a.get("lte:pci")
            if pci is None or (args.pci and pci not in args.pci):
                continue
            width = a.get("lte:width_mhz")
            n_prb = RB_OF_MHZ.get(width, {7.68e6: 25, 15.36e6: 50, 23.04e6: 75, 30.72e6: 100}.get(rate, 50))
            native = L.NATIVE_RATE[n_prb]
            seg = L.derotate(x[int(a["core:sample_start"]):], float(a.get("lte:frequency_offset_hz") or 0.0), rate)
            k = rate / native
            if k > 1.01 and abs(k - round(k)) < 1e-6:
                from scipy.signal import resample_poly
                seg = resample_poly(seg, 1, int(round(k))).astype(np.complex64)
                r_out = native
            else:
                r_out = rate
            ident = ""
            if a.get("lte:cell_id") is not None:
                ident = "%s Cell ID %s (eNB %s, sector %s) TAC %s" % (a.get("lte:operator") or "?", a["lte:cell_id"],
                                                                      a.get("lte:enb"), a.get("lte:sector"),
                                                                      a.get("lte:tac"))
            ports = a.get("lte:ports") or 2
            fn = "%s_pci%d.cf32" % (name, pci)
            write_sigmf(os.path.join(out, fn), seg, r_out, carrier, when,
                        "LTE EARFCN %s PCI %d, starts at subframe 0, %.2f Msps. %s" % (earfcn, pci, r_out / 1e6, ident),
                        {"lte:earfcn": earfcn, "lte:pci": pci, "lte:prb": n_prb, "lte:ports": ports})
            readme.append("  %-40s PCI %d, %d PRB, %.2f Msps %s" % (fn, pci, n_prb, r_out / 1e6,
                                                                    ("- " + ident) if ident else "(identity not known)"))
            readme.append("      pdsch_ue -i %s -p %d -c %d -P %d -Q -d" % (fn, n_prb, pci, ports))
            print("  %s  PCI %d, %d PRB, %.2f Msps  %s" % (fn, pci, n_prb, r_out / 1e6, ident))
        readme.append("")
    with open(os.path.join(out, "README.txt"), "w") as fh:
        fh.write("\n".join(readme) + "\n")
    print("written to %s (see README.txt)" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
