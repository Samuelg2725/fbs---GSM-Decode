#!/usr/bin/env python3
"""Join the IQ recordings of a tracking run (lte_track.py --save-iq) into one SigMF recording per sample rate.

  python3 sigmf_join.py ~/field/iq/earfcn76
  python3 sigmf_join.py ~/field/iq/earfcn76 --out ~/field/iq/earfcn76_joined

The recordings are NOT continuous: there is a gap of a few seconds between them while each one is processed.
The joined file keeps that information - SigMF "captures" mark where each recording starts (sample index) and
its real start time, and every PCI annotation is moved to its position in the joined file. Viewers such as
inspectrum show the recordings one after another.

Recordings at different sample rates (e.g. 15.36 Msps before the tracker learned the carrier was 15 MHz wide,
23.04 Msps after) cannot share a file: one joined file is written per sample rate.
Needs only Python 3 (no extra packages); works on the Pi or on a Mac.
"""
import argparse
import glob
import json
import os
import shutil
import sys

SAMPLE_BYTES = {"ci16_le": 4, "cf32_le": 8, "ci8": 2}


def main():
    ap = argparse.ArgumentParser(description="Join the SigMF recordings in a folder into one file per sample rate.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("folder", help="folder with the .sigmf-data / .sigmf-meta pairs")
    ap.add_argument("--out", default=None, help="output name without extension (default: <folder>_joined)")
    args = ap.parse_args()
    folder = os.path.expanduser(args.folder.rstrip("/"))
    metas = sorted(glob.glob(os.path.join(folder, "*.sigmf-meta")))
    if not metas:
        sys.exit("no .sigmf-meta files in %s" % folder)
    groups = {}
    for m in metas:
        with open(m) as fh:
            meta = json.load(fh)
        data = m[:-len(".sigmf-meta")] + ".sigmf-data"
        if not os.path.exists(data):
            print("skipping %s (no .sigmf-data)" % os.path.basename(m))
            continue
        g = meta["global"]
        key = (g["core:sample_rate"], g["core:datatype"], g.get("core:num_channels", 1))
        groups.setdefault(key, []).append((meta, data))
    base = os.path.expanduser(args.out) if args.out else folder + "_joined"
    for (rate, dtype, nch), items in sorted(groups.items()):
        items.sort(key=lambda it: it[0]["captures"][0].get("core:datetime", ""))
        name = base if len(groups) == 1 else "%s_%gMsps" % (base, rate / 1e6)
        per_sample = SAMPLE_BYTES.get(dtype, 4) * int(nch)
        joined = {"global": dict(items[0][0]["global"]), "captures": [], "annotations": []}
        joined["global"]["core:description"] = (
            joined["global"].get("core:description", "") +
            " Joined from %d recordings; they are not continuous - each capture segment is one recording, with "
            "its real start time." % len(items))
        offset = 0
        with open(name + ".sigmf-data", "wb") as out:
            for meta, data in items:
                cap = dict(meta["captures"][0])
                cap["core:sample_start"] = offset
                cap["lte:source_file"] = os.path.basename(data)
                joined["captures"].append(cap)
                for a in meta.get("annotations", []):
                    a = dict(a)
                    a["core:sample_start"] = a["core:sample_start"] + offset
                    joined["annotations"].append(a)
                with open(data, "rb") as fh:
                    shutil.copyfileobj(fh, out, 16 << 20)
                offset += os.path.getsize(data) // per_sample
        joined["annotations"].sort(key=lambda a: a["core:sample_start"])
        with open(name + ".sigmf-meta", "w") as fh:
            json.dump(joined, fh, indent=1)
        print("%s.sigmf-data: %d recordings, %.2f Msps, %.1f s of signal, %.0f MB" % (
            name, len(items), rate / 1e6, offset / rate, os.path.getsize(name + ".sigmf-data") / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
