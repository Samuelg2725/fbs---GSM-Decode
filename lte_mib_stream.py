#!/usr/bin/env python3
"""Every MIB of every cell on one LTE carrier, each with its exact time.

  python3 lte_mib_stream.py --earfcn 1363                 # all PCIs on EARFCN 1363, 1 s recording
  python3 lte_mib_stream.py --earfcn 1363 --secs 2 --loop # 2 s recordings, again and again (Ctrl-C to stop)
  python3 lte_mib_stream.py --earfcn 1363 --dual          # both antennas, combined per cell
  python3 lte_mib_stream.py --earfcn 1363 --pci 21 199    # only these PCIs (each searched for specifically,
                                                          # also when too weak for the blind search)

A cell sends its MIB (PBCH) in subframe 0 of every 10 ms radio frame. For each recording:
  1. every cell on the carrier is found and identified with lte_decode.py's decoder (cell search, weak-cell
     search, other sectors of each mast, interference cancellation; Cell ID from SIB1);
  2. each cell's signal is lined up on its own frame timing and the MIB of EVERY frame is demodulated on its own
     (srsRAN's PBCH decoder via the mib_stream helper; frames are not combined), so each line is one MIB that was
     actually received;
  3. each MIB gets its time: the moment the recording started + the sample position of that frame's subframe 0.
Timing: between MIBs, and between cells, exact to one sample (65 ns at 15.36 Msps). The absolute time comes from
the Pi's clock when the recording started (accurate to its NTP sync, typically a few ms); the recordings are not
continuous with --loop (the gaps are the processing time).

Output: one line per MIB (time to the microsecond, PCI, SFN, MIB bytes), the frame timing of each cell relative to
the first, ~/field/lte_mib_stream.csv, and every MIB in ~/field/pcap/lte_mib_<time>.pcapng with its exact time.
Run with LD_LIBRARY_PATH=/usr/local/lib.
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", str(max(1, (os.cpu_count() or 4) - 1)))
os.environ.setdefault("OMP_NUM_THREADS", os.environ["OPENBLAS_NUM_THREADS"])
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                                   # noqa: E402
import lte_decode as L                               # noqa: E402
import pcap_out                                      # noqa: E402

CSV = os.path.expanduser("~/field/lte_mib_stream.csv")


def helper():
    path = L.pdsch_ue().parent / "mib_stream"
    if not path.is_file():
        sys.exit("The mib_stream helper is not built (%s).\nRe-run the srsRAN patch and build it:\n"
                 "  python3 ~/fbs-upload/\"osprey 2\"/srsran_watch/print_sib1.py && cd ~/srsRAN_4G/build && "
                 "cmake .. && make -j3 mib_stream" % path)
    return str(path)


def decoder_args(args):
    """The settings lte_decode.py's decoder expects (its own defaults), plus this tool's."""
    return argparse.Namespace(
        gain=args.gain, gain_now=args.gain, gain_now2=args.gain, no_auto_gain=False, secs=args.secs, dual=args.dual,
        low_if_mhz=0.6, file=None, file_rate=None, search_secs=0.08, no_weak_search=False, pss_min=None,
        no_siblings=False, show_failed_siblings=False, only=bool(args.pci), no_cancel=False, retries=0,
        srsran_timeout=30, srsran_log=None, debug_search=False, always_decode=False, full_every=10, csv="",
        targets={args.earfcn: list(args.pci)} if args.pci else {}, _targets_now=list(args.pci or []),
        _earfcn_now=args.earfcn, _done_pcis=set(), loop=False)


def cell_signal(x1, x2, rate, cell):
    """The cell's signal from its subframe 0 on, frequency offset removed, both antennas combined, at its own
    standard rate -> (complex64 array, rate)."""
    off = float(cell.get("offset") or 0.0)
    if x2 is not None:
        w1, w2, _share = L.cell_weights(x1, x2, rate, cell["pci"], off)
        x = (w1 * x1 + w2 * x2).astype(np.complex64)
    else:
        x = x1
    sig = L.derotate(x, off, rate)[int(cell["start"]):]
    native = L.NATIVE_RATE.get(cell.get("mib") or 0, rate)
    k = rate / native
    if k > 1.01 and abs(k - round(k)) < 1e-6:
        from scipy.signal import resample_poly
        sig = resample_poly(sig, 1, int(round(k))).astype(np.complex64)
        return sig, native
    return sig.astype(np.complex64), rate


def without_others(x1, x2, rate, cell, cells):
    """The recording with every OTHER cell removed (interference cancellation, lte_cancel.py: their reference
    signals, sync signals and MIB channel rebuilt and subtracted, strongest first), so a cell under a stronger one
    gets its MIBs read without it."""
    from lte_cancel import cancel_crs
    own_rb = {7680000: 25, 15360000: 50, 23040000: 75, 30720000: 100}.get(int(round(rate)), 50)
    others = sorted([o for o in cells if o is not cell], key=lambda o: -(o.get("strength") or 0))
    for o in others:
        kw = dict(rate=rate, start=int(o["start"]), pci=o["pci"], n_rb=min(o.get("mib") or own_rb, own_rb),
                  ports=o.get("ports") or 2, cfo_hz=float(o.get("offset") or 0.0), mib_hex=o.get("mib_hex"))
        x1 = cancel_crs(x1, **kw)
        if x2 is not None:
            x2 = cancel_crs(x2, **kw)
    return x1, x2


def stream_cell(sig, rate, cell):
    """MIB of every frame -> [(frame index, sfn or None, hex, snr)]."""
    nframes = len(sig) // int(round(rate / 100))
    fh = tempfile.NamedTemporaryFile(prefix="mib_", suffix=".cf32", delete=False, dir="/dev/shm")
    fh.write(sig.tobytes())
    fh.close()
    try:
        out = subprocess.run([helper(), fh.name, str(cell["pci"]), str(cell["mib"]), str(cell.get("ports") or 2),
                              str(nframes)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=120,
                             preexec_fn=L._srsran_child).stdout.decode(errors="replace")
    finally:
        os.unlink(fh.name)
    rows = []
    for line in out.splitlines():
        p = dict(kv.split("=", 1) for kv in line.split()[2:] if "=" in kv)
        if not line.startswith("FRAME"):
            continue
        k = int(line.split()[1])
        if p.get("found") == "1":
            rows.append((k, int(p["sfn"]), p["hex"], float(p.get("snr", "nan"))))
        else:
            rows.append((k, None, None, float(p.get("snr", "nan"))))
    return rows


def utc_us(t):
    return time.strftime("%H:%M:%S", time.gmtime(t)) + ".%06d" % int(round((t % 1) * 1e6) % 1000000)


def one_recording(args, dargs, freq_mhz, pcap, csvfh, n_rec):
    rate = 15.36e6
    known = L.carrier_memory().get(str(args.earfcn), {})
    if known.get("width"):
        rate = L.carrier_rate(known["width"])
    if known.get("gain") is not None:
        dargs.gain_now = float(known["gain"])
        dargs.gain_now2 = float(known.get("gain2", known["gain"]))
    x1, x2 = L.auto_gain(dargs, freq_mhz, rate, *L.record(dargs, freq_mhz, rate))
    if x1 is None:
        C.warn("recording failed")
        return 0
    t0 = getattr(dargs, "_capture_t", time.time()) + 0.010   # the radio's first 10 ms are dropped
    cells = [c for c in L._decode_capture(dargs, freq_mhz, rate, x1, x2)
             if (c.get("mib") or c.get("ident")) and c.get("start") is not None]
    if args.pci:
        cells = [c for c in cells if c["pci"] in args.pci]
    if not cells:
        C.warn("no cell whose MIB decodes on EARFCN %d in this recording" % args.earfcn)
        return 0
    for c in cells:
        if not c.get("mib"):
            c["mib"] = known.get("width") or 50
        L.remember_carrier(args.earfcn, c["mib"], dargs.gain_now, dargs.gain_now2 if args.dual else None)
    C.heading("Recording %d   %s UTC   %.2f s   %d cell%s" % (n_rec, utc_us(t0), len(x1) / rate, len(cells),
                                                            "" if len(cells) == 1 else "s"))
    ref = min(cells, key=lambda c: c["start"])
    all_rows, stats = [], []
    for c in sorted(cells, key=lambda c: c["pci"]):
        C.status("  PCI %d: demodulating every MIB..." % c["pci"])
        sig, crate = cell_signal(x1, x2, rate, c)
        rows = stream_cell(sig, crate, c)
        n_ok = sum(1 for r in rows if r[1] is not None)
        if len(cells) > 1 and not args.no_cancel and rows and n_ok < 0.9 * len(rows):
            # MIBs missed: read them again with the other cells removed, keep whichever got more
            C.status("  PCI %d: %d of %d MIBs - again with the other %d cell%s removed..." % (
                c["pci"], n_ok, len(rows), len(cells) - 1, "" if len(cells) == 2 else "s"))
            y1, y2 = without_others(x1, x2, rate, c, cells)
            sig, crate = cell_signal(y1, y2, rate, c)
            rows2 = stream_cell(sig, crate, c)
            if sum(1 for r in rows2 if r[1] is not None) > n_ok:
                rows = rows2
        C.clear_status()
        x = c.get("ident") or {}
        label = "PCI %d%s" % (c["pci"], (" (Cell ID %d, %s)" % (x["cell_id"], x.get("operator") or "?")) if x else "")
        frame_s = (c["start"] - ref["start"]) / rate % 0.010
        ok = [r for r in rows if r[1] is not None]
        stats.append((c["pci"], label, len(ok), len(rows), frame_s, np.nanmean([r[3] for r in rows]) if rows else
                      float("nan")))
        for k, sfn, hx, snr in rows:
            t = t0 + c["start"] / rate + k * 0.010
            all_rows.append((t, c, k, sfn, hx, snr, label))
    all_rows.sort(key=lambda r: r[0])
    for t, c, k, sfn, hx, snr, label in all_rows:
        if sfn is None:
            if args.show_missed:
                print(C.c("  %s  PCI %-3d  frame %3d  no MIB (SNR %.1f dB)" % (utc_us(t), c["pci"], k, snr), "grey"))
            continue
        if not args.quiet:
            print("  %s  PCI %-3d  SFN %4d  MIB %s  SNR %5.1f dB" % (utc_us(t), c["pci"], sfn, hx, snr))
        if pcap is not None:
            pcap.lte_rrc(pcap_out.LTE_RRC_BCCH_BCH, bytes.fromhex(hx), earfcn=args.earfcn, pci=c["pci"], ts=t,
                         sfn=sfn, note="frame %d of the recording" % k)
        if csvfh is not None:
            csvfh.write("%s,%.6f,%d,%d,%d,%s,%.1f,%d\n" % (time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t)) +
                                                          ".%06d" % int(round((t % 1) * 1e6) % 1000000), t,
                                                          args.earfcn, c["pci"], sfn, hx, snr, n_rec))
    if pcap is not None:
        pcap.flush()
    if csvfh is not None:
        csvfh.flush()
    t = C.Table(["PCI", "Cell", "MIBs decoded", "Frame timing vs first cell", "Average SNR"],
                align={"PCI": "right", "MIBs decoded": "right", "Frame timing vs first cell": "right",
                       "Average SNR": "right"})
    for pci, label, ok, n, frame_s, snr in stats:
        t.add([pci, label.split(" ", 2)[2] if "(" in label else "-",
               ("%d of %d frames" % (ok, n), "green" if ok == n else ("yellow" if ok else "red")),
               "%+.3f us" % (frame_s * 1e6 if frame_s < 0.005 else (frame_s - 0.010) * 1e6), "%.1f dB" % snr])
    t.show("EARFCN %d: MIBs per cell (every frame demodulated on its own)" % args.earfcn)
    return sum(s[2] for s in stats)


def main():
    ap = argparse.ArgumentParser(description="Every MIB of every cell on one LTE carrier, with exact times.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--earfcn", type=int, required=True, help="the LTE carrier (EARFCN)")
    ap.add_argument("--pci", type=int, nargs="+", help="only these PCIs (default: every cell found)")
    ap.add_argument("--secs", type=float, default=1.0, help="seconds per recording (default 1 = 100 MIBs per cell)")
    ap.add_argument("--loop", action="store_true", help="record again and again (Ctrl-C to stop)")
    ap.add_argument("--dual", action="store_true", help="use both antennas, combined per cell")
    ap.add_argument("--gain", type=float, default=20.0, help="starting gain (adjusted automatically)")
    ap.add_argument("--quiet", action="store_true", help="no line per MIB, only the table per recording")
    ap.add_argument("--no-cancel", action="store_true",
                    help="read each cell's MIBs without first removing the other cells (faster, fewer MIBs)")
    ap.add_argument("--show-missed", action="store_true", help="also a line for frames whose MIB did not decode")
    pcap_out.add_args(ap)
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    if args.no_color:
        C.set_colour(False)
    helper()
    freq_mhz = L.earfcn_to_mhz(args.earfcn)
    if freq_mhz is None:
        sys.exit("EARFCN %d is not an LTE downlink channel" % args.earfcn)
    L.clean_leftovers()
    C.be_nice()
    import lte_fastsync
    lte_fastsync.install()
    try:
        L.RADIO[0] = L.LibRadio(args.dual)
    except Exception as e:
        C.warn("could not keep the radio open (%s) - using bladeRF-cli" % e)
    dargs = decoder_args(args)
    pcap = pcap_out.open_from_args(args, "lte_mib_%d" % args.earfcn)
    os.makedirs(os.path.dirname(CSV), exist_ok=True)
    new = not os.path.exists(CSV)
    csvfh = open(CSV, "a")
    if new:
        csvfh.write("utc,unix_time,earfcn,pci,sfn,mib_hex,snr_db,recording\n")
    C.info("EARFCN %d (%.1f MHz): every MIB of every cell, timed from the recording start (Pi clock) + sample position"
           % (args.earfcn, freq_mhz))
    C.stop_on_first_ctrl_c()
    n_rec, total = 0, 0
    try:
        while True:
            n_rec += 1
            total += one_recording(args, dargs, freq_mhz, pcap, csvfh, n_rec)
            if not args.loop:
                break
    except KeyboardInterrupt:
        print("")
        C.warn("stopped")
    finally:
        C.ignore_ctrl_c()
        csvfh.close()
    C.info("%d MIBs in %d recording%s; added to %s" % (total, n_rec, "" if n_rec == 1 else "s", CSV))
    pcap_out.report(pcap)
    return 0


if __name__ == "__main__":
    sys.exit(main())
