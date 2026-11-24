#!/usr/bin/env python3
import gc
import importlib.util
import os
import struct
import subprocess
import sys
import threading
import time
import re
import pmt
from gnuradio import blocks, filter, gr, gsm
from gnuradio.filter import firdes
import osmosdr

spec = importlib.util.spec_from_file_location("gsm_clean", "/home/revector/gsm_arfcn_clean.py")
clean = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clean)

RATE = 40e6
CENTER = 942.5e6
SECS = 2.0
GAIN = 50.0
PATH = os.path.expanduser("~/gsm_wide.cfile")
RADIO = "lime"
FMT = "cfile"            # cfile = complex float32 (Lime/osmosdr), sc16 = int16 I/Q (bladeRF-cli)
BLADE_LIBS = "/usr/local/lib"
SKIP_EMPTY = False       # --skip-empty: short sync probe first, skip channels with no GSM sync
STOP_EARLY = False       # --stop-early: stop a channel as soon as its cell ID (SI3) is decoded
PROBE_SECS = 2.0         # --probe-secs: give up on a channel if no sync within this much capture
RETRIES = 0              # --retries: re-decode channels that synced but gave no cell ID
HAVE_PROGRESS = True     # set False if this GNU Radio can't report samples processed

def channels():
    """Every E-GSM900 channel inside the capture (Owen's list, clipped to what was recorded)."""
    lo = max(925e6, CENTER - RATE / 2 + 0.3e6)
    hi = min(960e6, CENTER + RATE / 2 - 0.3e6)
    found = []
    for n in range(975, 1024):
        freq = 935e6 + 0.2e6 * (n - 1024)
        if lo <= freq <= hi:
            found.append((n, freq))
    for n in range(1, 125):
        freq = 935e6 + 0.2e6 * n
        if lo <= freq <= hi:
            found.append((n, freq))
    return found

def py_to_bytes(py):
    if isinstance(py, (bytes, bytearray)):
        return bytes(py)
    if isinstance(py, tuple):
        for part in py:
            got = py_to_bytes(part)
            if got and len(got) >= 12:
                return got
    if isinstance(py, list) and py and isinstance(py[0], int):
        try:
            return bytes(bytearray(py))
        except Exception:
            return None
    return None

def blob_of(msg):
    objs = []
    try:
        objs.append(pmt.to_python(msg))
    except Exception:
        pass
    try:
        objs.append(pmt.to_python(pmt.cdr(msg)))
    except Exception:
        pass
    try:
        objs.append(pmt.blob_data(pmt.cdr(msg)))
    except Exception:
        pass
    for obj in objs:
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

def si3_body(raw):
    if not raw:
        return None
    starts = []
    if len(raw) > 1 and 1 <= raw[1] <= 12:
        starts.append(raw[1] * 4)
    starts.extend((16, 0))
    for off in starts:
        if off + 12 <= len(raw) and raw[off + 2] == 0x1b:
            return raw[off:]
    return None

def si3_text(raw, arfcn, freq):
    body = si3_body(raw)
    if body is None:
        return None
    cid = (body[3] << 8) + body[4]
    if cid == 0:
        return None
    mcc = ((body[5] & 0xF) * 100) + (((body[5] & 0xF0) >> 4) * 10) + (body[6] & 0xF)
    mnc = ((body[7] & 0xF) * 10) + (body[7] >> 4)
    if ((body[6] & 0xF0) >> 4) < 10:
        mnc = mnc * 10 + ((body[6] & 0xF0) >> 4)
    lac = (body[8] << 8) + body[9]
    try:
        pwr = struct.unpack_from("b", raw, 6)[0]
    except Exception:
        pwr = 0
    return clean.cell_text(arfcn, freq, cid, lac, mcc, mnc, pwr)

def capture_bladerf():
    """BladeRF: capture natively with bladeRF-cli (osmosdr drops samples at high rates).
    Writes int16 I/Q (sc16) - to /dev/shm by default, the SD card is too slow."""
    n = int(RATE * SECS)
    if os.path.exists(PATH):
        os.remove(PATH)
    cmd = ["bladeRF-cli", "-e", "set agc rx off",
           "-e", "set frequency rx %d" % int(CENTER),
           "-e", "set samplerate rx %d" % int(RATE),
           "-e", "set bandwidth rx %d" % int(min(RATE * 0.9, 56e6)),
           "-e", "set gain rx %d" % int(GAIN),
           "-e", "rx config file=%s format=bin n=%d samples=16384 buffers=64 xfers=16 timeout=20000" % (PATH, n),
           "-e", "rx start", "-e", "rx wait %d" % int(SECS * 3000 + 15000)]
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = BLADE_LIBS + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    print("capturing %.1fs at %.1f MS/s, centre %.1f MHz, gain %g (bladeRF-cli)"
          % (SECS, RATE / 1e6, CENTER / 1e6, GAIN), flush=True)
    started = time.time()
    p = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    elapsed = time.time() - started
    size = os.path.getsize(PATH) if os.path.exists(PATH) else 0
    print("wrote %d bytes in %.1fs" % (size, elapsed), flush=True)
    if size < n * 4 * 0.9:
        print("capture is short - bladeRF-cli said:", flush=True)
        print(p.stdout.decode(errors="replace")[-800:], flush=True)
        return 0
    return size // 4

def capture():
    if RADIO == "bladerf":
        return capture_bladerf()
    src = osmosdr.source(args="numchan=1 driver=lime,soapy=0,nchan=1")
    src.set_sample_rate(RATE)
    got = float(src.get_sample_rate())
    if got <= 0:
        got = RATE
    if abs(got - RATE) / RATE > 0.02:
        print("lime stayed at %.1f MS/s, wanted 40" % (got / 1e6), flush=True)
        return 0
    src.set_center_freq(CENTER, 0)
    src.set_freq_corr(0, 0)
    src.set_gain(GAIN, 0)
    src.set_bandwidth(RATE, 0)
    src.set_antenna("LNAL", 0)
    n = int(got * SECS)
    tb = gr.top_block()
    head = blocks.head(gr.sizeof_gr_complex, n)
    sink = blocks.file_sink(gr.sizeof_gr_complex, PATH, False)
    tb.connect(src, head, sink)
    print("capturing %.1fs at %.1f MS/s, centre %.1f MHz" % (SECS, got / 1e6, CENTER / 1e6), flush=True)
    started = time.time()
    tb.start()
    tb.wait()
    elapsed = time.time() - started
    size = os.path.getsize(PATH)
    print("wrote %d bytes in %.1fs" % (size, elapsed), flush=True)
    del tb
    if size < n * 8 * 0.9:
        print("capture is short", flush=True)
        return 0
    if elapsed > SECS * 1.8:
        print("lime did not keep up, the recording has gaps", flush=True)
    return size // 8

class one(gr.top_block):
    def __init__(self, n_in, offset, fc, count_bursts=False):
        gr.top_block.__init__(self, "wide")
        if FMT == "sc16":
            src = blocks.file_source(gr.sizeof_short, PATH, False)
            head = blocks.head(gr.sizeof_short, 2 * int(n_in))
            self.items_per_sample = 2
            to_c = blocks.interleaved_short_to_complex(False, False)
            scale = blocks.multiply_const_cc(complex(1.0 / 2048.0, 0.0))
            front = [src, head, to_c, scale]
        else:
            src = blocks.file_source(gr.sizeof_gr_complex, PATH, False)
            head = blocks.head(gr.sizeof_gr_complex, int(n_in))
            self.items_per_sample = 1
            front = [src, head]
        self.head = head
        # Owen's two stages: /10 then /2 (40 MS/s -> 2 MS/s). Other rates scale stage 1.
        d1 = 10 if abs(RATE - 40e6) < 1 else max(1, int(RATE // 4e6))
        taps = firdes.low_pass(1.0, RATE, 500e3, 700e3)
        xlate = filter.freq_xlating_fir_filter_ccc(d1, taps, offset, RATE)
        taps2 = firdes.low_pass(1.0, RATE / d1, 100e3, 50e3)
        low = filter.fir_filter_ccc(2, taps2)
        out_rate = RATE / d1 / 2.0
        self.gsm_input = gsm.gsm_input(ppm=0, osr=4, fc=fc, samp_rate_in=out_rate)
        self.gsm_receiver = gsm.receiver(4, ([0]), ([]))
        self.gsm_clock = gsm.clock_offset_control(fc, out_rate, osr=4)
        self.demapper = gsm.gsm_bcch_ccch_demapper(0)
        self.decoder = gsm.control_channels_decoder()
        self.dbg = blocks.message_debug()
        self.connect(*(front + [xlate, low, self.gsm_input, self.gsm_receiver]))
        self.msg_connect(self.gsm_clock, "ctrl", self.gsm_input, "ctrl_in")
        self.msg_connect(self.gsm_receiver, "C0", self.demapper, "bursts")
        self.msg_connect(self.gsm_receiver, "measurements", self.gsm_clock, "measurements")
        self.msg_connect(self.demapper, "bursts", self.decoder, "bursts")
        self.msg_connect(self.decoder, "msgs", self.dbg, "store")
        if count_bursts:
            # the receiver only emits bursts once it has locked to a cell's FCCH+SCH,
            # so any burst at all = a real GSM carrier is there (even if too weak to decode)
            self.bursts = blocks.message_debug()
            self.msg_connect(self.gsm_receiver, "C0", self.bursts, "store")

def run_one(n_use, arfcn, freq, stop_early, probe=None):
    """Run Owen's decoder over the first n_use samples.
    probe = N samples: give up (no sync) if the receiver has not locked after N samples.
    Returns (si3 lines, synced or None, status)."""
    global HAVE_PROGRESS
    tb = one(n_use, freq - CENTER, freq, count_bursts=probe is not None)
    done = threading.Event()

    def waiter():
        tb.wait()
        done.set()

    tb.start()
    th = threading.Thread(target=waiter, daemon=True)
    th.start()
    lines, shown, seen = [], set(), 0
    synced = None if probe is None else False
    status = None

    def drain():
        nonlocal seen
        n = tb.dbg.num_messages()
        while seen < n:
            text = si3_text(blob_of(tb.dbg.get_message(seen)), arfcn, freq)
            seen += 1
            if text and text not in shown:
                shown.add(text)
                lines.append(text)

    while True:
        finished = done.is_set()
        drain()
        if probe is not None and not synced:
            # the receiver only emits C0 bursts after locking to FCCH+SCH = real GSM carrier
            if tb.bursts.num_messages() > 0 or lines:
                synced = True
            elif not finished:
                try:
                    done_samples = tb.head.nitems_written(0) / tb.items_per_sample
                except Exception:
                    HAVE_PROGRESS = False
                    done_samples = 0
                if done_samples >= probe:
                    status = "no sync"
                    tb.stop()
                    th.join(10)
                    break
        if stop_early and lines and not finished:
            tb.stop()
            th.join(10)
            break
        if finished:
            time.sleep(0.1)          # Owen's settle time for the last messages
            drain()
            if probe is not None and not synced and tb.bursts.num_messages() > 0:
                synced = True
            break
        done.wait(0.02)
    del tb
    gc.collect()
    if status is None:
        if lines:
            status = "decoded"
        elif probe is not None:
            status = "synced, no cell ID" if synced else "no sync"
        else:
            status = "nothing"
    return lines, synced, status

def decode(n_in, arfcn, freq):
    """Returns (lines, status). status: decoded / decoded on retry / no sync /
    synced, no cell ID / nothing."""
    if SKIP_EMPTY:
        probe = min(n_in, int(RATE * PROBE_SECS))
        if HAVE_PROGRESS:
            # single pass: abort only if no sync within the first PROBE_SECS of capture
            lines, synced, status = run_one(n_in, arfcn, freq, STOP_EARLY, probe=probe)
        else:
            # fallback (no sample counter): short probe pass, then a full pass if it synced
            lines, synced, status = run_one(probe, arfcn, freq, STOP_EARLY, probe=probe + 1)
            if synced and not (lines and STOP_EARLY) and probe < n_in:
                lines, _, status = run_one(n_in, arfcn, freq, STOP_EARLY)
                status = "decoded" if lines else "synced, no cell ID"
    else:
        lines, synced, status = run_one(n_in, arfcn, freq, STOP_EARLY)
    # gr-gsm's clock correction runs asynchronously, so the same file can decode differently
    # from run to run - retry carriers that synced (or anything, without --skip-empty)
    tries = 0
    while not lines and tries < RETRIES and status in ("synced, no cell ID", "nothing"):
        tries += 1
        lines, _, _ = run_one(n_in, arfcn, freq, STOP_EARLY)
        if lines:
            status = "decoded on retry"
    return lines, status

def main():
    global PATH, RADIO, FMT, RATE, CENTER, SECS, GAIN
    import argparse
    t_start = time.time()
    cap_secs = None
    ap = argparse.ArgumentParser(description="Wideband GSM: one capture, decode every channel in it.")
    ap.add_argument("file", nargs="?", help="decode this existing capture instead of recording")
    ap.add_argument("--radio", choices=["lime", "bladerf"], default="lime",
                    help="radio to record with (default lime, Owen's rig)")
    ap.add_argument("--rate", type=float, default=None, help="sample rate (default 40e6)")
    ap.add_argument("--center", type=float, default=None, help="centre MHz (default 942.5)")
    ap.add_argument("--secs", type=float, default=None, help="capture seconds (default 2)")
    ap.add_argument("--gain", type=float, default=None, help="default 50 lime / 40 bladerf")
    ap.add_argument("--format", choices=["auto", "cfile", "sc16"], default="auto",
                    help="file format; auto = .sc16 extension -> sc16, else cfile")
    ap.add_argument("--skip-empty", action="store_true",
                    help="probe each channel first; skip it if the receiver never syncs")
    ap.add_argument("--stop-early", action="store_true",
                    help="stop each channel at its first decoded cell ID")
    ap.add_argument("--retries", type=int, default=None,
                    help="re-decode channels that synced but gave no cell ID (default 0, --fast 1)")
    ap.add_argument("--fast", action="store_true", help="= --skip-empty --stop-early --retries 1")
    ap.add_argument("--probe-secs", type=float, default=2.0,
                    help="--skip-empty gives up if no sync within this much capture (default 2.0)")
    args = ap.parse_args()
    global SKIP_EMPTY, STOP_EARLY, PROBE_SECS, RETRIES
    SKIP_EMPTY = args.skip_empty or args.fast
    STOP_EARLY = args.stop_early or args.fast
    PROBE_SECS = args.probe_secs
    RETRIES = args.retries if args.retries is not None else (1 if args.fast else 0)
    RADIO = args.radio
    if args.rate is not None:
        RATE = args.rate
    if args.center is not None:
        CENTER = args.center * 1e6
    if args.secs is not None:
        SECS = args.secs
    GAIN = args.gain if args.gain is not None else (40.0 if RADIO == "bladerf" else 50.0)
    if args.file:
        PATH = os.path.abspath(os.path.expanduser(args.file))
    elif RADIO == "bladerf":
        PATH = "/dev/shm/gsm_wide.sc16"
    FMT = args.format if args.format != "auto" else ("sc16" if PATH.endswith(".sc16") else "cfile")
    if args.file:
        if not os.path.isfile(PATH):
            print("no file", PATH)
            return 1
        n_in = os.path.getsize(PATH) // (4 if FMT == "sc16" else 8)
        print("using %s (%s) samples %d = %.2fs at %.1f MS/s, centre %.1f MHz"
              % (PATH, FMT, n_in, n_in / RATE, RATE / 1e6, CENTER / 1e6), flush=True)
    else:
        t_cap = time.time()
        n_in = capture()
        cap_secs = time.time() - t_cap
        if not n_in:
            return 1
    chans = channels()
    mode = []
    if SKIP_EMPTY:
        mode.append("skip-empty (probe %.1fs)" % PROBE_SECS)
    if STOP_EARLY:
        mode.append("stop-early")
    if RETRIES:
        mode.append("retries %d" % RETRIES)
    print("decode mode: %s" % (" + ".join(mode) if mode else "full (Owen's original)"), flush=True)
    print("%d channels in the %.0f MS/s capture, one after another on this Pi" % (len(chans), RATE / 1e6), flush=True)
    cells = 0
    unique = {}
    status_count, status_time, synced_only = {}, {}, []
    t_dec = time.time()
    for index, (arfcn, freq) in enumerate(chans, 1):
        print("%d/%d  ARFCN %d  %.1f MHz" % (index, len(chans), arfcn, freq / 1e6), flush=True)
        t0 = time.time()
        lines, status = decode(n_in, arfcn, freq)
        took = time.time() - t0
        status_count[status] = status_count.get(status, 0) + 1
        status_time[status] = status_time.get(status, 0.0) + took
        if not lines:
            if status == "synced, no cell ID":
                synced_only.append((arfcn, freq))
                print("  SYNCED but no cell ID decoded  (%.1fs)" % took, flush=True)
            elif status == "no sync":
                print("  nothing on this channel  (no sync, %.1fs)" % took, flush=True)
            else:
                print("  nothing on this channel  (%.1fs)" % took, flush=True)
            continue
        for text in lines:
            cells += 1
            print(text, flush=True)
            plain = re.sub("\x1b\\[[0-9;]*m", "", text)
            cid = plain.split("CID:")[1].split(",")[0].strip() if "CID:" in plain else "?"
            unique.setdefault((arfcn, cid), plain.split(", Pwr")[0])
        print("  (%.1fs)" % took, flush=True)
    dec_secs = time.time() - t_dec
    print("cells", cells, flush=True)
    print("")
    print("===== %d unique cell(s) =====" % len(unique), flush=True)
    for key in sorted(unique, key=lambda k: (k[0] < 900, k[0])):
        print("  " + unique[key], flush=True)
    if synced_only:
        print("")
        print("GSM carrier present but no cell ID decoded (try longer --secs / other gain):", flush=True)
        print("  " + "  ".join("ARFCN %d (%.1f)" % (a, f / 1e6) for a, f in synced_only), flush=True)
    print("")
    print("----- timing -----", flush=True)
    if cap_secs is not None:
        print("  capture  %7.1fs" % cap_secs, flush=True)
    print("  decode   %7.1fs  (%d channels, avg %.2fs each)"
          % (dec_secs, len(chans), dec_secs / max(1, len(chans))), flush=True)
    for st in sorted(status_count):
        print("    %-22s %3d ch  %7.1fs" % (st, status_count[st], status_time[st]), flush=True)
    print("  total    %7.1fs" % (time.time() - t_start), flush=True)
    if SKIP_EMPTY and not HAVE_PROGRESS:
        print("  note: no sample counter in this GNU Radio - used two-pass probe instead", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
