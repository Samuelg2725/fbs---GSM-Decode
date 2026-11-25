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
PATH2 = "/dev/shm/gsm_wide_rx2.sc16"
RADIO = "lime"
FMT = "cfile"            # cfile = complex float32 (Lime/osmosdr), sc16 = int16 I/Q (bladeRF-cli)
BLADE_LIBS = "/usr/local/lib"
SKIP_EMPTY = False       # --skip-empty: short sync probe first, skip channels with no GSM sync
STOP_EARLY = False       # --stop-early: stop a channel as soon as its cell ID (SI3) is decoded
PROBE_SECS = 2.0         # --probe-secs: give up on a channel if no sync within this much capture
RETRIES = 0              # --retries: re-decode channels that synced but gave no cell ID
HAVE_PROGRESS = True
GAIN2 = None             # --gain2: bladerf RX2 at this gain, captured together with RX1 (two antennas)     # set False if this GNU Radio can't report samples processed

def channels():
    """Every UK GSM downlink channel inside the capture: E-GSM900 (Owen's list) + DCS1800."""
    lo = CENTER - RATE / 2 + 0.3e6
    hi = CENTER + RATE / 2 - 0.3e6
    found = []
    for n in range(975, 1024):
        freq = 935e6 + 0.2e6 * (n - 1024)
        if lo <= freq <= hi:
            found.append((n, freq))
    for n in range(1, 125):
        freq = 935e6 + 0.2e6 * n
        if lo <= freq <= hi:
            found.append((n, freq))
    for n in range(512, 886):                     # DCS1800: 1805.2-1879.8 MHz
        freq = 1805.2e6 + 0.2e6 * (n - 512)
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
    if GAIN2 is not None:
        return capture_bladerf_dual(n)
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

def capture_bladerf_dual(n):
    """RX1 + RX2 together (shared tuner, separate gains). bladeRF-cli writes them interleaved;
    split into PATH (RX1) and PATH2 (RX2). Gains are set after 'rx config' so both channels get them."""
    raw = PATH + ".dual"
    # 'set agc off' (both channels) - 'set agc rx off' only covers RX1 and RX2 then refuses a manual gain
    cmd = ["bladeRF-cli", "-e", "set agc off",
           "-e", "set frequency rx %d" % int(CENTER),
           "-e", "set samplerate rx %d" % int(RATE),
           "-e", "set bandwidth rx %d" % int(min(RATE * 0.9, 56e6)),
           "-e", "rx config file=%s format=bin n=%d channel=1,2 samples=16384 buffers=64 xfers=16 timeout=20000"
           % (raw, n),
           "-e", "set gain rx1 %d" % int(GAIN), "-e", "set gain rx2 %d" % int(GAIN2),
           "-e", "rx start", "-e", "rx wait %d" % int(SECS * 3000 + 15000)]
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = BLADE_LIBS + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    print("capturing %.1fs at %.1f MS/s, centre %.1f MHz, RX1 gain %g + RX2 gain %g (bladeRF-cli, two antennas)"
          % (SECS, RATE / 1e6, CENTER / 1e6, GAIN, GAIN2), flush=True)
    started = time.time()
    p = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    size = os.path.getsize(raw) if os.path.exists(raw) else 0
    print("wrote %d bytes in %.1fs" % (size, time.time() - started), flush=True)
    if size < n * 8 * 0.9:
        print("capture is short - bladeRF-cli said:", flush=True)
        print(p.stdout.decode(errors="replace")[-800:], flush=True)
        return 0
    import numpy as np
    d = np.memmap(raw, dtype=np.int16, mode="r")
    d = d[:(len(d) // 4) * 4].reshape(-1, 4)
    step = 1 << 20
    with open(PATH, "wb") as f1, open(PATH2, "wb") as f2:
        for i in range(0, len(d), step):
            f1.write(np.ascontiguousarray(d[i:i + step, 0:2]).tobytes())
            f2.write(np.ascontiguousarray(d[i:i + step, 2:4]).tobytes())
    del d
    os.remove(raw)
    for name, path in (("RX1", PATH), ("RX2", PATH2)):
        x = np.memmap(path, dtype=np.int16, mode="r")[:2000000].astype(np.float32)
        pk = np.abs(x).max() / 2048.0
        print("  %s level %.1f dBFS, peak %.2f%s" % (name, 10 * np.log10(np.mean(x * x) / 2048.0 ** 2 * 2 + 1e-20),
              pk, "  CLIPPING" if pk >= 0.99 else ""), flush=True)
    return size // 8


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

def decode_all(n_in, chans, tag=None):
    """Decode every channel in chans from PATH. Returns (unique, cells, status_count, status_time, synced_only)."""
    cells = 0
    unique = {}
    status_count, status_time, synced_only = {}, {}, []
    for index, (arfcn, freq) in enumerate(chans, 1):
        print("%s%d/%d  ARFCN %d  %.1f MHz" % ("[%s] " % tag if tag else "", index, len(chans), arfcn, freq / 1e6),
              flush=True)
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
    return unique, cells, status_count, status_time, synced_only


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
    ap.add_argument("--gain", type=float, default=None, help="default 50 lime / 20 bladerf (40 overloaded strong cells at home)")
    ap.add_argument("--gain2", type=float, default=None,
                    help="bladerf: ALSO record RX2 (second antenna) at this gain; RX2 is decoded only on channels "
                         "where RX1 found no cell (e.g. --gain 20 --gain2 45)")
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
    global SKIP_EMPTY, STOP_EARLY, PROBE_SECS, RETRIES, GAIN2
    GAIN2 = args.gain2 if (args.gain2 is not None and args.radio == "bladerf" and not args.file) else None
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
    GAIN = args.gain if args.gain is not None else (20.0 if RADIO == "bladerf" else 50.0)   # bladerf 40 overloads
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
    t_dec = time.time()
    unique, cells, status_count, status_time, synced_only = decode_all(n_in, chans, "RX1" if GAIN2 is not None else None)
    rx2_added = []
    if GAIN2 is not None:
        done = {a for (a, _cid) in unique}
        todo = [(a, f) for (a, f) in chans if a not in done]
        print("", flush=True)
        print("===== RX2 (gain %g): %d channel(s) RX1 found nothing on =====" % (GAIN2, len(todo)), flush=True)
        keep_path = PATH
        PATH = PATH2
        u2, c2, sc2, st2, so2 = decode_all(n_in, todo, "RX2")
        PATH = keep_path
        cells += c2
        for k, v in u2.items():
            if k not in unique:
                unique[k] = v
                rx2_added.append(k)
        for k in sc2:
            status_count["RX2 " + k] = sc2[k]
            status_time["RX2 " + k] = st2[k]
        synced_only = [x for x in synced_only if x[0] not in {a for a, _ in u2}] + \
            [x for x in so2 if x[0] not in {a for a, _ in synced_only}]
    dec_secs = time.time() - t_dec
    print("cells", cells, flush=True)
    print("")
    print("===== %d unique cell(s) =====" % len(unique), flush=True)
    for key in sorted(unique, key=lambda k: (k[0] < 900, k[0])):
        print("  " + unique[key] + ("   <- RX2 only" if key in rx2_added else ""), flush=True)
    if GAIN2 is not None:
        print("RX1 (gain %g) found %d, RX2 (gain %g) added %d" % (GAIN, len(unique) - len(rx2_added), GAIN2,
              len(rx2_added)), flush=True)
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
