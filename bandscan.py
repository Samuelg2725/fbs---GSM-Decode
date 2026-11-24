#!/usr/bin/env python3
"""bandscan.py - GSM cell scanner for the BladeRF 2.0: any band, wide OR narrow mode, with timings.

MODES (pick one per run)
  --mode wide     one wide capture per chunk (default 15.36 Msps), FFT finds the live
                  channels, each one is filtered out digitally and decoded (in parallel).
  --mode narrow   quick FFT survey, then ONE live radio stream (osmosdr, same as
                  gsm_arfcn_fast.py) tunes to each live channel, prints its cell as soon as
                  the SI3 arrives (~0.5-1 s per channel), then moves on.
  --compare       run wide then narrow on the same band and print a timing table.

BANDS (--band, comma list or repeat the flag)
  gsm450 gsm480 gsm850 pgsm900 egsm900 (= gsm900) rgsm900 dcs1800 pcs1900
  groups:  uk / eu = egsm900+dcs1800   us = gsm850+pcs1900   all = every band above
  custom:  --range 1805-1830 --range 900-960   (MHz, any span the BladeRF can tune, 47-6000)
           Points on a known GSM raster get their real ARFCN; anything else is scanned
           on a 200 kHz raster and labelled "custom".

EXAMPLES
  python3 bandscan.py                                   # E-GSM900, wide mode
  python3 bandscan.py --band all                        # every GSM band
  python3 bandscan.py --band uk --mode narrow
  python3 bandscan.py --range 1805-1830 --gain-sweep
  python3 bandscan.py --band gsm900 --compare --timing-csv ~/timings.csv
  python3 bandscan.py --band all --dry-run              # print the plan only, no radio

Every stage is timed (capture / survey / detect / decode / wait) and printed at the end of
each run. --timing-csv appends the numbers to a CSV so runs can be compared over time.
"""
import argparse
import csv
import math
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import threading
import time
from collections import namedtuple

import numpy as np

SHM = "/dev/shm"
CHANNEL_SPACING = 200e3

# band -> list of (first ARFCN, last ARFCN, downlink MHz of first ARFCN), 200 kHz steps
BANDS = {
    "gsm450": [(259, 293, 460.6)],
    "gsm480": [(306, 340, 489.0)],
    "gsm850": [(128, 251, 869.2)],
    "pgsm900": [(1, 124, 935.2)],
    "egsm900": [(975, 1023, 925.2), (0, 124, 935.0)],
    "rgsm900": [(955, 1023, 921.2), (0, 124, 935.0)],
    "dcs1800": [(512, 885, 1805.2)],
    "pcs1900": [(512, 810, 1930.2)],
}
# when two bands share a frequency, the first in this list names it
LABEL_ORDER = ["pgsm900", "egsm900", "rgsm900", "gsm850", "dcs1800", "pcs1900", "gsm450", "gsm480"]
GROUPS = {
    "gsm900": ["egsm900"], "900": ["egsm900"], "1800": ["dcs1800"], "dcs": ["dcs1800"],
    "850": ["gsm850"], "1900": ["pcs1900"], "pcs": ["pcs1900"],
    "uk": ["egsm900", "dcs1800"], "eu": ["egsm900", "dcs1800"],
    "us": ["gsm850", "pcs1900"],
    "all": ["gsm450", "gsm480", "gsm850", "rgsm900", "dcs1800", "pcs1900"],
}
OPERATORS = {(234, 2): "O2", (234, 10): "O2", (234, 11): "O2", (234, 15): "Vodafone",
             (234, 20): "Three", (234, 30): "EE", (234, 33): "EE", (234, 34): "EE",
             (235, 0): "Mundio", (235, 2): "EE"}

Chan = namedtuple("Chan", "band arfcn freq")   # freq in Hz (int)


# ----------------------------------------------------------------------------- band planning

def band_channels(name):
    out = []
    for first, last, f0 in BANDS[name]:
        for n in range(first, last + 1):
            out.append(Chan(name, n, int(round(f0 * 1e3 + 200 * (n - first))) * 1000))
    return out


def known_by_freq():
    table = {}
    for name in reversed(LABEL_ORDER):          # earlier names overwrite later ones
        for c in band_channels(name):
            table[c.freq] = c
    return table


def expand_bands(spec_list):
    names = []
    for spec in spec_list:
        for tok in spec.replace(" ", "").lower().split(","):
            if not tok:
                continue
            if tok in GROUPS:
                names.extend(GROUPS[tok])
            elif tok in BANDS:
                names.append(tok)
            else:
                raise SystemExit("unknown band '%s'. Known: %s, groups: %s"
                                 % (tok, " ".join(BANDS), " ".join(GROUPS)))
    seen, out = set(), []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def parse_range(text):
    for sep in ("-", ":", ","):
        if sep in text:
            a, b = text.split(sep, 1)
            lo, hi = sorted((float(a), float(b)))
            return lo * 1e6, hi * 1e6
    raise SystemExit("bad --range '%s' (use e.g. 1805-1830)" % text)


def build_channels(band_names, ranges):
    chans = {}
    for name in band_names:
        for c in band_channels(name):
            chans.setdefault(c.freq, c)
    if ranges:
        known = known_by_freq()
        for lo, hi in ranges:
            k0 = int(math.ceil(lo / CHANNEL_SPACING - 1e-9))
            k1 = int(math.floor(hi / CHANNEL_SPACING + 1e-9))
            for k in range(k0, k1 + 1):
                f = k * int(CHANNEL_SPACING)
                chans.setdefault(f, known.get(f, Chan("custom", None, f)))
    out = sorted(chans.values(), key=lambda c: c.freq)
    bad = [c for c in out if c.freq < 47e6 or c.freq > 6000e6]
    if bad:
        print("warning: %d channel(s) outside BladeRF 2.0 range (47-6000 MHz) dropped" % len(bad))
        out = [c for c in out if 47e6 <= c.freq <= 6000e6]
    return out


def usable_half(rate, edge_frac):
    """Half-width (Hz) of the part of a capture we trust (analog filter rolls off at the edges)."""
    return rate / 2 - max(0.3e6, edge_frac * rate)


def plan_chunks(chans, rate, edge_frac, manual_centres=None):
    """Group channels into wide-capture chunks. Returns list of (centre_hz, [chans])."""
    half = usable_half(rate, edge_frac)
    if manual_centres:
        plan = [(c, []) for c in manual_centres]
        for ch in chans:
            for centre, members in plan:
                if abs(ch.freq - centre) <= half:
                    members.append(ch)
                    break
        return [(c, m) for c, m in plan if m]
    # split into contiguous segments (gaps wider than one capture), then cut each segment
    # into the fewest equal-width chunks that fit - avoids a tiny leftover chunk at the end
    segments, cur = [], [chans[0]]
    for c in chans[1:]:
        if c.freq - cur[-1].freq > 2 * half:
            segments.append(cur)
            cur = [c]
        else:
            cur.append(c)
    segments.append(cur)
    plan = []
    for seg in segments:
        lo, hi = seg[0].freq, seg[-1].freq
        n = max(1, int(math.ceil((hi - lo) / (2 * half) - 1e-9)))
        width = (hi - lo) / n
        bins = [[] for _ in range(n)]
        for c in seg:
            bins[min(n - 1, int((c.freq - lo) / width)) if width > 0 else 0].append(c)
        for members in bins:
            if members:
                plan.append(((members[0].freq + members[-1].freq) / 2.0, members))
    return plan


def chan_name(c):
    return ("%s ARFCN %d" % (c.band, c.arfcn)) if c.arfcn is not None else "custom"


# ----------------------------------------------------------------------------- SI3 parsing

def py_to_bytes(py):
    if isinstance(py, (bytes, bytearray)):
        return bytes(py)
    if isinstance(py, tuple):
        for part in py:
            g = py_to_bytes(part)
            if g and len(g) >= 12:
                return g
    if isinstance(py, list) and py and isinstance(py[0], int):
        try:
            return bytes(bytearray(py))
        except Exception:
            return None
    return None


def blob_of(msg):
    import pmt
    for getter in (lambda: pmt.to_python(msg),
                   lambda: pmt.to_python(pmt.cdr(msg)),
                   lambda: pmt.blob_data(pmt.cdr(msg))):
        try:
            obj = getter()
        except Exception:
            continue
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


def si3(raw):
    if not raw:
        return None
    starts = []
    if len(raw) > 1 and 1 <= raw[1] <= 12:
        starts.append(raw[1] * 4)
    starts.extend((16, 0))
    for off in starts:
        if off + 12 <= len(raw) and raw[off + 2] == 0x1b:
            b = raw[off:]
            cid = (b[3] << 8) + b[4]
            if cid == 0:
                return None
            mcc = ((b[5] & 0xF) * 100) + (((b[5] & 0xF0) >> 4) * 10) + (b[6] & 0xF)
            mnc = ((b[7] & 0xF) * 10) + (b[7] >> 4)
            if ((b[6] & 0xF0) >> 4) < 10:
                mnc = mnc * 10 + ((b[6] & 0xF0) >> 4)
            lac = (b[8] << 8) + b[9]
            return (cid, lac, mcc, mnc)
    return None


# ----------------------------------------------------------------------------- decode worker

def _worker_init():
    """Pool workers: leave Ctrl-C to the main process and run at low priority so the
    Pi (ssh, USB) stays responsive while decoding."""
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        os.nice(10)
    except OSError:
        pass


def _warm(_):
    """Import GNU Radio in a pool worker up front so its start-up cost is timed separately."""
    from gnuradio import blocks, gr, gsm  # noqa: F401
    return os.getpid()


def decode_job(job):
    """Decode one channel straight from a raw sc16 capture. Runs in a pool worker.
    Stops as soon as an SI3 is decoded, or after job['decode_secs'] of samples."""
    t0 = time.time()
    info, nmsg, err = None, 0, None
    try:
        from gnuradio import blocks, gr, gsm
        from gnuradio import filter as grfilter
        from gnuradio.filter import firdes
        rate, f, centre, decim = job["rate"], job["freq"], job["centre"], job["decim"]
        chanrate = rate / decim
        tb = gr.top_block()
        chain = [blocks.file_source(gr.sizeof_short, job["path"], False)]
        if job["decode_secs"] > 0:
            chain.append(blocks.head(gr.sizeof_short, 2 * int(rate * job["decode_secs"])))
        chain.append(blocks.interleaved_short_to_complex(False, False))
        chain.append(blocks.multiply_const_cc(complex(1.0 / 2048.0, 0.0)))
        if decim > 1 or f != centre:
            taps = firdes.low_pass(1.0, rate, job["cutoff"], job["tw"])
            chain.append(grfilter.freq_xlating_fir_filter_ccf(decim, taps, f - centre, rate))
        gi = gsm.gsm_input(ppm=0, osr=4, fc=f, samp_rate_in=chanrate)
        rcv = gsm.receiver(4, ([0]), ([]))
        clk = gsm.clock_offset_control(f, chanrate, osr=4)
        dem = gsm.gsm_bcch_ccch_demapper(0)
        dcd = gsm.control_channels_decoder()
        dbg = blocks.message_debug()
        tb.connect(*(chain + [gi, rcv]))
        tb.msg_connect(clk, "ctrl", gi, "ctrl_in")
        tb.msg_connect(rcv, "C0", dem, "bursts")
        tb.msg_connect(rcv, "measurements", clk, "measurements")
        tb.msg_connect(dem, "bursts", dcd, "bursts")
        tb.msg_connect(dcd, "msgs", dbg, "store")
        tb.start()
        waiter = threading.Thread(target=tb.wait, daemon=True)
        waiter.start()
        seen = 0
        while True:
            finished = not waiter.is_alive()
            n = dbg.num_messages()
            while seen < n and info is None:
                info = si3(blob_of(dbg.get_message(seen)))
                seen += 1
            nmsg = n
            if info is not None or finished:
                break
            time.sleep(0.02)
        tb.stop()
        waiter.join(timeout=10)
    except Exception as exc:
        err = "%s: %s" % (type(exc).__name__, exc)
    out = dict(job)
    out.update(info=info, nmsg=nmsg, err=err, secs=time.time() - t0)
    return out


class _Done(object):
    """Looks like an AsyncResult for jobs run inline (--jobs 1)."""
    def __init__(self, value):
        self.value = value

    def ready(self):
        return True

    def get(self):
        return self.value


# ----------------------------------------------------------------------------- timing

class Timer(object):
    def __init__(self):
        self.stages = {}
        self.order = []

    def add(self, stage, secs, count=1):
        if stage not in self.stages:
            self.stages[stage] = [0, 0.0]
            self.order.append(stage)
        self.stages[stage][0] += count
        self.stages[stage][1] += secs

    def get(self, stage):
        return self.stages.get(stage, [0, 0.0])


STAGE_HELP = {
    "startup": "start decoder worker processes",
    "survey": "short wide captures to find live channels (narrow mode)",
    "capture": "wide radio captures (bladeRF-cli, ~2 s of it is device start-up)",
    "open": "narrow: open the live osmosdr stream",
    "live": "narrow: tune + wait for SI3, per channel (live osmosdr stream)",
    "detect": "FFT channel detection",
    "decode": "per-channel decode work (sum over all jobs, overlaps in parallel)",
    "wait": "time spent blocked waiting for decodes to finish",
}


# ----------------------------------------------------------------------------- the scanner

class Scanner(object):
    def __init__(self, args, chans, say, pool):
        self.a = args
        self.chans = chans
        self.say = say
        self.pool = pool
        self.env = dict(os.environ)
        if args.lib:
            self.env["LD_LIBRARY_PATH"] = args.lib + (":" + self.env["LD_LIBRARY_PATH"]
                                                      if self.env.get("LD_LIBRARY_PATH") else "")

    # -- radio --------------------------------------------------------------------------------
    def _radio_setup(self, rate, gain):
        bw = min(rate * 0.9, 56e6)
        return ["bladeRF-cli", "-e", "set agc rx off",
                "-e", "set samplerate rx %d" % int(rate),
                "-e", "set bandwidth rx %d" % int(bw),
                "-e", "set gain rx %d" % int(gain)]

    def _rx(self, freq, path, nsamp, secs):
        return ["-e", "set frequency rx %d" % int(freq),
                "-e", "rx config file=%s format=bin n=%d samples=16384 buffers=64 xfers=16 timeout=20000"
                % (path, nsamp),
                "-e", "rx start", "-e", "rx wait %d" % int(secs * 3000 + 15000)]

    def _run_cli(self, cmd, timeout):
        try:
            p = subprocess.run(cmd, env=self.env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                               timeout=timeout)
            return p.returncode, p.stderr.decode(errors="replace")
        except Exception as exc:
            return -1, str(exc)

    @staticmethod
    def _good(path, nsamp):
        return os.path.exists(path) and os.path.getsize(path) >= 0.9 * nsamp * 4

    def capture(self, centre, gain, rate, secs, path):
        """One wide capture. Returns (ok, seconds)."""
        nsamp = int(rate * secs)
        if os.path.exists(path):
            os.remove(path)
        t0 = time.time()
        rc, err = self._run_cli(self._radio_setup(rate, gain) + self._rx(centre, path, nsamp, secs),
                                timeout=secs * 5 + 60)
        dt = time.time() - t0
        ok = self._good(path, nsamp)
        if not ok:
            tail = (err.strip().splitlines() or ["no output"])[-1]
            self.say("  capture failed (rc %s): %s" % (rc, tail))
        return ok, dt

    # -- detection ----------------------------------------------------------------------------
    def detect(self, path, rate, centre, members):
        """FFT the start of a raw sc16 capture and return [(chan, snr_dB)] above --snr."""
        nfft = 4096
        raw = np.memmap(path, dtype=np.int16, mode="r")
        nsamp = raw.size // 2
        skip = min(int(0.25 * rate), max(0, nsamp - nfft))   # skip tuning/settling transient
        nblk = min((nsamp - skip) // nfft, max(1, int(self.a.detect_secs * rate) // nfft))
        if nblk < 1:
            return []
        win = np.hanning(nfft).astype(np.float32)
        acc = np.zeros(nfft)
        for b0 in range(0, nblk, 256):
            b1 = min(nblk, b0 + 256)
            s = np.asarray(raw[2 * (skip + b0 * nfft):2 * (skip + b1 * nfft)], dtype=np.float32)
            x = (s[0::2] + 1j * s[1::2]).reshape(b1 - b0, nfft)
            acc += (np.abs(np.fft.fft(x * win, axis=1)) ** 2).sum(axis=0)
        del raw
        psd = np.fft.fftshift(acc / nblk)
        binhz = rate / nfft
        axis = centre + (np.arange(nfft) - nfft / 2) * binhz
        half = usable_half(rate, self.a.edge)
        inband = np.abs(axis - centre) <= half
        # lower quartile, not median: 3G/4G carriers in the band pushed the median up and
        # real GSM cells (ARFCN 13/18/22) fell under the threshold at gain 40
        noise = np.percentile(psd[inband] if inband.any() else psd, self.a.noise_pct)
        w = max(1, int(100e3 / binhz))
        out = []
        for c in members:
            if abs(c.freq - centre) > half:
                continue
            k = int(np.argmin(np.abs(axis - c.freq)))
            p = psd[max(0, k - w):k + w].mean()
            snr = 10 * np.log10(p / noise + 1e-12)
            if snr >= self.a.snr:
                out.append((c, snr))
        out.sort(key=lambda x: -x[1])
        return out

    # -- decode bookkeeping -------------------------------------------------------------------
    def _submit(self, jobs, path):
        self.refs[path] = self.refs.get(path, 0) + len(jobs)
        if not jobs:
            self._release(path, 0)
        for job in jobs:
            if self.pool is not None:
                self.pending.append(self.pool.apply_async(decode_job, (job,)))
            else:
                t0 = time.time()
                self.pending.append(_Done(decode_job(job)))
                self.timer.add("wait", time.time() - t0, 0)

    def _release(self, path, n=1):
        self.refs[path] = self.refs.get(path, 0) - n
        if self.refs[path] <= 0 and not self.a.keep:
            try:
                os.remove(path)
            except OSError:
                pass

    def _drain(self, block):
        """Print results in the order they FINISH (not the order submitted)."""
        t0 = time.time()
        while True:
            still = []
            for ar in self.pending:
                if ar.ready():
                    self._handle(ar.get())
                else:
                    still.append(ar)
            self.pending = still
            if not block or not self.pending:
                break
            time.sleep(0.05)
        if block and self.pool is not None:
            self.timer.add("wait", time.time() - t0, 0)

    def _handle(self, r):
        c = r["chan"]
        if r.get("path"):
            self.timer.add("decode", r["secs"])
            self._release(r["path"])
        snr = r.get("snr")
        snr_txt = ("%+.0f dB" % snr) if snr is not None else "  -   "
        if r["err"]:
            self.say("    %-18s %9.1f MHz  ERROR %s" % (chan_name(c), c.freq / 1e6, r["err"]))
            return
        if r["info"] is None:
            if self.a.verbose:
                self.say("    %-18s %9.1f MHz  %s  no SI3 (%d msgs)  %.1fs"
                         % (chan_name(c), c.freq / 1e6, snr_txt, r["nmsg"], r["secs"]))
            return
        cid, lac, mcc, mnc = r["info"]
        op = OPERATORS.get((mcc, mnc), "?")
        new = c.freq not in self.found
        if new:
            self.found[c.freq] = dict(chan=c, cid=cid, lac=lac, mcc=mcc, mnc=mnc, op=op, snr=snr,
                                      gain=r["gain"], secs=r["secs"], t=time.time() - self.t_start)
        self.say("  %s %-18s %9.1f MHz  CID %-6d LAC %-5d PLMN %d-%02d [%s]  %s  took %.1fs  (%.0fs into run)"
                 % ("DECODED" if new else "again  ", chan_name(c), c.freq / 1e6, cid, lac, mcc, mnc, op,
                    snr_txt, r["secs"], time.time() - self.t_start))

    def _job(self, c, path, rate, centre, decim, gain, snr):
        a = self.a
        chanrate = rate / decim
        cutoff = a.chan_cutoff * 1e3
        tw = min(a.chan_tw * 1e3, max(50e3, chanrate - 300e3 - cutoff))
        return dict(chan=c, freq=float(c.freq), path=path, rate=float(rate), centre=float(centre),
                    decim=decim, cutoff=cutoff, tw=tw, decode_secs=a.decode_secs, gain=gain, snr=snr)

    # -- modes --------------------------------------------------------------------------------
    def run(self, mode):
        self.timer = Timer()
        self.found = {}
        self.pending = []
        self.refs = {}
        self.t_start = time.time()
        a = self.a
        plan = plan_chunks(self.chans, a.rate, a.edge, a.chunk_centres)
        self.say("")
        self.say("######## %s mode: %d channels, %d chunk(s), gains %s ########"
                 % (mode.upper(), len(self.chans), len(plan), ",".join("%g" % g for g in a.gain_list)))
        try:
            if mode == "wide":
                self._run_wide(plan)
            else:
                self._run_narrow(plan)
            self._drain(True)
        except KeyboardInterrupt:
            self.say("\ninterrupted - reporting what was found so far")
        wall = time.time() - self.t_start
        return self.found, self.timer, wall

    def _run_wide(self, plan):
        a = self.a
        decim = max(1, int(a.rate // 1.92e6))
        for gain in a.gain_list:
            for ci, (centre, members) in enumerate(plan):
                todo_members = [c for c in members if c.freq not in self.found]
                if not todo_members:
                    continue
                path = "%s/_bs_wide_%d_%g.sc16" % (SHM, ci, gain)
                self.say("=== chunk %d/%d  %.1f MHz (%.1f-%.1f)  gain %g ==="
                         % (ci + 1, len(plan), centre / 1e6, members[0].freq / 1e6,
                            members[-1].freq / 1e6, gain))
                ok, dt = self.capture(centre, gain, a.rate, a.secs, path)
                self.timer.add("capture", dt)
                self.say("  capture %.1fs" % dt)
                if not ok:
                    continue
                t0 = time.time()
                if a.detect == "fft":
                    live = self.detect(path, a.rate, centre, todo_members)
                else:
                    live = [(c, None) for c in todo_members]
                self.timer.add("detect", time.time() - t0)
                self.say("  %d live channel(s)  detect %.2fs" % (len(live), time.time() - t0))
                jobs = [self._job(c, path, a.rate, centre, decim, gain, s) for c, s in live]
                self._submit(jobs, path)
                self._drain(not a.pipeline)

    def _run_narrow(self, plan):
        a = self.a
        for gain in a.gain_list:
            live = []
            if a.detect == "fft":
                for ci, (centre, members) in enumerate(plan):
                    todo_members = [c for c in members if c.freq not in self.found]
                    if not todo_members:
                        continue
                    path = "%s/_bs_survey_%d.sc16" % (SHM, ci)
                    ok, dt = self.capture(centre, gain, a.rate, a.survey_secs, path)
                    self.timer.add("survey", dt)
                    if ok:
                        t0 = time.time()
                        got = self.detect(path, a.rate, centre, todo_members)
                        self.timer.add("detect", time.time() - t0)
                        live.extend(got)
                        self.say("  survey %.1f MHz gain %g: %d live  (%.1fs)"
                                 % (centre / 1e6, gain, len(got), dt))
                    if not a.keep and os.path.exists(path):
                        os.remove(path)
            else:
                live = [(c, None) for c in self.chans if c.freq not in self.found]
            live = [(c, s) for c, s in live if c.freq not in self.found]
            if not live:
                continue
            self.say("=== narrow (live stream): %d channel(s) at gain %g ===" % (len(live), gain))
            self._live_decode(live, gain)

    def _live_decode(self, live, gain):
        """One osmosdr stream, retuned channel to channel - the gsm_arfcn_fast.py way."""
        a = self.a
        from gnuradio import blocks, gr, gsm
        import osmosdr
        rate = a.narrow_rate
        first = float(live[0][0].freq)

        class LiveRx(gr.top_block):
            def __init__(self):
                gr.top_block.__init__(self, "bandscan_live")
                self.src = osmosdr.source(args=a.osmo_args)
                self.src.set_sample_rate(rate)
                self.src.set_center_freq(first, 0)
                self.src.set_freq_corr(0, 0)
                self.src.set_gain_mode(False, 0)
                self.src.set_gain(gain, 0)
                self.src.set_bandwidth(rate, 0)
                self.gi = gsm.gsm_input(ppm=0, osr=4, fc=first, samp_rate_in=rate)
                self.rcv = gsm.receiver(4, ([0]), ([]))
                self.clk = gsm.clock_offset_control(first, rate, osr=4)
                self.dem = gsm.gsm_bcch_ccch_demapper(0)
                self.dcd = gsm.control_channels_decoder()
                self.dbg = blocks.message_debug()
                self.connect(self.src, self.gi, self.rcv)
                self.msg_connect(self.clk, "ctrl", self.gi, "ctrl_in")
                self.msg_connect(self.rcv, "C0", self.dem, "bursts")
                self.msg_connect(self.rcv, "measurements", self.clk, "measurements")
                self.msg_connect(self.dem, "bursts", self.dcd, "bursts")
                self.msg_connect(self.dcd, "msgs", self.dbg, "store")

            def tune(self, f):
                self.stop()
                self.wait()
                self.src.set_center_freq(f, 0)
                try:
                    self.gi.set_fc(f)
                except Exception:
                    pass
                try:
                    self.clk.set_fc(f)
                except Exception:
                    pass
                self.start()

        t_open = time.time()
        tb = LiveRx()
        tb.start()
        self.timer.add("open", time.time() - t_open)
        try:
            for i, (c, snr) in enumerate(live):
                t0 = time.time()
                if i:
                    tb.tune(float(c.freq))
                seen = tb.dbg.num_messages()
                info, nmsg = None, 0
                settle = t0 + 0.1          # ignore anything left over from the last channel
                while time.time() - t0 < a.narrow_secs:
                    n = tb.dbg.num_messages()
                    while seen < n and info is None:
                        msg = tb.dbg.get_message(seen)
                        seen += 1
                        nmsg += 1
                        if time.time() >= settle:
                            info = si3(blob_of(msg))
                    if info is not None:
                        break
                    time.sleep(0.02)
                dt = time.time() - t0
                self.timer.add("live", dt)
                self._handle(dict(chan=c, info=info, nmsg=nmsg, err=None, secs=dt, gain=gain,
                                  snr=snr, path=None))
        finally:
            tb.stop()
            tb.wait()
            del tb
            time.sleep(0.5)        # let osmosdr release the BladeRF before the next bladeRF-cli


# ----------------------------------------------------------------------------- reporting

def print_cells(say, title, found):
    say("")
    say("===== %s: %d unique cell(s) =====" % (title, len(found)))
    for f in sorted(found):
        d = found[f]
        c = d["chan"]
        say("  %-18s %9.1f MHz  CID %-6d LAC %-5d PLMN %d-%02d  %-9s gain %-4g found @%.0fs"
            % (chan_name(c), f / 1e6, d["cid"], d["lac"], d["mcc"], d["mnc"], "[%s]" % d["op"],
               d["gain"], d["t"]))


def print_timing(say, title, timer, wall, ncells, startup):
    say("")
    say("----- timing: %s -----" % title)
    if startup:
        say("  %-8s %32s %8.1fs   %s" % ("startup", "", startup, STAGE_HELP["startup"]))
    for st in ["survey", "open", "capture", "detect", "live", "decode", "wait"]:
        n, s = timer.get(st)
        if s == 0 and n == 0:
            continue
        per = ("%4d x avg %6.2fs" % (n, s / n)) if n else ""
        say("  %-8s %32s %8.1fs   %s" % (st, per, s, STAGE_HELP[st]))
    rate = ncells / (wall / 60.0) if wall > 0 else 0
    say("  %-8s %32s %8.1fs   %d cell(s), %.1f cells/min" % ("TOTAL", "wall clock", wall, ncells, rate))


def print_compare(say, runs):
    say("")
    say("================ WIDE vs NARROW ================")
    names = [r[0] for r in runs]
    say("  %-22s %12s %12s" % ("", names[0], names[1]))
    for st in ["survey", "open", "capture", "detect", "live", "decode", "wait"]:
        vals = [r[2].get(st)[1] for r in runs]
        if any(vals):
            say("  %-22s %11.1fs %11.1fs" % (st, vals[0], vals[1]))
    say("  %-22s %11.1fs %11.1fs" % ("TOTAL wall", runs[0][3], runs[1][3]))
    say("  %-22s %12d %12d" % ("cells found", len(runs[0][1]), len(runs[1][1])))
    w, n = runs[0][3], runs[1][3]
    if w > 0 and n > 0:
        faster = names[0] if w < n else names[1]
        say("  -> %s was %.1fx faster" % (faster, max(w, n) / min(w, n)))
    a, b = set(runs[0][1]), set(runs[1][1])
    for label, only, src in ((names[0], a - b, runs[0][1]), (names[1], b - a, runs[1][1])):
        for f in sorted(only):
            d = src[f]
            say("  only in %-6s %-18s %9.1f MHz  CID %d" % (label, chan_name(d["chan"]), f / 1e6, d["cid"]))


def write_csv(path, mode, args, timer, wall, ncells, startup):
    path = os.path.expanduser(path)
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["when", "mode", "bands", "gains", "rate_msps", "secs", "jobs", "pipeline",
                        "startup_s", "survey_s", "open_s", "capture_s", "detect_s", "live_s", "decode_s", "wait_s",
                        "total_s", "cells"])
        w.writerow([time.strftime("%Y-%m-%d %H:%M:%S"), mode, args.band_desc,
                    "/".join("%g" % g for g in args.gain_list), "%.2f" % (args.rate / 1e6),
                    args.secs if mode == "wide" else args.narrow_secs, args.jobs, int(args.pipeline),
                    "%.2f" % startup] +
                   ["%.2f" % timer.get(s)[1] for s in ("survey", "open", "capture", "detect", "live", "decode", "wait")] +
                   ["%.2f" % wall, ncells])


# ----------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description="GSM band scanner for BladeRF (any band, wide or narrow).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    g = ap.add_argument_group("what to scan")
    g.add_argument("--band", action="append", default=None,
                   help="band(s) or group: %s | %s (default gsm900)"
                   % (",".join(BANDS), ",".join(GROUPS)))
    g.add_argument("--range", action="append", default=[], metavar="LO-HI",
                   help="custom MHz range, repeatable, e.g. --range 1805-1830")
    g.add_argument("--chunks", default=None, help="override wide chunk centres, MHz comma list")
    m = ap.add_argument_group("mode")
    m.add_argument("--mode", choices=["wide", "narrow"], default="wide")
    m.add_argument("--wide", dest="mode", action="store_const", const="wide", help="= --mode wide")
    m.add_argument("--narrow", dest="mode", action="store_const", const="narrow", help="= --mode narrow")
    m.add_argument("--compare", action="store_true", help="run wide then narrow and compare timings")
    m.add_argument("--detect", choices=["fft", "all"], default="fft",
                   help="fft = only decode channels with energy (default); all = try every channel")
    r = ap.add_argument_group("radio")
    r.add_argument("--gains", default="40", help="comma list, e.g. 35,40,45")
    r.add_argument("--gain-sweep", action="store_true", help="shortcut for --gains 35,40,45")
    r.add_argument("--rate", type=float, default=15.36e6, help="wide capture rate (default 15.36e6)")
    r.add_argument("--secs", type=float, default=4.0, help="wide capture length per chunk (default 4)")
    r.add_argument("--edge", type=float, default=0.1,
                   help="fraction of rate ignored at each capture edge (default 0.1)")
    r.add_argument("--narrow-rate", type=float, default=2e6, help="narrow stream rate (default 2e6)")
    r.add_argument("--narrow-secs", type=float, default=2.0,
                   help="narrow: max wait for SI3 per channel; SI3 repeats every ~0.94 s (default 2)")
    r.add_argument("--osmo-args", default="bladerf=0", help="narrow: osmosdr device args (default bladerf=0)")
    r.add_argument("--survey-secs", type=float, default=0.5,
                   help="narrow: length of the wide FFT survey per chunk (default 0.5)")
    r.add_argument("--lib", default="/usr/local/lib", help="LD_LIBRARY_PATH for bladeRF-cli ('' = none)")
    d = ap.add_argument_group("detection / decoding")
    d.add_argument("--snr", type=float, default=6.0, help="FFT detection threshold dB (default 6)")
    d.add_argument("--noise-pct", type=float, default=25.0,
                   help="percentile of the spectrum used as noise floor (default 25; old = 50)")
    d.add_argument("--detect-secs", type=float, default=1.0, help="seconds of capture to FFT (default 1)")
    d.add_argument("--decode-secs", type=float, default=2.0,
                   help="wide: max seconds of samples to decode per channel, 0 = all (default 2)")
    d.add_argument("--jobs", type=int, default=2,
                   help="wide: parallel decoders (default 2; 4 dropped ssh on the Pi; 1 = none)")
    d.add_argument("--pipeline", action="store_true",
                   help="capture the next chunk while the previous one decodes (watch for overflows)")
    d.add_argument("--chan-cutoff", type=float, default=250.0,
                   help="wide: channelizer low-pass cutoff kHz (default 250; old script used 100)")
    d.add_argument("--chan-tw", type=float, default=500.0,
                   help="wide: channelizer transition kHz (default 500; old script used 30 = ~1700 taps)")
    o = ap.add_argument_group("output")
    o.add_argument("--log", default=None, help="append all output to this file")
    o.add_argument("--timing-csv", default=None, help="append one timing row per run to this CSV")
    o.add_argument("--keep", action="store_true", help="keep capture files in /dev/shm")
    o.add_argument("-v", "--verbose", action="store_true", help="also print channels that did not decode")
    o.add_argument("--dry-run", action="store_true", help="print the scan plan and exit (no radio)")
    args = ap.parse_args()

    band_specs = args.band if args.band else ([] if args.range else ["gsm900"])
    band_names = expand_bands(band_specs)
    ranges = [parse_range(x) for x in args.range]
    args.band_desc = "+".join(band_names + ["%g-%g" % (lo / 1e6, hi / 1e6) for lo, hi in ranges])
    args.gain_list = [35.0, 40.0, 45.0] if args.gain_sweep else \
        [float(x) for x in args.gains.split(",") if x.strip()]
    args.chunk_centres = [float(x) * 1e6 for x in args.chunks.split(",") if x.strip()] if args.chunks else None
    args.jobs = max(1, args.jobs)
    chans = build_channels(band_names, ranges)
    if not chans:
        raise SystemExit("no channels to scan")

    logf = open(os.path.expanduser(args.log), "a") if args.log else None

    def say(s):
        print(s, flush=True)
        if logf:
            logf.write(s + "\n")
            logf.flush()

    modes = ["wide", "narrow"] if args.compare else [args.mode]
    plan = plan_chunks(chans, args.rate, args.edge, args.chunk_centres)
    say("bandscan %s | %s | %d channels %.1f-%.1f MHz | modes %s | gains %s | jobs %d%s"
        % (time.strftime("%Y-%m-%d %H:%M:%S"), args.band_desc, len(chans), chans[0].freq / 1e6,
           chans[-1].freq / 1e6, "+".join(modes), ",".join("%g" % x for x in args.gain_list),
           args.jobs, " | pipeline" if args.pipeline else ""))
    for i, (centre, members) in enumerate(plan):
        say("  chunk %d: centre %8.2f MHz  covers %8.1f-%8.1f  (%d ch)"
            % (i + 1, centre / 1e6, members[0].freq / 1e6, members[-1].freq / 1e6, len(members)))
    wide_mb = args.rate * args.secs * 4 / 1e6
    say("  wide: %d capture(s) x %.1fs @ %.2f Msps (%.0f MB each)  |  narrow: %d survey(s) x %.1fs, then "
        "live stream <=%.1fs per live channel @ %.1f Msps"
        % (len(plan) * len(args.gain_list), args.secs, args.rate / 1e6, wide_mb,
           len(plan) * len(args.gain_list), args.survey_secs, args.narrow_secs, args.narrow_rate / 1e6))
    if args.dry_run:
        return 0
    try:
        free = shutil.disk_usage(SHM).free / 1e6
        need = wide_mb * (len(plan) if args.pipeline else 1)
        if free < need * 1.2:
            say("  WARNING: /dev/shm has %.0f MB free, a wide capture needs ~%.0f MB (lower --secs/--rate)"
                % (free, need))
    except OSError:
        pass

    pool, startup = None, 0.0
    if args.jobs > 1 and "wide" in modes:
        t0 = time.time()
        pool = mp.get_context("spawn").Pool(args.jobs, initializer=_worker_init)
        pool.map(_warm, range(args.jobs))
        startup = time.time() - t0
        say("  %d decoder workers ready in %.1fs" % (args.jobs, startup))

    runs = []
    try:
        scanner = Scanner(args, chans, say, pool)
        for mode in modes:
            found, timer, wall = scanner.run(mode)
            print_cells(say, "%s mode, %s" % (mode.upper(), args.band_desc), found)
            print_timing(say, "%s mode" % mode, timer, wall, len(found), startup)
            if args.timing_csv:
                write_csv(args.timing_csv, mode, args, timer, wall, len(found), startup)
            runs.append((mode, found, timer, wall))
        if len(runs) == 2:
            print_compare(say, runs)
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if logf:
            logf.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
