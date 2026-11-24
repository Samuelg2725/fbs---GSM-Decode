#!/usr/bin/env python3
import importlib.util
import os
import signal
import struct
import sys
import time

import numpy
import pmt
from gnuradio import blocks, gsm

spec = importlib.util.spec_from_file_location("gsm_clean", "/home/revector/gsm_arfcn_clean.py")
clean = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clean)

def other_scanners():
    me = os.getpid()
    found = []
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        pid = int(name)
        if pid == me:
            continue
        try:
            cmd = open("/proc/%d/cmdline" % pid, "rb").read().replace(b"\x00", b" ").decode(errors="replace")
        except OSError:
            continue
        if "python" in cmd and ("gsm_arfcn" in cmd or "grgsm_scanner" in cmd or "gsm_sweep.py" in cmd):
            found.append(pid)
    return found

def kill_other_scanners():
    pids = other_scanners()
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    if not pids:
        return []
    time.sleep(1)
    for pid in pids:
        try:
            os.kill(pid, 0)
        except OSError:
            continue
        try:
            os.kill(pid, signal.SIGKILL)
        except OSError:
            pass
    time.sleep(1)
    return pids

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

def cell_from_body(body, raw, arfcn):
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
    freq = gsm.arfcn.arfcn2downlink(arfcn)
    return clean.cell_text(arfcn, freq, cid, lac, mcc, mnc, pwr), pwr

def freq_offsets(channels_num):
    bins = range(int(-numpy.floor(channels_num / 2)), int(numpy.floor((channels_num + 1) / 2)))
    return numpy.fft.ifftshift(numpy.array(list(bins)) * 2e5)

def heard_arfcn(raw, center_freq, offsets):
    if not raw or len(raw) < 6:
        return None
    field = struct.unpack_from(">H", raw, 4)[0] & 0x3FFF
    if field < len(offsets):
        return int(gsm.arfcn.downlink2arfcn(float(center_freq + float(offsets[field]))))
    if 0 <= field <= 1023 or 512 <= field <= 885:
        return field
    return None

def centers_for(lo, hi):
    out = []
    c = lo + 4
    while c <= hi:
        out.append(c)
        if c + 4 >= hi:
            break
        nxt = c + 8
        c = hi if nxt > hi else nxt
    return out

def all_centers():
    plan = []
    for name, lo, hi in (("900", 1, 124), ("900", 975, 1023), ("1800", 512, 885)):
        for c in centers_for(lo, hi):
            plan.append((name, c))
    return plan

def tune(tb, freq):
    src = getattr(tb, "rtlsdr_source", None) or getattr(tb, "sdr_source", None)
    if src is not None:
        src.set_center_freq(freq - 0.1e6, 0)
    rx = getattr(tb, "wideband_receiver", None)
    decoders = getattr(rx, "receivers_with_decoders", None) if rx else None
    if not decoders:
        return
    for item in decoders.values():
        if hasattr(item, "set_fc"):
            item.set_fc(freq)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Sweep every UK GSM channel (900 + 1800), loop until Ctrl-C")
    parser.add_argument("--radio", choices=["lime", "bladerf"], default="lime",
                        help="lime (default, Owen's rig) or bladerf")
    parser.add_argument("--gain", type=float, default=None, help="default 50 lime / 40 bladerf")
    args = parser.parse_args()
    gain = args.gain if args.gain is not None else (40.0 if args.radio == "bladerf" else 50.0)
    cleared = kill_other_scanners()
    if cleared:
        print("cleared leftover scan %s" % " ".join(str(pid) for pid in cleared), flush=True)
    plan = all_centers()
    samp_rate = 2e6
    channels_num = int(samp_rate / 0.2e6)
    offsets = freq_offsets(channels_num)
    first = gsm.arfcn.arfcn2downlink(plan[0][1])
    print("opening sweep, %d steps  radio %s  gain %g" % (len(plan), args.radio, gain), flush=True)
    mod = clean.load_scanner()
    real_source = None
    if args.radio == "bladerf":
        # grgsm_scanner builds osmosdr.source("numchan=1 " + args); numchan BREAKS the BladeRF,
        # so hand it a plain bladerf=0 source instead
        real_source = mod.osmosdr.source

        def blade_source(*a, **k):
            return real_source(args="bladerf=0")
        mod.osmosdr.source = blade_source
    muted = clean.silence()
    tb = None
    try:
        tb = mod.wideband_scanner(
            rec_len=1.2,
            sample_rate=samp_rate,
            carrier_frequency=first,
            ppm=0,
            gain=gain,
            args="driver=lime,soapy=0,nchan=1",
        )
        if args.radio == "bladerf":
            # grgsm_scanner turns AGC on; the BladeRF works with manual gain (40 clean)
            tb.rtlsdr_source.set_gain_mode(False, 0)
            tb.rtlsdr_source.set_gain(gain, 0)
        tb.dbg = blocks.message_debug()
        tb.msg_connect(tb.wideband_receiver, "msgs", tb.dbg, "store")
    finally:
        clean.restore(*muted)
        if real_source is not None:
            mod.osmosdr.source = real_source
    found = {}
    seen = 0
    label_freq = None
    label_center = None
    label_band = None
    try:
        while True:
            for n, (band, center) in enumerate(plan, 1):
                freq = gsm.arfcn.arfcn2downlink(center)
                if not freq:
                    print("skip channel %d" % center, flush=True)
                    continue
                tune(tb, freq)
                try:
                    tb.gsm_extract_system_info.reset()
                except Exception:
                    pass
                try:
                    tb.head.reset()
                except Exception:
                    pass
                try:
                    seen = tb.dbg.num_messages()
                except Exception:
                    pass
                time.sleep(0.3)
                tb.start()
                tb.wait()
                try:
                    tb.stop()
                except Exception:
                    pass
                time.sleep(0.4)
                try:
                    total = tb.dbg.num_messages()
                except Exception:
                    total = seen
                use_freq = freq if label_freq is None else label_freq
                use_center = center if label_center is None else label_center
                use_band = band if label_band is None else label_band
                print("step %d/%d  %s  channel %d" % (n, len(plan), use_band, use_center), flush=True)
                while seen < total:
                    raw = blob_of(tb.dbg.get_message(seen))
                    seen += 1
                    body = si3_body(raw)
                    if body is None:
                        continue
                    heard = heard_arfcn(raw, use_freq, offsets)
                    if heard is None or heard < 1 or abs(heard - use_center) > 5:
                        continue
                    parsed = cell_from_body(body, raw, heard)
                    if not parsed:
                        continue
                    text, pwr = parsed
                    key = (heard, int(text.split("CID:")[1].split(",")[0]))
                    prev = found.get(key)
                    if prev is None or pwr > prev[0]:
                        found[key] = (pwr, text)
                        if prev is None:
                            print(text, flush=True)
                label_freq = freq
                label_center = center
                label_band = band
            print("", flush=True)
            print("Cells found: %d" % len(found), flush=True)
            for key in sorted(found):
                print(found[key][1], flush=True)
    except KeyboardInterrupt:
        print("", flush=True)
        print("Cells found: %d" % len(found), flush=True)
        for key in sorted(found):
            print(found[key][1], flush=True)
    finally:
        if tb is not None:
            try:
                tb.stop()
            except Exception:
                pass
    return 0

if __name__ == "__main__":
    sys.exit(main())
