#!/usr/bin/env python3
import gc
import importlib.machinery
import importlib.util
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy
from gnuradio import gsm

OPS = {
    (234, 10): "O2",
    (234, 15): "Vodafone",
    (234, 20): "Three",
    (234, 30): "EE",
    (234, 33): "EE",
}
LOG = Path.home() / "gsm_cells.txt"

def load_scanner():
    loader = importlib.machinery.SourceFileLoader("grgsm_scanner_mod", "/usr/bin/grgsm_scanner")
    spec = importlib.util.spec_from_file_location("grgsm_scanner_mod", "/usr/bin/grgsm_scanner", loader=loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod

def now_stamp():
    t = datetime.now(timezone.utc)
    return t.strftime("%Y-%m-%dT%H:%M:%S.") + "%06dZ" % t.microsecond

def cell_text(arfcn, freq, cid, lac, mcc, mnc, pwr):
    name = OPS.get((mcc, mnc), "")
    label = (" " + name) if name else ""
    return "ARFCN: %4u, Freq: %6.1fM, CID: %5u, LAC: %5u, MCC: %3u, MNC: %3u%s, Pwr: %3i" % (
        arfcn, freq / 1e6, cid, lac, mcc, mnc, label, pwr)

def silence():
    sys.stdout.flush()
    sys.stderr.flush()
    null_fds = [os.open(os.devnull, os.O_RDWR) for _ in range(2)]
    save = os.dup(1), os.dup(2)
    os.dup2(null_fds[0], 1)
    os.dup2(null_fds[1], 2)
    return null_fds, save

def restore(null_fds, save):
    sys.stdout.flush()
    sys.stderr.flush()
    os.dup2(save[0], 1)
    os.dup2(save[1], 2)
    for fd in list(null_fds) + list(save):
        os.close(fd)

def listen(mod, arfcn, log, gain=50.0, seconds=3.0):
    freq = gsm.arfcn.arfcn2downlink(arfcn)
    if freq is None or freq < 100e6:
        return []
    samp_rate = 2e6
    channels_num = int(samp_rate / 0.2e6)
    muted = silence()
    tb = None
    try:
        tb = mod.wideband_scanner(
            rec_len=seconds,
            sample_rate=samp_rate,
            carrier_frequency=freq,
            ppm=0,
            gain=gain,
            args="driver=lime,soapy=0,nchan=1",
        )
        tb.start()
        tb.wait()
        tb.stop()
        stamp = now_stamp()
        info = tb.gsm_extract_system_info
        chans = list(info.get_chans() or [])
        if not chans:
            return []
        freq_offsets = numpy.fft.ifftshift(
            numpy.array(range(int(-numpy.floor(channels_num / 2)),
                              int(numpy.floor((channels_num + 1) / 2)))) * 2e5)
        chans_a = numpy.array(chans)
        found_freqs = freq + freq_offsets[chans_a]
        cids = numpy.array(info.get_cell_id())
        lacs = numpy.array(info.get_lac())
        mccs = numpy.array(info.get_mcc())
        mncs = numpy.array(info.get_mnc())
        pwrs = numpy.array(info.get_pwrs())
        found = []
        for i in range(len(chans_a)):
            got = int(gsm.arfcn.downlink2arfcn(float(found_freqs[i])))
            if got != arfcn:
                continue
            cid = int(cids[i])
            if cid == 0:
                continue
            text = "%s  %s" % (stamp, cell_text(
                arfcn, freq, cid, int(lacs[i]), int(mccs[i]), int(mncs[i]), int(pwrs[i])))
            log.write(text + "\n")
            log.flush()
            found.append(text)
        return found
    finally:
        if tb is not None:
            try:
                tb.stop()
            except Exception:
                pass
            del tb
        gc.collect()
        restore(*muted)

def main():
    if len(sys.argv) < 2:
        print("usage: gsm_arfcn_clean.py ARFCN [ARFCN ...]")
        return 1
    mod = load_scanner()
    rows = []
    seen = set()
    with LOG.open("a", buffering=1) as log:
        for n, arg in enumerate(sys.argv[1:]):
            if n:
                time.sleep(1)
            for text in listen(mod, int(arg), log):
                cid = int(text.split("CID:")[1].split(",")[0])
                if cid in seen:
                    continue
                seen.add(cid)
                rows.append(text)
                print(text)
    print("")
    print("Cells found: %d" % len(rows))
    for text in rows:
        print(text)
    return 0

if __name__ == "__main__":
    sys.exit(main())
