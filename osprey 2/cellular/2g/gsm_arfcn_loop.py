#!/usr/bin/env python3
import argparse
import importlib.util
import os
import signal
import sys
import time

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
        if "gsm_arfcn" in cmd or "grgsm_scanner" in cmd:
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

def cid_of(text):
    return int(text.split("CID:")[1].split(",")[0])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arfcn", type=int, required=True)
    parser.add_argument("--cellid", type=int)
    args = parser.parse_args()
    cleared = kill_other_scanners()
    if cleared:
        print("cleared leftover scan %s" % " ".join(str(pid) for pid in cleared))
    if args.cellid is None:
        print("watching ARFCN %d" % args.arfcn)
    else:
        print("watching ARFCN %d for cell %d" % (args.arfcn, args.cellid))
    mod = clean.load_scanner()
    seen_expected = False
    with clean.LOG.open("a", buffering=1) as log:
        while True:
            try:
                found = clean.listen(mod, args.arfcn, log)
            except KeyboardInterrupt:
                raise
            except Exception as exc:
                cleared = kill_other_scanners()
                if cleared:
                    print("cleared leftover scan %s" % " ".join(str(pid) for pid in cleared))
                    time.sleep(1)
                    continue
                print("%s  ARFCN %d  radio error: %s" % (clean.now_stamp(), args.arfcn, exc))
                print("Lime still held. Unplug it, plug it back in.")
                time.sleep(5)
                continue
            if not found:
                if args.cellid is None:
                    print("%s  ARFCN %d  no cell" % (clean.now_stamp(), args.arfcn))
                time.sleep(1)
                continue
            for text in found:
                cid = cid_of(text)
                if args.cellid is None:
                    print(text)
                    continue
                if cid == args.cellid:
                    if not seen_expected:
                        print(text)
                        seen_expected = True
                    continue
                print(text)
                print("ARFCN %d expected %d, got %d" % (args.arfcn, args.cellid, cid))
            time.sleep(1)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("")
        sys.exit(0)
