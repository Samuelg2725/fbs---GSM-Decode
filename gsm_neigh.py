#!/usr/bin/env python3
"""GSM neighbour lists: read which frequencies a cell tells phones to measure (SI2/SI2bis/SI2ter),
and use them to find cells without sweeping the whole band.

  LD_LIBRARY_PATH=/usr/local/lib python3 ~/gsm_neigh.py --arfcn 13            # one cell: ID + neighbours
  LD_LIBRARY_PATH=/usr/local/lib python3 ~/gsm_neigh.py --arfcn 13 119 --follow
        # start at 13 and 119, then visit every neighbour they list, then those cells' neighbours ...
  (no starting cell known? run the full sweep - gsm_arfcn_fast --fresh --sync-secs 1.2 - it is the
   reliable way; a shortcut that skipped channels missed a whole operator in testing, so it was removed)

Each channel gets a fresh decoder (same as gsm_arfcn_fast --fresh), gain 20 by default.
A channel is left after --sync-secs without GSM sync; with sync it waits up to --dwell for the cell ID
and the neighbour list(s). Output: per cell its ID and neighbour ARFCNs; at the end every cell found
and every listed neighbour that did NOT decode (weak / far / other operator band).
"""
import argparse
import os
import re
import signal
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.expanduser("~"))
import gsm_arfcn_fast as F          # noqa: E402
import gsm_freqlist as FL           # noqa: E402
from gnuradio import gsm            # noqa: E402

BAND900 = list(range(975, 1024)) + list(range(1, 125))
# The previous (stopped) flowgraph is kept alive until the next one exists: gr-osmosdr then re-uses the
# open BladeRF instead of closing + reopening it (~0.3 s per channel instead of ~2.8 s). This is what
# gsm_arfcn_fast --fresh does implicitly.
_PREV = [None]


def quietly(fn):
    """Run fn() with the radio driver's console chatter hidden. Ctrl-C is held back while output is
    redirected, so it can never leave the terminal silenced; it is delivered right after."""
    signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGINT})
    try:
        muted = F.clean.silence()
        try:
            return fn()
        finally:
            F.clean.restore(*muted)
    finally:
        signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGINT})


def stop_tb(tb):
    try:
        tb.stop()
        tb.wait()
    except Exception:
        pass


def visit(args, arfcn):
    """Tune, decode. Returns dict(cell=text or None, si2={name: arfcns}, synced=bool, secs)."""
    freq = gsm.arfcn.arfcn2downlink(arfcn)
    out = {"arfcn": arfcn, "cell": None, "lists": {}, "formats": {}, "synced": False, "secs": 0.0,
           "si2quater": False}
    if freq is None or freq < 100e6:
        return out
    def open_new():
        new = F.top(freq, args.gain, "bladerf", False, None)
        new.start()
        return new
    tb = quietly(open_new)
    quietly(lambda: _PREV.__setitem__(0, None))   # new radio flowgraph is up: the old one can go now
    t0 = time.time()
    seen = 0
    want_ext = False
    try:
        while time.time() - t0 < args.dwell:
            if not out["synced"]:
                out["synced"] = tb.synced_since(0)
                if not out["synced"] and time.time() - t0 > args.sync_secs:
                    break
            n = tb.dbg.num_messages()
            while seen < n:
                raw = F.blob_of(tb.dbg.get_message(seen))
                seen += 1
                if out["cell"] is None:
                    t = F.si3_text(raw, arfcn, freq)
                    if t:
                        out["cell"] = re.sub("\x1b\\[[0-9;]*m", "", t)
                nb = FL.neighbours(raw)
                if nb:
                    name, fmt, lst, ext = nb
                    out["lists"][name] = lst
                    out["formats"][name] = fmt
                    if name == "SI2" and ext:
                        want_ext = True
                t, _b = FL.rr_body(raw, (0x07,))
                if t:
                    out["si2quater"] = True
            have_lists = "SI2" in out["lists"] and (not want_ext or "SI2bis" in out["lists"])
            # SI2ter / SI2quater come on their own schedule: give them a little extra time once SI2 is in
            if out["cell"] and have_lists and time.time() - t0 > args.min_secs:
                break
            time.sleep(0.02)
    finally:
        out["secs"] = time.time() - t0
        quietly(lambda: stop_tb(tb))
        _PREV[0] = tb
    return out


def fmt_list(lst):
    return " ".join(str(a) for a in lst) if lst else "-"


def main():
    ap = argparse.ArgumentParser(description="GSM neighbour lists (SI2/2bis/2ter) + neighbour-guided cell search.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--arfcn", type=int, nargs="+", help="cell(s) to read / start from")
    ap.add_argument("--follow", action="store_true", help="visit every listed neighbour, then theirs, ...")
    ap.add_argument("--gain", type=float, default=20.0)
    ap.add_argument("--dwell", type=float, default=5.0,
                    help="max seconds per channel once synced (cell ID + neighbour lists; default 5)")
    ap.add_argument("--min-secs", type=float, default=2.5,
                    help="stay at least this long on a cell to catch SI2bis/2ter (default 2.5)")
    ap.add_argument("--sync-secs", type=float, default=1.2, help="leave a channel without GSM sync after this")
    ap.add_argument("--retry-dwell", type=float, default=8.0,
                    help="channels that synced without a cell ID get one more try this long (0 = off)")
    args = ap.parse_args()
    if not args.arfcn:
        sys.exit("give --arfcn N [N ...] (cells to read / start from)")

    F.kill_other_scanners()
    signal.signal(signal.SIGINT, F._first_ctrl_c_only)
    t_all = time.time()
    queue = list(args.arfcn or [])
    done, cells, listed_by = {}, {}, {}

    def handle(a, r, retry=False):
        allnb = sorted(set(x for lst in r["lists"].values() for x in lst))
        tag = " (retry)" if retry else ""
        if r["cell"]:
            cells[a] = r
            print("%s  %4.1fs  %s%s" % (F.clean.now_stamp(), r["secs"], r["cell"], tag), flush=True)
            for name in ("SI2", "SI2bis", "SI2ter"):
                if name in r["lists"]:
                    print("        %-6s (%s): %s" % (name, r["formats"][name], fmt_list(r["lists"][name])),
                          flush=True)
            if not r["lists"]:
                print("        no neighbour list decoded in %.1fs (try --dwell 8)" % r["secs"], flush=True)
            if r["si2quater"]:
                print("        SI2quater present (3G/LTE neighbours - not decoded yet)", flush=True)
            for x in allnb:
                listed_by.setdefault(x, set()).add(a)
            if args.follow:
                new = [x for x in allnb if x not in done and x not in queue]
                if new:
                    print("        -> queue %d new: %s" % (len(new), fmt_list(new)), flush=True)
                queue.extend(new)
        else:
            why = "synced, no cell ID" if r["synced"] else "no GSM sync"
            print("%s  %4.1fs  ARFCN %d  %s%s" % (F.clean.now_stamp(), r["secs"], a, why, tag), flush=True)

    try:
        while queue:
            a = queue.pop(0)
            if a in done:
                continue
            r = visit(args, a)
            done[a] = r
            handle(a, r)
        # second chance for channels that locked but gave no cell ID (weak cells, e.g. ARFCN 1)
        # (only channels a cell lists as a neighbour, and not right next to a decoded cell - those are spill-over)
        again = [a for a, r in done.items() if r["synced"] and not r["cell"] and a in listed_by
                 and not any(abs(a - c) == 1 for c in cells)]
        if again and args.retry_dwell > 0:
            print("retrying %d channel(s) that synced without a cell ID (dwell %.0fs): %s"
                  % (len(again), args.retry_dwell, fmt_list(again)), flush=True)
            keep = args.dwell
            args.dwell = args.retry_dwell
            for a in again:
                r = visit(args, a)
                if r["cell"]:
                    done[a] = r
                handle(a, r, retry=True)
            args.dwell = keep
    except KeyboardInterrupt:
        pass
    finally:
        try:
            signal.signal(signal.SIGINT, signal.SIG_IGN)
        except KeyboardInterrupt:
            pass
    quietly(lambda: _PREV.__setitem__(0, None))
    el = time.time() - t_all
    print("\n===== %d cell(s) in %.0fs, %d channel(s) visited =====" % (len(cells), el, len(done)))
    for a in sorted(cells, key=lambda k: (k < 900 and k > 500, k < 900, k)):
        print("  " + cells[a]["cell"].split(", Pwr")[0])
    missing = sorted(x for x in listed_by if x not in cells)
    if missing:
        print("listed as neighbours but not decoded here: %s" % fmt_list(missing))
    return 0


if __name__ == "__main__":
    sys.exit(main())
