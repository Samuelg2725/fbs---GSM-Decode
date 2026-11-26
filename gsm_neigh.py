#!/usr/bin/env python3
"""GSM neighbour lists: which channels each cell tells phones to measure (System Information 2, 2bis and
2ter), and a cell search that follows those lists from cell to cell.

  Read one cell's identity and neighbour list:
      python3 gsm_neigh.py --arfcn 13
  Start from known cells (e.g. one per operator from a scan) and follow every neighbour list:
      python3 gsm_neigh.py --arfcn 13 119 --follow

Operators usually list their whole channel block (e.g. Vodafone 1-22, O2 113-123), so "listed but not
decoded" channels are normally just frequencies with no cell nearby. To find cells with no starting point,
run a band scan first (gsm_arfcn_fast.py --band 900). Run with LD_LIBRARY_PATH=/usr/local/lib.
"""
import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C                     # noqa: E402
import gsm_arfcn_fast as G             # noqa: E402
import pcap_out                        # noqa: E402
import gsm_freqlist as FL              # noqa: E402


def ranges(nums):
    """[1,2,3,5,7,8] -> '1-3 5 7-8'."""
    nums = sorted(set(nums))
    if not nums:
        return "-"
    out, start, prev = [], nums[0], nums[0]
    for n in nums[1:] + [None]:
        if n is not None and n == prev + 1:
            prev = n
            continue
        out.append(str(start) if start == prev else "%d-%d" % (start, prev))
        if n is not None:
            start = prev = n
    return " ".join(out)


def visit(radio, args, arfcn):
    """Tune to one channel; return its cell identity (or None) and neighbour lists."""
    freq = G.downlink(arfcn)
    out = {"arfcn": arfcn, "cell": None, "lists": {}, "formats": {}, "synced": False, "secs": 0.0}
    if freq is None:
        return out
    rx = radio.tune(freq)
    t0, seen, want_ext = time.time(), 0, False
    while time.time() - t0 < args.dwell:
        if not out["synced"]:
            out["synced"] = rx.synced()
            if not out["synced"] and time.time() - t0 > args.sync_secs:
                break
        seen, msgs = G.new_messages(rx, seen)
        for m in msgs:
            raw = G.blob_of(m)
            if out["cell"] is None:
                out["cell"] = G.si3_cell(raw, arfcn, freq)
            nb = FL.neighbours(raw)
            if nb:
                name, fmt, lst, ext = nb
                out["lists"][name], out["formats"][name] = lst, fmt
                if name == "SI2" and ext:
                    want_ext = True           # the list continues in SI2bis
        have_lists = "SI2" in out["lists"] and (not want_ext or "SI2bis" in out["lists"])
        if out["cell"] and have_lists and time.time() - t0 > args.min_secs:
            break                              # SI2ter has its own schedule: give it --min-secs
        time.sleep(0.02)
    out["secs"] = time.time() - t0
    radio.stop()
    return out


def main():
    ap = argparse.ArgumentParser(description="GSM neighbour lists and neighbour-guided cell search.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--arfcn", type=int, nargs="+", required=True, help="cell(s) to read / start from")
    ap.add_argument("--follow", action="store_true", help="visit every listed neighbour, then theirs, and so on")
    ap.add_argument("--gain", type=float, default=20.0, help="receive gain in dB (default 20)")
    ap.add_argument("--dwell", type=float, default=5.0,
                    help="longest wait per channel for the cell ID and neighbour lists, s (default 5)")
    ap.add_argument("--min-secs", type=float, default=2.5,
                    help="minimum time on a cell so SI2bis / SI2ter can arrive, s (default 2.5)")
    ap.add_argument("--sync-secs", type=float, default=1.2, help="leave a channel without GSM sync after this")
    ap.add_argument("--retry-dwell", type=float, default=8.0,
                    help="listed channels that showed GSM signal but no cell ID get one more try this long (0 = off)")
    pcap_out.add_args(ap)
    ap.add_argument("--no-color", action="store_true", help="plain text output")
    args = ap.parse_args()
    args.agc, args.bw = False, None
    G.PCAP[0] = pcap_out.open_from_args(args, "gsm_neigh")
    if args.no_color:
        C.set_colour(False)
    C.stop_on_first_ctrl_c()
    G.stop_other_scanners()
    radio = G.Radio(args)
    queue, done, cells, listed_by = list(args.arfcn), {}, {}, {}
    t_all = time.time()

    def report(a, r, retry=False):
        tag = C.c(" (second try)", "grey") if retry else ""
        if r["cell"]:
            cells[a] = r
            print("%s%s" % (G.cell_line(r["cell"], r["secs"]), tag), flush=True)
            lists = sorted(set(x for lst in r["lists"].values() for x in lst))
            for name in ("SI2", "SI2bis", "SI2ter"):
                if name in r["lists"]:
                    print("        %s neighbour list: %s" % (name, ranges(r["lists"][name])), flush=True)
            if not r["lists"]:
                C.warn("        no neighbour list received in %.1f s (try --dwell 8)" % r["secs"])
            for x in lists:
                listed_by.setdefault(x, set()).add(a)
            if args.follow:
                new = [x for x in lists if x not in done and x not in queue]
                if new:
                    C.info("        following %d new channel(s): %s" % (len(new), ranges(new)))
                queue.extend(new)
        else:
            text = "GSM signal, no cell ID decoded" if r["synced"] else "no GSM signal"
            print("%s          ARFCN %4d  %s%s" % (C.c(C.now_stamp(), "grey"), a,
                                                  C.c(text, "yellow" if r["synced"] else "grey"), tag), flush=True)
    try:
        while queue:
            a = queue.pop(0)
            if a in done:
                continue
            done[a] = visit(radio, args, a)
            report(a, done[a])
        # second chance for listed channels that had GSM signal but no cell ID, skipping channels right next
        # to a decoded cell (that signal is usually spill-over from the neighbour)
        again = [a for a, r in done.items() if r["synced"] and not r["cell"] and a in listed_by
                 and not any(abs(a - c) == 1 for c in cells)]
        if again and args.retry_dwell > 0:
            C.info("second try (%.0f s each) for %s" % (args.retry_dwell, ranges(again)))
            keep, args.dwell = args.dwell, args.retry_dwell
            for a in again:
                r = visit(radio, args, a)
                if r["cell"]:
                    done[a] = r
                report(a, r, retry=True)
            args.dwell = keep
    except KeyboardInterrupt:
        print("")
    finally:
        C.ignore_ctrl_c()
        radio.close()
    t = C.Table(["ARFCN", "MHz", "Operator", "MCC-MNC", "LAC", "Cell ID", "Power", "Neighbour list"],
                align={k: "right" for k in ("ARFCN", "MHz", "LAC", "Cell ID", "Power")})
    for a in sorted(cells, key=lambda k: (512 <= k <= 885, k < 900, k)):
        x = cells[a]["cell"]
        lst = sorted(set(v for lst in cells[a]["lists"].values() for v in lst))
        t.add([a, "%.1f" % x["mhz"], C.op_cell(x["operator"]), "%03d-%02d" % (x["mcc"], x["mnc"]), x["lac"],
               (x["cid"], "green+bold"), x["power"] if x["power"] is not None else "-", ranges(lst)])
    t.show("GSM cells found: %d   (%d channel(s) visited, %.0f s)" % (len(cells), len(done), time.time() - t_all))
    missing = sorted(x for x in listed_by if x not in cells)
    if missing:
        C.info("Listed as neighbours but no cell decoded here: %s" % ranges(missing))
    pcap_out.report(G.PCAP[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
