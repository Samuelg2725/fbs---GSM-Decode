#!/usr/bin/env python3
"""Read every cell srsRAN can separate on one EARFCN.

This file does not decode. It starts the Pluto pdsch_ue three times, once
for each sync group, then repeats. A sync group is the last digit of the
physical cell ID modulo 3. pdsch_ue returns the strongest cell in that
group, then this script reads its 28-bit identity and moves to the next
group. A known identity is marked MATCH.
"""

from __future__ import annotations

import argparse
import os
import re
import select
import signal
import subprocess
import sys
import time
from datetime import datetime

# band, downlink start in tenths of a MHz, EARFCN offset, first, last.
# F_tenths = start + (earfcn - offset).
_BANDS = (
    (1, 21100, 0, 0, 599),
    (3, 18050, 1200, 1200, 1949),
    (7, 26200, 2750, 2750, 3449),
    (8, 9250, 3450, 3450, 3799),
    (20, 7910, 6150, 6150, 6449),
    (28, 7580, 9210, 9210, 9659),
    (38, 25700, 37750, 37750, 38249),
    (40, 23000, 38650, 38650, 39649),
)

_PDSCH = "/home/revector/srsRAN_4G/build/lib/examples/pdsch_ue"
_IDENTITY = re.compile(
    r"SIB1 cell identity (\d+)\s+eNB (\d+)\s+sector (\d+)\s+TAC (\d+)\s+PLMN ([\d-]+)"
)
_PCI = re.compile(r"\*?Found Cell_id:\s+(\d+)\s+(FDD|TDD).*PSR=([0-9.]+)")
_MIB = re.compile(r"Decoded MIB")
_RATE = re.compile(r"Setting rx sampling rate ([0-9.]+) MHz")
_MISS = re.compile(r"Could not find any cell")


def downlink_hz(earfcn):
    """Hz of one downlink EARFCN, or None when the table above does not cover it."""
    for _band, start, offset, first, last in _BANDS:
        if first <= earfcn <= last:
            return (start + (earfcn - offset)) * 100_000
    return None


def say(text):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S  ") + text, flush=True)


def take_lines(pending, chunk):
    """Split a byte stream on newlines. Return complete text lines and the remainder."""
    pending += chunk.replace(b"\r", b"\n")
    lines = []
    while b"\n" in pending:
        raw, pending = pending.split(b"\n", 1)
        text = raw.decode("utf-8", "replace")
        text = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", text).strip()
        if text:
            lines.append(text)
    return lines, pending


def read_group(earfcn, hz, gain, group, known):
    """Read the strongest cell in one sync group. Return its 28-bit identity, or None."""
    say("group %d  listening" % group)
    command = [
        "stdbuf", "-oL", "-eL",
        _PDSCH,
        "-I", "soapy",
        "-a", "driver=plutosdr",
        "-f", str(hz),
        "-g", str(gain),
        "-Q",
        "-F",
        "-d",
        "-l", str(group),
    ]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=0,
    )
    streams = [process.stdout, process.stderr]
    pending = {process.stdout: b"", process.stderr: b""}
    # 20 s to find a tower. Once one is found, 15 s for the cell ID.
    # A real ID has been arriving inside that window. A longer wait only
    # stretches the ones that never decode.
    started = time.monotonic()
    deadline = started + 20
    limit = started + 40
    last_note = started
    saw_pci = False
    pci_num = None
    try:
        while time.monotonic() < deadline:
            if time.monotonic() - last_note >= 5:
                left = int(deadline - time.monotonic())
                say("group %d  still searching, %d s left" % (group, max(left, 0)))
                last_note = time.monotonic()
            remaining = max(0, deadline - time.monotonic())
            ready, _, _ = select.select(streams, [], [], min(0.5, remaining))
            if not ready:
                if process.poll() is not None:
                    break
                continue
            for stream in ready:
                chunk = stream.read(4096)
                if not chunk:
                    if stream in streams:
                        streams.remove(stream)
                    continue
                lines, pending[stream] = take_lines(pending[stream], chunk)
                for text in lines:
                    if _MISS.search(text):
                        continue
                    pci = _PCI.search(text)
                    if pci and float(pci.group(3)) > 0:
                        pci_num = int(pci.group(1))
                        say("group %d  PCI %s %s  PSR %s" % (group, pci.group(1), pci.group(2), pci.group(3)))
                        last_note = time.monotonic()
                        if not saw_pci:
                            deadline = min(limit, time.monotonic() + 15)
                            say("group %d  tower found, waiting up to 15 s for the cell ID" % group)
                        saw_pci = True
                    if _RATE.search(text):
                        say("group %d  sample rate %s MHz" % (group, _RATE.search(text).group(1)))
                    if _MIB.search(text):
                        say("group %d  MIB decoded, reading SIB1" % group)
                        last_note = time.monotonic()
                    found = _IDENTITY.search(text)
                    if found:
                        if pci_num is None and process.stdout in streams:
                            ready_out, _, _ = select.select([process.stdout], [], [], 0.4)
                            if ready_out:
                                extra = process.stdout.read(4096)
                                if extra:
                                    extra_lines, pending[process.stdout] = take_lines(pending[process.stdout], extra)
                                    for extra_text in extra_lines:
                                        late = _PCI.search(extra_text)
                                        if late and float(late.group(3)) > 0:
                                            pci_num = int(late.group(1))
                        cell = int(found.group(1))
                        pci_text = ("PCI %d  " % pci_num) if pci_num is not None else ""
                        line = (
                            "group %d  %s28-bit %s  (0x%07X)  eNB %s  sector %s  TAC %s  PLMN %s"
                            % (
                                group,
                                pci_text,
                                found.group(1),
                                cell,
                                found.group(2),
                                found.group(3),
                                found.group(4),
                                found.group(5),
                            )
                        )
                        if known:
                            line += "  MATCH" if cell in known else "  not a known cell"
                        say(line)
                        return cell
            if not streams:
                break
        if saw_pci:
            say("group %d  saw a cell, SIB1 did not pass" % group)
        else:
            say("group %d  no cell" % group)
        return None
    finally:
        if process.poll() is None:
            process.send_signal(signal.SIGINT)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        for stream in (process.stdout, process.stderr):
            if stream:
                stream.close()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Read the cells on one EARFCN. This file does not decode."
    )
    parser.add_argument("--earfcn", type=int, required=True, help="LTE downlink channel number")
    parser.add_argument("--gain", type=float, default=70.0, help="Pluto gain passed to pdsch_ue")
    parser.add_argument(
        "--cell-id",
        type=int,
        action="append",
        default=[],
        help="known 28-bit cell identity. Repeat for more than one. A match is marked MATCH",
    )
    parser.add_argument(
        "--group",
        type=int,
        nargs="+",
        default=[0, 1, 2],
        help="which groups to search. Example: --group 0 1   or   --group 2",
    )
    args = parser.parse_args(argv)
    groups = []
    for group in args.group:
        if group not in (0, 1, 2):
            sys.exit("group must be 0, 1 or 2")
        if group not in groups:
            groups.append(group)

    hz = downlink_hz(args.earfcn)
    if hz is None:
        sys.exit("EARFCN %d is outside the band table in this launcher." % args.earfcn)
    if not os.path.isfile(_PDSCH) or not os.access(_PDSCH, os.X_OK):
        sys.exit("pdsch_ue was not found at %s" % _PDSCH)

    say("%d Hz  EARFCN %d" % (hz, args.earfcn))
    say("groups %s, strongest cell in each, then repeat. Ctrl+C stops." % " ".join(str(group) for group in groups))
    if args.cell_id:
        say("known cells: %s" % ", ".join(str(cell) for cell in args.cell_id))
    known = set(args.cell_id)
    try:
        while True:
            for group in groups:
                read_group(args.earfcn, hz, args.gain, group, known)
            say("pass done")
    except KeyboardInterrupt:
        return 0


if __name__ == "__main__":
    sys.exit(main())
