"""Capture one LTE carrier on a HackRF One and print the cell identity.

    python hackrf_lte_decode.py --self-test
    python hackrf_lte_decode.py --earfcn 6300
    python hackrf_lte_decode.py --freq-mhz 806.0 --save vodafone.cf32
    python hackrf_lte_decode.py --file vodafone.cf32 --freq-mhz 806.0
"""

from __future__ import annotations

import argparse
import datetime
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from cellular_resolver import lte_carrier
from hackrf_capture import DEFAULT_RATE, DEFAULT_SECONDS, LOW_IF_HZ, capture, load_iq, notch_tone, save_capture
from lte_cellid import _derotate, _locate_wide
from lte_pci_search import search_samples
from lte_sib1 import decode_sib1_subframe, identity_from_transport_block, synthesize_sib1_subframe

_BANDWIDTH = {6: "1.4", 15: "3", 25: "5", 50: "10", 75: "15", 100: "20"}
_RESET = "\033[0m"
_BOLD_GREEN = "\033[1;32m"
_BOLD_RED = "\033[1;31m"
_CYAN = "\033[36m"


def _homes() -> list:
    """sudo changes home to /root. The decoder and the bit file stay in the real user's home."""
    homes = []
    sudo_user = os.environ.get("SUDO_USER")
    if sudo_user and sudo_user != "root":
        homes.append(Path("/home") / sudo_user)
    home = Path.home()
    if home not in homes:
        homes.append(home)
    fixed = Path("/home/revector")
    if fixed not in homes:
        homes.append(fixed)
    return homes


def _pdsch() -> Path:
    name = Path("srsRAN_4G") / "build" / "lib" / "examples" / "pdsch_ue"
    for home in _homes():
        path = home / name
        if path.is_file():
            return path
    return _homes()[0] / name


def _sib_files() -> list:
    return [home / "sib1_bits.bin" for home in _homes()]


def _miss(freq_mhz: float, report: dict, capture_timestamp=None) -> str:
    text = (
        "No LTE cell at %.3f MHz\nBest PSS %.1fx at %+.1f kHz, SSS %.2f (needs %.2f)"
        % (
            freq_mhz,
            report["best_ratio"],
            report["best_offset_hz"] / 1000.0,
            report["best_sss"],
            report["sss_min"],
        )
    )
    if capture_timestamp is not None:
        text += "\nCapture timestamp: %r" % (capture_timestamp,)
    return text


def _align(samples, sample_rate: float, pci: int, cfo_hz: float):
    """Put subframe 0 at the start of a new array, with the tuner offset removed."""
    corrected = _derotate(np.asarray(samples), cfo_hz, sample_rate)
    start = _locate_wide(corrected, pci, 1024)
    if start is None:
        return None
    return np.asarray(corrected[start:], dtype=np.complex64)


def _read_mib(text: str) -> dict:
    def field(name):
        match = re.search(r"%s:\s+(\S+)" % name, text)
        return None if match is None else match.group(1)
    prb = field("PRB")
    ports = field("Nof ports")
    sfn = field("SFN")
    if prb is None or "Decoded MIB" not in text:
        return {}
    prb = int(prb)
    return {
        "prb": prb,
        "mhz": _BANDWIDTH.get(prb, "?"),
        "ports": None if ports is None else int(ports),
        "sfn": None if sfn is None else int(sfn.split(",")[0]),
        "phich": field("PHICH Length"),
        "ng": field("PHICH Resources"),
    }


def _clear_sib_bits() -> None:
    for path in _sib_files():
        if path.is_file():
            path.unlink()


def _saved_sib_bits() -> Path | None:
    found = [path for path in _sib_files() if path.is_file()]
    if not found:
        return None
    return max(found, key=lambda path: path.stat().st_mtime)


def _run_broadcast(aligned_path: Path, pci: int, n_prb: int, n_subframes: int):
    """Read the whole aligned recording. The identity is only sent every 20 ms."""
    _clear_sib_bits()
    completed = subprocess.run(
        [
            str(_pdsch()), "-d", "-i", str(aligned_path),
            "-p", str(n_prb), "-c", str(pci), "-P", "2", "-Q", "-n", str(n_subframes),
        ],
        capture_output=True,
        text=True,
    )
    text = completed.stdout + completed.stderr
    sib = _saved_sib_bits()
    payload = sib.read_bytes() if sib is not None else None
    return text, payload


def _broadcast(samples, sample_rate: float, pci: int, cfo_hz: float):
    """srsRAN reads the broadcast block. Returns the MIB fields and the identity."""
    reader = _pdsch()
    if not reader.is_file():
        return {}, None, "The broadcast reader is not built at %s" % reader
    aligned = _align(samples, sample_rate, pci, cfo_hz)
    if aligned is None:
        return {}, None, "The frame start was not found for PCI %d" % pci
    with tempfile.NamedTemporaryFile(suffix=".cf32", delete=False) as handle:
        path = Path(handle.name)
    try:
        aligned.tofile(path)
        subframe = max(1, int(round(sample_rate / 1000.0)))
        n_subframes = max(40, len(aligned) // subframe)
        text = ""
        mib = {}
        payload = None
        for n_prb in (50, 25, 75, 100):
            text, got = _run_broadcast(path, pci, n_prb, n_subframes)
            if got and payload is None:
                payload = got
            found = _read_mib(text)
            if found:
                mib = found
                if payload:
                    break
        if not mib:
            return {}, None, "The broadcast block did not decode"
        if not payload:
            return mib, None, "The broadcast block decoded, but the identity bits were not saved"
        identity = identity_from_transport_block(payload)
        if identity is None:
            return mib, None, "The identity bits did not contain a UK network code"
        return mib, identity, ""
    finally:
        path.unlink(missing_ok=True)


def _full_report(freq_mhz: float, earfcn, hit: dict, mib: dict, identity: dict, capture_timestamp=None, target_cell_id=None) -> str:
    if earfcn is None:
        lines = ["Frequency: %.3f MHz" % freq_mhz]
    else:
        lines = ["Frequency: %.3f MHz  EARFCN %d" % (freq_mhz, earfcn)]
    lines.append(
        "PCI %d  SSS %.2f  offset %+.1f kHz"
        % (hit["pci"], hit["sss_score"], hit["offset_hz"] / 1000.0)
    )
    if mib:
        lines.append(
            "MIB: %d RB (%s MHz), PHICH %s Ng=%s, SFN %s, antennas %s"
            % (mib["prb"], mib["mhz"], mib["phich"], mib["ng"], mib["sfn"], mib["ports"])
        )
    if identity:
        name = identity["operator"] or "unknown network"
        lines.append("Network: %d-%s %s" % (identity["mcc"], identity["mnc_digits"], name))
        lines.append("TAC: %d" % identity["tac"])
        cell_line = "Cell ID: %s%d%s (0x%08X)" % (_CYAN, identity["cell_id"], _RESET, identity["cell_id"])
        if target_cell_id is not None and identity["cell_id"] == target_cell_id:
            cell_line += " %s[TARGET MATCH]%s\a" % (_BOLD_GREEN, _RESET)
        elif target_cell_id is not None:
            cell_line += " %s[NO MATCH]%s" % (_BOLD_RED, _RESET)
        lines.append(cell_line + _RESET)
        lines.append("eNB ID: %s%d%s" % (_CYAN, identity["enb_id"], _RESET))
        lines.append("Sector: %d" % identity["sector_id"])
    if capture_timestamp is not None:
        local = datetime.datetime.fromtimestamp(capture_timestamp)
        clock = "%s.%03d" % (time.strftime("%H:%M:%S", time.localtime(capture_timestamp)), local.microsecond // 1000)
        lines.append("Capture timestamp: %r (%s)" % (capture_timestamp, clock))
    return "\n".join(lines)


def decode_capture(samples, sample_rate: float, freq_mhz: float, earfcn=None, pci=None, capture_timestamp=None, target_cell_id=None) -> str:
    report = search_samples(np.asarray(samples), sample_rate)
    hits = report["hits"]
    if not hits and pci is None:
        return _miss(freq_mhz, report, capture_timestamp)
    if hits:
        hit = hits[0]
    else:
        hit = {"pci": int(pci), "sss_score": report["best_sss"], "offset_hz": 0.0}
    print("Reading the broadcast block for PCI %d..." % hit["pci"], flush=True)
    mib, identity, problem = _broadcast(
        samples, sample_rate, int(hit["pci"]), float(hit["offset_hz"]),
    )
    text = _full_report(freq_mhz, earfcn, hit, mib, identity, capture_timestamp, target_cell_id)
    if problem:
        text += "\n" + problem
    if target_cell_id is not None and not identity:
        text += "\nTarget cell ID %d was not found" % target_cell_id
    extras = [int(item["pci"]) for item in hits[1:]]
    if extras:
        text += "\nAlso heard PCI %s" % ", ".join(str(item) for item in extras)
    return text


def self_test() -> bool:
    """A flat 10 MHz subframe 5 must return the cell identity that was written."""
    cell_id = 7159822
    samples, rate, _nfft = synthesize_sib1_subframe(430, 50, cell_id, mcc=(2, 3, 4), mnc=(1, 5), tac=1234)
    found = decode_sib1_subframe(samples, 430, 50, ng="1/6", sfn=0, sample_rate=rate)
    return bool(
        found
        and found["cell_id"] == cell_id
        and found["mcc"] == 234
        and found["mnc"] == 15
        and found["tac"] == 1234
        and found["enb_id"] == (cell_id >> 8)
        and found["operator"] == "Vodafone"
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Decode PCI, MIB and SIB1 from a HackRF capture.")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--earfcn", type=int)
    parser.add_argument("--freq-mhz", type=float)
    parser.add_argument("--file", help="complex64 or HackRF .cs8 capture")
    parser.add_argument("--save", help="write the baseband complex64 capture here")
    parser.add_argument("--rate", type=float, default=DEFAULT_RATE)
    parser.add_argument("--seconds", type=float, default=DEFAULT_SECONDS)
    parser.add_argument("--lna-gain", type=int, default=32)
    parser.add_argument("--vga-gain", type=int, default=28)
    parser.add_argument("--amp", action="store_true")
    parser.add_argument("--low-if-mhz", type=float, default=LOW_IF_HZ / 1e6)
    parser.add_argument("--notch-mhz", type=float, help="remove a tone at this offset, for an older capture")
    parser.add_argument("--pci", type=int, help="use this physical cell ID when the sync check is short")
    parser.add_argument("--cell-id", type=int, help="known 28-bit cell identity to match against this capture")
    args = parser.parse_args(argv)
    if args.self_test:
        ok = self_test()
        print("self-test %s" % ("OK" if ok else "FAILED"))
        return 0 if ok else 1
    earfcn = args.earfcn
    if earfcn is not None:
        carrier = lte_carrier(earfcn)
        if carrier is None:
            print("EARFCN %d is not a downlink channel" % earfcn, file=sys.stderr)
            return 1
        freq_mhz = carrier.downlink_mhz
    elif args.freq_mhz is not None:
        freq_mhz = float(args.freq_mhz)
    else:
        print("Pass --earfcn, --freq-mhz, or --self-test.", file=sys.stderr)
        return 1
    capture_timestamp = None
    if args.file:
        samples = load_iq(args.file, args.rate, low_if_hz=0.0)
        if args.notch_mhz:
            samples = notch_tone(samples, -abs(args.notch_mhz) * 1e6, args.rate)
            samples = notch_tone(samples, abs(args.notch_mhz) * 1e6, args.rate)
    else:
        samples, capture_timestamp = capture(
            freq_mhz * 1e6,
            sample_rate=args.rate,
            seconds=args.seconds,
            lna_gain=args.lna_gain,
            vga_gain=args.vga_gain,
            amp=args.amp,
            low_if_hz=args.low_if_mhz * 1e6,
        )
        if args.save:
            save_capture(args.save, samples)
            print("Saved %d samples to %s" % (len(samples), args.save))
    print(decode_capture(samples, args.rate, freq_mhz, earfcn, args.pci, capture_timestamp, args.cell_id))
    return 0


if __name__ == "__main__":
    sys.exit(main())
