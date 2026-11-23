"""Read Kraken channel 0 and print the LTE physical cell ID.

CellSearch cannot open the Kraken tuners. They are already owned by the
direction finder. This program copies one channel, resamples 2.4 MHz to the
1.92 MHz rate that LTE cell search expects, and runs the PSS/SSS search.
"""

from __future__ import annotations

import argparse
import json
import socket
import sys
import time
from pathlib import Path

import numpy as np
from scipy import signal

from lte_sync import FS, search_cell

TARGET_FILE = Path("active_target.json")
CAPTURE_S = 0.08  # CellSearch uses 153600 samples, which is 80 ms at 1.92 MHz
DECODER_UDP = "127.0.0.1:5005"
CELLID_SSS_MIN = 0.75
# (freq_mhz, pci) -> bandwidth, cell id, status, attempts. Not cleared between loops.
_CELLID_CACHE = {}


def resample_to_lte(samples: np.ndarray, input_hz: float) -> np.ndarray:
    """2.4e6 -> 1.92e6 is 4/5. Any other rate uses the same rational resample."""
    gcd = np.gcd(int(input_hz), FS)
    up = FS // gcd
    down = int(input_hz) // gcd
    if up == 1 and down == 1:
        return np.asarray(samples, dtype=np.complex64)
    return signal.resample_poly(samples, up, down).astype(np.complex64)


def load_target(path: Path):
    """PCI and CellID from active_target.json. Missing file means no target yet."""
    if not path.exists():
        return None, None
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    pci = data.get("target_pci")
    cell_id = data.get("target_cell_id")
    return (None if pci is None else int(pci), None if cell_id is None else int(cell_id))


def load_target_pci(path: Path):
    pci, _cell_id = load_target(path)
    return pci


def select_cells(hits, target_pci):
    """Prefer the target PCI when it cleared the sync check. Otherwise the strongest.

    hits are strongest-first. secondary_pcis keeps the rest in that same order.
    """
    valid = [hit for hit in hits if hit.get("pci") is not None]
    if not valid:
        return None, [], False
    chosen = valid[0]
    prioritized = False
    if target_pci is not None:
        for hit in valid:
            if int(hit["pci"]) == int(target_pci):
                prioritized = hit is not valid[0]
                chosen = hit
                break
    secondary = []
    chosen_pci = int(chosen["pci"])
    for hit in valid:
        pci = int(hit["pci"])
        if pci != chosen_pci and pci not in secondary:
            secondary.append(pci)
    return chosen, secondary, prioritized


def decoder_packet(centre_hz, hit, target_pci, target_cell_id, now=None, secondary_pcis=None, target_prioritized=False, cellid=None):
    """One UDP report. A miss is pci null so the direction finder sees the attempt."""
    pci = None if hit is None else int(hit["pci"])
    match = target_pci is not None and pci is not None and pci == int(target_pci)
    bandwidth = None if not cellid else cellid.get("bandwidth_mhz")
    decoded = None if not cellid else cellid.get("cell_id")
    cell_match = (
        target_cell_id is not None
        and decoded is not None
        and int(decoded) == int(target_cell_id)
    )
    return {
        "source": "lte_decoder",
        "timestamp": time.time() if now is None else now,
        "freq_mhz": None if centre_hz is None else round(float(centre_hz) / 1e6, 3),
        "pci": pci,
        "sss": None if hit is None else round(float(hit["sss_score"]), 2),
        "offset_khz": None if hit is None else round(float(hit["offset_hz"]) / 1000.0, 1),
        "target_pci": None if target_pci is None else int(target_pci),
        "target_cell_id": None if target_cell_id is None else int(target_cell_id),
        "is_target_match": match,
        "secondary_pcis": [int(item) for item in (secondary_pcis or [])],
        "target_prioritized": bool(target_prioritized and match),
        "mib_bw_mhz": None if bandwidth is None else float(bandwidth),
        "decoded_cell_id": None if decoded is None else int(decoded),
        "cell_id_match": cell_match,
    }


def remember_cellid(cache, freq_mhz, hit, decode):
    """Run the cell-identity decoder at most twice per frequency and PCI.

    A strong SSS lock is required. A stored cell id is never decoded again.
    decode is called with no arguments and returns the decode_cell_identity dict.
    """
    if hit is None or freq_mhz is None or float(hit.get("sss_score") or 0.0) < CELLID_SSS_MIN:
        if hit is None or freq_mhz is None:
            return None
        return cache.get((round(float(freq_mhz), 3), int(hit["pci"])))
    key = (round(float(freq_mhz), 3), int(hit["pci"]))
    entry = cache.get(key)
    if entry is not None and (entry.get("cell_id") is not None or entry.get("attempts", 0) >= 2):
        return entry
    result = decode()
    entry = {
        "bandwidth_mhz": result.get("bandwidth_mhz"),
        "cell_id": result.get("cell_id"),
        "cell_id_status": result.get("cell_id_status"),
        "attempts": 1 if entry is None else entry["attempts"] + 1,
    }
    cache[key] = entry
    return entry


def _udp_endpoint(text: str):
    host, separator, port = text.rpartition(":")
    if not separator or not host:
        raise ValueError("UDP address must be host:port")
    return host, int(port)


def _print_packet(packet: dict) -> None:
    if packet["pci"] is None:
        where = "unknown" if packet["freq_mhz"] is None else "%.3f MHz" % packet["freq_mhz"]
        print("%s no PCI" % where, flush=True)
        return
    where = "unknown" if packet["freq_mhz"] is None else "%.3f MHz" % packet["freq_mhz"]
    tail = "offset %+.1f kHz SSS %.2f" % (packet["offset_khz"], packet["sss"])
    others = packet.get("secondary_pcis") or []
    if not others:
        mark = " MATCH" if packet["is_target_match"] else ""
        print("%s PCI %d%s %s" % (where, packet["pci"], mark, tail), flush=True)
        return
    labels = []
    if packet.get("target_prioritized"):
        labels.append("PCI %d (Strongest)" % others[0])
        labels.append("PCI %d [TARGET MATCH PRIORITIZED]" % packet["pci"])
        labels.extend("PCI %d" % pci for pci in others[1:])
    elif packet["is_target_match"]:
        labels.append("PCI %d [TARGET MATCH]" % packet["pci"])
        labels.extend("PCI %d" % pci for pci in others)
    else:
        labels.append("PCI %d (Strongest)" % packet["pci"])
        labels.extend("PCI %d" % pci for pci in others)
    count = 1 + len(others)
    print("%s Found %d cells: %s %s" % (where, count, ", ".join(labels), tail), flush=True)


def _note_cellid(samples, rate: float, freq_mhz, hit) -> None:
    def _decode():
        from lte_cellid import decode_cell_identity
        return decode_cell_identity(
            samples,
            int(hit["pci"]),
            cfo_hz=float(hit.get("offset_hz") or 0.0),
            sample_rate=rate,
        )

    remember_cellid(_CELLID_CACHE, freq_mhz, hit, _decode)


def _one_search(target_path: Path) -> dict:
    try:
        target_pci, target_cell_id = load_target(target_path)
    except (OSError, ValueError) as exc:
        print("Could not read %s: %s" % (target_path, exc), flush=True)
        target_pci, target_cell_id = None, None
    try:
        samples, rate, centre_hz = capture_channel0()
        report = search_samples(samples, rate, target_pci)
        chosen, secondary, prioritized = select_cells(report["hits"], target_pci)
        freq_mhz = None if centre_hz is None else centre_hz / 1e6
        for hit in report["hits"]:
            _note_cellid(samples, rate, freq_mhz, hit)
        cached = None if chosen is None or freq_mhz is None else _CELLID_CACHE.get((round(freq_mhz, 3), int(chosen["pci"])))
        return decoder_packet(
            centre_hz,
            chosen,
            target_pci,
            target_cell_id,
            secondary_pcis=secondary,
            target_prioritized=prioritized,
            cellid=cached,
        )
    except Exception as exc:
        print("Search failed: %s" % exc, flush=True)
        return decoder_packet(None, None, target_pci, target_cell_id)


def run_loop(target_path: Path, interval: float, udp_text: str) -> int:
    """Keep decoding and send every result, including a miss, to the direction finder."""
    host, port = _udp_endpoint(udp_text)
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    print(
        "LTE decode loop. UDP %s:%d, %.1fs after each search." % (host, port, interval),
        flush=True,
    )
    try:
        while True:
            packet = _one_search(target_path)
            sock.sendto(json.dumps(packet).encode("utf-8"), (host, port))
            _print_packet(packet)
            time.sleep(interval)
    except KeyboardInterrupt:
        print("Stopped.")
        return 0
    finally:
        sock.close()


def capture_channel0() -> tuple[np.ndarray, float, float]:
    """One channel from the live Kraken buffers. Requires sudo.

    Returns samples, sample rate, and the tuned centre in Hz. Calibration
    frames are skipped: the noise source is not an LTE downlink.
    """
    from kraken_shm_reader import FRAME_TYPE_DATA, KrakenShmReader

    reader = KrakenShmReader()
    reader.open()
    try:
        deadline = time.time() + 3.0
        while time.time() < deadline:
            header, _mapped = reader._newest_header()
            if header is None:
                time.sleep(0.05)
                continue
            if header.noise_source_state or header.frame_type != FRAME_TYPE_DATA:
                time.sleep(0.05)
                continue
            try:
                mapped = _mapped_for(reader, header)
            except RuntimeError:
                time.sleep(0.02)
                continue
            channel = reader._copy_channel0(mapped, header)
            if channel is None:
                time.sleep(0.02)
                continue
            rate = float(header.sampling_freq)
            need = int(round(CAPTURE_S * rate))
            if len(channel) < need:
                raise RuntimeError("Kraken frame is shorter than 80 ms")
            return channel[:need], rate, float(header.rf_center_freq)
        raise RuntimeError("no Kraken data frame is ready (it may be calibrating)")
    finally:
        reader.close()


def _mapped_for(reader, header):
    for mapped in reader._maps:
        again = reader._header(mapped)
        if again is not None and again.cpi_index == header.cpi_index:
            return mapped
    raise RuntimeError("Kraken frame changed while it was being read")


def search_samples(samples: np.ndarray, input_hz: float, target_pci=None) -> dict:
    report = search_cell(resample_to_lte(samples, input_hz))
    for item in report["hits"]:
        item["matches_target"] = target_pci is not None and item["pci"] == target_pci
    return report


def _band_note(centre_hz: float) -> str:
    mhz = centre_hz / 1e6
    if 832.0 <= mhz <= 862.0:
        return "This is Band 20 uplink. The base station downlink is 791-821 MHz."
    if 880.0 <= mhz < 925.0:
        return "This is 900 MHz uplink. LTE downlink here is 925-960 MHz, or Band 20 at 791-821 MHz."
    if 791.0 <= mhz <= 821.0:
        return "Band 20 downlink. EE is centred at 798.5, Vodafone at 806.0, O2 at 816.0."
    if 925.0 <= mhz <= 960.0:
        return "Band 8 downlink. The centre has to be the middle of a cell, not the edge of the block."
    return ""


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Find an LTE physical cell ID from Kraken channel 0.")
    parser.add_argument("--self-test", action="store_true", help="decode a synthetic PCI 320 and exit")
    parser.add_argument("--file", help="complex64 file recorded at --rate, instead of the live Kraken")
    parser.add_argument("--rate", type=float, default=2_400_000.0)
    parser.add_argument("--target", default=str(TARGET_FILE), help="active_target.json written by kraken_commander")
    parser.add_argument("--loop", action="store_true", help="keep decoding and send each result on UDP")
    parser.add_argument("--interval", type=float, default=1.0, help="seconds to wait after each search in --loop")
    parser.add_argument("--udp", default=DECODER_UDP, help="host:port for decoder reports (default 127.0.0.1:5005)")
    args = parser.parse_args(argv)

    if args.loop and not args.self_test and not args.file:
        return run_loop(Path(args.target), args.interval, args.udp)

    target_pci = load_target_pci(Path(args.target))
    centre_hz = None
    if args.self_test:
        from lte_sync import synthesize_pair

        samples = synthesize_pair(320, 4000, 20000)
        report = search_cell(samples)
    elif args.file:
        samples = np.fromfile(args.file, dtype=np.complex64)
        print("Searching the file", flush=True)
        report = search_samples(samples, args.rate, target_pci)
    else:
        print("Reading Kraken channel 0", flush=True)
        samples, rate, centre_hz = capture_channel0()
        print("Searching %.3f MHz" % (centre_hz / 1e6), flush=True)
        report = search_samples(samples, rate, target_pci)

    hits = report["hits"]
    if not hits:
        print("No LTE cell ID found")
        if centre_hz is not None:
            print("Centre %.3f MHz" % (centre_hz / 1e6))
            note = _band_note(centre_hz)
            if note:
                print(note)
        print(
            "Best PSS %.1fx (needs %.1fx) at %+.1f kHz"
            % (report["best_ratio"], report["min_ratio"], report["best_offset_hz"] / 1000.0)
        )
        if report["best_ratio"] >= report["min_ratio"]:
            print(
                "PSS sequence %.2f, SSS %.2f (needs %.2f)"
                % (report["best_pss_seq"], report["best_sss"], report["sss_min"])
            )
        return 1
    if args.self_test:
        item = hits[0]
        print("PCI %d" % item["pci"])
        print("offset %+.1f kHz, SSS %.2f" % (item["offset_hz"], item["sss_score"]))
        return 0
    _cell_id = None
    if not args.self_test:
        try:
            _pci, _cell_id = load_target(Path(args.target))
        except (OSError, ValueError):
            _cell_id = None
    chosen, secondary, prioritized = select_cells(hits, target_pci)
    _print_packet(decoder_packet(
        centre_hz,
        chosen,
        target_pci,
        _cell_id,
        secondary_pcis=secondary,
        target_prioritized=prioritized,
    ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
