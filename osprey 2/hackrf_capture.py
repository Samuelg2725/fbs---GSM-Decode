"""Record one LTE capture from a HackRF One.

The radio is tuned 600 kHz above the cell. That moves the HackRF DC spike
off the centre 72 subcarriers. 1.92 MHz was a bad offset: cutting the
recording down to 1.92 MS/s folds that spike straight back onto the sync
signals. The spike is notched out and the samples are returned at baseband.

hackrf_transfer is the tool shipped with the HackRF host software. SoapySDR
is used when it is installed and hackrf_transfer is not on PATH.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

DEFAULT_RATE = 15_360_000
DEFAULT_SECONDS = 0.120
# 40 bins at 15 kHz. Outside the centre 72 subcarriers, and not a multiple
# of 1.92 MHz, so decimation does not fold it onto the PSS.
LOW_IF_HZ = 600_000


def load_iq(path, sample_rate: float = DEFAULT_RATE, low_if_hz: float = 0.0) -> np.ndarray:
    """Read a .cs8 HackRF file or a complex64 file and return baseband complex64."""
    file_path = Path(path)
    if file_path.suffix.lower() == ".cs8":
        raw = np.fromfile(file_path, dtype=np.int8)
        raw = raw[: raw.size // 2 * 2].astype(np.float32).reshape(-1, 2)
        samples = ((raw[:, 0] + 1j * raw[:, 1]) / 128.0).astype(np.complex64)
    else:
        samples = np.fromfile(file_path, dtype=np.complex64)
    return condition(samples, sample_rate, low_if_hz)


def notch_tone(samples, freq_hz: float, sample_rate: float) -> np.ndarray:
    """Remove one complex tone. Used for the HackRF DC spike after the shift."""
    samples = np.asarray(samples, dtype=np.complex128)
    index = np.arange(len(samples))
    tone = np.exp(1j * 2 * np.pi * float(freq_hz) * index / float(sample_rate))
    coef = np.vdot(tone, samples) / len(samples)
    return np.asarray(samples - coef * tone, dtype=np.complex64)


def condition(samples, sample_rate: float, low_if_hz: float = 0.0) -> np.ndarray:
    """Drop the DC mean, shift a low-IF capture back, and notch the spike."""
    samples = np.asarray(samples, dtype=np.complex64)
    samples = samples - np.mean(samples)
    if low_if_hz:
        # The radio is tuned above the cell, so the cell sits at -low_if.
        # Shift it up to baseband, which moves the DC spike to +low_if.
        index = np.arange(len(samples))
        samples = samples * np.exp(1j * 2 * np.pi * float(low_if_hz) * index / float(sample_rate))
        samples = notch_tone(samples, float(low_if_hz), sample_rate)
    return np.asarray(samples, dtype=np.complex64)


def capture(
    centre_hz: float,
    sample_rate: float = DEFAULT_RATE,
    seconds: float = DEFAULT_SECONDS,
    lna_gain: int = 32,
    vga_gain: int = 28,
    amp: bool = False,
    low_if_hz: float = LOW_IF_HZ,
):
    """One baseband capture. centre_hz is the LTE downlink, not the radio tune.

    Returns the baseband samples and the Unix time taken immediately before
    the radio is started.
    """
    count = int(round(float(sample_rate) * float(seconds)))
    tune = int(round(float(centre_hz) + float(low_if_hz)))
    if shutil.which("hackrf_transfer"):
        raw, capture_timestamp = _transfer(tune, int(sample_rate), count, int(lna_gain), int(vga_gain), bool(amp))
    else:
        raw, capture_timestamp = _soapy(tune, float(sample_rate), count, int(lna_gain), int(vga_gain), bool(amp))
    return condition(raw, sample_rate, low_if_hz), capture_timestamp


def save_capture(path, samples) -> None:
    """Write complex64 baseband samples."""
    np.asarray(samples, dtype=np.complex64).tofile(path)


def _transfer(tune_hz: int, sample_rate: int, count: int, lna_gain: int, vga_gain: int, amp: bool):
    with tempfile.TemporaryDirectory() as folder:
        target = Path(folder) / "capture.cs8"
        command = [
            "hackrf_transfer",
            "-r", str(target),
            "-f", str(tune_hz),
            "-s", str(sample_rate),
            "-n", str(count),
            "-l", str(lna_gain),
            "-g", str(vga_gain),
            "-a", "1" if amp else "0",
        ]
        capture_timestamp = time.time()
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "").strip()
            raise RuntimeError("hackrf_transfer failed: %s" % detail)
        raw = np.fromfile(target, dtype=np.int8)
    raw = raw[: raw.size // 2 * 2]
    if raw.size < count:
        raise RuntimeError("HackRF returned %d samples, wanted %d" % (raw.size // 2, count))
    pairs = raw[: count * 2].astype(np.float32).reshape(-1, 2)
    return ((pairs[:, 0] + 1j * pairs[:, 1]) / 128.0).astype(np.complex64), capture_timestamp


def _soapy(tune_hz: int, sample_rate: float, count: int, lna_gain: int, vga_gain: int, amp: bool):
    import SoapySDR
    from SoapySDR import SOAPY_SDR_CF32, SOAPY_SDR_RX

    device = SoapySDR.Device({"driver": "hackrf"})
    device.setSampleRate(SOAPY_SDR_RX, 0, float(sample_rate))
    device.setFrequency(SOAPY_SDR_RX, 0, float(tune_hz))
    device.setGain(SOAPY_SDR_RX, 0, "LNA", float(lna_gain))
    device.setGain(SOAPY_SDR_RX, 0, "VGA", float(vga_gain))
    device.setGain(SOAPY_SDR_RX, 0, "AMP", 14.0 if amp else 0.0)
    stream = device.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CF32)
    capture_timestamp = time.time()
    device.activateStream(stream)
    output = np.zeros(count, dtype=np.complex64)
    filled = 0
    try:
        while filled < count:
            chunk = np.zeros(min(1 << 16, count - filled), dtype=np.complex64)
            read = device.readStream(stream, [chunk], len(chunk))
            if read.ret <= 0:
                raise RuntimeError("HackRF stream stopped (%s)" % read.ret)
            output[filled:filled + read.ret] = chunk[: read.ret]
            filled += read.ret
    finally:
        device.deactivateStream(stream)
        device.closeStream(stream)
    return output, capture_timestamp
