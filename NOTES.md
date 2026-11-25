# Technical notes and test findings

Test location: home, November 2026. Reference: a phone network-monitor app at the same place.

## Hardware facts

- **BladeRF 2.0 micro (AD9361):** RX1 and RX2 share one tuner (same frequency, sample rate and analog filter)
  but have independent gains. Verified: RX1 at 10 dB and RX2 at 30 dB recorded together.
- With two channels, `bladeRF-cli rx config n=` counts samples of both channels; gains must be set after
  `rx config`, and `set agc off` (not `set agc rx off`) turns automatic gain off on both.
- Two-channel streaming into GNU Radio works through **gr-soapy** (`driver=bladerf`); gr-osmosdr with
  `numchan=2` does not work with the BladeRF.

## GSM

- **Gain:** 40 dB overloaded the radio at home (cells caught 21–32 % of the time, one never decoded);
  10 dB, 20 dB and AGC all caught 100 %. Default is 20 dB.
- **Cell ID rate:** System Information 3 (Cell ID) is broadcast every 0.94 s; some O2 cells send it three
  times per 1.9 s. Tracking now catches 100 % of broadcasts, so ~1 s is the limit.
- **Fresh decoder per channel:** a decoder kept across many empty channels drifts in frequency and then
  misses strong cells. Building a new decoder per channel fixed it; keeping the previous flowgraph alive until
  the new one runs lets gr-osmosdr re-use the open BladeRF (~0.3 s per channel instead of ~2.8 s).
- **Full 900 MHz scan:** 7/7 phone cells (incl. −90 and −97 dBm cells) in 257 s with one antenna.
  The two-antenna scan also finds 7/7 but takes ~680 s (the radio is reopened per channel).
- **Neighbour lists:** operators list their whole block (Vodafone 1–22, O2 113–123 here), so they do not
  speed up finding cells; a band scan is the reliable way.
- **Wideband (15 MHz recordings):** one gain must suit every cell in the slice and LTE carriers in the GSM
  band confuse a power-based channel filter; the narrow scan was kept.
- **1800 MHz:** the EE cell on ARFCN 649 (−92 dBm on the phone) was not received — likely the antennas.

## LTE

- **Decoder:** srsRAN `pdsch_ue` in file mode reads the MIB and SIB1. The capture is aligned to subframe 0
  first and srsRAN is told the PCI (`-c`), so any PCI that the cell search finds can be decoded.
  (srsRAN's live mode only locks to the strongest cell, and did not lock at all with the BladeRF.)
- **Low IF:** the radio is tuned 0.6 MHz above the carrier and the DC spike is notched out.
- **Capture rate:** 15.36 Msps files decode correctly for 10 MHz cells (resampling to 11.52 Msps lost SIB1).
  15/20 MHz cells are recorded again at 23.04/30.72 Msps when the MIB reports their width.
- **Weak-cell search:** the original search chooses one frequency offset per PSS group; a strong cell's false
  peak can win. A second pass at the strongest cell's offset finds weaker cells on the same carrier.
- **Sector siblings:** the sectors of one mast share its frame timing and are numbered 3·G+0/1/2. Decoding
  the other two PCIs with a decoded cell's timing found O2 PCI 0 on 816 and 952.5 MHz (sync too weak to
  detect on its own).
- **Two antennas:** per-cell maximum-ratio combining from each cell's PSS (simulated: +3 dB with equal
  antennas, never worse than the better antenna).
- **Speed:** the cell search was rewritten for speed (`lte_fastsync.py`, identical results, 5× faster on
  the Pi); a 10 MHz carrier takes ~7 s.
- **Result at home:** 7 of the 13 cells the phone lists decoded with full identity. Not received at all:
  EE 1846.7/1866.5 MHz and Three 2117.6 MHz (no sync even at 54 dB gain — likely the antennas);
  793.5 MHz cells (−93 dBm on the phone) are at the detection limit; O2 2144.7 MHz syncs but its MIB does
  not decode (cause not found).

## Tried and removed

Wideband GSM with rescue pass, keep-open GSM retune, srsRAN live mode, the all-Python LTE MIB decoder,
an LTE "deep search" (false detections on air) and resampling srsRAN's input to 11.52 Msps (lost SIB1).
The code for these is in the git history.
