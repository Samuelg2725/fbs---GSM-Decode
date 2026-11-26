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

- **Tracking one carrier (`lte_track.py`):** identities are decoded once per PCI; later cycles only run the
  cell search (~1.5 s) to see which PCIs are present, so a cycle is ~2–4 s. Raw IQ is kept as SigMF with each
  PCI's subframe-0 position, for re-decoding or direction finding (RX1/RX2 share one tuner, so `--dual`
  recordings keep the phase between antennas).

- **SIB1 combining (srsRAN patch, `sib1_print.cc`):** stock `pdsch_ue` clears the soft buffer before every SIB1
  transmission and gives up on a subframe whose PDCCH does not decode. The patch adds up PDCCH soft bits over the
  occasions with the same redundancy version (one sum per CFI; scaled by 1/sqrt(n) because srsRAN skips candidates
  whose mean soft value is below 0.3), keeps the PDSCH soft buffer over up to 8 transmissions, and reuses the last
  SIB1 allocation when a PDCCH is missed. Soft bits are divided by 4 first: srsRAN's 16-bit turbo decoder lost
  accuracy with the summed values (8 transmissions: 33 decodes with /4, 17 without; 16 transmissions was worse).
  Test with srsRAN's own transmitter (`sib1_enodeb`, SIB1 every 20 ms, RV per 36.321) and 1 s recordings: stock
  srsRAN 0/3 recordings at -18 and -19 dB (noise vs average power), patched 3/3; patched 2/3 at -20 dB.
  `SIB1_NO_COMBINE=1` switches it off, `SIB1_DEBUG=1` prints each SIB1 occasion.
- **Interference cancellation (`lte_cancel.py`):** CRS of a decoded cell rebuilt from its PCI and the channel
  measured on its own pilots (frequency smoothing after removing the timing phase slope, 0.25/0.5/0.25 across
  slots), subtracted with the cyclic prefix; PSS and SSS likewise (without that a 15 dB stronger cell's sync pulled
  the weak cell's timing). Test with two `sib1_enodeb` cells, different timing: weak cell 15/18/21 dB below the
  strong one - nothing without, MIB + SIB1 with; srsRAN's SNR of the weak cell +3.5 -> +17.8 dB at -15 dB. Synthetic
  load test: +12 dB signal-to-interference at 0% load, +3 dB at 10%, +1.5 dB at 30%.
- **PBCH and SIB1 cancellation:** the PBCH of a cell is rebuilt from its decoded MIB for every SFN (found by matching
  all 1024 in the first frames), SIB1 subframes (PCFICH + SIB1 PDCCH + PDSCH) by `sib1_regen` with srsRAN's encoder
  from the decoder's `SIB1dci` line and the SIB1 bytes (correlation with the transmitted subframe 1.000, all four
  RVs). The library's `conv_rate_match` (osprey 2/lte_coding.py) does not follow 36.212 5.1.4.2 (60% bit agreement
  with srsRAN's PBCH); lte_cancel.py has its own, which matches 100%. Two synchronised `sib1_enodeb` cells with real
  SIB1s on the same resource blocks: weak cell decoded at -10/-15/-20 dB with cancellation (nothing without), MIB
  only at -25. Office: PCI 76 (eNB 2264, not in the phone's list) decoded after removing PCI 199.
- Office measurements showed frame-synchronised cells: 1363 PCI 21 with 199 (7.054/7.055 ms), 3725 PCI 80 with 144.
- When the SSS cannot tell subframe 0 from 5 (weak cell), both frame positions are tried; the MIB CRC decides.

## Tried and removed

Wideband GSM with rescue pass, keep-open GSM retune, srsRAN live mode, the all-Python LTE MIB decoder,
an LTE "deep search" (false detections on air) and resampling srsRAN's input to 11.52 Msps (lost SIB1).
The code for these is in the git history.
