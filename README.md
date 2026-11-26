# GSM and LTE cell identity decoder (BladeRF 2.0 + Raspberry Pi 5)

Finds the mobile-network cells around you and decodes their identity:

- **GSM (2G):** channel (ARFCN), operator, LAC and Cell ID — scan a whole band, or track one cell continuously.
- **LTE (4G):** carrier, PCI, operator, TAC, eNB, sector and Cell ID — every cell the radio can hear on each
  carrier, not only the strongest.

Results are printed as colour-coded tables, and saved to CSV files in `~/field/`. Every decoded message is also
saved to a **pcap file** in `~/field/pcap/` that opens in Wireshark (see "pcap files" below).

---

## Quick start

Run everything from the repository folder, with the BladeRF library on the path:

```bash
cd ~/fbs-upload
export LD_LIBRARY_PATH=/usr/local/lib
```

| I want to… | Command |
|---|---|
| Find every GSM cell in the 900 MHz band (~4 min) | `python3 gsm_arfcn_fast.py --band 900` |
| Find every GSM cell in the 1800 MHz band | `python3 gsm_arfcn_fast.py --band 1800` |
| Check some GSM channels | `python3 gsm_arfcn_fast.py --arfcn 1 6 13 22 114 116 119` |
| Track one GSM cell continuously | `python3 gsm_arfcn_fast.py --arfcn 13` |
| Track one GSM cell with both antennas | `python3 gsm_dual.py --arfcn 13` |
| Read a GSM cell's neighbour list | `python3 gsm_neigh.py --arfcn 13` |
| Decode LTE carriers by frequency | `python3 lte_decode.py --freq-mhz 806 816 952.5` |
| Decode LTE carriers by channel number | `python3 lte_decode.py --earfcn 6300 6400 3725 223` |
| Decode LTE with both antennas | `python3 lte_decode.py --freq-mhz 816 --dual` |
| Check how much signal the antennas bring in per band (no decoding) | `python3 spectrum_check.py` |
| Compare antennas one after another (same gains) | `python3 spectrum_check.py --antennas whip panel` |
| Look for the cells a phone lists (EARFCN:PCI,PCI) | `python3 lte_decode.py --cells 6175:216,272 1617:389,461` |
| Decode only one cell | `python3 lte_decode.py --cells 6225:320 --only` |
| …and keep repeating until Ctrl-C | `python3 lte_decode.py --cells 6225:320 --only --loop` |
| Watch every LTE cell on one carrier, saving raw IQ | `python3 lte_track.py --earfcn 6400` |
| …the same with both antennas (for direction finding) | `python3 lte_track.py --earfcn 6400 --dual` |

Press **Ctrl-C** at any time: the tool stops cleanly and prints its summary table.
Add `--help` to any command to see all options.

Only one tool can use the BladeRF at a time.

---

## GSM

### Scan (`gsm_arfcn_fast.py --band …` or several `--arfcn`)
Visits each channel in turn (1-124 first, then 975-1023); a live line shows which ARFCN it is checking. Channels with no GSM signal are left after 1.2 s; a channel with a cell is held
until its Cell ID is decoded (up to 2.5 s). Every channel gets a fresh decoder, which keeps the scan reliable.

Useful options: `--passes 0` (repeat until Ctrl-C), `--verbose` (also list empty channels),
`--csv ~/field/scan.csv` (save the results).

### Track one cell (`gsm_arfcn_fast.py --arfcn N`)
Prints every Cell-ID broadcast (a cell sends it about once a second) and, on exit, a summary: broadcasts
decoded and missed, percentage caught, longest gap and signal level. The summary is also added to
`~/field/track_summary.csv`. Use `--label NAME` to name the run.

### Two antennas (`gsm_dual.py`)
RX1 and RX2 listen to the same channel at different gains (default 20 dB and 40 dB). A strong cell that
overloads the high-gain antenna is still decoded on the low-gain one, and a weak cell is picked up by the
high-gain one. Each decode is tagged **[RX1+RX2]**, **[RX1 only]** or **[RX2 only]**.
`--levels` shows both antennas' signal levels (and warns about clipping).

### Neighbour lists (`gsm_neigh.py`)
Shows the channels each cell tells phones to measure. `--follow` visits every listed channel in turn.
Operators usually list their whole channel block, so most listed channels have no cell.

### Gain
The default gain is **20 dB**. Higher gain (e.g. 40) overloads the radio near strong cells and loses decodes.
If a strong neighbouring channel causes trouble, `--bw 300e3` narrows the radio's filter.

---

## LTE

`lte_decode.py` records each carrier, searches for every cell (PCI) on it, and decodes each one with srsRAN.
It automatically:

- adjusts the gain for each band (800 MHz and 2100 MHz need very different gains);
- searches again for weaker cells on the same carrier;
- tries the other sectors of every mast it decodes (they share the mast's timing, so they can be decoded
  even when their own sync signal is too weak to detect);
- records again at a higher sample rate for 15 MHz and 20 MHz wide cells.

Two techniques get weak cells that srsRAN alone cannot:

- **SIB1 combining** (in the srsRAN patch): a cell repeats the same SIB1 every 20 ms in four redundancy versions.
  Stock srsRAN decodes each transmission on its own; the patched decoder adds up the control channel and the
  SIB1 data over the repeats. Measured with srsRAN's own test signals: stock srsRAN stops decoding about 3-4 dB
  earlier (1 s recordings: stock 0 of 3 at a noise level where the patched decoder gets 3 of 3).
- **Interference cancellation** (`lte_cancel.py`): every cell whose signals can be rebuilt is removed from the
  recording - reference signals and sync signals (from the PCI), the MIB channel (from its decoded MIB) and, for
  decoded cells, the whole SIB1 subframe (rebuilt by srsRAN's own transmitter code, `sib1_regen`). The cleaned
  recording is searched again and the remaining cells retried, in rounds. This is what decodes cells hidden under
  a stronger cell - including frame-synchronised cells, whose MIB and SIB1 are sent on exactly the same resource
  elements. Measured with srsRAN test cells carrying real SIB1s: a synchronised cell 10-20 dB weaker than the
  decoded one is fully decoded after cancellation and not at all before; cells with different timing likewise.
  User traffic of the strong cell cannot be rebuilt, so on a busy carrier the gain is smaller. `--no-cancel` turns
  it off. The line `frame timing (ms into the 10 ms frame)` shows which cells are synchronised.

With `--cells EARFCN:PCI,PCI …` (for example copied from a phone's network-monitor app) every listed PCI is
looked for directly, even when the search does not find it: knowing the PCI, its sync signal is added up over
the whole recording, which finds cells several dB weaker than the blind search. A "Cell list check" table at
the end shows which of the listed cells were decoded.

Results go to `~/field/lte_cells.csv`. Useful options: `--dual` (both antennas), `--secs 2` (longer recording
for weak sectors), `--no-auto-gain --gain N`.

### Tracking one carrier (`lte_track.py`)
Stays on one EARFCN and repeats: record → find every PCI → decode any PCI not identified yet. In a city one
EARFCN carries many cells (different PCIs); this shows them all as they come and go:

```
14:02:11.120  NEW PCI 117 (strength 9, weak-cell search)
           PCI 117 decoded  Cell ID 4968730 (eNB 19409, sector 26)  TAC 57379  Vodafone 234-15   [weak-cell search]
14:02:11.400  #12   2.6s  gain 20 -21 dBFS  PCI 1 O2 138630264 [22]  PCI 117 Vodafone 4968730 [9]  PCI 0 O2 138630254 [sector]   IQ saved
14:02:20.700  PCI 117 (Vodafone 4968730) not heard for 3 cycles
```

- **Raw IQ** is saved in `~/field/iq/` (SigMF: `.sigmf-data` 16-bit I/Q + `.sigmf-meta` JSON) whenever a new
  PCI appears (`--save-iq all` saves every cycle, `none` saves nothing). The `.sigmf-meta` file lists every PCI
  in the recording with its subframe-0 sample position, frame length, frequency offset and identity, so each
  cell can be picked out of the shared channel later. With `--dual` each sample holds RX1 and RX2 together
  (same tuner, so the phase between antennas is kept). One second is 61 MB (123 MB with two antennas);
  the oldest recordings are deleted above `--iq-max-gb` (default 20).
- Cycles that only check which PCIs are present use a short 0.2 s recording (`--watch-secs`); a full 1 s
  recording is taken when something must be decoded or saved. Each line shows the record and search times.
  The radio is kept open between recordings (libbladeRF), so a watch cycle takes about 2 s; `--radio cli`
  starts bladeRF-cli for each recording instead (slower, used automatically if the radio cannot be kept open).
- 5 MHz cells (e.g. Three on 800 MHz) are decoded from a recording reduced to 7.68 Msps, the rate srsRAN
  expects for that width.
- On Ctrl-C: a table of every PCI (identity, how often heard, strength, first/last seen), plus
  `~/field/lte_track_summary.csv` and a per-cycle log `~/field/lte_track_log.csv`.

UK carriers seen during testing (EARFCN → MHz): 6175 → 793.5, 6225 → 798.5, 6300 → 806, 6400 → 816,
3501 → 930.1, 3725 → 952.5, 1363 → 1821.3, 1617 → 1846.7, 1815 → 1866.5, 76 → 2117.6, 223 → 2132.3,
347 → 2144.7. A phone network-monitor app shows which ones are used where you are.

---

## pcap files

Every tool writes the messages it decodes to `~/field/pcap/<tool>_<date>_<time>.pcap` (`--pcap FILE` to choose
the name, `--no-pcap` to turn it off). Messages are stored as GSMTAP (UDP port 4729), which Wireshark decodes:

- **GSM:** every message gr-gsm decodes - System Information 1-4, 2bis/2ter/2quater, paging and other BCCH/CCCH
  messages. With two antennas the GSMTAP antenna field says which receiver decoded it (0 = RX1, 1 = RX2).
- **LTE:** the MIB (LTE RRC BCCH-BCH) and SIB1 (BCCH-DL-SCH) of every decoded cell, with the EARFCN in the
  GSMTAP ARFCN field and the PCI in the frame-number field.

```bash
wireshark ~/field/pcap/lte_20261126_103300.pcap
tshark -r ~/field/pcap/gsm_scan_20261126_103300.pcap -V | less
```

The LTE bytes come from the srsRAN patch: after updating the repository, re-run the patch and rebuild srsRAN
once (setup step 3).

---

## Example output

```
== GSM cells found: 7   (1 pass(es), 257 s) ==========================
+-------+-------+----------+---------+-------+---------+-------+------+
| ARFCN |   MHz | Operator | MCC-MNC |   LAC | Cell ID | Power | Seen |
+-------+-------+----------+---------+-------+---------+-------+------+
|     1 | 935.2 | Vodafone | 234-15  |    98 |   16631 |   -66 |  1/1 |
|    13 | 937.6 | Vodafone | 234-15  |    98 |     942 |   -29 |  1/1 |
|   119 | 958.8 | O2       | 234-10  | 22860 |   22679 |   -33 |  1/1 |
...
```

---

## Setup (once)

1. **BladeRF library** (libbladeRF 2.x and `bladeRF-cli`) installed in `/usr/local/lib`.
2. **GNU Radio 3.10, gr-gsm, gr-osmosdr and gr-soapy** with the SoapySDR BladeRF module
   (`SoapySDRUtil --find` should list the BladeRF).
3. **srsRAN 4G** (for LTE), built with the small patch in `osprey 2/srsran_watch` that prints the SIB1 cell
   identity:
   ```bash
   sudo apt install -y build-essential cmake libfftw3-dev libmbedtls-dev libboost-program-options-dev libconfig++-dev libsctp-dev
   cd ~ && git clone --depth 1 https://github.com/srsran/srsRAN_4G.git
   python3 ~/fbs-upload/"osprey 2"/srsran_watch/print_sib1.py
   mkdir -p ~/srsRAN_4G/build && cd ~/srsRAN_4G/build && cmake .. && make -j3 pdsch_ue sib1_regen
   ```
   The patch makes srsRAN print the decoded cell identity and the raw MIB / SIB bytes (for the pcap files), adds
   SIB1 combining and the PHICH fix, and builds the `sib1_regen` helper used by interference cancellation.
   Re-run these whenever anything in `osprey 2/srsran_watch` changes:
   ```bash
   python3 ~/fbs-upload/"osprey 2"/srsran_watch/print_sib1.py && cd ~/srsRAN_4G/build && cmake .. && make -j3 pdsch_ue sib1_regen
   ```
4. Optional check that the fast LTE search gives the same answers as the original code:
   `python3 lte_fastsync.py` (ends with `OK`).

---

## Files

| File | Purpose |
|---|---|
| `gsm_arfcn_fast.py` | GSM band / channel scan and single-cell tracking |
| `gsm_dual.py` | GSM with both antennas (tracking and scan) |
| `gsm_neigh.py` | GSM neighbour lists |
| `gsm_freqlist.py` | Decoder for GSM frequency lists (used by `gsm_neigh.py`) |
| `lte_decode.py` | LTE cell identity decoder |
| `spectrum_check.py` | Antenna check: signal power per frequency, with and without the antennas |
| `lte_track.py` | Tracks every LTE cell on one carrier and saves raw IQ |
| `sigmf_join.py` | Joins the IQ recordings of a tracking run into one file (per sample rate) |
| `iq_export.py` | Exports IQ recordings as decoder-ready files (centred, float32, one file per cell aligned to its subframe 0) with decode commands |
| `lte_fastsync.py` | Faster LTE cell search (identical results) |
| `pcap_out.py` | Writes decoded messages to pcap files (GSMTAP) for Wireshark |
| `common.py` | Shared tables, colours and operator names |
| `osprey 2/` | LTE signal-processing library and the srsRAN patch |
| `NOTES.md` | Technical notes and test findings |

Set `NO_COLOR=1` (or pass `--no-color`) for plain text output.

## Reference scan with the Teltonika CALYX modem (ALA440)

The CALYX's modem decodes cells like a phone does. `modem_scan.py` makes it scan (`AT^NETSCAN`, no SIM needed)
and lists every cell it decoded: EARFCN, PCI, operator, TAC and Cell ID. Each cell is marked with whether the
SDR (`lte_decode.py` / `lte_track.py`, via `~/field/lte_cells.csv`) has decoded it with the same Cell ID.
At the end it prints the `lte_decode.py --cells ...` command for the cells the SDR has not decoded yet.

    sudo modprobe option && echo 1d12 0101 | sudo tee /sys/bus/usb-serial/drivers/option1/new_id   # once per boot
    sudo python3 modem_scan.py                 # LTE, all bands
    sudo python3 modem_scan.py --band 20 3     # only bands 20 and 3 (faster)
    sudo python3 modem_scan.py --rat all       # LTE + GSM
    sudo python3 modem_scan.py --sdr           # then the SDR decodes the same cells, and the two are compared

`modem_probe.py --cmd 'AT...'` sends any single AT command. Other useful ALA440 commands:
`AT^CELLINFO=1` (serving cell), `AT^CELLINFO=2` (neighbour cells with PCI and RSRP),
`AT^CELLLOCK=1,3,1,<earfcn>,<pci>` (lock onto one cell; set it in flight mode `AT+CFUN=0`, then `AT+CFUN=1`;
`AT^CELLLOCK=0` unlocks, which needs a modem restart to take effect).

### Modem and SDR at the same time: `compare_scan.py`

    sudo python3 compare_scan.py                       # add --dual for both antennas
    python3 lte_decode.py --uk                         # SDR alone: every UK LTE carrier (20 EARFCNs), blind

The modem keeps scanning all LTE bands over and over while the SDR goes once through every carrier
(the UK list plus every EARFCN the modem reports). At the end, one table lists every cell either of them found,
with both Cell IDs: both / modem only / SDR only / different ID. Only this run's results count.
It is saved as ~/field/compare_<time>.csv.
