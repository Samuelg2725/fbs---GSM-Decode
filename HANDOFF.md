# fbs — GSM Cell Scanner + Direction Finding (BladeRF / LimeSDR + KrakenSDR)

## What this is
Scans the GSM (2G) band, decodes each cell's identity (CID, LAC, MCC, MNC = operator),
and with a KrakenSDR gives a compass bearing to each tower. Runs on a Raspberry Pi 5.
Two radios work as the decoder: LimeSDR (long working) and BladeRF 2.0 micro (focus of
recent work — now working).

## Hardware
- Raspberry Pi 5 (user `revector`, `revector.local`, 10.0.0.9)
- BladeRF 2.0 micro (xA4), USB 3, antenna on RX1. Use a powered USB 3 hub.
- LimeSDR + KrakenSDR = Owen's direction-finding rig (decoder Pi 192.168.2.10, Kraken 192.168.2.20)

## *** THE CRITICAL FIX — read first ***
For days the BladeRF captured clean-looking IQ but decoded ZERO GSM. Root cause:
the device FIRMWARE (v2.4.0) was TOO OLD for its FPGA (v0.16.0) — they used a
"legacy message size" that silently SCRAMBLED the sample stream (spectrum looked
perfect, bits were corrupt). The stock apt libbladeRF (2.5.0) didn't detect it.
FIX:
  1. Built libbladeRF + bladeRF-cli 2.6.1 from source (github.com/Nuand/bladeRF, host/).
  2. Flashed firmware v2.6.0 (https://www.nuand.com/fx3/bladeRF_fw_latest.img).
After that the BladeRF decodes GSM identically to the Lime.
New lib lives in /usr/local/lib — prefix radio commands with
`LD_LIBRARY_PATH=/usr/local/lib` (the apt lib 2.5.0 otherwise shadows it). (TODO: make permanent.)

Other hard-won lessons:
- `numchan`/`nchan` in osmosdr args BREAKS the BladeRF. The stock grgsm_scanner/livemon/
  decode binaries hardcode numchan=1 AND crash on GNU Radio 3.10 (pdu.socket_pdu removed),
  so they are unusable here — we decode with our own gr-gsm flowgraphs instead.
- BladeRF has no auto DC correction (tune off-centre + shift if ever needed; not needed post-firmware).
- Gain window on GSM900 here ~35-47 dB (20=too weak, 45≈clipping, 40=clean). Captured peak
  should be a few hundred to ~1500 of 2047.
- Streaming needs big USB buffers: samples=16384 buffers=64 xfers=16 (131072 fails NO_MEM).
  Overflow ("0000" flood) = dropped samples = broken GSM timing. Capture natively to /dev/shm.

## WORKING pipeline (BladeRF, capture -> decode)
Single channel:
  LD_LIBRARY_PATH=/usr/local/lib bladeRF-cli -e "set agc rx off" -e "set frequency rx 937600000" \
    -e "set samplerate rx 2000000" -e "set bandwidth rx 1500000" -e "set gain rx 40" \
    -e "rx config file=/dev/shm/cap.sc16 format=bin n=16000000 samples=16384 buffers=64 xfers=16 timeout=10000" \
    -e "rx start" -e "rx wait 30000"
  python3 -c "import numpy as np; r=np.fromfile('/dev/shm/cap.sc16',np.int16).astype(np.float32)/2048.0; (r[0::2]+1j*r[1::2]).astype(np.complex64).tofile('/dev/shm/cap.cfile')"
  python3 ~/filedec.py /dev/shm/cap.cfile      # -> CID 51478 LAC 98 MCC 234 MNC 15 (Vodafone)

Wideband (ONE capture -> many channels):
  LD_LIBRARY_PATH=/usr/local/lib bladeRF-cli -e "set agc rx off" -e "set frequency rx 942500000" \
    -e "set samplerate rx 15360000" -e "set bandwidth rx 14000000" -e "set gain rx 40" \
    -e "rx config file=/dev/shm/wide.sc16 format=bin n=61440000 samples=16384 buffers=64 xfers=16 timeout=20000" \
    -e "rx start" -e "rx wait 40000"
  python3 ~/widescan.py /dev/shm/wide.sc16 942.5    # FFTs, finds live channels, decodes each
  # Proven: 5 cells from one 4s capture (ARFCN 4,5,9,13,22).

Lime (live, Owen's rig): python3 gsm_arfcn_fast.py --arfcn 13 --df   (osmosdr driver=lime + Kraken bearing)

## Files
- experimental_lte/lte_scan.py  (PARKED 15:xx - user reverted to basic lte_decode.py) LTE MULTI-BAND SCAN (2026-11-25): per band survey (100 kHz power bins) -> candidate
                     centres (+ known B20 793.5/798.5/806/816) -> 0.2 s capture, low-IF 0.6, Python PSS/SSS
                     probe with +-0.3 MHz digital nudges -> srsRAN pdsch_ue -p 50 (MIB) then re-read at the
                     cell's own rate (decimate for 5 MHz, re-capture 23.04/30.72 for 15/20 MHz) -> table +
                     ~/field/lte_scan.csv. --bands, --freqs/--earfcns, --discover, --max-pci, --loop, --label.
                     FIRST PI RUN 12:4x: sync found the phone's PCIs (806 PCI 21, 816 PCI 1, 952.5 PCI 1) but
                     every row said 'MIB ok, no SIB1': our pdsch_ue patch PRINTS the SIB1 line (sib1_print.cc)
                     instead of saving sib1_bits.bin -> now parsed from the text. Survey centres were 1-2 MHz
                     off -> default is now the phone's carrier list (B20 6175/6225/6300/6400, B8 3501/3725,
                     B3 1363/1617, B1 223/347); --discover adds the survey. All PCIs per carrier decoded.
- gsm_dual.py        TWO-ANTENNA GSM TRACKING (2026-11-25): RX1 + RX2 on the same ARFCN, separate gains
                     (--gain1 20 / --gain2 45, or --agc1/--agc2), a full gr-gsm decoder per channel, merged
                     output tagged [RX1+RX2]/[RX1 only]/[RX2 only]; summary = caught % for RX1, RX2, COMBINED
                     + dBFS level/peak per channel -> ~/field/dual_summary.csv. --src auto tries gr-soapy,
                     osmosdr via soapy, osmosdr numchan=2, bladerf nchan=2 until both channels stream.
                     --levels = only print both levels (gain proof). Merge logic tested here with stubs;
                     radio path untested (needs Pi). Pi has SoapySDR + bladerf driver (SoapySDRUtil --find).
- gsm_neigh.py       NEIGHBOUR-GUIDED GSM SEARCH (2026-11-25): per channel (fresh decoder, g20, sync-skip)
                     decodes SI3 (cell ID) + SI2/SI2bis/SI2ter neighbour lists (gsm_freqlist.py: bitmap0,
                     range 1024/512/256/128, variable bitmap; range decoder checked against libosmocore's
                     explicit formulas). --follow visits neighbours of neighbours; --seed-sweep finds a first
                     cell itself. Flags SI2quater (3G/LTE neighbours) as present - decoding it = next step.
                     PI 14:30: --arfcn 13 119 --follow = 7/7 in 147 s, 33 ch (Vodafone SI2 lists 1-22, O2
                     113-123 = whole operator blocks, bitmap0; every cell has SI2quater). Seed sweep that
                     stopped at the first cell found only Vodafone (lists never cross operators) -> now:
                     every-4th-channel seed sweep -> 14:42 missed ALL Vodafone (1/6/13/22 off the grid).
                     REMOVED (seed sweep + stride + near-check). gsm_neigh now = start cells given with
                     --arfcn, --follow, retry weak listed channels once (8 s). No start cell -> full sweep.
                     Re-uses the open BladeRF (old flowgraph kept until the new one runs).
- experimental_lte/lte_track.py (PARKED) LTE TRACKING. Default = CAPTURE MODE: lock once (probe PCI/CFO, MIB -> RB), then loop
                     0.1 s capture -> align (first 30 ms only) -> pdsch_ue stopped at first SIB1; re-lock
                     after 3 misses. --live = srsRAN live receiver: FAILED at home 12:5x (endless cell search,
                     'buffer size 1024' warnings; likely DC spike on centre, live srsRAN has no low-IF).
                     (--live) wraps live pdsch_ue (-f, -g) on one carrier, parses the
                     SIB1 line from print_sib1.py: first-decode time, SIB1 per second, identity changes,
                     summary -> ~/field/lte_track_summary.csv. LTE equivalent of gsm_arfcn_fast tracking.
  lte_scan ALL-BANDS RUN 12:5x: 4 cells OK (806 V PCI21 s20, 930.1 V s22, 952.5 O2 s122, 2132.3 V 15MHz s24),
  all match phone; ~7.5 s per carrier, 92 s for 10. Fails: 816 'no LTE' (OK alone) -> probe now retries on
  the whole capture; 1821.3 Three + 2144.7 O2 'MIB not decoded' -> --debug-dir saves pdsch_ue output;
  793.5/798.5/1846.7 EE no sync (weak, -84..-92). Speed: pdsch_ue now stopped as soon as SIB1 prints
  (--no-early-stop = old), --timing shows capture/probe/align/srsran split.
- lte_decode.py      *** THE LTE TOOL (basic, unchanged since 10:58) *** LTE/4G cell identity on the BladeRF. bladeRF-cli capture -> pure-Python decoder in
                     "osprey 2" (lte_pci_search: PCI; lte_cellid + lte_sib1: MIB + SIB1 -> CellID, eNB,
                     sector, PLMN, TAC). No srsRAN. --freq-mhz/--earfcn, --rate auto (15.36 Msps,
                     re-captures at 23.04/30.72 for 15/20 MHz cells), --loop, --survey LO-HI (power map
                     to find carriers), --file/--save, --self-test. Needs python3-scipy.
                     Self-test + synthetic file decode pass (2026-11-25); NOT yet tried on air.
- scan_uk.sh         field scan: ./scan_uk.sh [900|1800|uk] [label]. 900 = 3 x 15.36 Msps
                     (930.9/942.5/954.1), 1800 = 7 x 15.36 Msps (1810.5..1874.5), 4 s each, --fast.
                     Logs + merged unique-cell summary in ~/field/<date>_<label>_<band>/.
                     Env: GAIN, SECS, EXTRA (decode flags), KEEP=1 (keep captures).
- gsm_wide.py channels() now = E-GSM900 + DCS1800 inside the capture span (UK 2G bands).
- gsm_arfcn_fast.py  (repo root) Owen's live narrow decoder incl. --df Kraken, now with
                     --radio lime|bladerf (default lime = Owen's exact behaviour) and --gain.
                     bladerf: args bladerf=0 (no numchan), manual gain 40.
- gsm_wide.py        (repo root) Owen's wideband: capture once, decode every channel in it.
                     --radio bladerf records with bladeRF-cli to /dev/shm/gsm_wide.sc16 and
                     decodes sc16 directly. Also decodes an existing file: gsm_wide.py FILE
                     [--rate --center --format]. Lime path unchanged.
- gsm_arfcn_fast.py  live single/multi-ARFCN decoder (osmosdr). Owen's copy has --df (Kraken) + colour.
                     BladeRF copy: bladerf=0, manual gain, no numchan.
- gsm_arfcn_clean.py helper (cell_text/now_stamp/silence); scripts import it from /home/revector/.
- filedec.py         decode one complex-float32 .cfile offline. Reliable workhorse.
- widescan.py        wideband scanner: one wide capture -> find live channels -> decode all.
- bandscan.py        MAIN SCANNER. Any GSM band or custom MHz range, wide OR narrow mode, timings.
                     See "bandscan.py" section below.
- blade_capture.py / capture_ram.py   Soapy IQ capture (superseded by native bladeRF-cli capture).
- iqview.py          ASCII spectrum of a capture.  compare.py  diff two captures (offset/inversion).
- srsran_watch/      older separate LTE cell-ID work.

## bandscan.py (main scanner)
  python3 bandscan.py                         # E-GSM900 925-960, wide mode
  python3 bandscan.py --band all              # gsm450/480/850/R-GSM900/DCS1800/PCS1900
  python3 bandscan.py --band uk               # E-GSM900 + DCS1800
  python3 bandscan.py --range 1805-1830       # any MHz span (repeatable); 200 kHz raster
  python3 bandscan.py --mode narrow           # tune each live channel at 2 Msps instead
  python3 bandscan.py --compare --timing-csv ~/timings.csv   # wide vs narrow timing table
  python3 bandscan.py --band all --dry-run    # print the chunk plan only
Other flags: --gain-sweep / --gains, --secs, --rate, --jobs, --pipeline, --detect all, --decode-secs,
--snr, --log, -v. Defaults to /usr/local/lib for bladeRF-cli (--lib '' to disable).
Wide mode: band split into equal chunks (usable = rate*0.8), each captured once, FFT finds live
channels, each decoded from the raw sc16 file in parallel worker processes, stops at first SI3.
Narrow mode: 0.5 s wide FFT survey per chunk, then ONE live osmosdr stream (bladerf=0, same as
gsm_arfcn_fast.py) retunes to each live channel and prints the cell when its SI3 arrives.
NOTE: SI3 is broadcast every ~0.94 s, so ~0.9 s between lines in gsm_arfcn_fast.py is the
network's repeat rate, NOT slow decoding. A first decode in 0.5-1 s per channel is the floor.
Live BladeRF decode via osmosdr bladerf=0 CONFIRMED working (ARFCN 13, every SI3, 2026-11-24).
SSH DROPPED on the Pi every time 4 parallel decoders ran (12:13, 12:15) -> default --jobs 2,
workers niced. Suspect power/undervoltage: check `vcgencmd get_throttled` + dmesg. Run long
scans inside tmux so a dropped ssh doesn't kill them.
Why the first wide version was SLOWER than narrow (fixed):
  - channel filter low_pass(100k cutoff, 30k transition) at 15.36 Msps = ~1700 taps per channel.
    gsm_input already does the final channel filter, so the channelizer only has to stop
    aliasing: now 250k/500k = ~100 taps (--chan-cutoff/--chan-tw to change; 100/30 = old).
  - it decoded the WHOLE 8 s file per channel even after SI3 was found -> now early exit + 3 s cap.
  - channels decoded one after another -> now --jobs parallel (default = CPU cores).
  - numpy int16->complex64 conversion wrote a ~1 GB file to /dev/shm -> now reads sc16 directly.

## RESULTS 2026-11-24 - DECISION: use Owen's code (gsm_arfcn_fast.py, gsm_wide.py) with --radio bladerf
- gsm_arfcn_fast.py --arfcn 13 --radio bladerf: identical to the old BladeRF copy, every SI3 (~0.94 s).
- gsm_arfcn_fast.py --arfcn 4 9 13 22 117 --radio bladerf: all 5 decode every cycle, ~1.1 s per
  channel, ~5.6 s per full cycle.
- gsm_wide.py --radio bladerf --rate 15.36e6 (2 s, 935.2-949.8 MHz): 16 unique cells
  (ARFCN 4,6,7,8,9,10,13,14,15,16,17,18,19,20,21,22), down to -61 dBm. Trying EVERY channel
  finds far more than FFT-detect-then-decode (bandscan.py found 4-8): weak cells sit under the
  FFT threshold but still decode. 940-950 MHz: nothing decodes (likely 3G/4G, not GSM).
- A 15.36 Msps file decoded with --rate 40e6 gives 0 cells: --rate MUST match the capture.
- 40 Msps BladeRF capture WORKS: 320,000,000 bytes in 3.4 s (no drops). 16 cells in one capture.
- Full band = 24 unique cells (union of all runs): Vodafone 234-15 ARFCN 2,4,5,6,7,8,9,10,13,14,
  15,16,17,18,19,20,21,22; O2 234-10 ARFCN 113,115,117,119,120,124. Nothing at 925-935 or 940-957.
- gsm_wide.py timings (single process, Pi 5, gain 40):
    40 Msps  2 s : capture 3.4 s, decode 244 s (173 ch, 1.0->1.6 s each, rising = heat throttling?) 17 cells
    15.36    2 s : capture 3.9 s, decode  95 s ( 74 ch, avg 1.28 s)                                  13 cells
    15.36    4 s : capture 5.9 s, decode 184 s ( 74 ch, avg 2.49 s)                                  14 cells
  Per-channel cost scales with capture SECONDS, barely with sample rate -> the gr-gsm chain on each
  channel dominates, not the wide filter. Empty channels (~85%) cost as much as live ones.
  4 s found ARFCN 12 (CID 3010, new), 15, 16 but lost 2, 20, 21 -> weak cells (-55..-61) are
  hit-or-miss run to run (fading), not just capture length. Total seen so far: 25 cells.
- 3 x 15.36 Msps (centres 930.9 / 942.5 / 954.1) found all 24; 1 x 40 Msps found 16. The 40 Msps
  run missed weak cells and O2 113/115/119/120 at 957-959 MHz, 15-17 MHz off centre where the
  analog filter rolls off. Channels in the outer ~1 MHz of ANY capture get lost (e.g. ARFCN 15/16
  missed at centre 930.9, ARFCN 2/5 missed at 942.5) -> overlapping captures are more complete.

## DECISION 2026-11-25: main method = continuous single-channel tracking
  LD_LIBRARY_PATH=/usr/local/lib python3 ~/gsm_arfcn_fast.py --radio bladerf --arfcn 4
User found the wideband approaches less accurate. Improvements now focus on this live tracker
(--dwell only affects multi-channel mode; with one ARFCN it is ignored).
MEASURE EVERY CHANGE: --label NAME [--log FILE] -> live stats every 60 s and on Ctrl-C a row in
~/field/track_summary.csv: decoded / expected (SI3 every 0.9415 s from first decode), rate %, missed,
longest gap, restarts, time to first decode, pwr avg/min/max, CIDs seen. Goal = rate % up.
SWEEP BUG SUSPECT (2026-11-25): long multi-ARFCN sweeps starting on empty channels decoded NOTHING,
yet single-channel ARFCN 8 decoded instantly 1 min after the 15:12 field sweep called it "no cell".
Short lists starting on a real cell work. Hypothesis: clock_offset_control learns a bogus correction
on empty channels and retune() keeps it. Test: --fresh (rebuilds decoder per channel) vs default.
HOME REFERENCE (phone, 2026-11-25 11:40): GSM ARFCN 13 CID 942 (-55), 119 CID 22679 (-63),
116 CID 12679 (-79), 6 CID 941 (-82), 22 CID 943 (-83), 660 EE CID 51861 (-99), 9 CID 30669 (-106).
gsm_wide --fast found 13, 22, 116, 119 (all CIDs match). The narrow 900 sweep (old retune) found
NOTHING, not even ARFCN 13 (-19 on BladeRF) -> sweep retune bug confirmed; --fresh is the fix candidate.
RESULTS 2026-11-25 11:32-11:46 (home, gain 40, after self-kill fix 1bc323a):
- Label vs no label, ARFCN 13: 23 vs 24 decodes in ~100 s -> label changes nothing (as designed).
  Rate only ~27% of the SI3 maximum, 2 restarts in the unlabelled run; Pwr -16 = very hot -> test lower gain.
- Short sweep 975-980,6,9,13,22,116,119: old retune found 6,22 (MISSED 13 at -16); --fresh found 6,13,22.
  --fresh also slightly faster (27.4 s vs 29.7 s). Sweep bug confirmed; --fresh fixes it.
- Full 900 sweep --fresh (173 ch, 477 s): 6, 13, 22. Still misses 116/119 (gsm_wide finds them) and 9 (-106).
- Ctrl-C under timeout|tee sent SIGINT twice and killed tee -> summary lost. Fixed: CSV written first,
  2nd SIGINT ignored, BrokenPipe tolerated.
GAIN TEST 2026-11-25 11:51-12:03 (home, ARFCN 13 = strongest cell, 2 min each):
  g20 124/124 = 100%, 0 restarts | g30 47/123 = 38% | g40 39/120 = 32%.  -> GAIN 40 OVERLOADS.
  Pwr barely moves (-19.7 @20, -16.5 @40) for +20 dB gain = front end compressed at 30/40.
  ARFCN 119 @40: 21.5%, 8 restarts. ARFCN 116 @40: 0 (52 restarts) though gsm_wide decodes it.
  Sweep --fresh --dwell 6: 6, 13, 22, 119 found (119 newly found vs 2.5 s dwell); 116, 9 still missed.
  Added --agc (AD9361 hardware AGC on bladerf) to test against fixed gains.
  RX1/RX2: AD9361 shares LO + sample rate + bandwidth, but GAIN IS PER CHANNEL -> RX1 low gain for
  strong cells, RX2 high gain for weak ones, same frequency, is possible (needs 2-ch streaming).
GAIN RESULT 2026-11-25 12:07-12:28 (home, 2 min each, single-channel tracking):
  ARFCN 13: g10 100%, g20 100%, AGC 100%  (g30 38%, g40 32%)
  ARFCN 116: g20 100%, AGC 100%  (g40: NEVER decoded, 52 restarts)
  ARFCN 119: g20 100%, AGC 100%  (g40 21.5%, 8 restarts)
  ARFCN 6, 22: g20 100%, AGC 100%.  Zero misses, zero restarts, longest gap 0.9 s everywhere.
  -> gain 40 was the cause of nearly all lost decodes. bladerf default is now 20 (--gain 40 = old).
  O2 116/119 send the cell ID 3x per 1.9 s (gaps 0.5/0.5/0.9) -> old rate showed 150%; rate is now
  decoded/(decoded+missed). First Ctrl-C now stops cleanly, further ones ignored (h13_g10 summary was lost).
  Still to test: weak cells (9 at -106, DCS 660 at -99) at g20 vs AGC vs g40; full sweep at g20/AGC.
DUAL RESULT 2026-11-25 13:10 (home): gr-soapy streams both channels ("streaming via soapy").
  ARFCN 13: RX1 g20 100%, RX2 g45 0% (+2 dBFS, peak 1.42 = fully clipped). ARFCN 9: nothing on either,
  endless restarts. Cause: the 2 MHz analog filter also passes ARFCN 13 (0.8 MHz from 9) -> at high gain
  13 clips the ADC and buries 9; at g20 9 (-106) is below noise. Added --bw (analog filter, min 200 kHz)
  to gsm_dual.py and gsm_arfcn_fast.py to filter the neighbour out before the ADC. Merge fix: a 2nd decode
  from the same channel within 0.3 s is a new broadcast (COMBINED was 1 below RX1).
DUAL PROVEN 2026-11-25 13:18: both g20 -> RX1 100%, RX2 100%, every broadcast on both. Swapped (RX1 g45,
  RX2 g20) -> RX1 32% (clipping), RX2 100%, combined 100%. RX2 decode path works; earlier zeros = clipping.
  Added two-antenna DISCOVERY: gsm_dual.py --arfcn A B C ... = sweep (both antennas per channel, rebuilt per
  channel, found if either decodes, summary says RX1/RX2 only); gsm_wide.py --gain2 G = record RX1+RX2
  together (bladeRF-cli channel=1,2), decode RX1, then RX2 only on channels RX1 found nothing on.
  gsm_wide bladerf default gain 40 -> 20.
DUAL SWEEP 2026-11-25 13:26 (gsm_dual.py, g20/g40, 9 ch, 39 s): found 6,13,22,116,119 (all known except 9).
  RX1 g20 5/5, RX2 g40 4/5 (missed 13: clipping) -> RX1 covered it. Cells decode in 0.2-1.0 s; empty
  channels cost the full dwell. gsm_wide --gain2: bladeRF-cli n= counts BOTH channels -> fixed (2n).
WIDEBAND FIXES 2026-11-25 (gsm_wide.py, all optional; --best = all of them):
  --gain auto (0.05 s test at g20 -> gain for -6 dBFS peak), --edge-mhz 1.0 (default, was 0.3: filter
  roll-off; 11.6 MHz slice spacing still overlaps), power map of the capture (C0 beacons always on ->
  "signal above floor" list), --gate DB (skip channels below DB above floor), --rescue (channels with
  signal but no cell ID -> gsm_dual.py sweep, 300 kHz filter, RX1/RX2 g20/g45). --best = --fast
  --gain auto --gate 6 --rescue. Chris: wideband less accurate for WEAK cells (one gain per 15 MHz slice,
  strongest cell sets it) -> rescue does the weak ones narrow. Power map tested here on synthetic sc16.
WIDEBAND --best RESULT 2026-11-25 13:32-13:42: decode itself found 6,13,22 (942.5) and 116,119 (954.1) at
  0.2-0.4 s each with auto gain 18-25. The power gate was fooled by LTE: 925-935 MHz (977-1022) = Vodafone
  LTE 930.1, 64-72 = O2 LTE 952.5 -> rescue spent ~7 min on LTE. DECISION (user): stick with NARROW.
  Narrow sweep cost was ~5.6 s/channel for a 3 s dwell (2.6 s = reopening the radio). Added
  gsm_dual.py --quick: radio stays open (retune + disconnect_all + fresh decoder chains), and a channel is
  left after --sync-secs (1.2) if neither antenna has GSM sync (receiver C0 bursts). Pi 13:44: sync skip works
  (5/5, 33.5 s vs 39 s) but retune failed: disconnect_all() leaves gsm_input's inner FIR behind
  ('insufficient connected'); per-edge disconnect then failed on msg_disconnect ('Port does not exist').
  NOW: --quick keeps decoders connected (Owen-style retune: set_fc on gsm_input/clock, reset() where the
  block has it, message counts re-based). RESULT 13:48: 9 channels 5/5 in 12.5 s (was 39 s) -> the old
  'retune bug' was the gain-40 overload, not decoder state. --quick is the sweep to use.
  WRONG - FULL BAND 13:50 (173 ch, 220 s): only 6 and 119 found; 13, 22, 116 'no GSM sync'. After ~50
  empty channels the kept decoders' frequency correction has drifted (9-ch test had only 3 empties first).
  Decoder carry-over IS real. Now: gsm_dual --quick = fresh radio+decoders per channel + sync skip
  (--keep-open = the unreliable retune). gsm_arfcn_fast --fresh --sync-secs 1.2 = single-antenna sweep with
  fresh decoders (~0.3 s rebuild) + sync skip. Phone at 14:07: 1(-97) 6 13 22 114(-90) 116 119 + EE 649.
  RESULT 13:57 gsm_arfcn_fast --fresh --sync-secs 1.2 --gain 20, full 900: 7/7 phone cells in 257 s incl.
  weak 114 (-90) and 1 (-97), CIDs match. 649 (EE 1800, -92): no sync at g20/g40/dual - antenna? parked.
Also LTE: O2 816 MHz decoded = eNB 541524 sector 120 TAC 50496 PLMN 234-10 (matches phone).
Next decode-rate candidates: restart threshold (2 s fixed now), fixed --ppm, gain/overload,
unbounded message_debug memory on long runs.

## (superseded) UPGRADE PLAN (agreed order - accuracy first, then speed)
1. DONE (needs Pi validation): gsm_wide.py --skip-empty / --stop-early / --fast / --probe-secs.
   Probe = decode first PROBE_SECS; receiver C0 bursts only appear after FCCH+SCH lock, so 0 bursts
   = no GSM carrier -> skip. Also reports "SYNCED but no cell ID" channels (carrier there, ID not
   decoded = candidates for longer capture). Default is still Owen's full decode.
   VALIDATE on ~/test_942_15M_2s.sc16 (actually 4 s @15.36, centre 942.5): full vs --fast cell
   lists must match.
   RESULT v1 (same file): full 160.8 s / 14 cells; --skip-empty(1 s) 87.1 s / 12 (+19,21 synced);
   --fast 50.5 s / 13. KEY FINDINGS: (a) the SAME file decodes differently run to run (19 decoded in
   full + fast, not in skip-empty; 12 only in skip/fast) -> gr-gsm clock_offset_control feedback is
   asynchronous, so results depend on CPU timing. Part of the 'weak cells come and go' is the
   decoder, not the radio. (b) 1 s probe missed ARFCN 5 and 15 (late lock). (c) full-run per-channel
   time rose 1.6->2.5 s for identical work = thermal throttling.
   v2: single-pass probe (abort only if no sync within --probe-secs, default 2.0; needs
   head.nitems_written, falls back to two-pass), --retries N (default 1 with --fast) for carriers
   that synced without an SI3.
   RESULT v2 (--fast, same file): 77.5 s / 14 cells (full: 160.8 s / 14). Union over all runs of
   this file = 15 (full missed 12, fast missed 5) -> --fast is as accurate as full at ~2x speed.
   Time now: 58 empty channels = 64 s (83%); 14 cells = 5.9 s total. Empty-channel cost rose
   0.8->1.2 s over the run = heat throttling again.
   Probe only saves time when the capture is LONGER than --probe-secs: use --secs 4 with --fast.
   Hard channels (5 no-sync; 7, 21 synced-no-ID) all sit next to much stronger cells (ARFCN 4 is
   -30, 21 dB above 5) -> adjacent-channel interference suspected; Owen's 2nd filter passes to
   150 kHz. Test a tighter channel filter as an option.
2. Cooling (Pi 5 Active Cooler); then 2 decodes in parallel.
3. Repeat-and-merge: N captures, union of cells with seen-count (confidence).
4. Longer captures (cheap once 1 is in: only live channels pay for the extra seconds).
5. RX2 diversity (same centre, 2nd antenna >=16 cm away), merge both.
6. gsm_wide prints a ready gsm_arfcn_fast.py --arfcn ... command for found cells (discover -> track/DF).
7. Later: polyphase channelizer (one pass for all channels) if still too slow.

## PENDING / PLANNED (for the next AI)
0. First Pi results (gain 40, E-GSM900): wide found 4 cells in 28.6 s; old-style narrow found
   ARFCN 4,9,10,13,18,22,113,117 before ssh dropped. Wide at gain 40 MISSED 13/18/22 in FFT
   detection (median noise floor too high) -> now 25th percentile (--noise-pct). Re-test.
   Many +6..+11 dB clusters (ARFCN 977-983, 62-67) never decode: probably 3G/4G, not GSM.
1. Run `bandscan.py --compare --timing-csv` on the Pi and record real numbers here.
2. If parallel decodes cause capture overflows with --pipeline, leave pipeline off.
3. DF (Kraken) as an opt-in flag on the scanner.
4. Make libbladeRF 2.6.1 permanent (no LD_LIBRARY_PATH); rebuild SoapyBladeRF against it for live path.
5. Test BladeRF LIVE decode (gsm_arfcn_fast.py via osmosdr bladerf=0 + new lib). If clean, use like Lime.
6. Wire BladeRF into Owen's --df / Kraken DF system.
7. RX2 (2nd receiver): NOT two bands at once - the AD9361 has ONE shared RX LO, both RX channels
   tune to the same centre. Use it for antenna DIVERSITY (2nd antenna >=16 cm away, decode both,
   merge cells) or 2-element DF. Costs: 2x USB data (40 Msps x2 too much; 15.36 x2 fine),
   2x /dev/shm, 2x decode time. Try --secs 4 first (free) for weak cells.

## SETUP FROM SCRATCH
apt: gnuradio gr-osmosdr gr-gsm soapysdr-tools soapysdr-module-bladerf python3-numpy
     python3-soapysdr bladerf libbladerf-dev
udev: plugdev rule for BladeRF if "insufficient permissions".
Build libbladeRF 2.6.1 (deps: build-essential cmake git libusb-1.0-0-dev pkg-config
     libcurl4-openssl-dev libncurses-dev) -> cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. ; make ; make install.
Flash firmware v2.6.0. Copy gsm_arfcn_clean.py to /home/<user>/.


## LTE ALL CELLS (2026-11-25 afternoon, lte_decode.py options; default behaviour unchanged)
Baseline re-checked 14:48: 806 -> 4968724 (eNB 19409 s20), 816 -> 138630264 (eNB 541524 s120), ~18.6 s each.
"srsRAN only does the strongest" = live mode. In file mode we pass -c PCI, so any PCI the Python search
finds can be decoded. New options:
  --all-pci (srsran engine): srsRAN per PCI found, labelled per PCI.
  --fast-srsran: stop pdsch_ue at the first SIB1 line (default runs widths 50/25/75/100 regardless = most of
                 the 16 s); skip other widths once the MIB's own width was tried.
  --weak-search: stock search picks ONE freq offset per PSS root (global max) - a strong cell's false peak at
                 a wrong offset (-39 kHz) hid a -10 dB cell. Second pass at the strongest cell's offset.
  --pss-min N:   lower sync threshold for the second pass (stock 6).
  Synthetic: stock finds only strongest; --weak-search finds -10 dB cells; + --pss-min 4 finds -16 dB;
  -20 dB not. Detected != decoded: a weak co-channel cell's SIB1 overlaps the strong cell (may need
  interference cancellation later). Phone references: O2 816 PCI 0+1, O2 952.5 PCI 0+1, EE 1846.7 389+461.
  PI 15:1x: all three options -> 816 still only PCI 1 (phone: PCI 0 -64, PCI 1 -62 = NOT weak); search took
  12-19 s (whole 0.2 s capture). Added: capture level/clipping print (LTE default gain 40 = GSM overload
  setting!), --search-secs 0.08, --debug-search (every sync peak per PSS group + accept/reject reason),
  align_narrow(): weak cell 'frame start not found' because H._align correlates PSS over the full 10 MHz
  where the strong cell drowns it -> fallback aligns at 1.92 Msps, SSS averaged over all repeats; exact on
  synthetic (-10 dB where stock fails). Full phone LTE list 15:14 saved in the conversation (13 cells).
  PI 15:2x (816, debug): g40 CLIPS (peak 1.00) -> use --gain 20 for LTE too. Group 0 (PCI 0's PSS) has NO
  real peak at either gain (best 9.7x, SSS 0.35 = noise; PCI 1 30x SSS 0.98): PCI 0 (sector 110 of the same
  eNB 541524 as PCI 1/sector 120) is not reaching our antenna as a detectable sync signal.
  Added --siblings: for each decoded PCI p, try 3*(p//3)+0/1/2 with p's frame alignment (co-sited sectors
  share the eNB clock); srsRAN CRC rejects absent PCIs.
  PI 15:3x --siblings g20: 952.5 PCI 0 (sibling of 1) = MIB OK, no SIB1 -> a cell the sync search cannot see
  IS there (phone 12:46 listed 541524:112); SIB1 lost to co-channel overlap with PCI 1. 816 PCI 0: MIB fails.
  806/930.1: single cell (= phone). 2132.3 (15 MHz) / 2144.7 (20 MHz): capture only 15.36 Msps + level
  -41..-46 dBFS at g20. Added: --auto-gain (clip -> -10 dB, peak<0.15 -> up to +20 dB, recapture once);
  wider re-capture in the srsran --all-pci path when the MIB says 75/100 RB (23.04/30.72 Msps);
  align_wide() = Owen's _align with nfft from the rate (his was fixed 1024); first srsRAN width = the
  capture's own. Try --secs 1.0 for SIB1 of overlapped cells (50 SIB1 chances instead of 10).

## LTE srsRAN (for lte_decode.py --engine srsran)
*** LTE DECODE CONFIRMED 2026-11-25 10:57 (home): BladeRF capture at 806 MHz (15.36 Msps, low-IF 0.6)
-> srsRAN pdsch_ue file mode -> SIB1 cell identity 4968724 = eNB 19409 sector 20, TAC 57379,
PLMN 234-15 (Vodafone). eNB + TAC match the phone scanner exactly. 10 SIB1s in a 0.2 s capture
(SIB1 every 20 ms) = every one decoded. The all-Python MIB decoder still fails by 3-5 bits on air.
First on-air test 2026-11-25 (home, 800 MHz): survey found 10 MHz carriers at 801-811 (centre 806)
and 811.5-821 (centre 816). Python engine found PCI 21 @806 and PCI 1 @816 but MIB failed by 3-5
bits every time (also with --low-if-mhz 0.6). Owen's on-air HackRF path used srsRAN pdsch_ue for
MIB+SIB1, so build it (paths are hard-coded to /home/revector):
  sudo apt install -y build-essential cmake libfftw3-dev libmbedtls-dev libboost-program-options-dev \
       libconfig++-dev libsctp-dev
  cd ~ && git clone --depth 1 https://github.com/srsran/srsRAN_4G.git
  python3 ~/fbs-upload/"osprey 2"/srsran_watch/print_sib1.py      # adds the SIB1 cell-ID print
  mkdir -p ~/srsRAN_4G/build && cd ~/srsRAN_4G/build && cmake .. && make -j3 pdsch_ue
Then: python3 ~/fbs-upload/lte_decode.py --freq-mhz 806 --low-if-mhz 0.6 --engine srsran
Built OK on the Pi 2026-11-25 incl. srsRAN's own BladeRF driver (libsrsran_rf_blade.so). Live:
  LD_LIBRARY_PATH=/usr/local/lib:$HOME/srsRAN_4G/build/lib/src/phy/rf ~/srsRAN_4G/build/lib/examples/pdsch_ue -f 806000000 -g 40
  (NO -a "device_name=bladerf": srsRAN passes it to bladerf_open -> "Invalid operation or parameter").
  Cell ID lines go to ~/sib1_cell.txt. Phone scanner reference at home: Vodafone PCI 21 eNB 19409
  TAC 57379; O2 PCI 1 eNB 541524 TAC 50496; B20 EARFCN 6175 (793.5) PCI 272/216.
