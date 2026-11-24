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
7. RX2 (2nd receiver): two bands at once (speed) OR antenna diversity (sensitivity). Doubles USB/CPU load.

## SETUP FROM SCRATCH
apt: gnuradio gr-osmosdr gr-gsm soapysdr-tools soapysdr-module-bladerf python3-numpy
     python3-soapysdr bladerf libbladerf-dev
udev: plugdev rule for BladeRF if "insufficient permissions".
Build libbladeRF 2.6.1 (deps: build-essential cmake git libusb-1.0-0-dev pkg-config
     libcurl4-openssl-dev libncurses-dev) -> cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. ; make ; make install.
Flash firmware v2.6.0. Copy gsm_arfcn_clean.py to /home/<user>/.

