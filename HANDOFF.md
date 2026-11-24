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
- blade_capture.py / capture_ram.py   Soapy IQ capture (superseded by native bladeRF-cli capture).
- iqview.py          ASCII spectrum of a capture.  compare.py  diff two captures (offset/inversion).
- srsran_watch/      older separate LTE cell-ID work.

## PENDING / PLANNED (for the next AI)
1. Optional CLI flags so features (wideband, gain-sweep, DF) are opt-in per run.
2. Full 925-960 scan: 3 wide chunks (~930 / 942.5 / 955 MHz) x gain sweep [35,40,45], merged report.
   widescan.py already does ONE chunk; wrap it in a loop over chunks+gains and dedupe cells.
3. Dynamic gain sweep (Chris's idea): capture each chunk at several gains, merge, to catch strong+weak.
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
