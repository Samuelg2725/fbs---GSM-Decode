# Experimental LTE scripts (parked 2026-11-25)

Built in a rush on 2026-11-25 afternoon; set aside so the basic, proven LTE tool is the only one in use.

- lte_scan.py  - multi-band scan (known UK EARFCN list -> sync probe -> srsRAN). Did decode 6 cells at home,
                 but had several rushed changes (early stop, align window) layered on top.
- lte_track.py - repeated capture tracking; srsRAN live mode never locked on the BladeRF.

They import lte_decode.py from the folder above; to try one, run it from the repo root:
  python3 experimental_lte/lte_scan.py --earfcns 6300 6400

The basic tool is ../lte_decode.py (unchanged since the 10:57 Vodafone/O2 decode).
