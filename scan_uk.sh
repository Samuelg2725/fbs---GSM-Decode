#!/bin/bash
# Field scan of UK GSM bands with the BladeRF, using gsm_wide.py.
#   ./scan_uk.sh                 # all UK 2G: 900 (3 captures) + 1800 (7 captures)
#   ./scan_uk.sh 900  carpark    # just 900, label "carpark"
#   ./scan_uk.sh 1800 roof
# Settings via environment (all optional):
#   GAIN=45 ./scan_uk.sh 1800    SECS=4 (default)    EXTRA="--fast" (default; "" = Owen's full decode)
#   KEEP=1 ./scan_uk.sh          also keep every capture file for re-decoding later (246 MB each)
BAND=${1:-uk}
LABEL=${2:-$(date +%H%M)}
GAIN=${GAIN:-40}
SECS=${SECS:-4}
EXTRA=${EXTRA---fast}
OUT=~/field/$(date +%Y%m%d)_${LABEL}_${BAND}
C900="930.9 942.5 954.1"
C1800="1810.5 1821.2 1831.8 1842.5 1853.2 1863.8 1874.5"
case "$BAND" in
  900)  CS="$C900" ;;
  1800) CS="$C1800" ;;
  uk)   CS="$C900 $C1800" ;;
  *)    echo "band must be 900, 1800 or uk"; exit 1 ;;
esac
mkdir -p "$OUT"
echo "scan $BAND  label $LABEL  gain $GAIN  secs $SECS  opts '$EXTRA'  -> $OUT" | tee "$OUT/summary.txt"
START=$SECONDS
for c in $CS; do
  echo "===== centre $c MHz =====" | tee -a "$OUT/summary.txt"
  T0=$SECONDS
  python3 ~/gsm_wide.py --radio bladerf --rate 15.36e6 --center "$c" --secs "$SECS" --gain "$GAIN" $EXTRA \
    2>&1 | tee "$OUT/scan_$c.log"
  echo "  centre $c took $((SECONDS - T0)) s" | tee -a "$OUT/summary.txt"
  [ -n "$KEEP" ] && cp /dev/shm/gsm_wide.sc16 "$OUT/cap_${c}_15M_${SECS}s_g${GAIN}.sc16"
done
{
  echo
  echo "===== ALL UNIQUE CELLS ($BAND, label $LABEL) ====="
  grep -h "^  ARFCN:" "$OUT"/scan_*.log | sort -u | sort -t: -k2 -n
  echo "total: $(grep -h "^  ARFCN:" "$OUT"/scan_*.log | sort -u | wc -l) cell(s)"
  echo
  echo "carrier present but no cell ID:"
  grep -h -A1 "^GSM carrier present" "$OUT"/scan_*.log | grep -v "^GSM carrier\|^--"
  echo
  echo "whole scan took $((SECONDS - START)) s"
} | tee -a "$OUT/summary.txt"
