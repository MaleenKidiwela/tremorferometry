#!/bin/bash
# E3 step 1 sensitivity — Muenchmeyer's own control shape: per-station-PER-HOUR shifts,
# |shift| uniform in [60, 300] s (>> the ~30 s max inter-station travel-time difference across this
# network, so real coincidence is still destroyed, but the burst envelope stays aligned at the
# few-minute scale -> the CONSERVATIVE, higher-chance end of the null).
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair/results; L=lfe_catalog_repair/logs
A="python lfe_catalog_repair/scripts/e2_associate.py --threads 7 --n-stations 3 --min-stations 3 \
   --shift-mode station_hour --shift-min 60 --shift-max 300"
run() {
  local F=$R/e3_nullhr_$1_$2_t$3_s$4.csv
  [ -f "$F" ] && { echo "[skip] $F"; return; }
  $A --set $1 --model $2 --threshold $3 --shift-null $4 --out $F \
     > $L/e3_nullhr_$1_$2_t$3_s$4.log 2>&1 && echo "[ok] $F"
}
export -f run; export R L A
JOBS=()
for S in 1 2 3 4 5; do
  JOBS+=("vi2010 cascadia 0.2 $S"); JOBS+=("vi2010 cascadia 0.3 $S"); JOBS+=("vi2017 cascadia 0.3 $S")
done
printf '%s\n' "${JOBS[@]}" | xargs -P 4 -I{} bash -c 'run {}'
echo "[e3 nullhr ALL DONE]"
