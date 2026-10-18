#!/bin/bash
# E4 step 4 — association RE-NULL at tolerances matched to the measured S-pick residual.
#   pick_match_tolerance 1.2 s (was 2.5), edt_pick_std 0.5 s (was 1.5), zlim 20-50 km (was 0-70),
#   n_picks = 4 and 5 run DIRECTLY (E3 ran n_picks=3 and post-filtered on n_stations),
#   thresholds 0.3 and 0.5, 5 constant-shift seeds + 2 hourly-shift seeds per configuration.
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair/results; L=lfe_catalog_repair/logs
A="python lfe_catalog_repair/scripts/e2_associate.py --threads 7 --pick-tol 1.2 --pick-std 0.5 --zlim 20,50"
obs() {  # thr N
  local F=$R/e4_assoc_vi2010_cascadia_t$1_n$2.csv
  [ -f "$F" ] && { echo "[skip] $F"; return; }
  $A --set vi2010 --model cascadia --threshold $1 --n-stations $2 --min-stations $2 --out $F \
     > $L/e4_assoc_t$1_n$2.log 2>&1 && echo "[ok] $F"
}
nul() {  # thr N seed mode
  local F=$R/e4_assocnull_vi2010_cascadia_t$1_n$2_$4_s$3.csv
  [ -f "$F" ] && { echo "[skip] $F"; return; }
  $A --set vi2010 --model cascadia --threshold $1 --n-stations $2 --min-stations $2 \
     --shift-null $3 --shift-mode $4 $( [ "$4" = station_hour ] && echo "--shift-min 60 --shift-max 300" ) --out $F > $L/e4_assocnull_t$1_n$2_$4_s$3.log 2>&1 \
     && echo "[ok] $F"
}
export -f obs nul; export R L A
JOBS=()
for T in 0.3 0.5; do for N in 4 5; do
  JOBS+=("obs $T $N")
  for S in 1 2 3 4 5; do JOBS+=("nul $T $N $S station"); done
  for S in 1 2;         do JOBS+=("nul $T $N $S station_hour"); done
done; done
printf '%s\n' "${JOBS[@]}" | xargs -P 3 -I{} bash -c '{}'
echo "[e4 assoc null ALL DONE]"
