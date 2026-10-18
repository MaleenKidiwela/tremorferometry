#!/bin/bash
# E3 step 1 driver (2026-09-11) — chance-association null: per-station constant pick shifts.
# 25 PyOcto runs: vi2010 cascadia t0.2/0.3/0.5 x 5 seeds, vi2017 cascadia t0.3 x 5, vi2010 jcms t0.3 x 5.
# 4 concurrent runs x 7 PyOcto threads = 28 threads (the pod's cap is 32; nproc lies with 176).
cd /home/jovyan/tremorferometry
R=lfe_catalog_repair/results; L=lfe_catalog_repair/logs
A="python lfe_catalog_repair/scripts/e2_associate.py --threads 7 --n-stations 3 --min-stations 3"
run() {  # set model thr seed
  local F=$R/e3_nullevents_$1_$2_t$3_s$4.csv
  [ -f "$F" ] && { echo "[skip] $F"; return; }
  $A --set $1 --model $2 --threshold $3 --shift-null $4 --out $F \
     > $L/e3_null_$1_$2_t$3_s$4.log 2>&1 && echo "[ok] $F"
}
export -f run; export R L A
JOBS=()
for S in 1 2 3 4 5; do
  for T in 0.2 0.3 0.5; do JOBS+=("vi2010 cascadia $T $S"); done
  JOBS+=("vi2017 cascadia 0.3 $S")
  JOBS+=("vi2010 jcms 0.3 $S")
done
printf '%s\n' "${JOBS[@]}" | xargs -P 4 -I{} bash -c 'run {}'
echo "[e3 null ALL DONE]"
