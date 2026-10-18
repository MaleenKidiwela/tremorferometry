#!/bin/bash
# E2 step 3a driver (resumed run 2026-09-10): PyOcto association sweep + scoring.
cd /home/jovyan/tremorferometry
export OMP_NUM_THREADS=8
R=lfe_catalog_repair/results; L=lfe_catalog_repair/logs
A="python lfe_catalog_repair/scripts/e2_associate.py --threads 8 --n-stations 3 --min-stations 3"
S="python lfe_catalog_repair/scripts/e2_assoc_score.py"
rm -f $R/e2_assoc_scores.csv
for SET in vi2010 vi2017; do
  for M in cascadia jcms; do
    for T in 0.2 0.3 0.5; do
      F=$R/e2_events_${SET}_${M}_t${T}.csv
      [ -f "$F" ] || $A --set $SET --model $M --threshold $T
      for N in 3 4 5; do
        [ -f "$F" ] && $S --events $F --min-stations $N
      done
    done
  done
done
echo "[e2 assoc ALL DONE]"
