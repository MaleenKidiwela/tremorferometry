#!/bin/bash
# Pairwise (event-grade) dv/v rollout: natural series (no injection), all stations, 2-4 s.
# Produces data/pw_<STA>_B?_series.csv via dvv_pairwise.py --amp 0; top-30 coverage patches per station.
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python
LOG=logs/run_pw.log; : > $LOG
for NPZ in data/long_window_daily_*.npz; do
  S=$(basename $NPZ .npz); S=${S#long_window_daily_}
  case $S in *_e23|*_despike) continue;; esac
  echo "=== $S pairwise ===" | tee -a $LOG
  $PY scripts/dvv_pairwise.py --npz $NPZ --station $S --amp 0 --n-patches 30 --workers 6 \
      --out-prefix pw >> $LOG 2>&1
done
# GNW (e23 stacks) and PGC (despike stacks) use their canonical npz
$PY scripts/dvv_pairwise.py --npz data/long_window_daily_GNW_e23.npz --station GNW --amp 0 --n-patches 30 --workers 6 --out-prefix pw >> $LOG 2>&1
$PY scripts/dvv_pairwise.py --npz data/long_window_daily_PGC_despike.npz --station PGC --amp 0 --n-patches 30 --workers 6 --out-prefix pw >> $LOG 2>&1
echo "PW ALL DONE" | tee -a $LOG
