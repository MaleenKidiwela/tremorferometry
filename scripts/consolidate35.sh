#!/bin/bash
# Bring the 7 new stations up to full product parity (true-scale _calT + _calT_des + pairwise),
# then rerun the 35-station multi-window inversion on true-scale data.
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/consolidate35.log; : > $L
say(){ echo "=== $(date +%H:%M:%S) $* ===" | tee -a $L; }
NEW="B017 B001 B005 B003 B045 B932 B049"

say "A. true-scale _calT for the 7 new stations (origin-anchored, 3 windows)"
for S in $NEW; do
  for W in "1.0 3.0 1to3" "2.0 4.0 2to4" "3.0 5.0 3to5"; do set -- $W
    [ -f data/daily_dvv_${S}_${3}_calT.csv ] || $PY scripts/dvv_roll30cal.py --station $S \
      --npz data/long_window_daily_${S}.npz --window $1 $2 --origin-anchor \
      --out data/daily_dvv_${S}_${3}_calT.csv --workers 12 >> $L 2>&1
  done
done

say "B. deseason -> _calT_des for the 7 new stations"
for S in $NEW; do
  for W in 1to3 2to4 3to5; do
    $PY scripts/deseason_cal.py --glob "data/daily_dvv_${S}_${W}_calT.csv" >> $L 2>&1
  done
done

say "C. pairwise (event-grade) for the 5 new stations missing it"
for S in B005 B003 B045 B932 B049; do
  [ -f data/pw_${S}_series.csv ] || $PY scripts/dvv_pairwise.py --npz data/long_window_daily_${S}.npz \
    --station $S --amp 0 --n-patches 30 --workers 6 --out-prefix pw >> $L 2>&1
done

say "D. rerun 35-station multi-window inversion on TRUE-SCALE (deseasoned + raw)"
SFX=_calT_des DESEASON=1 OUT=fault_tomography/inversion/fault_4d_multiwin_calT35.npz \
  $PY fault_tomography/inversion/invert_multiwin.py >> $L 2>&1
SFX=_calT DESEASON=0 OUT=fault_tomography/inversion/fault_4d_multiwin_calT35_raw.npz \
  $PY fault_tomography/inversion/invert_multiwin.py >> $L 2>&1

say "CONSOLIDATE35 DONE"
