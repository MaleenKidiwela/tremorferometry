#!/bin/bash
# 2-4 s rolling-30 + SVD-Wiener dv/v for all finalized margin stations.
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python
LOG=logs/run_2to4_roll.log; : > $LOG

# boreholes (B004 already done in the test) + non-boreholes that use the base npz
STA="B011 B013 B941 B018 B023 B928 B026 B014 B204 B028 B030 B032 B033 B036 B039 B040 B927 B020 B022 B935 B201 HDW NLLB COR COLT"
for S in $STA; do
  NPZ=data/long_window_daily_$S.npz
  [ -f "$NPZ" ] || { echo "MISSING $NPZ" | tee -a $LOG; continue; }
  echo "=== $S ===" | tee -a $LOG
  $PY scripts/dvv_2to4_roll30svd.py --station $S --npz $NPZ \
      --out data/daily_dvv_${S}_2to4_roll.csv --workers 28 >> $LOG 2>&1
done
# GNW / PGC use their finalized cleaned daily-stack sets
echo "=== GNW (e23) ===" | tee -a $LOG
$PY scripts/dvv_2to4_roll30svd.py --station GNW --npz data/long_window_daily_GNW_e23.npz \
    --out data/daily_dvv_GNW_2to4_roll.csv --workers 28 >> $LOG 2>&1
echo "=== PGC (despike) ===" | tee -a $LOG
$PY scripts/dvv_2to4_roll30svd.py --station PGC --npz data/long_window_daily_PGC_despike.npz \
    --out data/daily_dvv_PGC_2to4_roll.csv --workers 28 >> $LOG 2>&1
echo "ALL DONE" | tee -a $LOG
grep -E 'DONE|MISSING' $LOG
