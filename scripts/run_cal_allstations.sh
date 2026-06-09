#!/bin/bash
# TRUE 30-calendar-day trailing rolling dv/v for windows 1-3, 2-4, 3-5, all 28 stations.
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python
LOG=logs/run_cal.log; : > $LOG
STA="B004 B011 B013 B941 B018 B023 B928 B026 B014 B204 B028 B030 B032 B033 B036 B039 B040 B927 B020 B022 B935 B201 HDW NLLB COR COLT"
for WIN in "1.0 3.0 1to3" "2.0 4.0 2to4" "3.0 5.0 3to5"; do
  set -- $WIN; W1=$1; W2=$2; TAG=$3
  for S in $STA; do
    NPZ=data/long_window_daily_$S.npz
    [ -f "$NPZ" ] || continue
    echo "=== $S $TAG ===" | tee -a $LOG
    $PY scripts/dvv_roll30cal.py --station $S --npz $NPZ --window $W1 $W2 \
        --out data/daily_dvv_${S}_${TAG}_cal.csv --workers 28 >> $LOG 2>&1
  done
  # GNW (e23) and PGC (despike)
  $PY scripts/dvv_roll30cal.py --station GNW --npz data/long_window_daily_GNW_e23.npz --window $W1 $W2 \
      --out data/daily_dvv_GNW_${TAG}_cal.csv --workers 28 >> $LOG 2>&1
  $PY scripts/dvv_roll30cal.py --station PGC --npz data/long_window_daily_PGC_despike.npz --window $W1 $W2 \
      --out data/daily_dvv_PGC_${TAG}_cal.csv --workers 28 >> $LOG 2>&1
  echo "### $TAG window done ###" | tee -a $LOG
done
echo "ALL DONE" | tee -a $LOG
