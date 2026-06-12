#!/bin/bash
# Pipeline runner - called by agent to process stations
# Usage: bash run_pipeline.sh <station> <step>
set -e
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python

STA=$1
sta=$(echo $STA | tr '[:upper:]' '[:lower:]')
STEP=$2

case $STEP in
  dvv)
    $PY scripts/dvv_coda_parallel.py --npz data/long_window_daily_${STA}.npz \
      --window 1.0 4.0 --station $STA \
      --out-csv data/daily_dvv_${STA}_coda_1to4.csv \
      --out-fig figures/smoke_dvv_${STA}_raw.png --workers 20
    ;;
  rolling)
    for W in "1.0 3.0 1to3" "2.0 4.0 2to4" "3.0 5.0 3to5"; do
      read w1 w2 wname <<< $W
      $PY scripts/dvv_roll30cal.py --station $STA \
        --npz data/long_window_daily_${STA}.npz \
        --window $w1 $w2 \
        --out data/daily_dvv_${STA}_${wname}_cal.csv --workers 16
    done
    ;;
  *)
    echo "Unknown step: $STEP"
    exit 1
    ;;
esac
