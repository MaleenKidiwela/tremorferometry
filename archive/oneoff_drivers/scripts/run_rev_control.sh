#!/bin/bash
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
for spec in "B926 Z" "B926 H2" "B926 H1" "B011 Z" "B011 H2"; do
  set -- $spec; S=$1; C=$2; TAG=${S}p90f40
  npz=data/long_window_daily_${TAG}rev_${C}.npz
  out=data/daily_dvv_${TAG}rev_${C}_2to4.csv
  [ -f "$out" ] && { echo "skip $out"; continue; }
  [ -f "$npz" ] || { echo "MISSING $npz"; continue; }
  echo "=== rev dv/v $S $C ==="
  $PY scripts/dvv_roll30cal.py --station $S --npz $npz --window 2 4 --origin-anchor --workers 6 --out $out
done
echo "=== CONTROL RESULT ==="
$PY scripts/rev_dvv_control.py
echo "REV_CONTROL_DONE"
