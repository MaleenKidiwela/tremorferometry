#!/bin/bash
cd /home/jovyan/tremorferometry
# CPU discovery for buffered stations: candidate-detect + picker-score (clustering=GPU, done later)
# STA lat lon
declare -a S=( "B004 48.20 -124.43" "B001 48.04 -123.13" "B003 48.06 -124.14" "B927 49.22 -124.81" "B928 48.83 -125.13" )
for row in "${S[@]}"; do
  set -- $row; STA=$1; LA=$2; LO=$3
  la0=$(python3 -c "print(round($LA-0.9,2))"); la1=$(python3 -c "print(round($LA+0.9,2))")
  lo0=$(python3 -c "print(round($LO-1.35,2))"); lo1=$(python3 -c "print(round($LO+1.35,2))")
  cand=data/${STA,,}_pnsn_candidates_100km.parquet
  echo "=== [disc] $STA candidate detection $(date +%H:%M) ==="
  if [ ! -f "$cand" ]; then
    PYTHONPATH=src /home/jovyan/envs/tremorferometry/bin/python scripts/discover_nllb_pnsn_driven.py \
      --station $STA --wfdir data/waveforms --bbox $la0 $la1 $lo0 $lo1 \
      --pnsn catalogs/pnsn_tremor_cascadia_full.csv --fs 100 --candidates-only \
      --candidates-out $cand --workers 12
  fi
  echo "=== [disc] $STA picker scoring $(date +%H:%M) ==="
  cp -f lfe_features/models/tremor_picker_b011.joblib lfe_features/models/tremor_picker_${STA,,}.joblib
  PYTHONPATH=src /home/jovyan/envs/tremorferometry/bin/python lfe_features/score_candidates.py \
    --net PB --sta $STA --cand $cand --y0 2010 --y1 2026 --thr 0.7 --workers 12
done
echo "=== [disc] ALL BUFFERED STATIONS PREPPED (candidates+scores; clustering pending GPU) ==="
