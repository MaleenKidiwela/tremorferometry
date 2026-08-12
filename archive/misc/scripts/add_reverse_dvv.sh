#!/bin/bash
# Add the REVERSED noise-reference products to an already-forward-finalized station: build reversed
# stacks from the existing reversed mf catalog, then reversed Z dv/v (for the mandatory dv/v noise
# correction; certification stays forward/causality). CPU-only (reversed densify already ran).
# Usage: scripts/add_reverse_dvv.sh <STA>
cd /home/jovyan/tremorferometry
STA=$(echo "$1" | tr '[:lower:]' '[:upper:]'); sta=$(echo "$STA" | tr '[:upper:]' '[:lower:]')
PY=/home/jovyan/envs/tremorferometry/bin/python
STATUS=logs/rollout_status.log; log(){ echo "$(date +%H:%M) $*" >> "$STATUS"; echo "$(date +%H:%M) $*"; }
TAG=${STA}p90f40

[ "$(ls data/mf_${sta}p90f40rev_*.csv 2>/dev/null | wc -l)" -gt 0 ] || { log "$STA rev-add ABORT: no reversed mf files"; exit 1; }

if [ ! -f data/long_window_daily_${TAG}rev_Z.npz ]; then
  log "$STA rev stacks start (noise-reference for dv/v correction)"
  $PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_${sta}p90f40rev_*.csv" \
     --network PB --station "$STA" --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 \
     --despike-mad 8 --workers 8 --out-prefix data/long_window_daily_${TAG}rev \
     > logs/stack_${sta}f40rev.log 2>&1 || { log "$STA rev-add ABORT: rev stacks failed"; exit 2; }
  log "$STA rev stacks done"
fi

OUT=data/daily_dvv_${TAG}rev_Z_2to4.csv
if [ ! -f "$OUT" ]; then
  $PY scripts/dvv_roll30cal.py --station "$STA" --npz data/long_window_daily_${TAG}rev_Z.npz \
     --window 2 4 --origin-anchor --workers 8 --out "$OUT" > logs/dvv_${sta}_Zrev.log 2>&1 \
     || { log "$STA rev-add ABORT: rev Z dvv failed"; exit 3; }
fi
log "$STA REVERSE dv/v ADDED (rev stacks + rev Z dv/v) — ready for noise correction"
