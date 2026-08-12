#!/bin/bash
# CPU tail (retired-reverse method): FORWARD 3-comp stacks -> forward Z dv/v -> finalize_causality.
# Reverse certification retired (causality = fwd coda/fwd mirror reproduces it at 99%). No rev stacks,
# no H1/H2 dv/v. Runs in background so the next station's GPU densify overlaps it.
# Usage: scripts/borehole_cpu_tail.sh <STA>
cd /home/jovyan/tremorferometry
STA=$(echo "$1" | tr '[:lower:]' '[:upper:]'); sta=$(echo "$STA" | tr '[:upper:]' '[:lower:]')
PY=/home/jovyan/envs/tremorferometry/bin/python
STATUS=logs/rollout_status.log; mkdir -p logs
log(){ echo "$(date +%H:%M) $*" >> "$STATUS"; echo "$(date +%H:%M) $*"; }
TAG=${STA}p90f40

# FORWARD stacks only (3-comp written together; only Z used by causality finalizer)
if [ ! -f data/long_window_daily_${TAG}_Z.npz ]; then
  log "$STA fwd stacks start"
  $PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_${sta}p90f40_*.csv" \
     --network PB --station "$STA" --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 \
     --despike-mad 8 --workers 10 --out-prefix data/long_window_daily_${TAG} \
     > logs/stack_${sta}f40.log 2>&1 || { log "$STA ABORT: fwd stacks failed"; exit 9; }
  log "$STA fwd stacks done"
fi

# FORWARD Z dv/v only
OUT=data/daily_dvv_${TAG}_Z_2to4.csv
if [ ! -f "$OUT" ]; then
  $PY scripts/dvv_roll30cal.py --station "$STA" --npz data/long_window_daily_${TAG}_Z.npz \
     --window 2 4 --origin-anchor --workers 10 --out "$OUT" > logs/dvv_${sta}_Z.log 2>&1 \
     || { log "$STA ABORT: dvv Z failed"; exit 11; }
fi
log "$STA fwd Z dv/v done"

# Causality finalize (Z-only)
$PY scripts/finalize_causality.py "$STA" "$TAG" > logs/finalize_${sta}.log 2>&1 \
   || { log "$STA ABORT: causality finalize failed"; exit 12; }
zc=$($PY -c "import json;print(json.load(open('data/${sta}_3comp_summary.json'))['z_certified'])" 2>/dev/null)
zs=$($PY -c "import json;print(json.load(open('data/${sta}_3comp_summary.json')).get('z_std_pct'))" 2>/dev/null)
log "$STA CAUSALITY FINALIZE done: Z ${zc} causality-cert fam, std ${zs} %"
if [ "${zc:-0}" -lt 20 ] 2>/dev/null; then log "$STA FLAG: only ${zc} causality-cert families (<20) — likely starved or genuinely weak"; fi
