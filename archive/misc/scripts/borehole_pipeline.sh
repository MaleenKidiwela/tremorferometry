#!/bin/bash
# Per-station borehole 3-comp dv/v pipeline (replicates B011's known-good provenance).
# Usage: scripts/borehole_pipeline.sh <STA>
# GPU steps (cluster, fwd densify, rev densify) run foreground, gated behind ANY other GPU user.
# All stages are idempotent (skip if output already present). Exits nonzero on a hard failure
# so the batch driver can skip the station and continue.
cd /home/jovyan/tremorferometry
STA=$(echo "$1" | tr '[:lower:]' '[:upper:]')
sta=$(echo "$STA" | tr '[:upper:]' '[:lower:]')
PY=/home/jovyan/envs/tremorferometry/bin/python
STATUS=logs/rollout_status.log
mkdir -p logs
log(){ echo "$(date +%H:%M) $*" >> "$STATUS"; echo "$(date +%H:%M) $*"; }
gpu_env(){
  export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
  export CUDA_PATH=/opt/conda/targets/x86_64-linux
  export PYTHONPATH=src
}
wait_gpu_free(){
  # block until no cluster/densify process is running AND GPU mem < 2 GB
  while true; do
    busy=$(ps -eo pid,cmd | grep -E 'densify_gnw_gpu|discover_gpu' | grep -v grep)
    mem=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1)
    [ -z "$busy" ] && [ "${mem:-99999}" -lt 2000 ] && break
    sleep 120
  done
}

DISC=data/${sta}_disc_p70_2010_2026_m3
CAND=data/${sta}_cand_filtered.parquet
PICK=data/family_picker_${sta}_p70_2010_2026_m3.csv
SEL=${DISC}_sel300.summary.csv
NPZ40=${DISC}_40hz.npz
NPZ40R=${DISC}_40hz_rev.npz
TAG=${STA}p90f40

[ -f "$CAND" ] || { log "$STA ABORT: missing $CAND"; exit 2; }

# 1. CLUSTER (GPU ~1 min)
if [ ! -f "${DISC}.summary.csv" ]; then
  wait_gpu_free; gpu_env
  log "$STA cluster start"
  $PY scripts/discover_gpu.py --station "$STA" --candidates "$CAND" --out "${DISC}.npz" \
     --fs 100 --cc-threshold 0.80 --min-family-members 3 --min-years 1 --workers 12 \
     > logs/cluster_${sta}.log 2>&1 || { log "$STA ABORT: cluster failed"; exit 3; }
  nf=$(wc -l < "${DISC}.summary.csv"); log "$STA cluster done: $((nf-1)) families"
fi

# 2. APP B (CPU) -> picker csv (keep pred==LFE label downstream)
if [ ! -f "$PICK" ]; then
  log "$STA appB start"
  PYTHONPATH=src $PY lfe_features/score_family_stacks_picker.py --net PB --sta "$STA" \
     --members "${DISC}.members.parquet" --summary "${DISC}.summary.csv" --out "$PICK" --workers 12 \
     > logs/appb_${sta}.log 2>&1 || { log "$STA ABORT: appB failed"; exit 4; }
  nl=$($PY -c "import pandas as pd;print((pd.read_csv('$PICK').pred=='LFE').sum())" 2>/dev/null)
  log "$STA appB done: ${nl:-?} LFE-label families"
fi

# 3. SELECT (CPU) coverage-balanced ~300
if [ ! -f "$SEL" ]; then
  $PY scripts/select_families_coverage.py "$DISC" "$PICK" 300 > logs/select_${sta}.log 2>&1 \
     || { log "$STA ABORT: select failed"; exit 5; }
  ns=$(wc -l < "$SEL"); log "$STA select done: $((ns-1)) families"
fi

# 4. RESAMPLE 100->40 Hz + time-reversed (CPU)
if [ ! -f "$NPZ40R" ]; then
  $PY scripts/make_40hz_npz.py "$DISC" > logs/resample_${sta}.log 2>&1 \
     || { log "$STA ABORT: resample failed"; exit 6; }
  log "$STA resample 40hz done"
fi

# 5. FORWARD densify (GPU ~1.5-2h)
if ! grep -q 'ALL YEARS DONE' logs/densify_${sta}p90f40.log 2>/dev/null; then
  wait_gpu_free; gpu_env
  log "$STA fwd densify start"
  $PY scripts/densify_gnw_gpu.py --templates-npz "$NPZ40" --summary-csv "$SEL" --min-snr 0 \
     --network PB --station "$STA" --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 \
     --top-n 0 --despike-mad 8 --workers 8 --start-year 2010 --end-year 2026 \
     --out-dir data --out-prefix mf_${sta}p90f40_ > logs/densify_${sta}p90f40.log 2>&1
  grep -q 'ALL YEARS DONE' logs/densify_${sta}p90f40.log || { log "$STA ABORT: fwd densify incomplete"; exit 7; }
  log "$STA fwd densify done"
fi

# 6. REVERSED densify (GPU) — noise-match floor
if ! grep -q 'ALL YEARS DONE' logs/densify_${sta}p90f40rev.log 2>/dev/null; then
  wait_gpu_free; gpu_env
  log "$STA rev densify start"
  $PY scripts/densify_gnw_gpu.py --templates-npz "$NPZ40R" --summary-csv "$SEL" --min-snr 0 \
     --network PB --station "$STA" --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 \
     --top-n 0 --despike-mad 8 --workers 8 --start-year 2010 --end-year 2026 \
     --out-dir data --out-prefix mf_${sta}p90f40rev_ > logs/densify_${sta}p90f40rev.log 2>&1
  grep -q 'ALL YEARS DONE' logs/densify_${sta}p90f40rev.log || { log "$STA ABORT: rev densify incomplete"; exit 8; }
  log "$STA rev densify done"
fi

# 7. STACKS fwd + rev (CPU, 3-comp Z/H1/H2)
if [ ! -f data/long_window_daily_${TAG}_Z.npz ]; then
  log "$STA fwd stacks start"
  $PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_${sta}p90f40_*.csv" \
     --network PB --station "$STA" --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 \
     --despike-mad 8 --workers 10 --out-prefix data/long_window_daily_${TAG} \
     > logs/stack_${sta}f40.log 2>&1 || { log "$STA ABORT: fwd stacks failed"; exit 9; }
  log "$STA fwd stacks done"
fi
if [ ! -f data/long_window_daily_${TAG}rev_Z.npz ]; then
  log "$STA rev stacks start"
  $PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_${sta}p90f40rev_*.csv" \
     --network PB --station "$STA" --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 \
     --despike-mad 8 --workers 10 --out-prefix data/long_window_daily_${TAG}rev \
     > logs/stack_${sta}f40rev.log 2>&1 || { log "$STA ABORT: rev stacks failed"; exit 10; }
  log "$STA rev stacks done"
fi

# 8. dv/v Z + H2 (CPU) — H1 not needed by finalizer (anti-causal at boreholes)
for C in Z H2; do
  OUT=data/daily_dvv_${TAG}_${C}_2to4.csv
  if [ ! -f "$OUT" ]; then
    $PY scripts/dvv_roll30cal.py --station "$STA" --npz data/long_window_daily_${TAG}_${C}.npz \
       --window 2 4 --origin-anchor --workers 10 --out "$OUT" > logs/dvv_${sta}_${C}.log 2>&1 \
       || { log "$STA ABORT: dvv $C failed"; exit 11; }
  fi
done
log "$STA dvv done (Z,H2)"

# 9. FINALIZE (CPU): Z certification + free H2 gate + gated dv/v + figure + summary json
$PY scripts/finalize_3comp_dvv.py "$STA" "$TAG" > logs/finalize_${sta}.log 2>&1 \
   || { log "$STA ABORT: finalize failed"; exit 12; }
summ=$($PY -c "import json;d=json.load(open('data/${sta}_3comp_summary.json'));print('Z',d['z_certified'],'fam std',d.get('z_std_pct'),'% | H2',d['h2_pass'],'fam real',d['h2_real'])" 2>/dev/null)
log "$STA FINALIZE done: $summ"
