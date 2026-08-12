#!/bin/bash
# Concurrent-densify variant of borehole_pipeline.sh (for B004 onward): runs FORWARD and REVERSED
# densify SIMULTANEOUSLY (GPU sits ~0% idle during a single densify; both ~15-19GB fit the 46GB L40S,
# ~halving per-station wall-clock). Same rigor: full reversed densify kept. Serial fallback on failure.
# Usage: scripts/borehole_pipeline_conc.sh <STA>
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
  while true; do
    busy=$(ps -eo pid,cmd | grep -E 'densify_gnw_gpu|discover_gpu' | grep -v grep)
    mem=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1)
    [ -z "$busy" ] && [ "${mem:-99999}" -lt 2000 ] && break
    sleep 120
  done
}
densify_run(){ # $1=templates-npz $2=out-prefix $3=logfile
  $PY scripts/densify_gnw_gpu.py --templates-npz "$1" --summary-csv "$SEL" --min-snr 0 \
     --network PB --station "$STA" --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 \
     --top-n 0 --despike-mad 8 --workers 8 --start-year 2010 --end-year 2026 \
     --out-dir data --out-prefix "$2" > "$3" 2>&1
}
verify_densify(){ # $1=logfile $2=out-prefix-tag(fwd/rev)  -> logs PASS/FAIL
  local lg="$1" tag="$2"
  local tc; tc=$(grep -oE '[0-9]+ templates,' "$lg" | head -1 | grep -oE '[0-9]+')
  local want; want=$(($(wc -l < "$SEL")-1))
  local hrs
  local mf; mf=$(ls data/mf_${sta}p90f40$([ "$tag" = rev ] && echo rev)_20*.csv 2>/dev/null | head -1)
  hrs=$($PY -c "import pandas as pd,sys; d=pd.read_csv('$mf',usecols=['time'],nrows=400000); h=pd.to_datetime(d.time).dt.hour; print(int(h.min()),int(h.max()))" 2>/dev/null)
  if [ "${tc:-0}" = "$want" ] && [ "$hrs" = "0 23" ]; then
    log "$STA $tag densify VERIFY PASS: $tc templates (=sel), mf hours $hrs"
  else
    log "$STA $tag densify VERIFY WARN: templates=$tc want=$want hours=[$hrs]"
  fi
}

DISC=data/${sta}_disc_p70_2010_2026_m3
CAND=data/${sta}_cand_filtered.parquet
PICK=data/family_picker_${sta}_p70_2010_2026_m3.csv
SEL=${DISC}_sel300.summary.csv
NPZ40=${DISC}_40hz.npz
NPZ40R=${DISC}_40hz_rev.npz
TAG=${STA}p90f40
FLOG=logs/densify_${sta}p90f40.log
RLOG=logs/densify_${sta}p90f40rev.log

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

# 2. APP B (CPU) -> picker csv
if [ ! -f "$PICK" ]; then
  log "$STA appB start"
  PYTHONPATH=src $PY lfe_features/score_family_stacks_picker.py --net PB --sta "$STA" \
     --members "${DISC}.members.parquet" --summary "${DISC}.summary.csv" --out "$PICK" --workers 12 \
     > logs/appb_${sta}.log 2>&1 || { log "$STA ABORT: appB failed"; exit 4; }
  nl=$($PY -c "import pandas as pd;print((pd.read_csv('$PICK').pred=='LFE').sum())" 2>/dev/null)
  log "$STA appB done: ${nl:-?} LFE-label families"
fi

# 3. SELECT (CPU)
if [ ! -f "$SEL" ]; then
  $PY scripts/select_families_coverage.py "$DISC" "$PICK" 300 > logs/select_${sta}.log 2>&1 \
     || { log "$STA ABORT: select failed"; exit 5; }
  ns=$(wc -l < "$SEL"); log "$STA select done: $((ns-1)) families"
fi

# 4. RESAMPLE 100->40 Hz + reversed (CPU)
if [ ! -f "$NPZ40R" ]; then
  $PY scripts/make_40hz_npz.py "$DISC" > logs/resample_${sta}.log 2>&1 \
     || { log "$STA ABORT: resample failed"; exit 6; }
  log "$STA resample 40hz done"
fi

# 5+6. CONCURRENT forward + reversed densify (GPU) -----------------------------
fdone(){ grep -q 'ALL YEARS DONE' "$FLOG" 2>/dev/null; }
rdone(){ grep -q 'ALL YEARS DONE' "$RLOG" 2>/dev/null; }
if ! fdone || ! rdone; then
  wait_gpu_free; gpu_env
  log "$STA CONCURRENT fwd+rev densify start (2x workers 8 = 16 CPU)"
  fdone || { densify_run "$NPZ40"  "mf_${sta}p90f40_"    "$FLOG" & FPID=$!; }
  rdone || { densify_run "$NPZ40R" "mf_${sta}p90f40rev_" "$RLOG" & RPID=$!; }
  [ -n "${FPID:-}" ] && wait "$FPID"
  [ -n "${RPID:-}" ] && wait "$RPID"
  # serial fallback for whichever did not complete (e.g. GPU OOM under concurrency)
  if ! fdone; then
    log "$STA fwd densify concurrent INCOMPLETE -> serial retry"
    wait_gpu_free; gpu_env; densify_run "$NPZ40" "mf_${sta}p90f40_" "$FLOG"
    fdone || { log "$STA ABORT: fwd densify failed twice"; exit 7; }
  fi
  if ! rdone; then
    log "$STA rev densify concurrent INCOMPLETE -> serial retry"
    wait_gpu_free; gpu_env; densify_run "$NPZ40R" "mf_${sta}p90f40rev_" "$RLOG"
    rdone || { log "$STA ABORT: rev densify failed twice"; exit 8; }
  fi
  log "$STA fwd+rev densify done"
  verify_densify "$FLOG" fwd
  verify_densify "$RLOG" rev
fi

# 7. STACKS fwd + rev (CPU)
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

# 8. dv/v Z + H2 (CPU)
for C in Z H2; do
  OUT=data/daily_dvv_${TAG}_${C}_2to4.csv
  if [ ! -f "$OUT" ]; then
    $PY scripts/dvv_roll30cal.py --station "$STA" --npz data/long_window_daily_${TAG}_${C}.npz \
       --window 2 4 --origin-anchor --workers 10 --out "$OUT" > logs/dvv_${sta}_${C}.log 2>&1 \
       || { log "$STA ABORT: dvv $C failed"; exit 11; }
  fi
done
log "$STA dvv done (Z,H2)"

# 9. FINALIZE (CPU)
$PY scripts/finalize_3comp_dvv.py "$STA" "$TAG" > logs/finalize_${sta}.log 2>&1 \
   || { log "$STA ABORT: finalize failed"; exit 12; }
summ=$($PY -c "import json;d=json.load(open('data/${sta}_3comp_summary.json'));print('Z',d['z_certified'],'fam std',d.get('z_std_pct'),'% | H2',d['h2_pass'],'fam real',d['h2_real'])" 2>/dev/null)
log "$STA FINALIZE done: $summ"
