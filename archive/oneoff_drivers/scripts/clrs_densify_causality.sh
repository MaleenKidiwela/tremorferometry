#!/bin/bash
# CLRS decisive gate (Merlin step 3): densify the 48 picker families -> 30-day stacks -> coda dv/v ->
# CAUSALITY certify. Forward-only (finalized pipeline; mirror does the reverse job from fwd stacks).
# GATE = >=20 certified AND >=15% survival (fleet inclusion rule). Usage: nohup bash scripts/clrs_densify_causality.sh &
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src
log(){ echo "$(date +%H:%M) [clrs-dc] $*" | tee -a logs/clrs_densify_causality.log; }
TAG=CLRS
log "densify (forward) 48 picker families @40Hz, 2014-2026 (top-n 0 = cap off)"
$PY scripts/densify_gnw_gpu.py --templates-npz data/clrs_picker_families.npz \
  --summary-csv data/clrs_picker_families.summary.csv --min-snr 0 \
  --network CN --station CLRS --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 \
  --top-n 0 --despike-mad 8 --workers 8 --start-year 2014 --end-year 2026 \
  --out-dir data --out-prefix mf_clrsp90f40_ > logs/densify_clrs.log 2>&1 || { log "densify FAILED"; exit 1; }
log "densify done: $(ls data/mf_clrsp90f40_*.csv 2>/dev/null | wc -l) year files"
log "stacks (long-window 3comp)"
$PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_clrsp90f40_*.csv" \
  --network CN --station CLRS --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 \
  --workers 10 --out-prefix data/long_window_daily_${TAG} > logs/stack_clrs.log 2>&1 || { log "stacks FAILED"; exit 2; }
log "coda dv/v (2-4s, origin-anchored)"
$PY scripts/dvv_roll30cal.py --station CLRS --npz data/long_window_daily_${TAG}_Z.npz \
  --window 2 4 --origin-anchor --workers 10 --out data/daily_dvv_${TAG}_Z_2to4.csv > logs/dvv_clrs.log 2>&1 || { log "dvv FAILED"; exit 3; }
log "causality certification"
$PY scripts/finalize_causality.py ${TAG} ${TAG} > logs/finalize_clrs.log 2>&1 || { log "finalize FAILED"; exit 4; }
zc=$($PY -c "import json;d=json.load(open('data/clrs_3comp_summary.json'));print(d['z_certified'])" 2>/dev/null)
surv=$($PY -c "import json;d=json.load(open('data/clrs_3comp_summary.json'));print(f\"{100*d['z_certified']/48:.0f}\")" 2>/dev/null)
log "CLRS CAUSALITY DONE: ${zc}/48 certified (${surv}% survival) — GATE >=20 & >=15%"
