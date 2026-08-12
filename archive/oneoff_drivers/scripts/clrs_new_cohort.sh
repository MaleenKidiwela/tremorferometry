#!/bin/bash
# Densify the 151 NEW (min-years-1) CLRS families separately -> causality; combine with the original 16/48.
# Per-cohort survival tracked (Merlin junk-guard: new-cohort survival <10% => extra yield is junk).
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src
log(){ echo "$(date +%H:%M) [clrs-new] $*" | tee -a logs/clrs_new_cohort.log; }
log "densify 151 new families @40Hz 2014-2026"
$PY scripts/densify_gnw_gpu.py --templates-npz data/clrs_new_families.npz \
  --summary-csv data/clrs_new_families.summary.csv --min-snr 0 --network CN --station CLRS --fs 40 \
  --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 --top-n 0 --despike-mad 8 --workers 8 \
  --start-year 2014 --end-year 2026 --out-dir data --out-prefix mf_clrsnew_ > logs/densify_clrsnew.log 2>&1 || { log "densify FAILED"; exit 1; }
log "stacks (new cohort)"
$PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_clrsnew_*.csv" --network CN --station CLRS \
  --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 --workers 10 \
  --out-prefix data/long_window_daily_CLRSnew > logs/stack_clrsnew.log 2>&1 || { log "stacks FAILED"; exit 2; }
log "coda dv/v + causality (new cohort)"
$PY scripts/dvv_roll30cal.py --station CLRS --npz data/long_window_daily_CLRSnew_Z.npz --window 2 4 \
  --origin-anchor --workers 10 --out data/daily_dvv_CLRSnew_Z_2to4.csv > logs/dvv_clrsnew.log 2>&1 || { log "dvv FAILED"; exit 3; }
$PY scripts/finalize_causality.py CLRSnew CLRSnew > logs/finalize_clrsnew.log 2>&1 || { log "finalize FAILED"; exit 4; }
nc=$($PY -c "import json;print(json.load(open('data/clrsnew_3comp_summary.json'))['z_certified'])" 2>/dev/null)
log "NEW COHORT: ${nc}/151 certified ($(($nc*100/151))% survival) | + original 16/48 = TOTAL $((16+nc)) certified (GATE >=20)"
