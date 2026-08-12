#!/bin/bash
# OPTIMIZED 2017 gap fill for B926 + B011 (2017 was truncated at day 59).
# Only densifies the RELIABLE (certified, fwd/rev>1.5) families, FORWARD only, for 2017 -- the ringing
# families and the reversed run are unused for dv/v and the certification is already frozen. Then splices
# 2017 into the forward stacks, recomputes dv/v for the reliable families, and re-plots with frozen cert.
cd /home/jovyan/tremorferometry
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src
PY=/home/jovyan/envs/tremorferometry/bin/python
wait_gpu_free () { while ps -eo cmd | grep -E 'densify_gnw_gpu|discover_gpu' | grep -qv grep; do sleep 60; done; }
log () { echo "$(date +%H:%M) [2017fix] $*" | tee -a logs/rollout_status.log; }

# gate GPU work behind full batch completion (no 3-densify OOM race)
until ! ps -eo cmd | grep -E 'run_batch_conc|borehole_pipeline' | grep -qv grep; do sleep 120; done
log "batch complete -> 2017 gap fix (reliable families, forward-only)"

for S in B926 B011; do
  s=$(echo $S | tr A-Z a-z); TAG=${S}p90f40
  NPZ=data/${s}_disc_p70_2010_2026_m3_40hz.npz
  SEL=data/${s}_disc_p70_2010_2026_m3_sel300.summary.csv
  RSEL=data/${s}_reliable2017.summary.csv
  # reliable-only summary = sel300 families with fwd/rev ratio>1.5
  $PY - "$SEL" "data/${s}_fwd_vs_rev_coda.csv" "$RSEL" <<'PY'
import sys, pandas as pd
sel, fvr, out = sys.argv[1:4]
s = pd.read_csv(sel); cert = set(pd.read_csv(fvr).query("ratio>1.5").fam)
s[s.family_id.isin(cert)].to_csv(out, index=False)
print(f"reliable summary: {s.family_id.isin(cert).sum()} families -> {out}")
PY
  # 1. densify 2017 FORWARD only, reliable families
  wait_gpu_free
  log "$S densify 2017 fwd (reliable) start"
  $PY scripts/densify_gnw_gpu.py --templates-npz $NPZ --summary-csv $RSEL --min-snr 0 \
    --network PB --station $S --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 --top-n 0 \
    --despike-mad 8 --workers 8 --start-year 2017 --end-year 2017 \
    --out-dir data --out-prefix mf_${s}p90f40_ >> logs/densify_${s}p90f40.log 2>&1
  log "$S densify 2017 fwd done ($(wc -l < data/mf_${s}p90f40_2017.csv) rows)"
  # 2. stack ONLY 2017 (fwd) then splice into the existing full forward npz
  $PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_${s}p90f40_2017.csv" \
    --network PB --station $S --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 \
    --workers 10 --out-prefix data/long_window_daily_${TAG}_2017only
  for C in Z H1 H2; do
    $PY scripts/merge_year_into_npz.py data/long_window_daily_${TAG}_${C}.npz \
      data/long_window_daily_${TAG}_2017only_${C}.npz 2017
  done
  # 3. recompute dv/v (Z,H2) on the spliced forward stacks + re-plot with frozen certification
  for C in Z H2; do
    $PY scripts/dvv_roll30cal.py --station $S --npz data/long_window_daily_${TAG}_${C}.npz \
      --window 2 4 --origin-anchor --workers 10 --out data/daily_dvv_${TAG}_${C}_2to4.csv
  done
  $PY scripts/replot_gated_dvv.py $S $TAG
  rm -f data/long_window_daily_${TAG}_2017only_*.npz
  log "$S 2017-FILLED (reliable/fwd-only) -> $(tr -d '\n ' < data/${s}_3comp_summary.json)"
done
log "2017 GAP FIX COMPLETE (reliable families only, forward-only, frozen cert)"
