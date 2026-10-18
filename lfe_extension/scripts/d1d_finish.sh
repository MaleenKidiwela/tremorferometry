#!/bin/bash
# Track D, D1(d) -- wait for the withheld-network passes, apply the frozen HF veto to each (Track D copy
# of a1_hf_veto.py that also looks in lfe_extension/traces), then run the seasonal stack test.
cd /home/jovyan/tremorferometry/lfe_extension
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=3
MONTHS="2015-02 2014-02 2014-12 2014-08 2015-08 2014-07"
for m in $MONTHS; do
  while [ ! -f results/e5r_detections_d1d_wh_$m.parquet ]; do sleep 20; done
done
echo "[$(date -u +%FT%TZ)] all withheld passes present" >> logs/d1d_withheld.log
run_hf () {
  local m=$1 y=${1%%-*}
  [ -f results/a1_hf_ratio_d1d_wh_$m.parquet ] && return 0
  /opt/conda/bin/python scripts/d1d_hf_veto.py --det results/e5r_detections_d1d_wh_$m.parquet \
    --tag d1d_wh_$m --L 12 --year $y > logs/d1d_hf_$m.log 2>&1
  echo "hf $m exit=$? $(date -u +%FT%TZ)" >> logs/d1d_withheld.log
}
set -- $MONTHS
while [ $# -gt 0 ]; do
  run_hf $1 & p1=$!
  if [ $# -gt 1 ]; then run_hf $2 & p2=$!; wait $p2; fi
  wait $p1
  shift 2 2>/dev/null || shift
done
echo "[$(date -u +%FT%TZ)] HF passes done" >> logs/d1d_withheld.log
/opt/conda/bin/python scripts/d1d_season_stack_test.py --tiers T_core:8.75,T_ext:8.25 \
  --out results/trackD/season_stack_test.csv > logs/d1d_stack_test.log 2>&1
echo "stack_test exit=$? $(date -u +%FT%TZ)" >> logs/d1d_withheld.log
