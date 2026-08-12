#!/bin/bash
cd /home/jovyan/tremorferometry
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src
PY=/home/jovyan/envs/tremorferometry/bin/python

# wait for forward dv/v (downstream chain) AND reversed densify to finish
until grep -q 'B011 3-COMPONENT dv/v COMPLETE' logs/chain_b011_downstream.log 2>/dev/null; do sleep 120; done
echo "=== B011 fwd dv/v done $(date) ==="
until grep -q 'ALL YEARS DONE' logs/densify_b011p90f40rev.log 2>/dev/null; do sleep 120; done
echo "=== B011 reversed densify done, building reversed 3-ch stacks $(date) ==="

$PY scripts/build_long_window_3comp.py \
  --mf-csv-glob 'data/mf_b011p90f40rev_*.csv' --network PB --station B011 --fs 40 \
  --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 --workers 10 \
  --out-prefix data/long_window_daily_B011p90f40rev
echo "=== B011 reversed stacks done, finalizing 3-comp $(date) ==="

$PY scripts/finalize_3comp_dvv.py B011 B011p90f40
echo "=== B011 3-COMP FINALIZE COMPLETE $(date) ==="
