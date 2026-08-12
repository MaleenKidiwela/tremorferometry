#!/bin/bash
cd /home/jovyan/tremorferometry
until grep -q 'B011 FORWARD DENSIFY DONE' logs/chain_b011.log 2>/dev/null || grep -q 'ALL YEARS DONE' logs/densify_b011p90f40.log 2>/dev/null; do sleep 120; done
echo "=== B011 densify done, building 3-ch stacks $(date) ==="
/home/jovyan/envs/tremorferometry/bin/python scripts/build_long_window_3comp.py \
  --mf-csv-glob 'data/mf_b011p90f40_*.csv' --network PB --station B011 --fs 40 \
  --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 --workers 10 \
  --out-prefix data/long_window_daily_B011p90f40
echo "=== B011 stacks done, 3-comp dv/v $(date) ==="
for C in Z H1 H2; do
  /home/jovyan/envs/tremorferometry/bin/python scripts/dvv_roll30cal.py \
    --station B011 --npz data/long_window_daily_B011p90f40_${C}.npz \
    --window 2 4 --origin-anchor --workers 10 --out data/daily_dvv_B011p90f40_${C}_2to4.csv
done
echo "=== B011 3-COMPONENT dv/v COMPLETE $(date) ==="
