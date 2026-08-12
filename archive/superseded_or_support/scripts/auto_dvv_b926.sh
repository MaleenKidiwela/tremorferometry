#!/bin/bash
cd /home/jovyan/tremorferometry
until ls data/long_window_daily_B926p90f40_H2.npz >/dev/null 2>&1; do sleep 60; done
echo "=== stacks landed, running 3-comp dv/v $(date) ==="
for C in Z H1 H2; do
  /home/jovyan/envs/tremorferometry/bin/python scripts/dvv_roll30cal.py \
    --station B926 --npz data/long_window_daily_B926p90f40_${C}.npz \
    --window 2 4 --origin-anchor --workers 10 \
    --out data/daily_dvv_B926p90f40_${C}_2to4.csv
  echo "=== dv/v $C done $(date) ==="
done
echo "=== B926 3-COMPONENT dv/v COMPLETE ==="
