#!/bin/bash
# E3 step 2 config C — synthetic S-only geometry test, 4 noise levels in parallel.
cd /home/jovyan/tremorferometry
L=lfe_catalog_repair/logs
G="python lfe_catalog_repair/scripts/e3_geom_synth.py --reps 2 --min-n 3 --max-n 6"
$G --sigma 0.05  --tag s005   > $L/e3_geom_s005.log   2>&1 &
$G --sigma 0.3   --tag s030   > $L/e3_geom_s030.log   2>&1 &
$G --sigma 0.5   --tag s050   > $L/e3_geom_s050.log   2>&1 &
$G --sigma heavy --tag sheavy > $L/e3_geom_sheavy.log 2>&1 &
wait
echo "[e3 geom ALL DONE]"
