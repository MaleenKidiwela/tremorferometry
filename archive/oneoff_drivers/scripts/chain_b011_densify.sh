#!/bin/bash
cd /home/jovyan/tremorferometry
until grep -q 'ALL YEARS DONE' logs/densify_b926rev.log 2>/dev/null; do sleep 120; done
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src
/home/jovyan/envs/tremorferometry/bin/python scripts/densify_gnw_gpu.py \
  --templates-npz data/b011_disc_p70_2010_2026_m3_40hz.npz \
  --summary-csv data/b011_disc_p70_2010_2026_m3_sel300.summary.csv \
  --min-snr 0 --network PB --station B011 --fs 40 --fmin 2 --fmax 8 --threshold 0.8 \
  --min-gap-s 6 --top-n 0 --despike-mad 8 --workers 10 \
  --start-year 2010 --end-year 2026 --out-dir data --out-prefix mf_b011p90f40_ \
  >> logs/densify_b011p90f40.log 2>&1
echo "B011 FORWARD DENSIFY DONE"
