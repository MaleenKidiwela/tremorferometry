#!/bin/bash
# E2 battery queue (PLAN_CAT §3/§6). Marker: catalogue/E2_DONE
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
export OMP_NUM_THREADS=1
SCAN_YEAR=2011 python scripts/cat_e2.py battery 2011 >> logs/e2.out 2>&1 || exit 1
SCAN_YEAR=2012 python scripts/cat_e2.py battery 2012 >> logs/e2.out 2>&1 || exit 1
python scripts/cat_e2.py summary >> logs/e2.out 2>&1 || exit 1
TZ=America/Los_Angeles date > E2_DONE
