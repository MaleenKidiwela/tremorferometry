#!/bin/bash
# E3/E4 queue (PLAN_CAT §8). Marker: catalogue/E3_DONE. Resume: re-run this script (finished MF jobs are skipped).
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
export OMP_NUM_THREADS=1
python scripts/cat_e3.py mf >> logs/e3.out 2>&1 || exit 1
python scripts/cat_e3.py detect >> logs/e3.out 2>&1 || exit 1
SCAN_YEAR=2011 python scripts/cat_e3.py battery 2011 >> logs/e3.out 2>&1 || exit 1
SCAN_YEAR=2012 python scripts/cat_e3.py battery 2012 >> logs/e3.out 2>&1 || exit 1
python scripts/cat_e3.py summary >> logs/e3.out 2>&1 || exit 1
TZ=America/Los_Angeles date > E3_DONE
