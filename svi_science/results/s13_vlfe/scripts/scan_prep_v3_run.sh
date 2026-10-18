#!/bin/bash
# VLFE scan pilot: prep v3 for SVI 2011 (PLAN_SCAN.md §2). Marker SCAN_PREP_V3_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1 S13_DISP=$PWD/data/disp1sps_v3
python scripts/scan_prep_v3.py CN.PGC,CN.VGZ,CN.SNB,CN.NLLB,CN.PFB,CN.YOUB,UW.LRIV,CN.LZB 2011-01-01 2011-12-31 --workers 8 > logs/scan_prep_v3.out 2>&1 && touch SCAN_PREP_V3_DONE
