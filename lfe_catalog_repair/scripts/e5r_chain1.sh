#!/bin/bash
# E5 step 5: offsets (run 2, consistency rule) -> split rebuilt on the 2012 offsets -> from-scratch templates (15 stations)
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=8
set -e
python scripts/e5r_offsets.py > logs/e5r_offsets_run2.log 2>&1
python scripts/e5r_groundtruth.py --offsets results/e5r_family_offsets_2012.csv > logs/e5r_groundtruth_final.log 2>&1
export OMP_NUM_THREADS=1
python scripts/e5r_templates.py --workers 15 > logs/e5r_templates.log 2>&1
echo "[$(date -u +%FT%TZ)] chain1 DONE"
