#!/bin/bash
# E4 step 1 — pass 1 (no moveout) then pass 2 (empirical moveout measured in pass 1).
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
set -e
python scripts/e4_template_test.py --tag ring_p1 --n-random 600 --doy-from 222 --doy-to 255
python scripts/e4_template_test.py --tag ring_p2 --n-random 600 --doy-from 222 --doy-to 255 \
       --moveout-in results/e4_moveout_ring_p1.csv
echo "[e4 step1 DONE]"
