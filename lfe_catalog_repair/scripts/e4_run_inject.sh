#!/bin/bash
cd /home/jovyan/tremorferometry
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
python lfe_catalog_repair/scripts/e4_inject.py --days 3 --per-day 140 --nstack 20,0
echo "[e4 inject DONE]"
