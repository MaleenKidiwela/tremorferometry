#!/bin/bash
# serialize the single GPU across the D3 stations
cd /home/jovyan/tremorferometry
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
export CUDA_PATH=/opt/conda/targets/x86_64-linux
PY=/home/jovyan/envs/tremorferometry/bin/python
$PY lfe_catalog_repair/scripts/d3_redetect.py "$@"
