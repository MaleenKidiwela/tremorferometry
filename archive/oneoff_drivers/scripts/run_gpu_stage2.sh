#!/bin/bash
# Wrapper to run discover_gpu.py with correct CUDA_PATH
export CUDA_PATH=/opt/conda/targets/x86_64-linux
export PYTHONPATH=/home/jovyan/tremorferometry/src
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYBIN="/home/jovyan/envs/tremorferometry/bin/python"
$PYBIN scripts/discover_gpu.py "$@"
