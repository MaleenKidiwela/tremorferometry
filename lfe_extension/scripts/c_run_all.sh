#!/bin/bash
# Track C — start (or restart) the production watcher.  Resumable and idempotent: it picks up at the
# oldest runnable month that is not yet in results/trackC/months_done.csv and never re-runs a finished one.
#   nohup bash scripts/c_run_all.sh > /dev/null 2>&1 &
#   tail -f logs/trackC_watch.log
#   kill $(cat logs/trackC_watch.pid)          # stop; safe at any point
cd /home/jovyan/tremorferometry/lfe_extension || exit 1
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
export CUDA_PATH=/opt/conda/targets/x86_64-linux
export OMP_NUM_THREADS=4
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
mkdir -p logs results/trackC
exec /opt/conda/bin/python scripts/c_watch.py "$@" >> logs/trackC_watch.log 2>&1
