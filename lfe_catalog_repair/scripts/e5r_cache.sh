#!/bin/bash
# E5 step 5 -- 20-Hz 1-8 Hz 3-C cache of the 2012 traces (scripts/e4_cache20.py, unchanged), one process per station.
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
for ns in "$@"; do
  ( python scripts/e4_cache20.py --stations $ns --year 2012 > logs/e5_cache_${ns#*.}.log 2>&1; echo "[$(date -u +%FT%TZ)] cache $ns: $(tail -1 logs/e5_cache_${ns#*.}.log)" ) &
done
wait
echo "[$(date -u +%FT%TZ)] cache batch done: $*"
