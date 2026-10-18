#!/bin/bash
# real run then the +1800 s null control (same code, same table)
cd /home/jovyan/tremorferometry/lfe_oregon
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R=results/chain_paper
python scripts/xs_event_locate.py --out $R/event_locs --workers 12 > logs/event_locs.log 2>&1 || { echo fail > $R/event_locs/EVLOC_FAILED; exit 1; }
python scripts/xs_event_locate.py --out $R/event_locs_null --shift 1800 --workers 12 > logs/event_locs_null.log 2>&1 || { echo fail > $R/event_locs_null/EVLOC_FAILED; exit 1; }
echo ok > $R/EVLOC_QUEUE_DONE
