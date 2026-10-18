#!/bin/bash
# POST-HOC positive control for R3 (flagged in NOTES/RESULTS): same pipeline on Table-S1 ETS-window detections.
# 4 event workers + aux + agg = 6 processes. Marker R3_ETS_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro/scripts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R3_TAG=ets
python hg_r3.py rows > ../logs/r3_ets_rows.log 2>&1 || exit 1
( python hg_r3.py aux > ../logs/r3_ets_aux.log 2>&1 ) &
( python hg_r3.py agg > ../logs/r3_ets_agg.log 2>&1 ) &
SCAN_YEAR=2011 python hg_r3.py events 2011 > ../logs/r3_ets_events2011.log 2>&1
SCAN_YEAR=2012 python hg_r3.py events 2012 > ../logs/r3_ets_events2012.log 2>&1
wait
touch ../R3_ETS_DONE
