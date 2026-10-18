#!/bin/bash
# Amendment R3 queue: events 2011, events 2012 (4 workers each, sequential; with aux + agg = 6 processes), then aux and agg (single process each, run in
# parallel). Markers R3_EVENTS_DONE R3_AUX_DONE R3_AGG_DONE R3_COMPUTE_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro/scripts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
M=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro
( python hg_r3.py aux > ../logs/r3_aux.log 2>&1 && touch $M/R3_AUX_DONE ) &
( python hg_r3.py agg > ../logs/r3_agg.log 2>&1 && touch $M/R3_AGG_DONE ) &
SCAN_YEAR=2011 python hg_r3.py events 2011 > ../logs/r3_events2011.log 2>&1
SCAN_YEAR=2012 python hg_r3.py events 2012 > ../logs/r3_events2012.log 2>&1 && touch $M/R3_EVENTS_DONE
wait
[ -f $M/R3_EVENTS_DONE ] && [ -f $M/R3_AUX_DONE ] && [ -f $M/R3_AGG_DONE ] && touch $M/R3_COMPUTE_DONE
