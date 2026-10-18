#!/bin/bash
# R3 events re-run after the pickling fix (4 workers; aux + agg run alongside = 6 processes). Marker R3_EVENTS_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro/scripts
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
SCAN_YEAR=2011 python hg_r3.py events 2011 > ../logs/r3_events2011.log 2>&1 && \
SCAN_YEAR=2012 python hg_r3.py events 2012 > ../logs/r3_events2012.log 2>&1 && touch ../R3_EVENTS_DONE
