#!/bin/bash
# Scan extension v2 validation on 2013 (Amendment R). Waits for PILOT2018_DONE (CPU cap). Marker SCANV2_2013_DONE.
# FAIL lines -> logs/q_main.log. Resume: re-run (STEP markers in data/pilot2013/).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 SCAN_YEAR=2013 NINJ=400 PILOT=2013
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') v2 $*" >> $L; }
st(){ local m=$1; shift; [ -f data/pilot2013/STEP_V2_${m}_DONE ] && return 0; log "start $m"
  "$@" >> logs/scanv2.out 2>&1 || { log "FAIL $m (see logs/scanv2.out)"; exit 1; }; touch data/pilot2013/STEP_V2_${m}_DONE; log "end $m"; }
until [ -f PILOT2018_DONE ]; do pgrep -f r_queue.sh > /dev/null || { log "FAIL r_queue not running and no PILOT2018_DONE"; exit 1; }; sleep 60; done
SD=$D/data/scan_2013
st setup python scripts/cat_scanv2.py setup 2013
st sources python scripts/cat_scanv2.py sources
export SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1
st detect bash -c "cd $D && python scripts/scan_detect.py real_v2 --bp $SD/bp_2013v2 --ndraw 0 --shift 0 && test -f $SD/runs/real_v2/DONE"
export SCAN_RUN=real_v2 SCAN_BP=$SD/bp_2013v2 PILOT_SFX=_v2
st battery python scripts/cat_pilot.py battery
st inj python scripts/cat_pilot.py inj
st summary python scripts/cat_pilot.py summary
st R python scripts/cat_r.py
TZ=America/Los_Angeles date > SCANV2_2013_DONE; log "scan v2 2013 done"
