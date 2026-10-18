#!/bin/bash
# Scan extension v2 for one pilot year (Amendments R, S). Marker SCANV2_${Y}_DONE.
Y=$1   # usage: v2_year.sh YEAR (pilot dir data/pilotYEAR must exist)
# FAIL lines -> logs/q_main.log. Resume: re-run (STEP markers in data/pilot$Y/).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 SCAN_YEAR=$Y NINJ=400 PILOT=$Y
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') v2 $*" >> $L; }
st(){ local m=$1; shift; [ -f data/pilot$Y/STEP_V2_${m}_DONE ] && return 0; log "start $m"
  "$@" >> logs/scanv2_$Y.out 2>&1 || { log "FAIL $m (see logs/scanv2_$Y.out)"; exit 1; }; touch data/pilot$Y/STEP_V2_${m}_DONE; log "end $m"; }
SD=$D/data/scan_$Y
st setup python scripts/cat_scanv2.py setup $Y
[ -f $D/catalogue/data/inj/sources_v1.npz ] || st sources python scripts/cat_scanv2.py sources
export SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1
st detect bash -c "cd $D && python scripts/scan_detect.py real_v2 --bp $SD/bp_${Y}v2 --ndraw 0 --shift 0 && test -f $SD/runs/real_v2/DONE"
export SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2 PILOT_SFX=_v2
st battery python scripts/cat_pilot.py battery
st inj python scripts/cat_pilot.py inj
st summary python scripts/cat_pilot.py summary
st R python scripts/cat_r.py
TZ=America/Los_Angeles date > SCANV2_${Y}_DONE; log "scan v2 $Y done"
