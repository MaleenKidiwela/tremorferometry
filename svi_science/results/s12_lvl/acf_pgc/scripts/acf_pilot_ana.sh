#!/bin/bash
# S12c pilot analysis (PLAN §4): waits for PILOT_PRE_DONE, then visibility (3 stations, 2 masks) and the epoch analysis
# (3 stations x 2 bands x 2 masks, 4 at a time). Markers PILOT_VIS_DONE, EPOCH_pilot2_*_DONE, PILOT_ANA_DONE.
# tag pilot2 = after the maxshift 0.5 -> 0.2 s cycle-skip fix (tag pilot = first PGC run, kept).
cd /home/jovyan/tremorferometry
WS=svi_science/results/s12_lvl/acf_pgc
while [ ! -f $WS/PILOT_PRE_DONE ]; do sleep 30; done
if [ ! -f $WS/PILOT_VIS_DONE ]; then
  for M in tremorfree none; do python $WS/scripts/acf02_visibility.py --stations CN.PGC,CN.VGZ,CN.SNB --years 2018-2021 --mask $M --tag pilot; done
  date > $WS/PILOT_VIS_DONE
fi
run() { K=$WS/EPOCH_pilot2_$1_$2_$3_DONE; [ -f $K ] && return
  python $WS/scripts/acf03_epoch.py --station $1 --years 2018-2021 --band $2 --mask $3 --tag pilot2 --ndraw 2000 > $WS/logs/epoch_pilot2_$1_$2_$3.log 2>&1 && date > $K; }
for S in CN.PGC CN.VGZ CN.SNB; do
  for B in b05_2 b1_4; do run $S $B tremorfree & run $S $B none & done
  wait
done
date > $WS/PILOT_ANA_DONE
