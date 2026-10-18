#!/bin/bash
# newvlfe stage 1: scrambled-control scan runs (GPU, 1 process), 2013 and 2018. Marker newvlfe/SCANCTRL_DONE.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; N=$D/newvlfe; cd $D
export OMP_NUM_THREADS=1 SCAN_GAIN_FALLBACK=1
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $N/logs/queue.log; }
for Y in 2013 2018; do
  export SCAN_YEAR=$Y SCAN_GRID=$D/data/scan_$Y/templates_grid_v2.npz
  BP=$D/data/scan_$Y/bp_${Y}v2
  for spec in "polB:4:0" "shA3:0:3" "shA7:0:7" "shA11:0:11"; do
    IFS=: read nm nd sh <<< "$spec"; O=$N/data/runs/$Y/$nm
    [ -f $O/DONE ] && continue
    log "start $Y $nm"
    python $N/scripts/nv_detect.py $nm --bp $BP --ndraw $nd --shift $sh --out $O >> $N/logs/nv_detect_$Y.out 2>&1 || { log "FAIL $Y $nm"; exit 1; }
    log "end $Y $nm"
  done
done
TZ=America/Los_Angeles date > $N/SCANCTRL_DONE; log "scan controls done"
