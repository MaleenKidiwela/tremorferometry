#!/bin/bash
# NV3 (PLAN Amendment NV3): W2b scan-v2 runs per year (polBW = real + 4 polarity draws; shA3W/7W/11W station-day shifts).
# One GPU process at a time; waits while the GPU has < 22 GB free (shares with the catalogue queues). Marker NV3_SCAN_DONE.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; N=$D/newvlfe; cd $D
export OMP_NUM_THREADS=1 SCAN_GAIN_FALLBACK=1 NV_W2B=1
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $N/logs/queue.log; }
for Y in 2016 2019 2011 2013; do
  SD=$D/data/scan_$Y; [ $Y = 2011 ] && SD=$D/data/scan
  export SCAN_YEAR=$Y SCAN_GRID=$SD/templates_grid_v2.npz
  for spec in "polBW:4:0" "shA3W:0:3" "shA7W:0:7" "shA11W:0:11"; do
    IFS=: read nm nd sh <<< "$spec"; O=$N/data/runs/$Y/$nm
    [ -f $O/DONE ] && continue
    while :; do fr=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1); [ "$fr" -ge 22000 ] && break; sleep 60; done
    mkdir -p $N/data/runs/$Y; log "NV3 scan start $Y $nm"
    python $N/scripts/nv_detect.py $nm --bp $SD/bp_${Y}v2 --ndraw $nd --shift $sh --out $O >> $N/logs/nv_detect_$Y.out 2>&1 || { log "FAIL NV3 scan $Y $nm"; exit 1; }
    log "NV3 scan end $Y $nm"
  done
done
TZ=America/Los_Angeles date > $N/NV3_SCAN_DONE; log "NV3 scans done"
