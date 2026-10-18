#!/bin/bash
# newvlfe W2b re-run 2018 (coordinator 2026-10-16 ~11:15 am PT): scan v2 with catalogue W2b station-component-month
# exclusions (qc_exclusions.csv, 2018 = UW.LRIV N/E/Z), 1 GPU pass with 4 polarity draws (variant 0 = real, run polBW),
# then batteries realW / polBW (2 workers). Marker newvlfe/W2B2018_DONE; failures -> FAIL line in logs/queue.log.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; N=$D/newvlfe; cd $N
export OMP_NUM_THREADS=1 SCAN_GAIN_FALLBACK=1 SCAN_YEAR=2018 SCAN_GRID=$D/data/scan_2018/templates_grid_v2.npz NV_W2B=1
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $N/logs/queue.log; }
O=$N/data/runs/2018/polBW
if [ ! -f $O/DONE ]; then
  log "W2b start 2018 polBW scan"
  python scripts/nv_detect.py polBW --bp $D/data/scan_2018/bp_2018v2 --ndraw 4 --shift 0 --out $O >> logs/nv_detect_2018W.out 2>&1 || { log "FAIL W2b scan 2018 polBW"; exit 1; }
  log "W2b end 2018 polBW scan"
fi
for R in realW polBW; do
  [ -f data/bat_2018_${R}.parquet ] && continue
  log "W2b battery start 2018 $R"
  python scripts/nv_battery.py 2018 $R --workers 2 >> logs/nv_battery_2018W.out 2>&1 || { log "FAIL W2b battery 2018 $R"; exit 1; }
  log "W2b battery end 2018 $R"
done
TZ=America/Los_Angeles date > $N/W2B2018_DONE; log "W2b 2018 done"
