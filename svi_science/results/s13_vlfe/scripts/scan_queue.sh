#!/bin/bash
# VLFE scan pilot year queue (PLAN_SCAN.md §4-6). GPU path: real (+50 polarity draws), 5 station-day-shift controls,
# injection. CPU path (<= 7 workers): injection list + injected band-pass, then columns. Marker SCAN_QUEUE_DONE.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe
cd $D
export OMP_NUM_THREADS=1
L=$D/logs/scan_queue.log
log(){ echo "$(date '+%F %T') $*" >> $L; }
log "queue start"
( # CPU path: injection data
  [ -s data/scan/inject_2011.parquet ] || python scripts/scan_inject.py >> $L 2>&1
  python scripts/scan_bp.py $D/data/scan/bp_2011inj --inject $D/data/scan/inject_2011.parquet --workers 6 >> $L 2>&1
  log "inj bp done" ) &
CPU=$!
for r in real; do
  [ -f data/scan/runs/$r/DONE ] || python scripts/scan_detect.py $r --bp data/scan/bp_2011g --ndraw 50 >> $L 2>&1; log "detect $r rc $?"
done
for s in 3 5 7 11 13; do
  r=ctrlA_s$s
  [ -f data/scan/runs/$r/DONE ] || python scripts/scan_detect.py $r --bp data/scan/bp_2011g --ndraw 0 --shift $s >> $L 2>&1; log "detect $r rc $?"
done
wait $CPU
[ -f data/scan/runs/inj/DONE ] || python scripts/scan_detect.py inj --bp data/scan/bp_2011inj --ndraw 0 >> $L 2>&1; log "detect inj rc $?"
for r in real inj ctrlA_s3 ctrlA_s5 ctrlA_s7 ctrlA_s11 ctrlA_s13; do
  bp=data/scan/bp_2011g; [ $r = inj ] && bp=data/scan/bp_2011inj
  [ -s data/scan/runs/$r/cand_cols.parquet ] || python scripts/scan_columns.py $r --bp $bp --workers 7 >> $L 2>&1; log "columns $r rc $?"
done
log "queue done"
touch $D/SCAN_QUEUE_DONE
