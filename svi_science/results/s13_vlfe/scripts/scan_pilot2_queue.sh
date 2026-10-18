#!/bin/bash
# VLFE scan pilot 2 queue (PLAN_SCAN.md Amendment S7). Idempotent (per-run DONE / cand2.parquet). Marker SCAN2_QUEUE_DONE.
#   nohup bash scripts/scan_pilot2_queue.sh > logs/scan2_queue.out 2>&1 &     log logs/scan2_queue.log, PID logs/scan2_queue.pid
# GPU path (sequential): cleaned detector runs real_c (+50 polarity draws), ctrlA_s{3,5,7,11,13}_c, inj2_c, inj3_c.
# CPU path (5 workers, parallel with the GPU): injected band-pass for inj3 (near) and for the references (inj2, inj3).
# Then: decay summary (needs STACK_scan2dec_svi_2011-08_DONE and _2012-08_DONE), columns (8 workers), validate2 if present.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe
cd $D
export OMP_NUM_THREADS=1
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=$D/logs/scan2_queue.log
echo $$ > logs/scan2_queue.pid
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $L; }
REFS=UW.GNW,UW.DOSE,CC.OBSR,UW.LEBA,UW.RATT,UW.ELW,TA.D04D,UW.WISH,UW.SP2,TA.B05D,US.NLWA,TA.D03D
log "queue start"
( [ -f data/scan/bp_2011inj3/UW.LRIV.npy ] || python scripts/scan_bp.py $D/data/scan/bp_2011inj3 --inject $D/data/scan/inject3_2011.parquet --workers 5 >> $L 2>&1
  python scripts/scan_bp_ref.py $D/data/scan/bp2_ref_inj2 --stations $REFS --inject $D/data/scan/inject_2011.parquet --workers 5 >> $L 2>&1
  python scripts/scan_bp_ref.py $D/data/scan/bp2_ref_inj3 --stations $REFS --inject $D/data/scan/inject3_2011.parquet --workers 5 >> $L 2>&1
  log "cpu bp done" ) &
CPU=$!
det(){ # run bp ndraw shift
  [ -f data/scan/runs/$1/DONE ] || python scripts/scan_detect2.py $1 --bp $2 --ndraw $3 --shift $4 >> $L 2>&1; log "detect2 $1 rc $?"; }
det real_c data/scan/bp_2011g 50 0
for s in 3 5 7 11 13; do det ctrlA_s${s}_c data/scan/bp_2011g 0 $s; done
det inj2_c data/scan/bp_2011inj2 0 0
wait $CPU
det inj3_c data/scan/bp_2011inj3 0 0
# decay summary (stacks launched separately; wait for their markers without polling the GPU path)
until [ -f STACK_scan2dec_svi_2011-08_DONE ] && [ -f STACK_scan2dec_svi_2012-08_DONE ]; do sleep 30; done
python scripts/scan_decay.py summary >> $L 2>&1; log "decay summary rc $?"
col(){ # run bp refbp [inj]
  [ -s data/scan/runs/$1/cand2.parquet ] || python scripts/scan_columns2.py $1 --bp $2 --refbp $3 ${4:+--inj $4} --workers 8 >> $L 2>&1; log "columns2 $1 rc $?"; }
col real_c data/scan/bp_2011g data/scan/bp2_ref
for s in 3 5 7 11 13; do col ctrlA_s${s}_c data/scan/bp_2011g data/scan/bp2_ref; done
col inj2_c data/scan/bp_2011inj2 data/scan/bp2_ref_inj2 data/scan/inject_2011.parquet
col inj3_c data/scan/bp_2011inj3 data/scan/bp2_ref_inj3 data/scan/inject3_2011.parquet
[ -f scripts/scan_validate2.py ] && { python scripts/scan_validate2.py > logs/scan2_validate.out 2>&1; log "validate2 rc $?"; }
log "queue done"
touch $D/SCAN2_QUEUE_DONE
