#!/bin/bash
# VLFE scan pilot 3 (PLAN_SCAN.md Amendment S8): SVI calendar 2012, ratio veto R1. Idempotent (per-step outputs / DONE files).
#   nohup bash scripts/scan3_queue.sh > logs/scan3_queue.out 2>&1 &      log logs/scan3_queue.log, PID logs/scan3_queue.pid
# <= 5 CPU workers (track B holds the other 3), OMP_NUM_THREADS=1. Marker SCAN3_QUEUE_DONE (validate3 runs separately).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe
cd $D
export OMP_NUM_THREADS=1 SCAN_YEAR=2012 S13_DISP=$D/data/disp1sps_v3 NWORK=5
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
SD=$D/data/scan_2012
mkdir -p $SD
L=$D/logs/scan3_queue.log
echo $$ > logs/scan3_queue.pid
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $L; }
PIL=CN.PGC,CN.VGZ,CN.SNB,CN.NLLB,CN.PFB,CN.YOUB,UW.LRIV,CN.LZB
REG=UW.GNW,UW.DOSE,CC.OBSR,UW.LEBA,UW.RATT,UW.ELW,TA.D04D,UW.WISH,UW.SP2,TA.B05D,US.NLWA,TA.D03D
FAR=IU.COR,US.NEW,US.HAWA,UW.OMAK,TA.I02D
REFS=$REG,$FAR
W=5
log "queue start"
# shared, year-independent inputs (same grid templates and reference templates as 2011; same ref table)
[ -e $SD/templates_grid.npz ] || ln -s $D/data/scan/templates_grid.npz $SD/templates_grid.npz
[ -e $SD/templates_ref.npz ] || ln -s $D/data/scan/templates_ref.npz $SD/templates_ref.npz
[ -e $SD/hfenv ] || ln -s $D/data/scan/hfenv $SD/hfenv
step(){ # name command...
  local n=$1; shift
  [ -f $SD/STEP_${n}_DONE ] && return 0
  log "start $n"; "$@" >> $L 2>&1; local rc=$?; log "end $n rc $rc"
  [ $rc -eq 0 ] && touch $SD/STEP_${n}_DONE
  return $rc
}
step prep_pilot python scripts/scan_prep_v3.py $PIL 2012-01-01 2012-12-31 --workers $W
step prep_reg python scripts/scan_prep_ref.py $REG 2012-01-01 2012-12-31 --workers $W
until [ -f REF2012_DOWNLOAD_DONE ]; do sleep 60; done
step prep_far python scripts/scan_prep_ref.py $FAR 2012-01-01 2012-12-31 --workers $W
step hfenv python scripts/scan_hfenv.py CN.PGC,CN.VGZ,CN.SNB,CN.NLLB,CN.PFB,CN.YOUB,UW.LRIV 2012-01-01 2012-12-31 --workers $W
step noise python scripts/scan3_noise.py
step masks python scripts/scan_masks2.py
step bp python scripts/scan_bp.py $SD/bp_2012g --workers $W
step bp_ref python scripts/scan_bp_ref.py $SD/bp2_ref --stations $REFS --workers $W
step refqc python scripts/scan_refqc.py
step inj2_list python scripts/scan_inject.py
step inj3_list python scripts/scan_inject3.py
# CPU (4 workers): injected band-pass sets; GPU (1 CPU) in parallel: real + controls
( step bp_inj2 python scripts/scan_bp.py $SD/bp_2012inj2 --inject $SD/inject_2012.parquet --workers 4
  step bp_ref_inj2 python scripts/scan_bp_ref.py $SD/bp2_ref_inj2 --stations $REFS --inject $SD/inject_2012.parquet --workers 4
  step bp_inj3 python scripts/scan_bp.py $SD/bp_2012inj3 --inject $SD/inject3_2012.parquet --workers 4
  step bp_ref_inj3 python scripts/scan_bp_ref.py $SD/bp2_ref_inj3 --stations $REFS --inject $SD/inject3_2012.parquet --workers 4 ) &
CPU=$!
det(){ # run bp ndraw shift
  [ -f $SD/runs/$1/DONE ] || { python scripts/scan_detect.py $1 --bp $2 --ndraw $3 --shift $4 >> $L 2>&1; log "detect $1 rc $?"; }; }
det real $SD/bp_2012g 50 0
for s in 3 5 7 11 13; do det ctrlA_s$s $SD/bp_2012g 0 $s; done
wait $CPU
det inj2 $SD/bp_2012inj2 0 0
det inj3 $SD/bp_2012inj3 0 0
col(){ # run bp refbp [inj]
  [ -s $SD/runs/$1/cand2.parquet ] || { python scripts/scan_columns2.py $1 --bp $2 --refbp $3 ${4:+--inj $4} --workers $W >> $L 2>&1; log "columns2 $1 rc $?"; }; }
col real $SD/bp_2012g $SD/bp2_ref
for s in 3 5 7 11 13; do col ctrlA_s$s $SD/bp_2012g $SD/bp2_ref; done
col inj2 $SD/bp_2012inj2 $SD/bp2_ref_inj2 $SD/inject_2012.parquet
col inj3 $SD/bp_2012inj3 $SD/bp2_ref_inj3 $SD/inject3_2012.parquet
n=0; for r in real ctrlA_s3 ctrlA_s5 ctrlA_s7 ctrlA_s11 ctrlA_s13 inj2 inj3; do [ -s $SD/runs/$r/cand2.parquet ] && n=$((n+1)); done
log "queue end: $n of 8 runs with columns"
[ $n -eq 8 ] && touch $D/SCAN3_QUEUE_DONE
