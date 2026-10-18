#!/bin/bash
# Amendment Q1 queue (2011-12 LOO re-test). Waits for the pilot-2014 re-run (P14_DONE). Marker Q1_DONE.
# FAIL lines -> logs/q_main.log. Resume: re-run (stages with STEP markers in data/q1/ are skipped).
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
export OMP_NUM_THREADS=1
L=logs/q_main.log; mkdir -p data/q1
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') q1 $*" >> $L; }
run(){ local n=$1; shift; [ -f data/q1/STEP_${n}_DONE ] && return 0
  log "start $n"; "$@" >> logs/q1.out 2>&1; local rc=$?
  if [ $rc -ne 0 ]; then log "FAIL $n rc $rc (see logs/q1.out)"; exit $rc; fi
  touch data/q1/STEP_${n}_DONE; log "end $n"; }
until [ -f P14_DONE ]; do pgrep -f p_queue.sh > /dev/null || { log "FAIL p_queue not running and no P14_DONE"; exit 1; }; sleep 30; done
run cont python scripts/cat_cont.py CN.LZB,UW.WISH,US.NLWA,CC.OBSR,TA.D04D,UW.ELW,UW.LEBA,UW.RATT
run bank python scripts/cat_q1.py bank
run mf python scripts/cat_q1.py mf
run detect python scripts/cat_q1.py detect
run bat2011 env SCAN_YEAR=2011 python scripts/cat_q1.py battery 2011
run bat2012 env SCAN_YEAR=2012 python scripts/cat_q1.py battery 2012
run summary python scripts/cat_q1.py summary
TZ=America/Los_Angeles date > Q1_DONE; log "done"
