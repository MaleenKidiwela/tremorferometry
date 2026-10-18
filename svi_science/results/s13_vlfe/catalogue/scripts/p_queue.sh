#!/bin/bash
# Amendment P queue: (B) early-era injection 2011-12, then (A) pilot 2014. Waits for the 2014 scan year.
# Markers: INJ_EARLY_DONE, P14_DONE (catalogue/). Every non-zero exit writes a FAIL line to logs/p_queue.log.
# Resume: re-run; finished stages are skipped (STEP markers in data/p14/, INJ_EARLY_DONE).
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
export OMP_NUM_THREADS=1
L=logs/p_queue.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $L; }
run(){ # name cmd...  (skips if data/p14/STEP_name_DONE)
  local n=$1; shift
  [ -f data/p14/STEP_${n}_DONE ] && return 0
  log "start $n"; "$@" >> logs/p_queue.out 2>&1; local rc=$?
  if [ $rc -ne 0 ]; then log "FAIL $n rc $rc (see logs/p_queue.out)"; exit $rc; fi
  touch data/p14/STEP_${n}_DONE; log "end $n"; }
mkdir -p data/p14
until [ -f ../data/scan_2014/SCAN2014_DONE ]; do
  pgrep -f scan2014_queue.sh > /dev/null || { log "FAIL scan2014 queue not running and no SCAN2014_DONE"; exit 1; }
  sleep 60; done
if [ ! -f INJ_EARLY_DONE ]; then
  run early2011 env SCAN_YEAR=2011 python scripts/cat_inj.py early 2011
  run early2012 env SCAN_YEAR=2012 python scripts/cat_inj.py early 2012
  run early_summary python scripts/cat_inj.py early_summary
  TZ=America/Los_Angeles date > INJ_EARLY_DONE; log "B done"
fi
run p14_cont python scripts/cat_p14.py cont
run p14_bank python scripts/cat_p14.py bank
run p14_mf python scripts/cat_p14.py mf
run p14_detect python scripts/cat_p14.py detect
run p14_battery env SCAN_YEAR=2014 python scripts/cat_p14.py battery
run p14_inj_list python scripts/cat_p14.py inj_list
run p14_inj env SCAN_YEAR=2014 python scripts/cat_p14.py inj
run p14_summary python scripts/cat_p14.py summary
TZ=America/Los_Angeles date > P14_DONE; log "A done"
