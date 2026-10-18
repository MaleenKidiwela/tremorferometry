#!/bin/bash
# Amendment R queue: 2013 injections at n=800 + R summary; then pilot 2018 (after scan year 2018) + R summary.
# Markers R2013_DONE, PILOT2018_DONE. FAIL lines -> logs/q_main.log. Resume: re-run (STEP markers skip).
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue
export OMP_NUM_THREADS=1 NINJ=400
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') R $*" >> $L; }
st(){ local Y=$1 m=$2; shift 2; [ -f data/pilot$Y/STEP_R_${m}_DONE ] && return 0; log "start $Y $m"
  "$@" >> logs/pilot$Y.out 2>&1 || { log "FAIL $Y $m (see logs/pilot$Y.out)"; exit 1; }; touch data/pilot$Y/STEP_R_${m}_DONE; log "end $Y $m"; }
until [ -f ../data/scan_2018/SCAN2018_DONE ]; do pgrep -f "scan_year_queue.sh 2018" > /dev/null || { log "FAIL scan2018 not running, no marker"; exit 1; }; sleep 60; done
# 2013 at n = 800
if [ ! -f R2013_DONE ]; then
  [ -f data/pilot2013/inj_list_n200.parquet ] || { mv data/pilot2013/inj_list.parquet data/pilot2013/inj_list_n200.parquet; mv data/pilot2013/inj.parquet data/pilot2013/inj_n200.parquet; }
  st 2013 inj_list env PILOT=2013 SCAN_YEAR=2013 python scripts/cat_pilot.py inj_list
  st 2013 inj env PILOT=2013 SCAN_YEAR=2013 python scripts/cat_pilot.py inj
  st 2013 summary env PILOT=2013 SCAN_YEAR=2013 python scripts/cat_pilot.py summary
  st 2013 R env PILOT=2013 SCAN_YEAR=2013 python scripts/cat_r.py
  TZ=America/Los_Angeles date > R2013_DONE
fi
mkdir -p data/pilot2018
for m in cont bank mf detect battery inj_list inj summary; do st 2018 $m env PILOT=2018 SCAN_YEAR=2018 python scripts/cat_pilot.py $m; done
st 2018 R env PILOT=2018 SCAN_YEAR=2018 python scripts/cat_r.py
TZ=America/Los_Angeles date > PILOT2018_DONE; log "pilot 2018 done"
