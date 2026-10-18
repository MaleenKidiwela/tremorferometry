#!/bin/bash
# rf_gosselin processing queue (PLAN.md §3-5). Idempotent: stages with a DONE marker are skipped.
# Waits for download markers; CPU <= 8 (rf03 Pool(8); rf04/rf05 single process).
W=/home/jovyan/tremorferometry/svi_science/results/s12_lvl/rf_gosselin
RAW=/home/jovyan/tremorferometry/data/rf_teleseis/raw
cd $W
PY="python -P"
log(){ echo "[$(TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M %Z')] $*" >> $W/logs/queue.log; }
CORE="CN.PGC CN.SNB CN.VGZ CN.LZB CN.PFB"
EXT_SVI="CN.NLLB CN.YOUB CN.MGB CN.GOBB"
EXT_WA="UW.GNW UW.SP2 UW.LON UW.OFR UW.WISH US.NLWA"
waitfor(){ for k in "$@"; do until [ -f $RAW/$k.DONE ]; do sleep 60; done; done; }
stage(){ name=$1; shift; if [ -f $W/${name}_DONE ]; then log "skip $name"; return; fi
  log "start $name: $*"; if "$@" >> $W/logs/$name.log 2>&1; then touch $W/${name}_DONE; log "done $name"; else log "FAILED $name"; fi; }
log "queue start (pid $$)"
waitfor $CORE
stage RF_CORE $PY scripts/rf03_rf.py $CORE
C=$(echo $CORE | tr ' ' ',')
stage REPRO   $PY scripts/rf04_analysis.py --tag repro   --t0 1993-01-01 --t1 2019-06-01 --ets SVI_ALBH --band 05 --stations $C
stage REPRO07 $PY scripts/rf04_analysis.py --tag repro07 --t0 1993-01-01 --t1 2019-06-01 --ets SVI_ALBH --band 07 --stations $C
stage OOS     $PY scripts/rf04_analysis.py --tag oos2019 --t0 2019-06-01 --t1 2026-12-31 --ets SVI_ALBH --band 05 --stations $C
stage FULL    $PY scripts/rf04_analysis.py --tag full    --t0 1993-01-01 --t1 2026-12-31 --ets SVI_ALBH --band 05 --stations $C
stage REC_REPRO $PY scripts/rf05_recovery.py --tag repro --stations $C
stage REC_FULL  $PY scripts/rf05_recovery.py --tag full --stations $C
waitfor $EXT_SVI $EXT_WA
stage RF_EXT $PY scripts/rf03_rf.py $EXT_SVI $EXT_WA
S=$(echo $EXT_SVI | tr ' ' ','); A=$(echo $EXT_WA | tr ' ' ',')
stage EXT_SVI $PY scripts/rf04_analysis.py --tag ext_svi --t0 1993-01-01 --t1 2026-12-31 --ets SVI_ALBH --band 05 --stations $S
stage EXT_WA  $PY scripts/rf04_analysis.py --tag ext_wa  --t0 2007-01-01 --t1 2026-12-31 --ets WA_SEAT --band 05 --stations $A
stage REC_EXT_SVI $PY scripts/rf05_recovery.py --tag ext_svi --stations $S
stage REC_EXT_WA  $PY scripts/rf05_recovery.py --tag ext_wa --stations $A
stage SLIP_SVI $PY scripts/rf05_recovery.py --slipfront --ets SVI_ALBH --stations $C,$S
stage SLIP_WA  $PY scripts/rf05_recovery.py --slipfront --ets WA_SEAT --stations $A
touch $W/QUEUE_DONE; log "queue finished"
