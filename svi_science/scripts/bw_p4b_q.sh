#!/usr/bin/env bash
# P4-2b step 5: completeness-ladder GPU queue (ONE GPU process at a time, one shard) + step-3 sidecar flags (CPU).
# Old station set, frozen detector, backward mode = the p5f configuration except --sset eff --excl-file (P4-2 mechanism).
# Each year of a ladder waits for P2_<Y>_DONE.  Markers p4b/LAD_<OLD>_DONE|FAILED; log p4b/logs/q.log.  Nothing deleted.
# Rules: results/backward/p4b/P4B.md.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
P=$O/p4b
echo $$ > $P/Q_PID
export SVI_REF_DIR=$P/mf
CD=$CB,$FWD
lg(){ echo "[q $(TZ=America/Los_Angeles date '+%F %T %Z')] $*" >> $P/logs/q.log; }
waitm(){ while [ ! -f $O/$1 ] && [ ! -f $O/${1/DONE/FAILED} ]; do sleep 600; done; [ -f $O/$1 ]; }
guard(){ python -c "import sys; sys.path.insert(0,'$LFE/scripts'); import svi_ref_common as R; R.guard('mf_$1','mf_$1_s*')" >> $P/logs/guard.log 2>&1; }
flags(){ python $S/bw_p4b_flags.py >> $P/logs/flags.log 2>&1 || lg "flags pass FAILED rc=$?"; }
lad(){ # OLD BOUNDARY(YYYYMMDD) YEAR...
  local old=$1 b=$2; shift 2; local run=lad_${old}_${b} ok=1
  [ -f $P/LAD_${old}_DONE ] && return 0
  for y in "$@"; do
    if ! waitm P2_${y}_DONE; then lg "$run: P2_${y} FAILED"; ok=0; continue; fi
    lg "$run $y start"
    python $S/bw_gate_mf.py --backward -- --run $run --kind prod --sset eff --excl-file $P/tables/${run}_excl.csv \
      --clean 0 --days-file $P/tables/${run}_days_$y.json --cache-dirs $CD --shard 0 --nshard 1 >> $P/logs/$run.log 2>&1 || ok=0
    lg "$run $y end ok=$ok"
  done
  guard $run
  [ $ok = 1 ] && date -u +%FT%TZ > $P/LAD_${old}_DONE || date -u +%FT%TZ > $P/LAD_${old}_FAILED
}
lg "queue start PID $$"
lad CN4 20050701 2005
lad B1 20051102 2005
flags
lad B2 20071113 2007 2008
flags
lad B3 20090827 2009
flags
waitm P5F_2008_DONE; flags
waitm P5F_2009_DONE; flags
lg "queue end"
date -u +%FT%TZ > $P/Q_DONE
