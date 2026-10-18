#!/usr/bin/env bash
# P2 v2 (replaces bw_prod_p2.sh from 2026-10-09): caches behind P1 queue A v2: KELB, then each ETS window, then each year
# (queue B is complete for every year).  Markers P2_KELB_DONE, P2_W<name>_DONE, P2_<Y>_DONE.  12 workers (~30 GB RAM).
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
echo $$ > $O/P2_PID
NS=CN.PGC,CN.VGZ,CN.NLLB,CN.SNB,UW.LRIV,PB.B001,PB.B003,PB.B004,PB.B005,PB.B006,PB.B007,PB.B009,PB.B010,PB.B011,PB.B926
NSE=C8.SHVB,C8.MGCB,C8.TWBB,C8.SHDB,C8.KHVB,PO.TWKB,PO.TSJB,PO.KLNB,PO.SILB,PO.SSIB,PO.TWGB,CN.LZB,PO.KELB
waitm(){ while [ ! -f $O/$1 ] && [ ! -f ${O}/${1/DONE/FAILED} ]; do sleep 300; done; }
cache(){ # marker start end
  [ -f $O/$1 ] && return 0
  python $S/bw_gate_cache.py --stations $NS,$NSE --start $2 --end $3 --workers 12 --status-dir $O/tables \
    > $O/logs/p2_${1%_DONE}.log 2>&1 && date -u +%FT%TZ > $O/$1 || date -u +%FT%TZ > ${O}/${1/DONE/FAILED}
}
waitm P1A_KELB_DONE; cache P2_KELB_DONE 2003-01-01 2004-06-03
for w in $WINDOWS; do IFS=: read n a b <<< "$w"; waitm P1A_${n}_DONE; cache P2_${n}_DONE $a $b; done
for Y in 2004 2005 2006 2007 2008 2009; do waitm P1A_${Y}_DONE; cache P2_${Y}_DONE $Y-01-01 $Y-12-31; done
date -u +%FT%TZ > $O/P2_DONE
