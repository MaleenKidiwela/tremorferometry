#!/usr/bin/env bash
# P2: 20 Hz / 1-8 Hz caches (forward chain, bw_gate_cache.py) behind P1, one year at a time, as soon as BOTH P1 queues
# marked that year (2010-2013 need only queue B).  12 workers (~30 GB RAM).  Marker P2_<Y>_DONE.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
echo $$ > $O/P2_PID
NS=CN.PGC,CN.VGZ,CN.NLLB,CN.SNB,UW.LRIV,PB.B001,PB.B003,PB.B004,PB.B005,PB.B006,PB.B007,PB.B009,PB.B010,PB.B011,PB.B926
NSE=C8.SHVB,C8.MGCB,C8.TWBB,C8.SHDB,C8.KHVB,PO.TWKB,PO.TSJB,PO.KLNB,PO.SILB,PO.SSIB,PO.TWGB,CN.LZB
for Y in 2010 2011 2012 2013 2003 2004 2005 2006 2007 2008 2009; do
  [ -f $O/P2_${Y}_DONE ] && continue
  while true; do
    a=1; [ $Y -le 2009 ] && { [ -f $O/P1A_${Y}_DONE ] || [ -f $O/P1A_${Y}_FAILED ] || a=0; }
    { [ -f $O/P1B_${Y}_DONE ] || [ -f $O/P1B_${Y}_FAILED ]; } && [ $a = 1 ] && break
    sleep 300
  done
  L=$NS; [ $Y -le 2006 ] && L=$NS,$NSE
  python $S/bw_gate_cache.py --stations $L --start $Y-01-01 --end $Y-12-31 --workers 12 --status-dir $O/tables \
    > $O/logs/p2_$Y.log 2>&1 && date -u +%FT%TZ > $O/P2_${Y}_DONE || date -u +%FT%TZ > $O/P2_${Y}_FAILED
done
date -u +%FT%TZ > $O/P2_DONE
