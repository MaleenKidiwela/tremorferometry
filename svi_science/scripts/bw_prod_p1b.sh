#!/usr/bin/env bash
# P1 queue B (EarthScope-bound): PB + LRIV every year, CN4 2010-2013 (EarthScope first), C8 POLARIS 2003-2006.
# 6 EarthScope lanes, 1 NRCan lane (fallback only).  Tier 1 (2010-2013) first.  Marker P1B_<Y>_DONE.  Ledger ledger_B.csv.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
export BW_LEDGER=$O/tables/ledger_B.csv
echo $$ > $O/P1B_PID
DL="python $S/bw_gate_download.py --lanes-nrcan 1 --lanes-es 6 --tries 20 --skip-cached-in $FWD"
for Y in 2010 2011 2012 2013 2003 2004 2005 2006 2007 2008 2009; do
  [ -f $O/P1B_${Y}_DONE ] && continue
  ST=$PB,LRIV; [ $Y -ge 2010 ] && ST=$ST,$CN4; [ $Y -le 2006 ] && ST=$ST,$C8
  ok=1
  retry $DL --stations $ST --start $Y-01-01 --end $Y-12-31 --tag p1b_$Y || ok=0
  if [ $Y -ge 2010 ]; then retry $DL --stations $CN4 --start $Y-01-01 --end $Y-12-31 --repair --tag p1b_${Y}_repair || ok=0; fi
  if [ $Y -le 2006 ]; then retry $DL --stations $C8 --start $Y-01-01 --end $Y-12-31 --repair --tag p1b_${Y}_repair || ok=0; fi
  [ $ok = 1 ] && date -u +%FT%TZ > $O/P1B_${Y}_DONE || date -u +%FT%TZ > $O/P1B_${Y}_FAILED
done
date -u +%FT%TZ > $O/P1B_DONE
