#!/usr/bin/env bash
# P1 queue A (NRCan-bound): CN4 2003-2009, POLARIS PO 2003-2006, LZB 2003-01..2006-06.  7 NRCan lanes (cap 8 incl. queue B's 1).
# Per year: download -> 2 repair passes (re-ask only missing channels) -> marker P1A_<Y>_DONE.  Ledger tables/ledger_A.csv.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
export BW_LEDGER=$O/tables/ledger_A.csv
echo $$ > $O/P1A_PID
DL="python $S/bw_gate_download.py --lanes-nrcan 7 --lanes-es 1 --tries 20 --skip-cached-in $FWD"
for Y in 2003 2004 2005 2006 2007 2008 2009; do
  [ -f $O/P1A_${Y}_DONE ] && continue
  if [ $Y -le 2005 ]; then ST=$CN4,$PO,LZB; E=$Y-12-31; elif [ $Y = 2006 ]; then ST=$CN4,$PO; E=$Y-12-31; else ST=$CN4; E=$Y-12-31; fi
  ok=1
  retry $DL --stations $ST --start $Y-01-01 --end $E --tag p1a_$Y || ok=0
  [ $Y = 2006 ] && { retry $DL --stations LZB --start 2006-01-01 --end 2006-06-30 --tag p1a_${Y}_lzb || ok=0; }
  for r in 1 2; do
    retry $DL --stations $ST --start $Y-01-01 --end $E --repair --tag p1a_${Y}_repair$r || ok=0
  done
  [ $ok = 1 ] && date -u +%FT%TZ > $O/P1A_${Y}_DONE || date -u +%FT%TZ > $O/P1A_${Y}_FAILED
done
date -u +%FT%TZ > $O/P1A_DONE
