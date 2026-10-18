#!/usr/bin/env bash
# P1 queue A v2 (NRCan; replaces bw_prod_p1a.sh from 2026-10-09): (1) PO.KELB 2003-01-01..2004-06-03; (2) the Bostock ETS
# windows +-2 months (WINDOWS in bw_prod_env.sh) for every needed NRCan station; (3) the rest, year by year 2004..2009.
# Each block: download -> 2 repair passes (only days on disk with a missing channel; days that failed every try of every
# channel are real gaps and are not re-asked).  Markers P1A_KELB_DONE, P1A_W<name>_DONE, P1A_<Y>_DONE.  Ledger ledger_A.csv.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
export BW_LEDGER=$O/tables/ledger_A.csv
echo $$ > $O/P1A_PID
DL="python $S/bw_gate_download.py --lanes-nrcan 7 --lanes-es 1 --tries 20 --skip-cached-in $FWD"
ALL=$CN4,$PO,LZB,KELB                       # stations outside their metadata epochs (or LZB after 2006-06) are skipped
block(){ # marker stations start end tag
  [ -f $O/$1 ] && return 0
  local ok=1
  retry $DL --stations $2 --start $3 --end $4 --tag $5 || ok=0
  for r in 1 2; do retry $DL --stations $2 --start $3 --end $4 --repair --tag ${5}_repair$r || ok=0; done
  [ $ok = 1 ] && date -u +%FT%TZ > $O/$1 || date -u +%FT%TZ > ${O}/${1/DONE/FAILED}
}
block P1A_KELB_DONE KELB 2003-01-01 2004-06-03 p1a_kelb
for w in $WINDOWS; do IFS=: read n a b <<< "$w"; block P1A_${n}_DONE $ALL $a $b p1a_$n; done
for Y in 2004 2005 2006 2007 2008 2009; do block P1A_${Y}_DONE $ALL $Y-01-01 $Y-12-31 p1a_$Y; done
date -u +%FT%TZ > $O/P1A_DONE
