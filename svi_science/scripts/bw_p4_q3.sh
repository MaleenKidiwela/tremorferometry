#!/usr/bin/env bash
# P4 CPU queue (no GPU): CN4 hourly energy + quiet flag year by year behind P2, the final P4-1 anchor (2010-13) and the
# quiet-flag validation behind P5F_2012/2013, then P4-3/P4-4 (era thresholds) behind P5 2003-2006.
# Markers in p4/: P4_E_<Y>_DONE (energy+flag of year Y), P4_1_DONE|P4_1_FAILED (final anchor), P4_QV_DONE (flag
# validation written, pass or fail inside qflag_validation.json), P4_3_DONE|P4_3_FAILED.  Nothing deleted.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
P=$O/p4
echo $$ > $P/Q3_PID
L=$P/logs/q3.log
w2(){ while [ ! -f $O/P2_$1_DONE ] && [ ! -f $O/P2_$1_FAILED ]; do sleep 300; done; }
wf(){ while [ ! -f $O/$1_DONE ] && [ ! -f $O/$1_FAILED ]; do sleep 300; done; }
for Y in 2012 2013 2003 2004 2005 2006 2007 2008 2009; do
  if [ ! -f $P/P4_E_${Y}_DONE ]; then
    w2 $Y
    echo "[q3] energy $Y $(date -u +%FT%TZ)" >> $L
    python $S/bw_p4_energy.py --years $Y --workers 8 >> $L 2>&1 && python $S/bw_p4_qflag.py >> $L 2>&1 \
      && date -u +%FT%TZ > $P/P4_E_${Y}_DONE || echo "[q3] energy/flag $Y FAILED" >> $L
  fi
  if [ $Y = 2013 ]; then
    wf P5F_2012; wf P5F_2013
    python $S/bw_p4_qflag.py >> $L 2>&1                       # runs the one-time 2012-13 validation
    [ -f $P/tables/qflag_validation.json ] && date -u +%FT%TZ > $P/P4_QV_DONE
    python $S/bw_p4_anchor.py > $P/logs/p4_1_final.log 2>&1
    python -c "import json,sys; d=json.load(open('$P/tables/p4_1_decision.json')); sys.exit(0 if d['final'] and d['pass_'] else 1)" \
      && date -u +%FT%TZ > $P/P4_1_DONE || date -u +%FT%TZ > $P/P4_1_FAILED
    python $S/bw_p4_era.py > $P/logs/p4_flagdiag.log 2>&1
  fi
done
for Y in 2003 2004 2005 2006; do wf P5F_$Y; wf P5E_$Y; done
python $S/bw_p4_era.py > $P/logs/p4_3.log 2>&1 && date -u +%FT%TZ > $P/P4_3_DONE || date -u +%FT%TZ > $P/P4_3_FAILED
