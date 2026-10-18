#!/bin/bash
# S13 scale-up queue (PLAN.md Amendment 1). Idempotent: stages with DONE markers are skipped; prep skips existing days.
# Run: nohup bash scripts/s13_scale_queue.sh > logs/scale_queue.out 2>&1 &     Log: logs/scale_queue.log
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1
export S13_DISP=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/disp1sps_v2   # Amendment 2 (short gaps joined)
LOG=logs/scale_queue.log
ts(){ TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M PT'; }
SVI_STA=CN.PGC,CN.VGZ,CN.SNB,CN.NLLB,CN.LZB,CN.PFB,CN.YOUB,CN.SYMB,CN.CLRS,CN.PTRF,C8.SOKB,C8.SHVB,C8.KHVB,C8.MGCB,C8.TWBB,C8.SHDB,C8.BPCB,UW.LRIV,UW.OPC,TA.B04A,UW.SQM,UW.ODUC,UW.OSQM,UW.SLDQ,UW.ODEER,UW.OTR,XU.BS11
WA_STA=UW.DOSE,UW.GNW,UW.LRIV,UW.SQM,UW.OSQM,UW.ODEER,UW.NATEM,UW.CHIMA,UW.EQUIL,UW.HDW,UW.HANS,UW.COYL,UW.KNGST,UW.KTSAP,UW.SKOKO,UW.MORSE,UW.WISH,UW.SP2,UW.OPC,UW.PRIN,UW.WINDI,UW.EAGLE,US.NLWA,TA.D03D,TA.D04E,TA.D04D,TA.B04A,YW.FACB,XU.W020,XU.W030,XU.W040,XU.C04A,XU.BS11,XU.PL11,XU.W060,XU.N040,XU.N050,XU.N060,CN.PGC,CN.VGZ
echo "$(ts) queue start" >> $LOG
for R in svi wa; do
  [ -f UNITS_${R}_DONE ] || { python scripts/s13_units.py $R >> $LOG 2>&1 && touch UNITS_${R}_DONE; }
  STA=$([ $R = svi ] && echo $SVI_STA || echo $WA_STA)
  if [ ! -f PREP_${R}_DONE ]; then
    echo "$(ts) prep $R" >> $LOG
    python scripts/s13_prep.py $STA x x --days-file tables/${R}_prep_days.csv --workers 8 > logs/prep_${R}.out 2>&1 && ! grep -q "ModuleNotFoundError" logs/prep_${R}.out && touch PREP_${R}_DONE
  fi
  # episodes, 2 at a time (each stacker is single-threaded)
  n=0
  for f in data/rows/${R}_*.parquet; do
    tag=scale_$(basename $f .parquet)
    [ -f STACK_${tag}_DONE ] && continue
    ( python scripts/s13_stack2.py $f $tag $STA > logs/${tag}.out 2>&1 && echo "$(ts) $tag done" >> $LOG || echo "$(ts) $tag FAILED" >> $LOG ) &
    n=$((n+1)); if [ $n -ge 4 ]; then wait -n; n=$((n-1)); fi
  done
  wait
done
python scripts/s13_scale_summary.py > logs/scale_summary.out 2>&1 && touch SCALE_SUMMARY_DONE
ls STACK_scale_*_DONE | wc -l >> $LOG
miss=0; for f in data/rows/*.parquet; do [ -f STACK_scale_$(basename $f .parquet)_DONE ] || miss=$((miss+1)); done   # guard added 2026-10-14
if [ -f SCALE_SUMMARY_DONE ] && [ $miss -eq 0 ]; then touch SCALE_QUEUE_DONE; else rm -f SCALE_SUMMARY_DONE; echo "$(ts) NOT done: $miss stacks missing" >> $LOG; fi
echo "$(ts) queue end" >> $LOG
