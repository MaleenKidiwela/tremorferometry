#!/bin/bash
# S11 ballistic round 3 queue (ballistic/round3/PLAN.md §5). Each stage writes a DONE marker; finished stages are skipped
# on restart (idempotent). Single-threaded BLAS per worker, 8 workers; one region in memory at a time (< 40 GB).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/jovyan/tremorferometry/svi_science/scripts
OUT=/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic/round3
mkdir -p $OUT/logs
echo $$ > $OUT/logs/queue.pid
echo "queue start $(TZ=America/Los_Angeles date)" >> $OUT/logs/queue.log
[ -f $OUT/MEASURE_SVI_DONE ] || python s11_ball3_measure.py --region SVI --workers 8 >> $OUT/logs/measure_svi_nohup.log 2>&1
[ -f $OUT/MEASURE_WA_NW_DONE ] || python s11_ball3_measure.py --region WA_NW --workers 8 >> $OUT/logs/measure_wa_nohup.log 2>&1
# analysis stages (separate file, each stage with its own DONE marker)
if [ -f $OUT/MEASURE_SVI_DONE ] && [ -f $OUT/MEASURE_WA_NW_DONE ] && [ -f s11_ball3_analysis.sh ]; then
  bash s11_ball3_analysis.sh >> $OUT/logs/analysis.log 2>&1
fi
[ -f $OUT/MEASURE_SVI_DONE ] && [ -f $OUT/MEASURE_WA_NW_DONE ] && [ -f $OUT/ANALYSIS_DONE ] && TZ=America/Los_Angeles date > $OUT/QUEUE_DONE
echo "queue end $(TZ=America/Los_Angeles date)" >> $OUT/logs/queue.log
