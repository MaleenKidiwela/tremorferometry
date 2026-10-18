#!/bin/bash
# S11 ballistic round 2 queue (round2/PLAN.md §6). Each stage writes a DONE marker; finished stages are skipped on restart.
# Single-threaded BLAS per worker (8 workers total, other queues share the machine).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/jovyan/tremorferometry/svi_science/scripts
OUT=/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic/round2
echo $$ > $OUT/logs/queue.pid
echo "queue start $(date)" >> $OUT/logs/queue.log
[ -f $OUT/MEASURE_DONE ] || python s11_ball2_measure.py --workers 8 >> $OUT/logs/measure_nohup.log 2>&1
[ -f $OUT/LOCUS_MEASURE_DONE ] || python s11_ball2_locus.py measure --workers 3 >> $OUT/logs/locus_nohup.log 2>&1
[ -f $OUT/MEASURE_DONE ] && [ -f $OUT/LOCUS_MEASURE_DONE ] && date > $OUT/QUEUE_DONE
echo "queue end $(date)" >> $OUT/logs/queue.log
