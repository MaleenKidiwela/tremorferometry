#!/bin/bash
# S11 ballistic round 3, Amendment 5 (decade-half DH + 4-year Q4 stacks). DONE markers; idempotent.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0 S11B3_SCHEMES=DH,Q4
cd /home/jovyan/tremorferometry/svi_science/scripts
OUT=/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic/round3
echo $$ > $OUT/logs/queue_dh.pid
[ -f $OUT/MEASURE_SVI_dh_DONE ] || python s11_ball3_measure.py --region SVI --workers 8 --tag _dh >> $OUT/logs/measure_dh_nohup.log 2>&1
[ -f $OUT/MEASURE_WA_NW_dh_DONE ] || python s11_ball3_measure.py --region WA_NW --workers 8 --tag _dh >> $OUT/logs/measure_dh_nohup.log 2>&1
if [ -f $OUT/MEASURE_SVI_dh_DONE ] && [ -f $OUT/MEASURE_WA_NW_dh_DONE ] && [ ! -f $OUT/A_DH_DONE ] && [ -f s11_ball3_dh.py ]; then
  python s11_ball3_dh.py >> $OUT/logs/analysis_dh.log 2>&1 && date > $OUT/A_DH_DONE
fi
[ -f $OUT/A_DH_DONE ] && date > $OUT/QUEUE_DH_DONE
