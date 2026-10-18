#!/bin/bash
# PLAN §10 Amendment C (C1) lead test: 3 lanes (<= 3 CPUs, OMP 1); DONE markers per run, LEAD_QUEUE_DONE at the end.
W=/home/jovyan/tremorferometry/svi_science/results/s12_lvl/acf_pgc; cd $W
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run(){ s=$1; b=$2; m=$3; k=LEAD_${s}_${b}_${m}_DONE; [ -f $k ] && return
  python -P scripts/acf08_lead_test.py --station $s --band $b --mask $m --ndraw 2000 > logs/lead_${s}_${b}_${m}.log 2>&1 && touch $k; }
( run CN.VGZ b05_2 tremorfree; run CN.VGZ b1_4 none; run CN.SNB b05_2 none ) &
( run CN.VGZ b05_2 none; run CN.SNB b05_2 tremorfree ) &
( run CN.VGZ b1_4 tremorfree; run CN.PGC b05_2 tremorfree ) &
wait; touch LEAD_QUEUE_DONE
