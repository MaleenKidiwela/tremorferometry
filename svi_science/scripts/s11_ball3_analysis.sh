#!/bin/bash
# S11 ballistic round 3 analysis stages (called by s11_ball3_queue.sh after both measurements). DONE markers per stage.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
cd /home/jovyan/tremorferometry/svi_science/scripts
OUT=/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic/round3
st() { local name=$1; shift; [ -f $OUT/${name}_DONE ] && return 0
       echo "== $name start $(TZ=America/Los_Angeles date)"; "$@" && date > $OUT/${name}_DONE; }
st A_CLOCK_SVI   python s11_ball3_clock.py --region SVI
st A_CLOCK_WA    python s11_ball3_clock.py --region WA_NW
st A_LONG_SVI    python s11_ball3_longterm.py --region SVI --workers 8
st A_LONG_WA     python s11_ball3_longterm.py --region WA_NW --workers 8
st A_PHASE       python s11_ball3_phase.py
st A_COMMON      python s11_ball3_common.py
st A_JOINTP      python s11_ball3_jointp.py
st A_FIGURES     python s11_ball3_figures.py
n=$(ls $OUT/A_*_DONE 2>/dev/null | wc -l)
[ "$n" -ge 8 ] && date > $OUT/ANALYSIS_DONE
echo "analysis end $(TZ=America/Los_Angeles date): $n/8 stages done"
