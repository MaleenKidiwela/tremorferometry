#!/bin/bash
# hg_repro main queue (PLAN.md §4). Markers in hg_repro/: PREP_DONE CONT_DONE MF_REAL_DONE MF_NOISE_DONE ANALYZE_DONE.
# usage: nohup bash scripts/hg_queue.sh FETCH_PID CNPREP_PID > logs/queue.log 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ts() { TZ=America/Los_Angeles date '+%F %T PT'; }
echo "$(ts) queue start"
if [ ! -f PREP_DONE ]; then
  tail --pid=$2 -f /dev/null           # CN prep (4 workers) finishes first
  python scripts/hg_prep.py TA.A04D,TA.B05D,TA.D03D 2011-06-01 2013-12-01 --workers 4 > logs/prep_ta.out 2>&1
  tail --pid=$1 -f /dev/null           # fetch finished
  [ -f FETCH_DONE ] || { echo "$(ts) fetch did not finish cleanly"; exit 1; }
  python scripts/hg_prep.py UW.SP2,UW.LRIV,UW.DOSE,UW.GNW 2011-06-01 2013-12-01 --workers 4 > logs/prep_uw.out 2>&1
  touch PREP_DONE; echo "$(ts) PREP_DONE"
fi
if [ ! -f CONT_DONE ]; then
  python scripts/hg_cont.py CN.NLLB,CN.VGZ,UW.LRIV,UW.DOSE,UW.GNW,UW.SP2,TA.A04D,TA.B05D,TA.D03D --workers 3 > logs/cont.out 2>&1 && touch CONT_DONE
  echo "$(ts) CONT_DONE"
fi
if [ ! -f MF_REAL_DONE ]; then
  for K in 1 2 3 4; do python scripts/hg_mf.py real $K > logs/mf_real_$K.out 2>&1 & done; wait
  python scripts/hg_mf.py real 5 > logs/mf_real_5.out 2>&1
  ls data/mf/real_T{1,2,3,4,5}.npz > /dev/null && touch MF_REAL_DONE; echo "$(ts) MF_REAL_DONE"
fi
if [ ! -f MF_NOISE_DONE ]; then
  [ -f data/mf/noise_times.csv ] || python scripts/hg_mf.py noisesel > logs/mf_noisesel.out 2>&1
  for c in "0 25" "25 50" "50 75" "75 100"; do python scripts/hg_mf.py noise $c > logs/mf_noise_${c// /_}.out 2>&1 & done; wait
  [ $(ls data/mf/noise_[0-9][0-9][0-9].npz | wc -l) -eq 100 ] && touch MF_NOISE_DONE; echo "$(ts) MF_NOISE_DONE"
fi
python scripts/hg_analyze.py > logs/analyze.out 2>&1 && touch ANALYZE_DONE
echo "$(ts) queue end"
