#!/bin/bash
# hg_repro queue 2: after the main queue -> GF caches -> stack-and-invert -> controls. Markers GF_DONE STACKINV_DONE CONTROLS_DONE.
# usage: nohup bash scripts/hg_queue2.sh MAINQUEUE_PID > logs/queue2.log 2>&1 &
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ts() { TZ=America/Los_Angeles date '+%F %T PT'; }
[ -n "$1" ] && tail --pid=$1 -f /dev/null
[ -f ANALYZE_DONE ] || { echo "$(ts) main queue did not finish (no ANALYZE_DONE)"; exit 1; }
echo "$(ts) queue2 start"
[ -f GF_DONE ] || { python scripts/hg_stackinv.py gf > logs/gf.out 2>&1 && touch GF_DONE; echo "$(ts) GF_DONE"; }
python scripts/hg_controls.py > logs/controls.out 2>&1 && touch CONTROLS_DONE; echo "$(ts) CONTROLS"
python scripts/hg_stackinv.py run > logs/stackinv.out 2>&1 && touch STACKINV_DONE; echo "$(ts) STACKINV"
echo "$(ts) queue2 end"
