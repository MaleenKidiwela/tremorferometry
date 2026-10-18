#!/bin/bash
# S13 pilot queue (PLAN.md §3). Skips stages with DONE markers. Run: nohup bash scripts/s13_pilot_queue.sh &
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1
[ -f STACK_pilot_2011-08_DONE ] || python scripts/s13_stack.py 2011-08 2011-08-12 2011-09-03 CN.PGC,CN.SNB,CN.VGZ,CN.NLLB,CN.LZB,CN.PFB,CN.YOUB,UW.LRIV --fams 56,10,43,0,3 --ndraw 500 --tag pilot > logs/stack_pilot_2011.out 2>&1 &
[ -f STACK_pilot_2012-08_DONE ] || python scripts/s13_stack.py 2012-08 2012-08-31 2012-09-28 CN.PGC,CN.SNB,CN.VGZ,CN.NLLB,CN.LZB,CN.YOUB,UW.LRIV --fams 56,10,43,0,3 --ndraw 500 --tag pilot > logs/stack_pilot_2012.out 2>&1 &
wait
[ -f STACK_pilot_2011-08_DONE ] && [ -f STACK_pilot_2012-08_DONE ] && python scripts/s13_pilot_figs.py > logs/pilot_figs.out 2>&1 && touch PILOT_QUEUE_DONE
