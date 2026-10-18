#!/bin/bash
# S10 multi queue (rules: svi_science/results/s10_dvv/multi/PLAN.md). Each step writes a DONE marker in multi/.
# Windows (s10_multi_windows.py) was launched separately; this waits for its WINDOWS_DONE marker.
# nohup bash svi_science/scripts/s10_multi_queue.sh > svi_science/results/s10_dvv/multi/logs/queue.out 2>&1 &
set -e
cd /home/jovyan/tremorferometry/svi_science/scripts
P=/home/jovyan/tremorferometry/svi_science/results/s10_dvv/multi
until [ -f $P/WINDOWS_DONE ]; do sleep 30; done
[ -f $P/MEASURE_DONE ] || python s10_multi_measure.py --workers 12
[ -f $P/ANALYSIS_DONE ] || python s10_multi_analysis.py
python s10_multi_figures.py
touch $P/QUEUE_DONE
