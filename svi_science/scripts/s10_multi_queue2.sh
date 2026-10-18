#!/bin/bash
# S10 multi queue, Revision 1 re-run (control episode stacks added; post-hoc battery). DONE markers in multi/.
set -e
cd /home/jovyan/tremorferometry/svi_science/scripts
P=/home/jovyan/tremorferometry/svi_science/results/s10_dvv/multi
[ -f $P/MEASURE_DONE ] || python s10_multi_measure.py --workers 12
python s10_multi_analysis.py
python s10_multi_figures.py
touch $P/QUEUE2_DONE
