#!/bin/bash
# S10 multi: analysis + figures after Revision 1 (measure already done). DONE markers in multi/.
set -e
cd /home/jovyan/tremorferometry/svi_science/scripts
P=/home/jovyan/tremorferometry/svi_science/results/s10_dvv/multi
python s10_multi_analysis.py
python s10_multi_figures.py
touch $P/QUEUE3_DONE
