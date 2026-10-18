#!/bin/bash
# S14 full run, tests only (stores done 2:13 pm PT). Markers TESTS_<REGION>_DONE, QUEUE_TESTS_DONE.
set -u
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /home/jovyan/tremorferometry/svi_science/scripts
L=/home/jovyan/tremorferometry/svi_science/results/s14_cc/logs/queue_tests.log
echo "$(TZ=America/Los_Angeles date) start" >> $L
S11_REGION=SVI   python s14_tests.py --nperm 1000 >> $L 2>&1; echo "tests SVI rc $?" >> $L
S11_REGION=WA_NW python s14_tests.py --nperm 1000 >> $L 2>&1; echo "tests WA rc $?" >> $L
python s14_figures.py SVI_pilot SVI WA_NW >> $L 2>&1; echo "figures rc $?" >> $L
echo "$(TZ=America/Los_Angeles date) end" >> $L
touch /home/jovyan/tremorferometry/svi_science/results/s14_cc/QUEUE_TESTS_DONE
