#!/bin/bash
# S14 full run: store columns + tests, SVI then WA-NW. <= 4 workers, OMP 1. Markers in svi_science/results/s14_cc/.
set -u
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /home/jovyan/tremorferometry/svi_science/scripts
L=/home/jovyan/tremorferometry/svi_science/results/s14_cc/logs/queue.log
echo "$(TZ=America/Los_Angeles date) start" >> $L
S11_REGION=SVI   python s14_store.py --workers 4 >> $L 2>&1; echo "store SVI rc $?" >> $L
S11_REGION=SVI   python s14_tests.py --nperm 1000 >> $L 2>&1; echo "tests SVI rc $?" >> $L
S11_REGION=WA_NW python s14_store.py --workers 4 >> $L 2>&1; echo "store WA rc $?" >> $L
S11_REGION=WA_NW python s14_tests.py --nperm 1000 >> $L 2>&1; echo "tests WA rc $?" >> $L
echo "$(TZ=America/Los_Angeles date) end" >> $L
touch /home/jovyan/tremorferometry/svi_science/results/s14_cc/QUEUE_DONE
