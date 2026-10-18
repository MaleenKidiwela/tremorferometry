#!/bin/bash
# Amendment S6: injection re-run with tau = D/2.33. Marker SCAN_INJ2_DONE.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D; export OMP_NUM_THREADS=1
python scripts/scan_inject.py > logs/scan_inj2.log 2>&1 && \
python scripts/scan_bp.py $D/data/scan/bp_2011inj2 --inject $D/data/scan/inject_2011.parquet --workers 8 >> logs/scan_inj2.log 2>&1 && \
python scripts/scan_detect.py inj2 --bp data/scan/bp_2011inj2 --ndraw 0 >> logs/scan_inj2.log 2>&1 && \
python scripts/scan_columns.py inj2 --bp data/scan/bp_2011inj2 --workers 8 >> logs/scan_inj2.log 2>&1 && touch SCAN_INJ2_DONE
