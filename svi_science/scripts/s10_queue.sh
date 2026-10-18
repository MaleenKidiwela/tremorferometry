#!/bin/bash
# S10 pilot queue (rules: svi_science/results/s10_dvv/pilot/PLAN.md). Each step writes a DONE marker in the pilot folder.
# nohup bash svi_science/scripts/s10_queue.sh > svi_science/results/s10_dvv/pilot/logs/queue.out 2>&1 &
set -e
cd /home/jovyan/tremorferometry/svi_science/scripts
P=/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot
python s10_inventory.py && touch $P/INVENTORY_DONE
python s10_windows.py verify
python s10_windows.py run --workers 12
python s10_windows.py run --c1f --workers 12
python s10_measure.py --workers 10
python s10_gates.py
python s10_kelb_merge.py
python s10_measure.py --stations PO.KLNBm --rs $P/tables/row_station_kelbm.parquet --tag _kelbm --workers 1
python s10_figures.py && touch $P/FIGURES_DONE
