#!/bin/bash
# S12c pilot preprocessing: daily ACFs 2018-2021, PGC first, then VGZ, SNB. Resumable (skips finished station-years).
cd /home/jovyan/tremorferometry
WS=svi_science/results/s12_lvl/acf_pgc
for S in CN.PGC CN.VGZ CN.SNB; do
  python $WS/scripts/acf01_daily.py --stations $S --years 2018-2021 --workers 8 || echo "FAIL $S"
done
date > $WS/PILOT_PRE_DONE
