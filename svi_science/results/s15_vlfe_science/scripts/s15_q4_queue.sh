#!/bin/bash
# S15 Q4 GF-source injections (PLAN Amendment 5). usage: s15_q4_queue.sh [sources|all]. Marker S15_Q4_DONE (all).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; S15=/home/jovyan/tremorferometry/svi_science/results/s15_vlfe_science
cd $D/catalogue
export OMP_NUM_THREADS=1 CAT_NW=4
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
LOG=$S15/tables/q4_queue.log; ts(){ TZ=America/Los_Angeles date '+%F %H:%M PT'; }
MODE=${1:-all}
env12(){ unset PROD_QC PROD_RELAY PROD_BANK_YEARS PROD_NAME RELAY2018_W2B RELAY_RUN_PREFIX SCAN_REG_EXTRA INJ_SOURCES
  SD=$D/data/scan_2012; export SCAN_YEAR=2012 PROD_YEAR=2012 PILOT_SFX=_v2 SCAN_REFBP=$SD/bp2_ref_v2 \
    SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_2012v2; }
env24(){ unset SCAN_REFBP
  SD=$D/data/scan_2024; export SCAN_YEAR=2024 PROD_YEAR=2024 PILOT_SFX=_v2 PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=$(cat data/fwdy_bank_years.txt) \
    PROD_NAME=prodfy RELAY2018_W2B=1 RELAY_RUN_PREFIX=prodfy SCAN_REG_EXTRA=,UW.OTR,UW.HDW SCAN_GRID=$SD/templates_grid_v2.npz \
    SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_2024v2 INJ_SOURCES=$D/catalogue/data/inj/sources_y12.npz; }
for Y in 2012 2024; do
  if [ $Y = 2012 ]; then env12; else env24; fi
  echo "$(ts) $Y sources" >> $LOG
  python $S15/scripts/s15_q4_inj.py sources $Y > $S15/tables/q4_sources_$Y.out 2>&1 || { echo "$(ts) $Y sources FAILED" >> $LOG; exit 1; }
  [ $MODE = sources ] && continue
  echo "$(ts) $Y inj start" >> $LOG
  python $S15/scripts/s15_q4_inj.py inj $Y > $S15/tables/q4_inj_$Y.out 2>&1 || { echo "$(ts) $Y inj FAILED" >> $LOG; exit 1; }
  echo "$(ts) $Y inj done" >> $LOG
done
[ $MODE = all ] && touch $S15/S15_Q4_DONE
