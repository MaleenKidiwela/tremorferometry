#!/bin/bash
# Amendment Y13 (<= 2 workers): after Y9B_DONE -> cut event windows, cluster, then per year 2010-2026 MF + battery, summary.
# 2021+ years wait for the prodfy run's MF (noise templates) of that year. Markers data/y13/STEP_*_DONE, Y13_DONE.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 CAT_NW=2 RELAY2018_W2B=1
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log; Q=data/y13; mkdir -p $Q
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y13 $*" >> $L; }
runenv(){ local Y=$1; local R
  case $Y in 2013|2016|2017|2018) R=prodq;; 2019|2020) R=prodf;; 2021|2022|2023|2024|2025|2026) R=prodfy;; *) R=prod;; esac
  unset PROD_QC PROD_RELAY PROD_BANK_YEARS PROD_NAME RELAY_RUN_PREFIX SCAN_REFBP SCAN_REG_EXTRA
  export PROD_YEAR=$Y PILOT_SFX=_v2
  [ $R != prod ] && export PROD_QC=1
  if [ $R = prodf ] || [ $R = prodfy ]; then export PROD_RELAY=1 PROD_BANK_YEARS=2018 PROD_NAME=$R SCAN_REG_EXTRA=,UW.OTR,UW.HDW; fi
  [ $R = prodfy ] && export RELAY_RUN_PREFIX=prodfy
  local SD=$D/data/scan_$Y; [ $Y = 2011 ] && SD=$D/data/scan
  { [ $Y = 2011 ] || [ $Y = 2012 ]; } && export SCAN_REFBP=$SD/bp2_ref_v2
  export SCAN_YEAR=$Y SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  RUNNAME=$R; }
st(){ local m=$1; shift; [ -f $Q/STEP_${m}_DONE ] && return 0; log "start $m"; "$@" >> logs/y13.out 2>&1 || { log "FAIL $m (logs/y13.out)"; exit 1; }; touch $Q/STEP_${m}_DONE; log "end $m"; }
until [ -f Y9B_DONE ]; do pgrep -f "scripts/yq_master.sh" > /dev/null || { [ -f Y9B_DONE ] || { log "FAIL master queue not running and no Y9B_DONE"; exit 1; }; }; sleep 120; done
for Y in 2014 2016 2017 2020; do ( runenv $Y; st cut$Y python scripts/cat_y13.py cut $Y ) || exit 1; done
st cluster bash -c "python scripts/cat_y13.py cluster > logs/y13_cluster.out 2>&1"
for Y in 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026; do
  if [ $Y -ge 2021 ]; then
    until [ -f data/prodfy$Y/STEP_battery_DONE ] && [ -f data/prodfy$Y/STEP_mf_DONE ]; do pgrep -f "scripts/prodfy_queue.sh" > /dev/null || { log "FAIL prodfy not running and prodfy$Y not ready"; exit 1; }; sleep 300; done
  fi
  ( runenv $Y; st mf$Y python scripts/cat_y13.py mf $Y; st bat$Y python scripts/cat_y13.py battery $Y ) || exit 1
done
st summary bash -c "python scripts/cat_y13.py summary > logs/y13_summary.out 2>&1"
TZ=America/Los_Angeles date > Y13_DONE; log "done"
