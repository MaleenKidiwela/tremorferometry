#!/bin/bash
# Amendment Y11 + Y9b compute queue (<= 2 workers). Waits for Y10B_DONE. Steps: (1) Mw 3.8 injections of prodf2021 with
# loso_dropped/max_share_ch -> Y11_INJ_DONE; (2) ladder extended to 2019-2021 -> Y11_LADDER_DONE; (3) Y9b matched noise banks
# 2014-2017 (scripts/y9b_nb.sh) -> Y9B_DONE. FAIL -> logs/q_main.log. Resume: re-run (done steps skipped).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 CAT_NW=2 NINJ=400 RELAY2018_W2B=1
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y11q $*" >> $L; }
until [ -f Y10B_DONE ]; do pgrep -f "scripts/y10b_queue.sh" > /dev/null || { log "FAIL y10b not running and no Y10B_DONE"; exit 1; }; sleep 60; done
if [ ! -f Y11_INJ_DONE ]; then
  ( Y=2021; SD=$D/data/scan_$Y
    export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=2018 PROD_NAME=prodf SCAN_REG_EXTRA=,UW.OTR,UW.HDW \
           SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
    log "start inj38 2021"; python scripts/cat_y11.py inj38 >> logs/y11.out 2>&1 ) || { log "FAIL inj38"; exit 1; }
  TZ=America/Los_Angeles date > Y11_INJ_DONE; log "inj38 done"
fi
[ "$1" = "--inj-only" ] && exit 0
if [ ! -f Y11_LADDER_DONE ]; then
  for Y in 2019 2020 2021; do
    [ -f data/ranked/R2BAT_${Y}_DONE ] && continue
    ( SD=$D/data/scan_$Y
      RN=prodf; [ $Y = 2021 ] && [ -f Y12_PILOT2021_DONE ] && RN=prodfy
      python scripts/cat_ranked.py year $Y $RN >> logs/y11_ladder.out 2>&1 || exit 1
      export SCAN_YEAR=$Y SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2 SCAN_REG_EXTRA=,UW.OTR,UW.HDW
      python scripts/cat_ranked2.py battery $Y >> logs/y11_ladder.out 2>&1 ) || { log "FAIL ladder battery $Y"; exit 1; }
    touch data/ranked/R2BAT_${Y}_DONE; log "ladder battery $Y done"
  done
  export LADDER_EXTRA=2019,2020,2021
  [ -f Y12_PILOT2021_DONE ] && export LADDER_RUNS=2021=prodfy
  python scripts/cat_ranked3.py score >> logs/y11_ladder.out 2>&1 || { log "FAIL ladder score"; exit 1; }
  for Y in 2019 2020 2021; do python scripts/cat_y9.py emerge $Y >> logs/y11_ladder.out 2>&1 || { log "FAIL emerge $Y"; exit 1; }; done
  python scripts/cat_y9.py ladder > logs/y11_ladder_table.out 2>&1 || { log "FAIL ladder"; exit 1; }
  unset LADDER_EXTRA
  TZ=America/Los_Angeles date > Y11_LADDER_DONE; log "ladder 2019-2021 done"
fi
[ -f Y9B_DONE ] || bash scripts/y9b_nb.sh || { log "FAIL y9b"; exit 1; }
log "queue done"
