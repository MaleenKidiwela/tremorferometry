#!/bin/bash
# Amendment Y10b: re-run inj + summary + R + psum for prodrq2018, prodf2019, prodf2020 with the W2b-masked relay sources.
# Waits for Y9_STEP4_DONE (CPU cap: forward chain 6 + this 2). Old outputs copied to tables/preY10b/ and
# data/<run><Y>/inj_v2_preY10b.parquet first. Marker Y10B_DONE; FAIL -> logs/q_main.log. Resume: re-run.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400 CAT_NW=2 RELAY2018_W2B=1
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y10b $*" >> $L; }
until [ -f Y9_STEP4_DONE ]; do pgrep -f "scripts/y9_nb.sh" > /dev/null || { log "FAIL y9_nb not running and no Y9_STEP4_DONE"; exit 1; }; sleep 60; done
mkdir -p tables/preY10b
for RY in prodrq2018 prodf2019 prodf2020; do
  R=${RY%????}; Y=${RY: -4}; PD=data/$RY; SD=$D/data/scan_$Y
  [ -f $PD/STEP_Y10b_DONE ] && continue
  cp -n tables/${RY}_* tables/preY10b/ 2>/dev/null; cp -n $PD/inj_v2.parquet $PD/inj_v2_preY10b.parquet
  unset PROD_NAME SCAN_REG_EXTRA
  export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=2018 \
         SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  [ $R = prodf ] && export PROD_NAME=prodf SCAN_REG_EXTRA=,UW.OTR,UW.HDW
  log "start $RY inj"
  python scripts/cat_pilot.py inj >> logs/y10b_$RY.out 2>&1 || { log "FAIL $RY inj"; exit 1; }
  python scripts/cat_pilot.py summary >> logs/y10b_$RY.out 2>&1 || { log "FAIL $RY summary"; exit 1; }
  python scripts/cat_r.py >> logs/y10b_$RY.out 2>&1 || { log "FAIL $RY R"; exit 1; }
  python scripts/cat_prod.py summary $Y >> logs/y10b_$RY.out 2>&1 || { log "FAIL $RY psum"; exit 1; }
  touch $PD/STEP_Y10b_DONE; log "end $RY"
done
TZ=America/Los_Angeles date > Y10B_DONE; log "done"
