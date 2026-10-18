#!/bin/bash
# Amendment Y12 pilot 2021: prep + QC the new stations, add QC-passing ones to the 2021 v2 scan, re-run battery/inj/summary/R
# as prodfy2021 (MF reused), then f vs f' re-score. Markers data/prodfy2021/STEP_*_DONE, Y12_PILOT2021_DONE. FAIL -> q_main.log.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 CAT_NW=2 NINJ=400 RELAY2018_W2B=1
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log; Y=2021; SD=$D/data/scan_$Y; PD=data/prodfy$Y; mkdir -p $PD
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y12p $*" >> $L; }
st(){ local m=$1; shift; [ -f $PD/STEP_${m}_DONE ] && return 0; log "start $m"; "$@" >> logs/y12_pilot.out 2>&1 || { log "FAIL $m (logs/y12_pilot.out)"; exit 1; }; touch $PD/STEP_${m}_DONE; log "end $m"; }
until [ -f FETCH_Y12_DONE ] || grep -q "CN.SYMB 2022-01-0[3-9]" logs/fetch_y12.log; do pgrep -f fetch_y12.py > /dev/null || { [ -f FETCH_Y12_DONE ] || { log "FAIL fetch not running"; exit 1; }; }; sleep 60; done
CAND=CN.SYMB,UW.OSQM,UW.RNWY,UW.SLDQ
st y12prep bash -c "cd $D && S13_DISP=$D/data/disp1sps_v3 python scripts/scan_prep_ref.py $CAND $Y-01-01 $Y-12-31 --workers 2"
st y12qc bash -c "python scripts/cat_y12.py qc $Y $CAND > logs/y12_qc_$Y.out 2>&1"
PASS=$(grep '^PASS' logs/y12_qc_$Y.out | awk '{print $2}'); log "QC pass: $PASS"
[ -n "$PASS" ] || { log "FAIL no station passes QC"; exit 1; }
st y12bpref bash -c "cd $D && SCAN_YEAR=$Y python scripts/scan_bp_ref.py $SD/bp2_ref --stations $PASS --workers 2"
[ -f $SD/templates_grid_v2.preY12.npz ] || cp $SD/templates_grid_v2.npz $SD/templates_grid_v2.preY12.npz
[ -d $SD/runs/real_v2_preY12 ] || mv $SD/runs/real_v2 $SD/runs/real_v2_preY12
st y12setup env SCAN_YEAR=$Y Y12_EXTRA=$PASS python scripts/cat_scanv2.py setup $Y
st y12detect bash -c "cd $D && SCAN_YEAR=$Y SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 python scripts/scan_detect.py real_v2 --bp $SD/bp_${Y}v2 --ndraw 0 --shift 0 && test -f $SD/runs/real_v2/DONE"
st y12sources python scripts/cat_y12.py sources
[ -f $PD/bank.npz ] || python scripts/cat_y12.py mkrun $Y >> logs/y12_pilot.out 2>&1
export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=2018 PROD_NAME=prodfy SCAN_REG_EXTRA=,UW.OTR,UW.HDW \
       SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2 INJ_SOURCES=$D/catalogue/data/inj/sources_y12.npz
for m in battery inj summary; do st $m python scripts/cat_pilot.py $m; done
st R python scripts/cat_r.py
st injx0p5 env INJ_SOURCES=$D/catalogue/data/inj/sources_y12_x0p5.npz python scripts/cat_y12.py injvar x0p5
st injx2 env INJ_SOURCES=$D/catalogue/data/inj/sources_y12_x2.npz python scripts/cat_y12.py injvar x2
st psum python scripts/cat_prod.py summary $Y
st rescore bash -c "python scripts/cat_y12.py rescore $Y > logs/y12_rescore_$Y.out 2>&1"
TZ=America/Los_Angeles date > Y12_PILOT2021_DONE; log "pilot 2021 done"
