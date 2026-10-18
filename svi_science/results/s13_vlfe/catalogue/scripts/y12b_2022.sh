#!/bin/bash
# Amendment Y12b: re-run prodfy2022 with CN.SYMB (teleseism ratio no longer a gate). Keeps the first run as the no-SYMB variant.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400 CAT_NW=2 RELAY2018_W2B=1 RELAY_RUN_PREFIX=prodfy
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log; Y=2022; SD=$D/data/scan_$Y; PD=data/prodfy$Y; T=prodfy${Y}_y12b
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y12b $*" >> $L; }
st(){ local d=$1 m=$2; shift 2; [ -f $d/STEP_${m}_DONE ] && return 0; log "start $m"; "$@" >> logs/$T.out 2>&1 || { log "FAIL $m (logs/$T.out)"; exit 1; }; touch $d/STEP_${m}_DONE; log "end $m"; }
if [ ! -d data/prodfy2022_noSYMB ]; then
  cp -a $PD data/prodfy2022_noSYMB; mkdir -p tables/noSYMB_2022; cp -a tables/prodfy2022_* tables/noSYMB_2022/
  cp -a $SD/templates_grid_v2.npz $SD/templates_grid_v2.y12noSYMB.npz; mv $SD/runs/real_v2 $SD/runs/real_v2_y12noSYMB
  cp -a logs/y12_qc_2022.out logs/y12_qc_2022.noSYMB.out
  for m in y12qc y12bpref v2setup v2detect y12sources battery inj summary R psum; do [ -f $PD/STEP_${m}_DONE ] && mv $PD/STEP_${m}_DONE $PD/STEP_${m}_DONE.noSYMB; done
  [ -f $PD/relay/STEP_rbank_DONE ] && mv $PD/relay/STEP_rbank_DONE $PD/relay/STEP_rbank_DONE.noSYMB
  log "no-SYMB variant kept"
fi
CAND=CN.SYMB,UW.OSQM,UW.RNWY,UW.SLDQ
st $PD y12qc bash -c "python scripts/cat_y12.py qc $Y $CAND > logs/y12_qc_$Y.out 2>&1"
PASS=$(grep '^PASS' logs/y12_qc_$Y.out | awk '{print $2}'); log "QC pass (Y12b): $PASS"
st $PD y12bpref bash -c "cd $D && SCAN_YEAR=$Y python scripts/scan_bp_ref.py $SD/bp2_ref --stations $PASS --workers $CAT_NW"
export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=2018 PROD_NAME=prodfy SCAN_REG_EXTRA=,UW.OTR,UW.HDW
st $PD v2setup env Y12_EXTRA=$PASS python scripts/cat_scanv2.py setup $Y
st $PD v2detect bash -c "cd $D && SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 python scripts/scan_detect.py real_v2 --bp $SD/bp_${Y}v2 --ndraw 0 --shift 0 && test -f $SD/runs/real_v2/DONE"
st $PD y12sources python scripts/cat_y12.py sources
export SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2 INJ_SOURCES=$D/catalogue/data/inj/sources_y12.npz
for m in battery inj summary; do st $PD $m python scripts/cat_pilot.py $m; done
st $PD R python scripts/cat_r.py
st $PD psum python scripts/cat_prod.py summary $Y
st $PD/relay rbank python scripts/cat_relay.py bank
TZ=America/Los_Angeles date > Y12B_2022_DONE; log "prodfy2022 with SYMB done"
