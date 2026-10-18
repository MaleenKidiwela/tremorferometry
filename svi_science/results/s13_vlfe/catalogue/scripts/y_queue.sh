#!/bin/bash
# Amendment Y1: 2018 under W2b (2011 bank), 2018 relay bank re-cut under W2b, re-run with 2011 + new relay if accepted.
# Marker Y_2018W2B_DONE (writes tables/y_relay2018_accepted.txt = 1/0). FAIL lines -> logs/q_main.log. CAT_NW=4.
# Resume: re-run (STEP markers in data/prodq2018/, data/prodq2018/relay/, data/prodrq2018/).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400 CAT_NW=${CAT_NW:-4}
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y $*" >> $L; }
st(){ local d=$1 m=$2 tag=$3; shift 3; [ -f $d/STEP_${m}_DONE ] && return 0; log "start $tag $m"
  "$@" >> logs/$tag.out 2>&1 || { log "FAIL $tag $m (see logs/$tag.out)"; exit 1; }; touch $d/STEP_${m}_DONE; log "end $tag $m"; }
Y=2018; SD=$D/data/scan_$Y
export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 PROD_QC=1 SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
unset PROD_RELAY PROD_NAME PROD_BANK_YEARS
PD=data/prodq$Y; mkdir -p $PD
for m in bank mf detect battery inj_list inj summary; do st $PD $m prodq$Y python scripts/cat_pilot.py $m; done
st $PD R prodq$Y python scripts/cat_r.py
st $PD psum prodq$Y python scripts/cat_prod.py summary $Y
RD=$PD/relay; mkdir -p $RD
st $RD rbank relayq$Y python scripts/cat_relay.py bank
ACC=0
if [ -f tables/prodq${Y}_relay_bank.csv ]; then
  for m in mf loo battery verdict; do st $RD r$m relayq$Y python scripts/cat_relay.py $m; done
  ACC=$(python -c "import json; print(int(json.load(open('tables/prodq${Y}_relay_verdict.json'))['accepted']))")
fi
echo $ACC > tables/y_relay2018_accepted.txt; log "2018 W2b relay accepted=$ACC"
if [ "$ACC" = "1" ]; then
  [ -f data/inj/sources_relay2018_preW2b.npz ] || { cp data/inj/sources.npz data/inj/sources_relay2018_preW2b.npz; cp tables/inj_sources.csv tables/inj_sources_relay2018_preW2b.csv; }
  st $RD rsources relayq$Y python scripts/cat_inj.py relay_sources 2018 prodq
  export PROD_RELAY=1 PROD_BANK_YEARS=2018 RELAY2018_W2B=1
  PR=data/prodrq$Y; mkdir -p $PR
  for f in $PR/STEP_*_DONE; do [ -f $f ] && mv $f $f.preY; done 2>/dev/null
  for m in bank mf detect battery inj_list inj summary; do st $PR $m prodrq$Y python scripts/cat_pilot.py $m; done
  st $PR R prodrq$Y python scripts/cat_r.py
  st $PR psum prodrq$Y python scripts/cat_prod.py summary $Y
fi
TZ=America/Los_Angeles date > Y_2018W2B_DONE; log "2018 W2b done"
