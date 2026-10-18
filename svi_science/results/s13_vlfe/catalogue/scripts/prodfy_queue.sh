#!/bin/bash
# Amendment Y12: forward chain 2022-2026 on the expanded scan network (run name prodfy; extended-bank re-run prodfyr).
# MF reused from prodf<Y> where the pre-Y12 chain finished it; Y12 stations prepped/QC'd/bp_ref'd per year and added to the v2
# scan (Y12_EXTRA); injections with sources_y12 (GF sources on new stations). Y7/Y9/Y10/Y10b/W2b rules unchanged; no stop rule.
# Markers data/prodfy<Y>/STEP_*_DONE, PRODFY_<Y>_DONE, FWDY_DONE. FAIL -> logs/q_main.log. Resume: re-run.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400 CAT_NW=${CAT_NW:-6} RELAY2018_W2B=1 RELAY_RUN_PREFIX=prodfy
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') fwdy $*" >> $L; }
st(){ local d=$1 m=$2 tag=$3; shift 3; [ -f $d/STEP_${m}_DONE ] && return 0; log "start $tag $m"
  "$@" >> logs/$tag.out 2>&1 || { log "FAIL $tag $m (see logs/$tag.out)"; exit 1; }; touch $d/STEP_${m}_DONE; log "end $tag $m"; }
[ -f data/fwdy_bank_years.txt ] || echo 2018 > data/fwdy_bank_years.txt
for Y in 2022 2023 2024 2025 2026; do
  [ -f PRODFY_${Y}_DONE ] && continue
  SD=$D/data/scan_$Y; PD=data/prodfy$Y; mkdir -p $PD; T=prodfy$Y
  unset SCAN_GRID SCAN_GAIN_FALLBACK SCAN_RUN SCAN_BP SCAN_REFBP PROD_NAME
  export SCAN_REG_EXTRA=,UW.OTR,UW.HDW
  if [ ! -f $SD/SCAN${Y}_DONE ]; then log "start scan year $Y"; bash $D/scripts/scan_year_queue.sh $Y || { log "FAIL scan year $Y"; exit 1; }; log "end scan year $Y"; fi
  if [ $Y = 2023 ]; then   # the pre-Y12 chain was stopped during prodf2023 inj_list: let that orphan step finish, then mark it
    while pgrep -f "cat_pilot.py inj_list" > /dev/null; do sleep 30; done
    [ -f data/prodf2023/inj_list.parquet ] && touch data/prodf2023/STEP_inj_list_DONE
  fi
  [ -f data/prodf$Y/STEP_detect_DONE ] && [ ! -f $PD/bank.npz ] && python scripts/cat_y12.py mkrun $Y >> logs/$T.out 2>&1
  BY=$(cat data/fwdy_bank_years.txt)
  export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=$BY PROD_NAME=prodfy
  MISS=$(python scripts/cat_prod.py missing $Y)
  if [ -n "$MISS" ] && [ ! -f $PD/STEP_prepmiss_DONE ] && [ ! -f data/prodf$Y/STEP_prepmiss_DONE ]; then
    log "start $Y prep missing $MISS"
    (cd $D && S13_DISP=$D/data/disp1sps_v3 python scripts/scan_prep_ref.py $MISS $Y-01-01 $Y-12-31 --workers $CAT_NW > /dev/null 2>&1) || { log "FAIL $Y prepmiss"; exit 1; }
    touch $PD/STEP_prepmiss_DONE; log "end $Y prepmiss"
  fi
  st $PD stations $T python scripts/cat_prod.py stations $Y
  st $PD cont $T python scripts/cat_pilot.py cont
  st $PD diag $T python scripts/cat_diag.py $Y
  st $PD qcx $T python scripts/qc_exclusions.py $Y
  # Y12 stations
  CAND=CN.SYMB,UW.OSQM,UW.RNWY,UW.SLDQ; [ $Y -ge 2023 ] && CAND=$CAND,UW.BUCKS; [ $Y -ge 2024 ] && CAND=$CAND,UW.ODUC
  st $PD y12prep $T bash -c "cd $D && S13_DISP=$D/data/disp1sps_v3 python scripts/scan_prep_ref.py $CAND $Y-01-01 $Y-12-31 --workers $CAT_NW"
  st $PD y12qc $T bash -c "python scripts/cat_y12.py qc $Y $CAND > logs/y12_qc_$Y.out 2>&1"
  PASS=$(grep '^PASS' logs/y12_qc_$Y.out | awk '{print $2}'); log "$Y Y12 QC pass: $PASS"
  [ -n "$PASS" ] && st $PD y12bpref $T bash -c "cd $D && SCAN_YEAR=$Y python scripts/scan_bp_ref.py $SD/bp2_ref --stations $PASS --workers $CAT_NW"
  [ -f $SD/templates_grid_v2.npz ] && [ ! -f $SD/templates_grid_v2.preY12.npz ] && cp $SD/templates_grid_v2.npz $SD/templates_grid_v2.preY12.npz
  [ -d $SD/runs/real_v2 ] && [ ! -d $SD/runs/real_v2_preY12 ] && [ ! -f $PD/STEP_v2detect_DONE ] && mv $SD/runs/real_v2 $SD/runs/real_v2_preY12
  st $PD v2setup $T env Y12_EXTRA=$PASS python scripts/cat_scanv2.py setup $Y
  st $PD v2detect $T bash -c "cd $D && SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 python scripts/scan_detect.py real_v2 --bp $SD/bp_${Y}v2 --ndraw 0 --shift 0 && test -f $SD/runs/real_v2/DONE"
  st $PD y12sources $T python scripts/cat_y12.py sources
  export SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2 INJ_SOURCES=$D/catalogue/data/inj/sources_y12.npz
  for m in bank mf detect battery inj_list inj summary; do st $PD $m $T python scripts/cat_pilot.py $m; done
  st $PD R $T python scripts/cat_r.py
  st $PD psum $T python scripts/cat_prod.py summary $Y
  RD=$PD/relay; mkdir -p $RD
  st $RD rbank relayfy$Y python scripts/cat_relay.py bank
  ACC=0
  if [ -f tables/prodfy${Y}_relay_bank.csv ]; then
    for m in mf loo battery verdict; do st $RD r$m relayfy$Y python scripts/cat_relay.py $m; done
    ACC=$(python -c "import json; print(int(json.load(open('tables/prodfy${Y}_relay_verdict.json'))['accepted']))")
  fi
  if [ "$ACC" = "1" ]; then
    echo $Y > data/fwdy_bank_years.txt
    st $RD rsources relayfy$Y python scripts/cat_inj.py relay_sources $Y prodfy
    st $RD y12sources2 relayfy$Y python scripts/cat_y12.py sources
    PR=data/prodfyr$Y; mkdir -p $PR
    export PROD_NAME=prodfyr PROD_BANK_YEARS=$Y
    for m in bank mf detect battery inj_list inj summary; do st $PR $m prodfyr$Y python scripts/cat_pilot.py $m; done
    st $PR R prodfyr$Y python scripts/cat_r.py
    st $PR psum prodfyr$Y python scripts/cat_prod.py summary $Y
  fi
  TZ=America/Los_Angeles date > PRODFY_${Y}_DONE; log "year $Y done: relay accepted=$ACC bank=$(cat data/fwdy_bank_years.txt) Y12 stations=$PASS"
done
TZ=America/Los_Angeles date > FWDY_DONE; log "prodfy chain done"
