#!/bin/bash
# Amendment X: forward relay chain 2019-2026. Bank starts at 2011 + 2018 relay (data/fwd_bank_years.txt), W2b QC on,
# scan v2, weighted, n = 800 per Mw per period (C0/C1 + all accepted relay clusters). Per year: production (prodf) ->
# relay bank (accept:Y7 LOO>=85%; scrambled col only) -> re-run with the extended bank (prodfr). Amendment Y10: no stop
# rule; low-completeness years are flagged and catalogued. Markers PRODF_<Y>_DONE, RELAYF_<Y>_DONE, FWD_DONE. FAIL lines -> logs/q_main.log. CAT_NW=4 (shares CPU).
# Resume: re-run (STEP markers in data/prodf<Y>/, data/prodf<Y>/relay/, data/prodfr<Y>/).
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 NINJ=400 CAT_NW=${CAT_NW:-4} SCAN_REG_EXTRA=,UW.OTR,UW.HDW RELAY2018_W2B=1   # Amendment Y1: W2b 2018 relay bank
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
L=logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') fwd $*" >> $L; }
st(){ local d=$1 m=$2 tag=$3; shift 3; [ -f $d/STEP_${m}_DONE ] && return 0; log "start $tag $m"
  "$@" >> logs/$tag.out 2>&1 || { log "FAIL $tag $m (see logs/$tag.out)"; exit 1; }; touch $d/STEP_${m}_DONE; log "end $tag $m"; }
for Y in 2019 2020 2021 2022 2023 2024 2025 2026; do
  [ -f FWD_DONE ] && break
  if [ -f RELAYF_${Y}_DONE ]; then continue; fi
  SD=$D/data/scan_$Y
  unset SCAN_GRID SCAN_GAIN_FALLBACK SCAN_RUN SCAN_BP SCAN_REFBP
  if [ ! -f $SD/SCAN${Y}_DONE ]; then log "start scan year $Y"; bash $D/scripts/scan_year_queue.sh $Y || { log "FAIL scan year $Y"; exit 1; }; log "end scan year $Y"; fi
  BY=$(cat data/fwd_bank_years.txt)
  PD=data/prodf$Y; mkdir -p $PD
  export SCAN_YEAR=$Y PROD_YEAR=$Y PILOT_SFX=_v2 PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=$BY PROD_NAME=prodf
  unset SCAN_GRID SCAN_GAIN_FALLBACK SCAN_RUN SCAN_BP SCAN_REFBP
  MISS=$(python scripts/cat_prod.py missing $Y)
  if [ -n "$MISS" ] && [ ! -f $PD/STEP_prepmiss_DONE ]; then
    log "start $Y prep missing $MISS"
    (cd $D && S13_DISP=$D/data/disp1sps_v3 python scripts/scan_prep_ref.py $MISS $Y-01-01 $Y-12-31 --workers $CAT_NW > /dev/null 2>&1) || { log "FAIL $Y prepmiss"; exit 1; }
    touch $PD/STEP_prepmiss_DONE; log "end $Y prepmiss"
  fi
  st $PD stations prodf$Y python scripts/cat_prod.py stations $Y
  st $PD cont prodf$Y python scripts/cat_pilot.py cont
  st $PD diag prodf$Y python scripts/cat_diag.py $Y
  st $PD qcx prodf$Y python scripts/qc_exclusions.py $Y
  st $PD v2setup prodf$Y python scripts/cat_scanv2.py setup $Y
  st $PD v2detect prodf$Y bash -c "cd $D && SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 python scripts/scan_detect.py real_v2 --bp $SD/bp_${Y}v2 --ndraw 0 --shift 0 && test -f $SD/runs/real_v2/DONE"
  export SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
  for m in bank mf detect battery inj_list inj summary; do st $PD $m prodf$Y python scripts/cat_pilot.py $m; done
  st $PD R prodf$Y python scripts/cat_r.py
  st $PD psum prodf$Y python scripts/cat_prod.py summary $Y
  [ -f PRODF_${Y}_DONE ] || { TZ=America/Los_Angeles date > PRODF_${Y}_DONE; log "prodf $Y done (bank years $BY)"; }
  RD=$PD/relay; mkdir -p $RD
  st $RD rbank relayf$Y python scripts/cat_relay.py bank
  ACC=0
  if [ -f tables/prodf${Y}_relay_bank.csv ]; then
    for m in mf loo battery verdict; do st $RD r$m relayf$Y python scripts/cat_relay.py $m; done
    ACC=$(python -c "import json; print(int(json.load(open('tables/prodf${Y}_relay_verdict.json'))['accepted']))")
  fi
  FINAL=prodf
  if [ "$ACC" = "1" ]; then
    BY=$Y; echo $BY > data/fwd_bank_years.txt     # Amendment Y1: bank capped at 2011 + latest accepted relay
    st $RD rsources relayf$Y python scripts/cat_inj.py relay_sources $Y prodf
    PR=data/prodfr$Y; mkdir -p $PR
    export PROD_NAME=prodfr PROD_BANK_YEARS=$BY
    for m in bank mf detect battery inj_list inj summary; do st $PR $m prodfr$Y python scripts/cat_pilot.py $m; done
    st $PR R prodfr$Y python scripts/cat_r.py
    st $PR psum prodfr$Y python scripts/cat_prod.py summary $Y
    FINAL=prodfr
  fi
  FAILED=$(python -c "
import json; s = json.load(open('tables/${FINAL}${Y}_summary_v2.json')); c = s['completeness']['at_3p8']
v = [c.get('ETS', 0), c.get('between', 0)]; print(int(all((x if isinstance(x, float) else 0) < 0.5 for x in v) and '$ACC' != '1'))")
  F=$(cat data/fwd_fails.txt); if [ "$FAILED" = "1" ]; then F=$((F+1)); else F=0; fi; echo $F > data/fwd_fails.txt
  TZ=America/Los_Angeles date > RELAYF_${Y}_DONE; log "year $Y done: relay accepted=$ACC final=$FINAL failed=$FAILED consecutive=$F bank=$(cat data/fwd_bank_years.txt)"
  [ $F -ge 2 ] && log "FLAG (Amendment Y10, no stop rule): $F consecutive low-completeness years (last $Y); continuing"
done
TZ=America/Los_Angeles date > FWD_DONE; log "forward chain done"
