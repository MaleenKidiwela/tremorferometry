#!/bin/bash
# Amendment Y14: re-score every final run 2010-2023 under f' (summary, R, psum); relay vetting counted afterwards from events_v2 (relay banks NOT rebuilt). No MF/battery/inj re-runs.
# f (pre-Y14) outputs copied to tables/preY14/ and data/<run>/rows_scored_v2_preY14.parquet first. Marker Y14_RESCORE_DONE.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D/catalogue
export OMP_NUM_THREADS=1 CAT_NW=1 RELAY2018_W2B=1
L=logs/q_main.log; log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') Y14 $*" >> $L; }
mkdir -p tables/preY14
for RY in prod2010 prod2011 prod2012 prodq2013 prod2014 prod2015 prodq2016 prodq2017 prodq2018 prodrq2018 prodf2019 prodf2020 prodfy2021 prodfy2022 prodfy2023; do
  R=${RY%????}; Y=${RY: -4}; PD=data/$RY; [ -f $PD/STEP_Y14_DONE ] && continue
  cp -n tables/${RY}_* tables/preY14/ 2>/dev/null; [ -f $PD/rows_scored_v2.parquet ] && cp -n $PD/rows_scored_v2.parquet $PD/rows_scored_v2_preY14.parquet
  ( unset PROD_QC PROD_RELAY PROD_BANK_YEARS PROD_NAME RELAY_RUN_PREFIX SCAN_REFBP SCAN_REG_EXTRA
    SD=$D/data/scan_$Y; [ $Y = 2011 ] && SD=$D/data/scan
    export PROD_YEAR=$Y PILOT_SFX=_v2 SCAN_YEAR=$Y SCAN_GRID=$SD/templates_grid_v2.npz SCAN_GAIN_FALLBACK=1 SCAN_RUN=real_v2 SCAN_BP=$SD/bp_${Y}v2
    { [ $Y = 2011 ] || [ $Y = 2012 ]; } && export SCAN_REFBP=$SD/bp2_ref_v2
    case $R in prodq) export PROD_QC=1;; prodrq) export PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=2018;;
      prodf|prodfy) export PROD_QC=1 PROD_RELAY=1 PROD_BANK_YEARS=2018 PROD_NAME=$R SCAN_REG_EXTRA=,UW.OTR,UW.HDW;; esac
    [ $R = prodfy ] && export RELAY_RUN_PREFIX=prodfy
    log "start $RY"
    python scripts/cat_pilot.py summary >> logs/y14_rescore.out 2>&1 || exit 1
    python scripts/cat_r.py >> logs/y14_rescore.out 2>&1 || exit 1
    python scripts/cat_prod.py summary $Y >> logs/y14_rescore.out 2>&1 || exit 1
    ) || { log "FAIL $RY"; exit 1; }
  touch $PD/STEP_Y14_DONE; log "end $RY"
done
cp -n tables/y9_ladder_table.csv tables/preY14/y9_ladder_table.csv; cp -n tables/catalogue_ladder_2010_2026.csv tables/preY14/catalogue_ladder_2010_2026.csv
cp -n tables/y9_ladder_completeness.csv tables/preY14/y9_ladder_completeness.csv; cp -n tables/y9_checks.json tables/preY14/y9_checks.json
LADDER_EXTRA=2019,2020,2021 LADDER_RUNS=2021=prodfy python scripts/cat_y9.py ladder > logs/y14_ladder.out 2>&1 || { log "FAIL ladder"; exit 1; }
TZ=America/Los_Angeles date > Y14_RESCORE_DONE; log "rescore done"
