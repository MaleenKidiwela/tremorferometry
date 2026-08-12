#!/bin/bash
# PGC dv/v — the co-located deliverable. PGC (~290 m from borehole B011) has no independent picker need:
# use B011's reliable-family detection times as the answer key to build PGC's own daily stacks, then
# causality-certify AT PGC and measure the coda stretch. Split PER ERA (BHZ 2010-2017 / HHZ 2018-2026,
# decimated 100->40) with PER-ERA stretch references (Merlin: never average the reference across the
# Aug-2017 sensor swap). Only causality-reliable families feed dv/v. Usage: nohup bash scripts/pgc_dvv_pipeline.sh &
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
mkdir -p data/pgc_era1_mf data/pgc_era2_mf logs
rm -f data/pgc_era1_mf/* data/pgc_era2_mf/*
for y in 2010 2011 2012 2013 2014 2015 2016 2017; do [ -f data/mf_b011p90f40_${y}.csv ] && ln -sf ../mf_b011p90f40_${y}.csv data/pgc_era1_mf/; done
for y in 2018 2019 2020 2021 2022 2023 2024 2025 2026; do [ -f data/mf_b011p90f40_${y}.csv ] && ln -sf ../mf_b011p90f40_${y}.csv data/pgc_era2_mf/; done
log(){ echo "$(date +%H:%M) [pgc-dvv] $*" | tee -a logs/pgc_dvv_pipeline.log; }

run_era(){
  local NAME=$1 GLOB=$2 TAG=$3 sta=$(echo "$3" | tr A-Z a-z)
  log "$NAME: PGC daily stacks at B011 reliable-family times (fs40, decimate HHZ)"
  $PY scripts/build_long_window_3comp.py --mf-csv-glob "$GLOB" --network CN --station PGC --fs 40 \
     --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 --workers 10 \
     --out-prefix data/long_window_daily_$TAG > logs/pgc_stack_$TAG.log 2>&1 || { log "$NAME stacks FAILED"; return 1; }
  log "$NAME: coda dv/v (2-4s, origin-anchored, PER-ERA reference)"
  $PY scripts/dvv_roll30cal.py --station PGC --npz data/long_window_daily_${TAG}_Z.npz \
     --window 2 4 --origin-anchor --workers 10 --out data/daily_dvv_${TAG}_Z_2to4.csv > logs/pgc_dvv_$TAG.log 2>&1 || { log "$NAME dvv FAILED"; return 1; }
  log "$NAME: causality certification (reliable families only)"
  $PY scripts/finalize_causality.py $TAG $TAG > logs/pgc_finalize_$TAG.log 2>&1 || { log "$NAME finalize FAILED"; return 1; }
  local zc=$($PY -c "import json;print(json.load(open('data/${sta}_3comp_summary.json'))['z_certified'])" 2>/dev/null)
  local zs=$($PY -c "import json;print(json.load(open('data/${sta}_3comp_summary.json'))['z_std_pct'])" 2>/dev/null)
  log "$NAME DONE: ${zc} PGC-certified reliable families, dv/v std ${zs}%"
}

run_era "BHZ 2010-2017" "data/pgc_era1_mf/*.csv" PGCbhz
run_era "HHZ 2018-2026" "data/pgc_era2_mf/*.csv" PGChhz
log "PGC dv/v COMPLETE (2 segments) — next: B011 instrumental-vs-real comparison across Aug-2017"
