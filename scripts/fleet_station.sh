#!/bin/bash
# FLEET BROADBAND/LAND dv/v driver — one station, the VALIDATED recipe (notes/FLEET_BROADBAND_PIPELINE.md).
# 40 Hz throughout; picker --target-fs 40; cluster --min-years 1; forward densify; causality gate >=20 & >=15%.
# Usage: nohup bash scripts/fleet_station.sh <NET> <STA> <START-YYYY-MM-DD> &
set -o pipefail
cd /home/jovyan/tremorferometry
NET=$1; STA=$2; START=${3:-2009-01-01}; CHAN=${4:-"BH?,HH?,EH?,SH?"}; s=$(echo "$STA" | tr A-Z a-z)
PY=/home/jovyan/envs/tremorferometry/bin/python
gpu_env(){ export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH"
           export CUDA_PATH=/opt/conda/targets/x86_64-linux; export PYTHONPATH=src; }
wait_gpu_free(){ while ps -eo cmd | grep -qE '[d]ensify_gnw_gpu|[d]iscover_gpu'; do sleep 60; done; }  # serialize the single GPU
log(){ echo "$(date +%H:%M) [fleet $STA] $*" | tee -a logs/fleet_${s}.log; }

# station coords from the fleet order (fallback FDSN)
read LA LO < <($PY -c "
import pandas as pd
d=pd.read_csv('data/broadband_fleet_order.csv'); r=d[d.sta=='$STA']
print(f'{r.lat.iloc[0]:.4f} {r.lon.iloc[0]:.4f}') if len(r) else print('')" 2>/dev/null)
[ -n "$LA" ] || { log "ABORT: no coords for $STA"; exit 1; }
if [ "$CHAN" = "auto" ]; then CHAN=$($PY scripts/pick_band.py $NET $STA 2>/dev/null); log "auto-band -> $CHAN"; fi
[ -n "$CHAN" ] || CHAN="BH?,HH?,EH?,SH?"
la0=$($PY -c "print(round($LA-0.9,2))"); la1=$($PY -c "print(round($LA+0.9,2))")
lo0=$($PY -c "print(round($LO-1.35,2))"); lo1=$($PY -c "print(round($LO+1.35,2))")
log "start: coords $LA $LO  bbox [$la0 $la1 $lo0 $lo1]"

# 1. DOWNLOAD — ALWAYS run (resumable: skips existing days, completes any partial). Never skip on a fragment.
log "download $NET.$STA from $START (chan $CHAN)"
$PY scripts/download_broadband.py --net $NET --sta $STA --start $START --end 2026-08-01 --workers 6 \
   --chan "$CHAN" > logs/dl_${s}.log 2>&1 || { log "download FAILED"; exit 2; }
log "download: $(find data/waveforms/$NET.$STA -name '*.mseed' | wc -l) day-files"

# 2. DETECT @40
CAND=data/${s}_pnsn_candidates.parquet
[ -f "$CAND" ] || { log "detect @40Hz"; PYTHONPATH=src $PY scripts/discover_nllb_pnsn_driven.py --station $STA \
   --wfdir data/waveforms --bbox $la0 $la1 $lo0 $lo1 --pnsn catalogs/pnsn_tremor_cascadia_full.csv --fs 40 \
   --candidates-only --candidates-out "$CAND" --workers 16 > logs/cand_${s}.log 2>&1 || { log "detect FAILED"; exit 3; }; }
NCAND=$($PY -c "import pandas as pd;print(len(pd.read_parquet('$CAND')))" 2>/dev/null)
log "detect: ${NCAND:-0} candidates"
if [ "${NCAND:-0}" -lt 100 ]; then log "GATE: ${NCAND:-0} candidates (<100) -> FLAG (no tremor/data)"; log "STATION DONE: $STA"; exit 0; fi

# 3. SCORE @40 (broadband picker, --target-fs 40)
cp -f lfe_features/models/picker_broadband_pgc.joblib lfe_features/models/tremor_picker_${s}.joblib
Y0=$($PY -c "print(int('$START'[:4]))");
[ -f data/${s}_cand_baseline.parquet ] || { log "score @40 (--target-fs 40)"; PYTHONPATH=src $PY lfe_features/score_candidates.py \
   --net $NET --sta $STA --cand "$CAND" --y0 $Y0 --y1 2026 --thr 0.5 --target-fs 40 --workers 16 \
   > logs/score_${s}.log 2>&1 || { log "score FAILED"; exit 4; }; }

# 4. SELECT adaptive top-30k
$PY scripts/adaptive_cand_threshold.py $STA 30000 > logs/adaptive_${s}.log 2>&1 || { log "adaptive FAILED"; exit 5; }
log "select: $(tail -1 logs/adaptive_${s}.log)"

# 5. CLUSTER @40 --min-years 1
wait_gpu_free; gpu_env
DISC=data/${s}_disc
[ -f ${DISC}.summary.csv ] || { log "cluster @40 --min-years 1"; $PY scripts/discover_gpu.py --station $STA \
   --wfdir data/waveforms --candidates data/${s}_cand_filtered.parquet --out ${DISC}.npz --fs 40 \
   --cc-threshold 0.80 --min-family-members 3 --min-years 1 --workers 24 > logs/cluster_${s}.log 2>&1 \
   || { log "cluster FAILED"; exit 6; }; }
NFAM=$(($(wc -l < ${DISC}.summary.csv)-1))
log "cluster: $NFAM families"
# early-exit: <50 families -> too few to reach >=20 certified even at best (~38%) survival; FLAG, skip densify.
# (was <60, but 53-59-family stations at high-survival regions can clear the bar -> gave GOBB/SYMB a re-check)
if [ "$NFAM" -lt 50 ]; then log "GATE: $NFAM families (<50) -> FLAG (skip densify)"; log "STATION DONE: $STA"; exit 0; fi

# 5b. Coverage-select the 300 densify budget by family SNR (validated predictive of causality survival, AUC 0.68;
#     the family-stack picker P(LFE) is ANTI-predictive AUC 0.40 -> dropped). Causality is the SOLE LFE verdict.
SEL=${DISC}_sel300.summary.csv
$PY scripts/select_families_coverage.py "$DISC" - 300 snr > logs/select_${s}.log 2>&1 || { log "select FAILED"; exit 7; }
NSEL=$(($(wc -l < "$SEL" 2>/dev/null)-1)); log "coverage-select (SNR rank, cap 300): $NSEL families"

# 6. DENSIFY (forward only) — coverage-selected, LFE-labeled families (<=300 budget)
log "densify (fwd) @40Hz: $NSEL families"
$PY scripts/densify_gnw_gpu.py --templates-npz ${DISC}.npz --summary-csv "$SEL" --min-snr 0 \
   --network $NET --station $STA --fs 40 --fmin 2 --fmax 8 --threshold 0.8 --min-gap-s 6 --top-n 0 \
   --despike-mad 8 --workers 8 --start-year $Y0 --end-year 2026 --out-dir data --out-prefix mf_${s}p90f40_ \
   > logs/densify_${s}.log 2>&1 || { log "densify FAILED"; exit 8; }
log "densify: $(ls data/mf_${s}p90f40_*.csv 2>/dev/null | wc -l) year files"

# 7-9. STACKS -> dv/v -> CAUSALITY
TAG=$(echo $STA | tr a-z A-Z)
$PY scripts/build_long_window_3comp.py --mf-csv-glob "data/mf_${s}p90f40_*.csv" --network $NET --station $STA \
   --fs 40 --fmin 2 --fmax 8 --cc-min 0.80 --min-det 20 --despike-mad 8 --workers 10 \
   --out-prefix data/long_window_daily_${TAG} > logs/stack_${s}.log 2>&1 || { log "stacks FAILED"; exit 8; }
$PY scripts/dvv_roll30cal.py --station $STA --npz data/long_window_daily_${TAG}_Z.npz --window 2 4 \
   --origin-anchor --workers 10 --out data/daily_dvv_${TAG}_Z_2to4.csv > logs/dvv_${s}.log 2>&1 || { log "dvv FAILED"; exit 9; }
$PY scripts/finalize_causality.py ${TAG} ${TAG} > logs/finalize_${s}.log 2>&1 || { log "finalize FAILED"; exit 10; }

# 10. GATE — survival = certified / DENSIFIED (NSEL), the fraction that passed causality
ZC=$($PY -c "import json;print(json.load(open('data/${s}_3comp_summary.json'))['z_certified'])" 2>/dev/null)
SURV=$($PY -c "print(f'{100*$ZC/$NSEL:.0f}')" 2>/dev/null)
VERDICT=$($PY -c "print('INCLUDE' if $ZC>=20 and $ZC/$NSEL>=0.15 else 'FLAG')" 2>/dev/null)
log "GATE: $ZC certified / $NSEL densified (${SURV}% survival) -> $VERDICT  [bar >=20 & >=15%]"
log "STATION DONE: $STA"
