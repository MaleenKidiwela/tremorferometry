#!/usr/bin/env bash
# SH_EXT2026 (2026-09-30): extend Shelly (2025)'s 27 Mendocino LFE families from 2024-11-30 to the latest day with the
# validated SH_V3B it0 detector (17 Text S1 stations, v2 cleaning, t11_mf_v2, Text S1 catalogue rules).
# Stages (all guarded, resumable): epochs -> days -> continuity check on 2024 held-out days (GATE) -> plan/fetch/cache
# (2023-12-01..LAST, unchanged t2_fetch/t3_cache) -> prep (17 stations) -> MF segments A (it0) / NK (it0 without KMR) ->
# KMR re-templates for its two new sensor epochs -> MF B / C -> catalogue -> tplcheck (non-fatal) -> report.
# Markers results/SH_EXT_{DONE,FAILED,HEARTBEAT}; pid results/SH_EXT_PID; lock locks/sh_ext2026.lock; log logs/sh_ext2026_queue.log
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p results/sh_ext2026 logs locks
exec 9>locks/sh_ext2026.lock
flock -n 9 || { echo "[ext] another instance holds the lock; exiting"; exit 0; }
echo $$ > results/SH_EXT_PID
STAGE="init"; FINISHED=0; HB=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/SH_EXT_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 300; do [ -f results/SH_EXT_HEARTBEAT ] && touch results/SH_EXT_HEARTBEAT; done ) & HB=$!; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; }
fail(){ echo "[ext] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/SH_EXT_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f results/SH_EXT_FAILED ] && return 0
          echo "[ext] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/SH_EXT_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[ext] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
soft(){ STAGE="$1"; shift; beat; echo "[ext] === $STAGE (non-fatal) $(date -u +%FT%TZ)"; "$@" || echo "[ext] WARN $STAGE rc=$? -- continuing"; }
rm -f results/SH_EXT_FAILED results/SH_EXT_DONE
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
V=results/sh_ext2026; DD=$V/days
export SH_PREP_DIR=$WS/$V/prep SH_V3_DIR=$WS/$V SH_EXT_LAST=${SH_EXT_LAST:-2026-09-28}
TEND=$(python -c "import pandas as pd;print((pd.Timestamp('$SH_EXT_LAST')+pd.Timedelta(days=1)).date())")
TPL0=results/sh_v3b/tpl_it0
MF="python scripts/t11_mf_v2.py"
PY="python scripts/sh_ext2026.py"
echo "[ext] start $(date -u +%FT%TZ) LAST=$SH_EXT_LAST TEND=$TEND"; df -h /home/jovyan | tail -1

# ---- 1. station epochs + day lists ---------------------------------------------------------------------------------
[ -f $V/epoch_decisions.csv ] || run epochs $PY epochs
[ -f $DD/C.json ] || run days $PY days

# ---- 2. continuity check (2024 test days inside his catalogue; cache already present) -- GATE -------------------------
if [ ! -f $V/CONTPREP_DONE ]; then run contprep $PY prep --days-file $DD/cont.json --workers 12; date -u +%FT%TZ > $V/CONTPREP_DONE; fi
[ -f $V/mf_cont.parquet ] || run mf_cont $MF --templates $TPL0 --days-file $DD/cont.json --truth shelly --t0 2024-01-01 --t1 2024-12-01 --out $V/mf_cont.parquet
grep -q PASS $V/CONT_GATE.txt 2>/dev/null || run cont $PY cont

# ---- 3. data: 17 stations, 2023-12-01 .. LAST (identical fetch / cache code) ------------------------------------------
[ -f $V/download_plan.csv ] || run plan $PY plan
if [ ! -f $V/FETCH_DONE ]; then
  run fetch python scripts/t2_fetch.py $V/download_plan.csv --workers 8 --manifest $V/download_manifest.csv
  run fetch_retry python scripts/t2_fetch.py $V/download_plan.csv --workers 4 --manifest $V/download_manifest.csv
  date -u +%FT%TZ > $V/FETCH_DONE
fi
if [ ! -f $V/CACHE_DONE ]; then
  run cache python scripts/t3_cache.py $V/download_plan.csv --workers 12 --out $V/cache_ext.parquet
  date -u +%FT%TZ > $V/CACHE_DONE
fi
soft datasum $PY datasum
if [ ! -f $V/PREP_DONE ]; then run prep $PY prep --days-file $DD/all.json --workers 12; date -u +%FT%TZ > $V/PREP_DONE; fi

# ---- 4. matched filter segments -------------------------------------------------------------------------------------
[ -f $V/mf_A.parquet ] || run mf_A $MF --templates $TPL0 --days-file $DD/A.json --truth shelly --t0 2023-12-01 --t1 2025-09-17 --out $V/mf_A.parquet
[ -f $V/tpl_noKMR.csv ] || run tpl_nokmr $PY tplset --stage nokmr
[ -f $V/mf_NK.parquet ] || run mf_NK $MF --templates $V/tpl_noKMR --days-file $DD/NK.json --t0 2025-09-17 --t1 $TEND --out $V/mf_NK.parquet
for ep in B C; do
  [ -f $V/cat_NK_$ep.parquet ] || run s1cat_$ep $PY s1cat --inp $V/mf_NK.parquet --days-file $DD/$ep.json --out $V/cat_NK_$ep.parquet
  [ -f $V/tpl_kmr${ep}_all.csv ] || run retpl_$ep python scripts/sh_ext2026_restack.py --stations NC.KMR -- --mode det \
      --det $V/cat_NK_$ep.parquet --picks $TPL0.csv --days-file $DD/$ep.json --thr 9 --workers 4 --out $V/tpl_kmr$ep.npz
done
[ -f $V/tpl_C.csv ] || run tplset $PY tplset --stage kmr
[ -f $V/mf_B.parquet ] || run mf_B $MF --templates $V/tpl_B --days-file $DD/B.json --t0 2025-09-18 --t1 2026-03-12 --out $V/mf_B.parquet
[ -f $V/mf_C.parquet ] || run mf_C $MF --templates $V/tpl_C --days-file $DD/C.json --t0 2026-03-13 --t1 $TEND --out $V/mf_C.parquet

# ---- 5. catalogue, template stability, report ---------------------------------------------------------------------------
run catalog $PY catalog
for W in W0 W1 W4; do
  [ -f $V/tplcheck_${W}_all.csv ] || soft tplcheck_$W python scripts/sh_ext2026_restack.py -- --mode det --det $V/cat_S1_9x_internal.parquet \
      --picks $TPL0.csv --days-file $DD/$W.json --thr 9 --nmax 150 --workers 4 --out $V/tplcheck_$W.npz
done
soft tplcheck $PY tplcheck
run report $PY report
date -u +%FT%TZ > results/SH_EXT_DONE
STAGE="done"; FINISHED=1; beat
echo "[ext] DONE $(date -u +%FT%TZ)"
