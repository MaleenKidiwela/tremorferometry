#!/usr/bin/env bash
# SH_BACK2026 (2026-09-30): backward extension of Shelly (2025)'s 27 Mendocino families before 2018.
#   BK0 epochs (full FDSN history, 31 stations) -> S6 = the S1 stations on the build sensor model before 2018, S2 = KMPB+KMR
#   BK1 known answer: v3b it0 detector + Text S1 rules with ONLY S6 (and S2) on the 385 held-out 2018-2024 test days,
#       scored exactly like v3b (sh_v3_rescore.py) -- runs immediately (cache already present), 4 prep workers
#   gate: BK2 only if BK1 S6 at 9xMAD (S1, source) any-family precision >= 0.5
#   BK2 (after the forward queue has finished: waits for results/SH_EXT_DONE or SH_EXT_FAILED): df -h + volume,
#       fetch/cache S6 2008-01-01..2017-12-31 (same-sensor periods), prep, MF, S1 catalogue, per-family daily counts, report
# Markers results/SH_BACK_{DONE,FAILED,HEARTBEAT}; pid results/SH_BACK_PID; lock locks/sh_back.lock; log logs/sh_back_queue.log
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p results/sh_back logs locks
exec 9>locks/sh_back.lock
flock -n 9 || { echo "[back] another instance holds the lock; exiting"; exit 0; }
echo $$ > results/SH_BACK_PID
STAGE="init"; FINISHED=0; HB=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/SH_BACK_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 300; do [ -f results/SH_BACK_HEARTBEAT ] && touch results/SH_BACK_HEARTBEAT; done ) & HB=$!; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; }
fail(){ echo "[back] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/SH_BACK_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f results/SH_BACK_FAILED ] && return 0
          echo "[back] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/SH_BACK_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[back] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
soft(){ STAGE="$1"; shift; beat; echo "[back] === $STAGE (non-fatal) $(date -u +%FT%TZ)"; "$@" || echo "[back] WARN $STAGE rc=$? -- continuing"; }
rm -f results/SH_BACK_FAILED results/SH_BACK_DONE
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
B=results/sh_back
PY="python scripts/sh_back2026.py"
TEST=results/aa_shelly/stage2_days_test.json
echo "[back] start $(date -u +%FT%TZ)"

# ---- BK0 / BK1 ------------------------------------------------------------------------------------------------------
[ -f $B/back_station_summary.csv ] || run epochs $PY epochs
for S in S6 S2; do
  [ -f $B/tpl_$S.csv ] || run tpl_$S $PY tplsub --set $S
  STAS=$(python -c "import sys;sys.path.insert(0,'scripts');import sh_back2026 as b;print(','.join(b.SETS['$S']))")
  if [ ! -f $B/PREP_${S}_TEST_DONE ]; then
    run prep_$S python scripts/sh_ext2026.py prep --days-file $TEST --workers 4 --stas $STAS --prep-dir $WS/$B/prep_$S --summary $WS/$B/prep_summary_$S.csv
    date -u +%FT%TZ > $B/PREP_${S}_TEST_DONE
  fi
  [ -f $B/mf_bk1_$S.parquet ] || run mf_bk1_$S env SH_PREP_DIR=$WS/$B/prep_$S python scripts/t11_mf_v2.py --templates $B/tpl_$S \
      --days-file $TEST --truth shelly --t0 2018-01-01 --t1 2024-12-01 --out $B/mf_bk1_$S.parquet
  [ -f $B/rescore_bk1_$S.csv ] || run rescore_$S python scripts/sh_v3_rescore.py --ours $B/mf_bk1_$S.parquet --tag bk1_$S --out $B/rescore_bk1_$S.csv
done
run gate $PY gate
soft report_bk1 $PY report

if ! grep -q PASS $B/BK1_GATE.txt; then
  echo "[back] gate FAILED -> BK2 not run ($(cat $B/BK1_GATE.txt))"
  date -u +%FT%TZ > results/SH_BACK_DONE; STAGE="done"; FINISHED=1; beat; echo "[back] DONE (BK1 only) $(date -u +%FT%TZ)"; exit 0
fi

# ---- BK2: only after the forward queue has finished (it must not block / compete with the forward work) -----------------
STAGE="wait_forward"; beat
while [ ! -f results/SH_EXT_DONE ] && [ ! -f results/SH_EXT_FAILED ]; do sleep 600; done
echo "[back] forward queue finished ($(ls results/SH_EXT_DONE results/SH_EXT_FAILED 2>/dev/null)); starting BK2 $(date -u +%FT%TZ)"
df -h /home/jovyan | tail -1
[ -f $B/days_bk2.json ] || run days2 $PY days2
[ -f $B/download_plan_bk2.csv ] || run plan2 $PY plan2
soft report_vol $PY report
if [ ! -f $B/FETCH_DONE ]; then
  run fetch python scripts/t2_fetch.py $B/download_plan_bk2.csv --workers 8 --manifest $B/download_manifest_bk2.csv
  run fetch_retry python scripts/t2_fetch.py $B/download_plan_bk2.csv --workers 4 --manifest $B/download_manifest_bk2.csv
  date -u +%FT%TZ > $B/FETCH_DONE
fi
if [ ! -f $B/CACHE_DONE ]; then
  run cache python scripts/t3_cache.py $B/download_plan_bk2.csv --workers 12 --out $B/cache_bk2.parquet
  date -u +%FT%TZ > $B/CACHE_DONE
fi
if [ ! -f $B/PREP_BK2_DONE ]; then
  run prep_bk2 python scripts/sh_ext2026.py prep --days-file $B/days_bk2.json --workers 12 --stas NC.KMPB,NC.KMR,PB.B049,PB.B932,PB.B933,PB.B934 \
      --prep-dir $WS/$B/prep_S6 --summary $WS/$B/prep_summary_S6.csv
  date -u +%FT%TZ > $B/PREP_BK2_DONE
fi
[ -f $B/mf_bk2.parquet ] || run mf_bk2 env SH_PREP_DIR=$WS/$B/prep_S6 python scripts/t11_mf_v2.py --templates $B/tpl_S6 \
    --days-file $B/days_bk2.json --t0 2008-01-01 --t1 2018-01-01 --out $B/mf_bk2.parquet
run catalog2 $PY catalog2
run report $PY report
date -u +%FT%TZ > results/SH_BACK_DONE
STAGE="done"; FINISHED=1; beat
echo "[back] DONE $(date -u +%FT%TZ)"
