#!/usr/bin/env bash
# DUC_FIX (2026-09-30) -- Ducellier & Creager (2022) 66 families, his templates, his stations, extended
# 2008-01 .. 2026-09 with the FIXED detector.  ONE detached queue, flock, guarded stages, markers:
#   results/DUC_FIX_HEARTBEAT  stage + pid + UTC (rewritten continuously)
#   results/DUC_FIX_PID        queue pid
#   results/DUC_FIX_DONE | results/DUC_FIX_FAILED
# Stage 1 (known answer) failing stops everything.
# Loops: fetch (8 threads) -> cache (8 workers, GPU FFT) -> prep (3 workers) -> MF pass 1 (GPU) + re-templating, then
# pass 2 re-scans of months touched by new-epoch templates, then the report.
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p results/duc_fix/{data,mf1,mf2,prep} logs locks
exec 9>locks/duc_fix.lock
flock -n 9 || { echo "[duc-fix] another instance holds the lock; exiting"; exit 0; }
FX=results/duc_fix
STAGE="init"; FINISHED=0; CHILDREN=""
echo $$ > results/DUC_FIX_PID
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/DUC_FIX_HEARTBEAT; }
killkids(){ for p in $CHILDREN; do kill "$p" 2>/dev/null; done; }
fail(){ echo "[duc-fix] FAILED at stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/DUC_FIX_FAILED; killkids; exit "$1"; }
onexit(){ rc=$?; killkids; [ "$FINISHED" = 1 ] && return 0; [ -f results/DUC_FIX_FAILED ] && return 0
          echo "[duc-fix] EXITED unexpectedly at stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/DUC_FIX_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[duc-fix] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
soft(){ STAGE="$1"; shift; beat; echo "[duc-fix] === $STAGE (non-fatal) $(date -u +%FT%TZ)"; "$@" || echo "[duc-fix] WARN $STAGE rc=$? - continuing"; }
rm -f results/DUC_FIX_FAILED
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
PY=python
MONTHS=$($PY - <<'PY'
import os, pandas as pd; print(" ".join(str(m) for m in pd.period_range("2008-01", os.environ.get("DUC_END", "2026-09"), freq="M")))
PY
)
nextm(){ $PY -c "import pandas as pd,sys; print(pd.Period(sys.argv[1],'M')+1)" "$1"; }

# ---------- 0. station epochs ----------
[ -f $FX/station_epochs.csv ] || run epochs $PY scripts/duc_epochs.py

# ---------- 1. KNOWN ANSWER: his 2008-04-01..05-10 reproduction on the fixed detector ----------
if [ ! -f $FX/detector_choice.json ]; then
  run s1_prep $PY scripts/duc_fix_prep.py --t0 2008-04-01 --t1 2008-05-12 --workers 6 --tag repro
  run stage1  $PY scripts/duc_fix_stage1.py
fi
[ -f $FX/detector_choice.json ] || fail 3
export DUC_MUTE=$($PY -c "import json; print(json.load(open('$FX/detector_choice.json'))['mute'])")
echo "[duc-fix] detector: DUC_MUTE=$DUC_MUTE"
soft report_s1 $PY scripts/duc_fix_report.py

# ---------- 2. fetch loop (background, I/O) -> fetched_<m> ; cache loop (background, 8 CPU) -> cached_<m> ----------
(
  for M in $MONTHS; do
    [ -f $FX/data/fetched_$M ] && continue
    for att in 1 2 3; do
      $PY scripts/duc_fix_data.py --month $M --stage fetch --fetch-workers 8 >> logs/duc_fix_fetch.log 2>&1 && break
      echo "[fetch] $M attempt $att failed $(date -u +%FT%TZ)" >> logs/duc_fix_fetch.log; sleep 120
    done
    date -u +%FT%TZ > $FX/data/fetched_$M
  done
) &
FETCH_PID=$!; CHILDREN="$CHILDREN $FETCH_PID"
(
  for M in $MONTHS; do
    [ -f $FX/data/cached_$M ] && continue
    until [ -f $FX/data/fetched_$M ]; do
      kill -0 $FETCH_PID 2>/dev/null || [ -f $FX/data/fetched_$M ] || { echo "[cache] fetch loop gone before $M" >> logs/duc_fix_cache.log; exit 4; }
      sleep 30
    done
    $PY scripts/duc_fix_data.py --month $M --stage cache --workers 8 >> logs/duc_fix_cache.log 2>&1 || echo "$M" >> $FX/data/cache_failed_months.txt
    date -u +%FT%TZ > $FX/data/cached_$M
  done
) &
DATA_PID=$!; CHILDREN="$CHILDREN $DATA_PID"
# ---------- 3. prep loop (background): per-day cleaning prep -> ready_<m> ----------
(
  for M in $MONTHS; do
    [ -f $FX/data/ready_$M ] && continue
    until [ -f $FX/data/cached_$M ]; do
      kill -0 $DATA_PID 2>/dev/null || [ -f $FX/data/cached_$M ] || { echo "[prep] data loop gone before $M" >> logs/duc_fix_prep.log; exit 4; }
      sleep 30
    done
    N=$(nextm $M)
    $PY scripts/duc_fix_prep.py --t0 $M-01 --t1 $N-01 --workers 3 --tag $M >> logs/duc_fix_prep.log 2>&1 || echo "$M" >> $FX/data/prep_failed_months.txt
    date -u +%FT%TZ > $FX/data/ready_$M
  done
) &
PREP_PID=$!; CHILDREN="$CHILDREN $PREP_PID"

# ---------- 4. MF pass 1 (chronological) + re-templating of new epochs ----------
LAST=${DUC_END:-2026-09}; n=0
for M in $MONTHS; do
  N=$(nextm $M)
  STAGE="mf1_wait_$M"; beat
  until [ -f $FX/data/ready_$M ] && { [ "$M" = "$LAST" ] || [ -f $FX/data/ready_$N ]; }; do
    kill -0 $PREP_PID 2>/dev/null || [ -f $FX/data/ready_$M ] || fail 5
    sleep 60; beat
  done
  if [ ! -f $FX/mf1/$M.parquet ]; then
    run mf1_$M $PY scripts/duc_fix_mf.py --month $M --out $FX/mf1/$M.parquet
  fi
  if [ ! -f $FX/data/retpl_$M ]; then
    run retpl_$M $PY scripts/duc_fix_retemplate.py --upto $M --workers 2
    date -u +%FT%TZ > $FX/data/retpl_$M
  fi
  n=$((n+1)); [ $((n % 12)) = 0 ] && soft report_$M $PY scripts/duc_fix_report.py
done
# DUC_DIAG (2026-09-30): --final only at the REAL end of the record; a partial scan (DUC_END earlier) must
# leave epochs whose 24-month window is still open PENDING (the 2012-07 partial run force-dropped 26 epochs).
[ "$LAST" = "2026-09" ] && [ ! -f $FX/data/retpl_final ] && { run retpl_final $PY scripts/duc_fix_retemplate.py --upto $LAST --final --workers 4; date -u +%FT%TZ > $FX/data/retpl_final; }

# ---------- 5. pass 2: re-scan months touched by new-epoch templates ----------
if [ -f $FX/rerun_months.txt ]; then
  for M in $(sort -u $FX/rerun_months.txt); do
    [ -f $FX/mf2/$M.parquet ] && continue
    run mf2_$M $PY scripts/duc_fix_mf.py --month $M --out $FX/mf2/$M.parquet
  done
fi

# ---------- 6. report ----------
wait $FETCH_PID 2>/dev/null; wait $DATA_PID 2>/dev/null; wait $PREP_PID 2>/dev/null
run report $PY scripts/duc_fix_report.py
date -u +%FT%TZ > results/DUC_FIX_DONE
STAGE="done"; FINISHED=1; beat
echo "[duc-fix] DONE $(date -u +%FT%TZ)"
