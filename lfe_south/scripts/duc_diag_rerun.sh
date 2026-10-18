#!/usr/bin/env bash
# DUC_DIAG re-run (2026-09-30) -- after the cache GPU-OOM diagnosis (results/duc_fix/DIAGNOSIS.md).
# 2008-01 .. 2012-07 only.  Launch-and-exit; markers results/DUC_DIAG_{HEARTBEAT,DONE,FAILED}; flock.
#   0. set aside (MOVE, never delete) the damaged prep days / mf1 months / re-templating outputs
#   A. cache loop (bg, 7 workers, fixed GPU FFT):  rebuild only the errored station-days -> diag_cached_<m>
#   B. prep (1 worker) + MF pass 1 (GPU) per month, 2009-06 first (test month) -> mf1/<m>.parquet
#   C. re-templating --upto 2012-07 WITHOUT --final (windows still open stay pending) + a forced diagnostic
#   D. pass 2 if any template was kept; report; summary into DIAGNOSIS.md
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p logs locks
exec 9>locks/duc_diag.lock
flock -n 9 || { echo "[duc-diag] another instance holds the lock; exiting"; exit 0; }
exec 8>locks/duc_fix.lock
flock -n 8 || { echo "[duc-diag] the DUC_FIX queue is running; exiting"; exit 0; }
FX=results/duc_fix; STAGE=init; FINISHED=0; CHILDREN=""
export DUC_END=2012-07 DUC_DIAG_END=2012-07
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/DUC_DIAG_HEARTBEAT; }
fail(){ echo "[duc-diag] FAILED at stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/DUC_DIAG_FAILED
        for p in $CHILDREN; do kill "$p" 2>/dev/null; done; exit "$1"; }
onexit(){ rc=$?; for p in $CHILDREN; do kill "$p" 2>/dev/null; done; [ "$FINISHED" = 1 ] && return 0; [ -f results/DUC_DIAG_FAILED ] && return 0
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/DUC_DIAG_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[duc-diag] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
rm -f results/DUC_DIAG_FAILED results/DUC_DIAG_DONE
echo $$ > results/DUC_DIAG_PID
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUC_MUTE=0
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
PY=python
nextm(){ $PY -c "import pandas as pd,sys; print(pd.Period(sys.argv[1],'M')+1)" "$1"; }

[ -f $FX/pre_diag/SETASIDE_DONE ] || { run setaside $PY scripts/duc_diag.py setaside; date -u +%FT%TZ > $FX/pre_diag/SETASIDE_DONE; }
AFF=$($PY scripts/duc_diag.py affected) || fail 2
ORDER="2009-06 $(echo $AFF | tr ' ' '\n' | grep -v '^2009-06$' | tr '\n' ' ')"
echo "[duc-diag] affected months: $(echo $AFF | wc -w); order: $ORDER"

# ---------- A. cache loop ----------
(
  for M in $ORDER; do
    [ -f $FX/data/diag_cached_$M ] && continue
    $PY scripts/duc_fix_data.py --month $M --stage cache --workers 7 >> logs/duc_diag_cache.log 2>&1 || echo "$M" >> $FX/data/diag_cache_failed.txt
    echo "[cache $M] $($PY scripts/duc_diag.py cachecheck $M)" >> logs/duc_diag_cache.log
    date -u +%FT%TZ > $FX/data/diag_cached_$M
  done
) &
CPID=$!; CHILDREN="$CHILDREN $CPID"

# ---------- B. prep + MF per month ----------
for M in $ORDER; do
  STAGE="wait_$M"; beat
  until [ -f $FX/data/diag_cached_$M ]; do
    kill -0 $CPID 2>/dev/null || [ -f $FX/data/diag_cached_$M ] || fail 4
    sleep 60; beat
  done
  N=$(nextm $M)
  run prep_$M $PY scripts/duc_fix_prep.py --t0 $M-01 --t1 $N-01 --workers 1 --tag diag_$M
  [ -f $FX/mf1/$M.parquet ] || run mf1_$M $PY scripts/duc_fix_mf.py --month $M --out $FX/mf1/$M.parquet
done
wait $CPID 2>/dev/null

# ---------- C. re-templating ----------
run retpl $PY scripts/duc_fix_retemplate.py --upto 2012-07 --workers 4
run retpl_forced $PY scripts/duc_diag.py retpl_forced

# ---------- D. pass 2 (only if templates were kept), report, summary ----------
if [ -f $FX/rerun_months.txt ]; then
  for M in $(sort -u $FX/rerun_months.txt); do
    [ -f $FX/mf2/$M.parquet ] && continue
    run mf2_$M $PY scripts/duc_fix_mf.py --month $M --out $FX/mf2/$M.parquet
  done
fi
run report $PY scripts/duc_fix_report.py
run summary $PY scripts/duc_diag.py summary
date -u +%FT%TZ > results/DUC_DIAG_DONE
STAGE=done; FINISHED=1; beat
echo "[duc-diag] DONE $(date -u +%FT%TZ)"
