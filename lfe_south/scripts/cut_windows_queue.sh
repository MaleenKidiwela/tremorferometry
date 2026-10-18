#!/usr/bin/env bash
# CUT REPEAT WINDOWS queue (2026-09-30).  Waits for an extension run to finish, then
#   cut (scripts/cut_windows.py, 6 workers) -> verify (20 random windows vs raw, sample by sample) ->
#   stackcheck (one family x good station x month stacked vs its v3b template) -> deletelist (LIST ONLY, nothing deleted)
# Usage: bash scripts/cut_windows_queue.sh [sh_ext|sh_back|duc_fix]      (default sh_ext = Shelly forward)
#   waits for results/<UP>_DONE (stops on <UP>_FAILED), UP = SH_EXT | SH_BACK | DUC_FIX
# Markers (sh_ext): results/CUT_{DONE,FAILED,HEARTBEAT}, pid results/CUT_PID; other runs: results/CUT_<RUN>_{...}
# Lock locks/cut_<run>.lock; log logs/cut_<run>_queue.log (when launched with nohup ... >> log).
set -u -o pipefail
RUN=${1:-sh_ext}
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
case "$RUN" in
  sh_ext)  UP=SH_EXT;  MK=CUT ;;
  sh_back) UP=SH_BACK; MK=CUT_SH_BACK ;;
  duc_fix) UP=DUC_FIX; MK=CUT_DUC_FIX ;;
  *) echo "unknown run $RUN"; exit 2 ;;
esac
mkdir -p results/cut logs locks
exec 9>locks/cut_$RUN.lock
flock -n 9 || { echo "[cut] another instance holds the lock; exiting"; exit 0; }
echo $$ > results/${MK}_PID
STAGE="init"; FINISHED=0; HB=""
beat(){ echo "run=$RUN stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/${MK}_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 300; do [ -f results/${MK}_HEARTBEAT ] && touch results/${MK}_HEARTBEAT; done ) & HB=$!; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; }
fail(){ echo "[cut] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "run=$RUN stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/${MK}_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f results/${MK}_FAILED ] && return 0
          echo "[cut] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "run=$RUN stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/${MK}_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[cut] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
rm -f results/${MK}_FAILED results/${MK}_DONE
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "[cut] start run=$RUN waiting for results/${UP}_DONE $(date -u +%FT%TZ)"

# ---- wait for the extension run (shell sleep loop; stop if it failed) --------------------------------------------------
STAGE="wait_${UP}"; beat
while [ ! -f results/${UP}_DONE ]; do
  if [ -f results/${UP}_FAILED ]; then
    STAGE="upstream_failed"; echo "[cut] results/${UP}_FAILED present ($(cat results/${UP}_FAILED)); not cutting"; fail 7
  fi
  sleep 600; beat
done
echo "[cut] ${UP}_DONE present ($(cat results/${UP}_DONE)); starting $(date -u +%FT%TZ)"; df -h /home/jovyan | tail -1

PY="python scripts/cut_windows.py"
run cut        $PY cut --run $RUN --workers 6
run verify     $PY verify --run $RUN --n 20
run stackcheck $PY stackcheck --run $RUN
run deletelist $PY deletelist --run $RUN
date -u +%FT%TZ > results/${MK}_DONE
STAGE="done"; FINISHED=1; beat
echo "[cut] DONE run=$RUN $(date -u +%FT%TZ) -- nothing deleted; list in results/cut/DELETE_CANDIDATES_*.csv"
