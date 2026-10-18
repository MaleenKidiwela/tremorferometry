#!/usr/bin/env bash
# =====================================================================================
# Royer & Bostock (2014) Washington reproduction -- agent AB's queue.
#
# Launch (from a FOREGROUND Bash call, never inside a background task):
#   cd /home/jovyan/tremorferometry/lfe_array_discovery && \
#     nohup setsid -f bash scripts/ab_royer_queue.sh >> logs/ab_royer_queue.log 2>&1 < /dev/null
#
# Markers written under results/ :
#   AB_ROYER_AUDIT_DONE     after stage 0 (the zero-fill / validity-mask quantification)
#   AB_ROYER_AUTOCORR_DONE  after both subarrays' network autocorrelation
#   AB_ROYER_DONE           the literal chain finished for ALL THREE subarrays (NW, SW, NVI)
#   AB_ROYER_EXT_DONE / AB_ROYER_EXT_FAILED   the extension after it
#   AB_ROYER_FAILED         any stage failed, or the queue was signalled
# Single instance is enforced with flock; every stage is guarded by its own output so a relaunch
# resumes instead of redoing hours of work.  Never use pkill -f to stop it: read the pid from
# results/AB_ROYER_PID and `kill -TERM` that pid.
# =====================================================================================
set -u -o pipefail

cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
STUDY=$PWD
RES=$STUDY/results
OUT=$RES/ab_royer
LOGS=$STUDY/logs
LOCKS=$STUDY/locks
mkdir -p "$OUT" "$LOGS" "$LOCKS"
PY=$(command -v python)
# GPU libs (project memory: cupy/torch CUDA-12 libpath)
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}

# ---- single instance -----------------------------------------------------------------
exec 9>"$LOCKS/ab_royer.lock"
if ! flock -n 9; then
  echo "[ab] another ab_royer_queue.sh holds the lock; exiting"
  exit 0
fi
echo $$ > "$RES/AB_ROYER_PID"

STAGE="init"
finish() {
  rc=$?
  if [ -f "$RES/AB_ROYER_DONE" ]; then
    if [ ! -f "$RES/AB_ROYER_EXT_DONE" ] && [ ! -f "$RES/AB_ROYER_EXT_FAILED" ]; then
      { echo "stage=$STAGE"; echo "rc=$rc"; date -u +%FT%TZ; } > "$RES/AB_ROYER_EXT_FAILED"
    fi
    exit 0
  fi
  { echo "stage=$STAGE"; echo "rc=$rc"; date -u +%FT%TZ; } > "$RES/AB_ROYER_FAILED"
  echo "[ab] FAILED in stage $STAGE (rc=$rc)"
  exit $rc
}
trap finish EXIT
trap 'STAGE="$STAGE(signalled)"; exit 143' TERM INT HUP

hb() { date -u +%FT%TZ > "$RES/AB_ROYER_HEARTBEAT"; echo "[ab] $(date -u +%FT%TZ) $*"; }
hb_loop() { while :; do date -u +%FT%TZ > "$RES/AB_ROYER_HEARTBEAT"; sleep 120; done; }
hb_loop & HB=$!
trap 'kill -TERM $HB 2>/dev/null; finish' EXIT

rm -f "$RES/AB_ROYER_FAILED"

# run <stage-name> <guard-file> <env-assignments...> -- <command...>
run() {
  STAGE=$1; guard=$2; shift 2
  if [ -n "$guard" ] && [ -e "$guard" ]; then hb "skip $STAGE (guard $guard exists)"; return 0; fi
  hb "start $STAGE"
  local envs=() ; while [ "$1" != "--" ]; do envs+=("$1"); shift; done; shift
  env "${envs[@]}" "$@" || return 1
  hb "done $STAGE"
}

# =====================================================================================
# STAGE 0 -- the zero-fill / validity-mask quantification on real Washington station-days.
#            Read-only: it opens cached station-days and measures gaps, constant runs, the
#            incidence of the degenerate `emed == 0` guard, and the per-hour MAD exposure.
# =====================================================================================
run stage0_gapscan "$OUT/gapscan_cache20.csv" AB_STAGE=0 -- \
  "$PY" scripts/ab_gapscan.py || finish
date -u +%FT%TZ > "$RES/AB_ROYER_AUDIT_DONE"

# =====================================================================================
# STAGE W -- pre-2009-08 tremor days: Wech & Creager (2008)-style envelope cross-correlation,
#            calibrated against the PNSN catalogue days 2009-08-06..2011-12-31 (ab_wecc.py)
# =====================================================================================
TR=$OUT/tremor
run stageW1_wecc_fetch "$TR/WECC_FETCH_DONE" AB_STAGE=W -- "$PY" scripts/ab_wecc.py fetch || finish
date -u +%FT%TZ > "$TR/WECC_FETCH_DONE"
run stageW2_wecc_env "$TR/WECC_ENV_DONE" AB_STAGE=W -- "$PY" scripts/ab_wecc.py env || finish
date -u +%FT%TZ > "$TR/WECC_ENV_DONE"
run stageW3_wecc_detect "$TR/wecc_windows_SW.parquet" AB_STAGE=W -- "$PY" scripts/ab_wecc.py detect || finish
run stageW4_wecc_calibrate "$TR/wecc_days_SW.csv" AB_STAGE=W -- "$PY" scripts/ab_wecc.py calibrate || finish

# =====================================================================================
# STAGE 0b -- sec.3.3 "all available tremor dates" (ab_tremor_days.py; uses the WECC days)
# =====================================================================================
run stage0b_tremor_days_v3 "$TR/SCAN_DAYS_V3_DONE" AB_STAGE=0 -- "$PY" scripts/ab_tremor_days.py || finish
date -u +%FT%TZ > "$TR/SCAN_DAYS_V3_DONE"

# =====================================================================================
# STAGE 1/2 -- fetch + 40 sps cache: autocorrelation stations on the sec.3.1 days and the
#              subarray pools (NW 25, SW 28, NVI 27) on every scan day.  CN from NRCan.
# =====================================================================================
run stage1_fetch_v3 "$OUT/FETCH_V3_DONE" AB_STAGE=1 -- "$PY" scripts/ab_fetch.py --workers 8 || finish
date -u +%FT%TZ > "$OUT/FETCH_V3_DONE"
run stage2_cache_v3 "$OUT/CACHE_V3_DONE" AB_STAGE=2 -- \
  "$PY" scripts/ab_cache40.py --sub NW,SW,NVI --workers 24 || finish
date -u +%FT%TZ > "$OUT/CACHE_V3_DONE"

SUBS="NW SW NVI"
# =====================================================================================
# STAGE 3 -- sec.3.1 network autocorrelation: 15 s / 0.5 s windows, per-WINDOW-record 8 x MAD
#            (v3 reading), pair-overlap cull; NVI: first 7 of the 10 named stations per hour.
# =====================================================================================
for SUB in $SUBS; do
  run "stage3_auto_$SUB" "$OUT/$SUB/auto_window/_params.json" AB_STAGE=3 -- \
    "$PY" scripts/ab_auto.py --sub "$SUB" --cols window || finish
done
date -u +%FT%TZ > "$RES/AB_ROYER_AUTOCORR_DONE"
# null-1 autocorrelation + cluster sweep: REPORTED columns (null count per linkage/cut/min size)
for SUB in $SUBS; do
  run "stage3n_null1_auto_$SUB" "$OUT/$SUB/null1/auto_window/_params.json" AB_STAGE=3 AB_NULL=1 -- \
    "$PY" scripts/ab_auto.py --sub "$SUB" --cols window --out "$OUT/$SUB/null1/auto_window" || finish
  run "stage4n_null1_sweep_$SUB" "$OUT/$SUB/null1/clusters_sweep/sweep_long.csv" AB_STAGE=4 AB_NULL=1 -- \
    "$PY" scripts/ab_cluster.py --sub "$SUB" --auto "$OUT/$SUB/null1/auto_window" \
      --out "$OUT/$SUB/null1/clusters_sweep" --sweep-only || finish
done

# =====================================================================================
# STAGE 4 -- sec.3.2 cluster analysis: Rowe et al. (2002) pair-group (UPGMA) dendrogram, cut at
#            median + 8 MAD similarity, every cluster (>= 2) stacked over the pool.
# =====================================================================================
for SUB in $SUBS; do
  run "stage4_cluster_$SUB" "$OUT/$SUB/clusters/summary.json" AB_STAGE=4 -- \
    "$PY" scripts/ab_cluster.py --sub "$SUB" \
      --null-sweep "$OUT/$SUB/null1/clusters_sweep/sweep_long.csv" || finish
done

# =====================================================================================
# v4 (2026-09-25, user): WASHINGTON FIRST, no scrambled-control reruns on the critical path.
# Per subarray in order NW, SW, NVI: STAGE 5 growth (sec.3.3) then STAGE 6 location (Hyp2000 +
# hypoDD, Table S1; GSC 1-D headline; L9/L12 secondary). AB_ROYER_WA_DONE after NW+SW.
# The stage-7 whole-chain scrambled controls are removed (run later on finished families if wanted).
# =====================================================================================
grow_locate() {
  local SUB=$1 LAST
  run "stage5_grow_$SUB" "$OUT/$SUB/grow/summary.json" AB_STAGE=5 -- \
    "$PY" scripts/ab_grow.py --sub "$SUB" --max-iters 8 --tol 0.01 || finish
  LAST=$(ls "$OUT/$SUB/grow"/det_it*.parquet 2>/dev/null \
         | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
  if [ -z "${LAST:-}" ]; then STAGE="stage6_locate_$SUB"; finish; fi
  run "stage6_hypo_$SUB" "$OUT/$SUB/hypo/loc_summary.csv" AB_STAGE=6 -- \
    "$PY" scripts/ab_hypo.py --sub "$SUB" --det "$OUT/$SUB/grow/det_it${LAST}.parquet" \
      --out "$OUT/$SUB/hypo" --boot 60 || finish
  run "stage6b_locate_secondary_$SUB" "$OUT/$SUB/locate/depth_summary.csv" AB_STAGE=6 -- \
    "$PY" scripts/ab_locate.py --sub "$SUB" \
      --det "$OUT/$SUB/grow/det_it${LAST}.parquet" \
      --tpl "$OUT/$SUB/grow/tpl_it${LAST}.npz" --boot 60 || true
}
grow_locate NW
grow_locate SW
STAGE=stage6_score_WA; "$PY" scripts/ab_score.py --sub NW,SW || true
date -u +%FT%TZ > "$RES/AB_ROYER_WA_DONE"
hb "WASHINGTON (NW+SW) GROWN + LOCATED"
grow_locate NVI
STAGE=stage8_score
"$PY" scripts/ab_score.py --sub NW,SW,NVI || finish
date -u +%FT%TZ > "$RES/AB_ROYER_DONE"
hb "LITERAL REPRODUCTION DONE (NW, SW, NVI)"

# =====================================================================================
# EXTENSION (user-approved; results/ab_royer/ext/, never overwrites literal outputs).
# Literal autocorrelation detections and initial templates unchanged; public stations within
# 75/100/150 km added ONLY at growth and location.  The null gets the identical treatment.
# Markers: AB_ROYER_EXT_DONE / AB_ROYER_EXT_FAILED.
# =====================================================================================
EXTF="$RES/AB_ROYER_EXT_FAILED"
ext_fail() { { echo "stage=$STAGE"; date -u +%FT%TZ; } > "$EXTF"; echo "[ab] EXT FAILED in $STAGE"; exit 1; }
rm -f "$EXTF"
EXT=$OUT/ext
run ext0_pool "$EXT/pool_NW.csv" AB_STAGE=E0 -- "$PY" scripts/ab_ext_pool.py || ext_fail
run ext1_fetch "$EXT/FETCH_DONE" AB_STAGE=E1 -- "$PY" scripts/ab_fetch.py --ext --workers 8 || ext_fail
date -u +%FT%TZ > "$EXT/FETCH_DONE"
run ext2_cache "$EXT/CACHE_DONE" AB_STAGE=E2 -- "$PY" scripts/ab_cache40.py --ext --sub NW,SW,NVI --workers 24 || ext_fail
date -u +%FT%TZ > "$EXT/CACHE_DONE"
ext_chain() {   # <SUB> <R> <K: 0 = real>
  local SUB=$1 R=$2 K=$3 SRC DST NV
  if [ "$K" = "0" ]; then SRC=$OUT/$SUB; DST=$EXT/R$R/$SUB; NV=0
  else SRC=$OUT/$SUB/null$K; DST=$EXT/R$R/$SUB/null$K; NV=$K; fi
  [ -f "$SRC/clusters/tpl_it0.npz" ] || return 0
  run "ext_grow_R${R}_${SUB}_k$K" "$DST/grow/summary.json" AB_STAGE=E3 AB_NULL=$NV -- \
    "$PY" scripts/ab_grow.py --sub "$SUB" --tpl "$SRC/clusters/tpl_it0.npz" --pool "ext:$R" \
      --out "$DST/grow" --max-iters 8 --tol 0.01 || return 1
  local LAST
  LAST=$(ls "$DST/grow"/det_it*.parquet 2>/dev/null | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
  [ -n "${LAST:-}" ] || return 0
  run "ext_hypo_R${R}_${SUB}_k$K" "$DST/hypo/loc_summary.csv" AB_STAGE=E4 AB_NULL=$NV -- \
    "$PY" scripts/ab_hypo.py --sub "$SUB" --det "$DST/grow/det_it${LAST}.parquet" \
      --pool "ext:$R" --out "$DST/hypo" --boot 60 || return 1
}
for K in 0; do   # v4: scrambled controls removed from the extension too
  for R in 75 100 150; do
    for SUB in NW SW NVI; do
      ext_chain "$SUB" "$R" "$K" || ext_fail
    done
  done
  STAGE="ext_score_k$K"; "$PY" scripts/ab_ext_score.py || ext_fail
done
date -u +%FT%TZ > "$RES/AB_ROYER_EXT_DONE"
kill -TERM $HB 2>/dev/null
hb "ALL STAGES DONE (literal + extension)"
exit 0
