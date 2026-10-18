#!/usr/bin/env bash
# MENDO (2026-09-29): blind LFE family DISCOVERY at the Mendocino triple junction (Shelly et al. 2025 region)
# with the Royer chain, recipe copied from scripts/lag15_nvi_sw_queue.sh (NW_lag15 / SW_lag15 / NVI_lag15):
#   40 sps 1-8 Hz cache from raw (ab_cache40.py; lfe_south/traces, 25 stations) + MENDO dead-channel / stuck-span
#   guard (ab_common.dead_guard, port of lfe_south SH_V2 step 1) -> burst mute + line notch (ab_mute.py K=15 qmax=2
#   qlong=1.3 pad 3) -> self-repeat station screen (mendo_station_screen.py; flagged stations go last in the
#   autocorrelation preference) -> ab_auto.py (15 s / 0.5 s / 8 MAD, per hour the first 7 three-component stations
#   present, 80-day contiguous autocorrelation block) -> ab_cluster.py --lag-max 15 --lag-mode zeropad
#   -> ab_grow_fast.py --merge most, 14 iters, tol 0.01, prep 6, restack 4 over all 572 scan days (runaway guard
#   15 M detections = NW's 5 M x 572/188 days) -> ab_hypo_mendo.py (Hyp2000 + hypoDD, USGS CVM v1.7 mendo_MTJ,
#   --boot 60) -> ab_restack28 (workers 4) -> review sheets -> ab_nomerge_score.py (auto_score.csv n_sites_3km)
#   -> mendo_compare.py (vs Shelly: location / waveform / shared times) -> mendo_summary.py.
# Then the SCRAMBLED CONTROL: the identical chain with AB_NULL=1 (each site's day series shifted by whole days within
# each year's scan-day list; stations < 3 km apart shift together), outputs results/ab_royer/MENDO_control/.
# First waits for lfe_south Shelly v3b (results/SH_V3B_DONE | _FAILED, or pid 3608906 gone).
# Markers results/MENDO_{DONE,FAILED,HEARTBEAT}, pid results/MENDO_PID, lock locks/mendo_discovery.lock,
# log logs/mendo_discovery_queue.log.  Stages are guarded: a relaunch skips finished stages.
#   nohup setsid bash scripts/mendo_discovery_queue.sh >> logs/mendo_discovery_queue.log 2>&1 < /dev/null &
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/mendo_discovery.lock; flock -n 9 || { echo "[queue] already running"; exit 0; }
R=results; B=results/ab_royer/MENDO; BC=results/ab_royer/MENDO_control
SOUTH=/home/jovyan/tremorferometry/lfe_south
mkdir -p $B $BC
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1
unset AB_NULL AB_MENDO_SCREEN AB_MUTE_DIR
( exec 9>&-; while true; do date -u +%FT%TZ > $R/MENDO_HEARTBEAT; sleep 120; done ) &
HB=$!
GU=
trap 'kill $HB $GU 2>/dev/null' EXIT
rm -f $R/MENDO_FAILED $R/MENDO_DONE
echo "$$" > $R/MENDO_PID
fail(){ echo "stage=$1 $(date -u +%FT%TZ)" > $R/MENDO_FAILED; echo "[queue] FAILED stage=$1 $(date -u +%FT%TZ)"; exit 1; }
echo "[queue] start $(date -u +%FT%TZ) pid $$"

# ---- wait for the Shelly v3b run to release the GPU ----
while true; do
  [ -f $SOUTH/results/SH_V3B_DONE ] && { echo "[queue] SH_V3B_DONE seen"; break; }
  [ -f $SOUTH/results/SH_V3B_FAILED ] && { echo "[queue] SH_V3B_FAILED seen"; break; }
  kill -0 3608906 2>/dev/null || { echo "[queue] v3b pid 3608906 gone"; break; }
  sleep 300
done
echo "[queue] GPU wait over $(date -u +%FT%TZ)"
df -h /home/jovyan | tail -1

# ---- shared data stages (real + control) ----
python scripts/mendo_setup.py 9>&- || fail setup
[ -f "$B/CACHE_DONE" ] || { python scripts/ab_cache40.py --sub MENDO --workers 24 9>&- && date -u +%FT%TZ > "$B/CACHE_DONE"; } || fail cache
echo "[queue] cache done $(date -u +%FT%TZ)"; df -h /home/jovyan | tail -1
python scripts/mendo_cvm_profile.py 9>&- || fail cvm_profile
for M in energy lines mask; do
  [ -f "$B/mute/${M}_DONE" ] || { python scripts/ab_mute.py $M --sub MENDO --out "$B/mute" --K 15 --qmax 2 --qlong 1.3 --pad 3 --workers 16 9>&- \
    && date -u +%FT%TZ > "$B/mute/${M}_DONE"; } || fail mute_$M
done
[ -f "$B/mute/params.json" ] || fail mask_missing
export AB_MUTE_DIR=$PWD/$B/mute
[ -f results/ab_royer/tremor/auto_screen_MENDO.csv ] || python scripts/mendo_station_screen.py --every 8 --hours 0,6,12,18 9>&- || fail screen
python -c "import sys; sys.path.insert(0,'scripts'); import ab_common as A; print('[queue] auto preference', A.SUBARRAY['MENDO']['stas']); print('[queue] pool', A.POOL['MENDO'])" 9>&-

run_chain(){
  local D=$1; local TAG=$2
  [ -f "$D/CHAIN_DONE" ] && { echo "[$TAG] already done"; return 0; }
  echo "[$TAG] start $(date -u +%FT%TZ) AB_NULL=${AB_NULL:-0}"
  [ -f "$D/auto_window/_params.json" ] || \
    python scripts/ab_auto.py --sub MENDO --cols window --out "$D/auto_window" 9>&- || fail ${TAG}_auto
  [ -f "$D/clusters/summary.json" ] || \
    python scripts/ab_cluster.py --sub MENDO --auto "$D/auto_window" --out "$D/clusters" \
      --lag-max 15 --lag-mode zeropad \
      --headline-from results/ab_royer/NW/clusters/summary.json 9>&- || fail ${TAG}_cluster
  echo "[$TAG] cluster done $(date -u +%FT%TZ)"
  if [ ! -f "$D/grow/GROW_DONE" ]; then
    python scripts/ab_grow_fast.py --sub MENDO --tpl "$D/clusters/tpl_it0.npz" --out "$D/grow" \
      --merge most --max-iters 14 --tol 0.01 --prep-workers 6 --restack-workers 4 9>&- &
    local GP=$!; echo "$GP" > $R/MENDO_${TAG}_GROW_PID
    ( exec 9>&-; while kill -0 $GP 2>/dev/null; do
        if [ -f "$D/grow/iters.csv" ] && python -c "import pandas as p,sys; d=p.read_csv('$D/grow/iters.csv'); sys.exit(0 if d.n_detections.max()>15e6 else 1)"; then
          echo "runaway $(date -u +%FT%TZ)" > $R/MENDO_${TAG}_RUNAWAY; kill $(pgrep -P $GP) $GP 2>/dev/null; break; fi
        sleep 60; done ) &
    GU=$!
    wait $GP; local RC=$?
    kill $GU 2>/dev/null; GU=
    [ -f $R/MENDO_${TAG}_RUNAWAY ] && fail ${TAG}_runaway
    [ $RC -eq 0 ] || fail ${TAG}_grow
    date -u +%FT%TZ > "$D/grow/GROW_DONE"
  fi
  local LAST=$(ls $D/grow/det_it*.parquet | sed 's/.*det_it//; s/\.parquet//' | sort -n | tail -1)
  echo "[$TAG] grow done, last iteration $LAST $(date -u +%FT%TZ)"
  [ -f "$D/hypo/loc_headline.csv" ] || \
    python scripts/ab_hypo_mendo.py --sub MENDO --det "$D/grow/det_it${LAST}.parquet" --out "$D/hypo" --boot 60 9>&- || true
  [ -f "$D/hypo/loc_headline.csv" ] || fail ${TAG}_hypo
  [ -f "$D/hypo_28/stacks_loc.npz" ] || \
    python scripts/ab_restack28.py --sub MENDO --det "$D/grow/det_it${LAST}.parquet" --out "$D/hypo_28/stacks_loc.npz" --workers 4 9>&- \
    || fail ${TAG}_restack28
  [ -f "$D/hypo_28/review/SHEETS_DONE" ] || { python scripts/ab_review_sheets.py --stacks "$D/hypo_28/stacks_loc.npz" \
    --det "$D/grow/det_it${LAST}.parquet" --out "$D/hypo_28/review" 9>&- && date -u +%FT%TZ > "$D/hypo_28/review/SHEETS_DONE"; } \
    || fail ${TAG}_sheets
  python scripts/ab_nomerge_score.py "$D/hypo_28/stacks_loc.npz" "$D/hypo_28/review/auto_score.csv" MENDO 9>&- \
    || fail ${TAG}_score
  date -u +%FT%TZ > "$D/CHAIN_DONE"
  echo "[$TAG] done $(date -u +%FT%TZ)"
}

run_chain $B REAL
python scripts/mendo_compare.py $B 9>&- || fail compare
python scripts/mendo_summary.py 9>&- || fail summary
export AB_NULL=1
run_chain $BC CONTROL
unset AB_NULL
python scripts/mendo_summary.py 9>&- || fail summary_final
date -u +%FT%TZ > $R/MENDO_DONE
echo "[queue] end $(date -u +%FT%TZ)"
