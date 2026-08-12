#!/bin/bash
# Build PGC mirror dv/v (both eras) = the artifact/circularity NULL for the B011-vs-PGC validation (Merlin
# step 2). Mirror = pre-arrival noise time-flipped into the coda slot, run through the IDENTICAL stretch, so
# it shares PGC's triggers/days/noise-field but carries NO LFE energy. PGC's coda dv/v must beat this null
# (and beat the same null correlated against B011) to count as real. Per-era npz -> inherent per-era ref.
# Usage: nohup bash scripts/pgc_mirror_build.sh &
cd /home/jovyan/tremorferometry
PY=/home/jovyan/envs/tremorferometry/bin/python
log(){ echo "$(date +%H:%M) [pgc-mirror] $*" | tee -a logs/pgc_mirror_build.log; }
for TAG in PGCbhz PGChhz; do
  log "$TAG: building mirror npz (pre-arrival [-3,-1]s -> coda slot, time-flipped)"
  $PY scripts/build_mirror_npz.py $TAG > logs/pgc_mirror_npz_$TAG.log 2>&1 || { log "$TAG mirror npz FAILED"; exit 1; }
  log "$TAG: mirror dv/v (2-4s, origin-anchor, per-era ref) — same params as the coda run"
  $PY scripts/dvv_roll30cal.py --station PGC --npz data/long_window_daily_${TAG}_MIRROR.npz \
     --window 2 4 --origin-anchor --workers 10 --out data/daily_dvv_${TAG}_MIRROR_2to4.csv \
     > logs/pgc_mirror_dvv_$TAG.log 2>&1 || { log "$TAG mirror dvv FAILED"; exit 1; }
  log "$TAG mirror DONE: $(($(wc -l < data/daily_dvv_${TAG}_MIRROR_2to4.csv)-1)) dv/v rows"
done
log "PGC mirror COMPLETE (both eras) — next: assemble 6 series + mirror-corrected residual corr vs mirror null"
