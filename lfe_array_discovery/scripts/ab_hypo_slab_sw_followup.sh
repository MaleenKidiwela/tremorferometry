#!/usr/bin/env bash
# Detached follow-up: wait for the ab_royer queue's SW free-depth headline
# (results/ab_royer/SW/hypo/loc_headline.csv), then run ab_hypo_slab.py for SW.
# Markers: results/ab_royer/SW/hypo_slab/DONE | FAILED; heartbeat HEARTBEAT; pid PID.
# Gives up (FAILED) if the queue has written AB_ROYER_FAILED or its pid is gone without the file.
# Launch: nohup setsid -f bash scripts/ab_hypo_slab_sw_followup.sh > logs/ab_hypo_slab_SW.log 2>&1
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
RES=$PWD/results
OUT=$RES/ab_royer/SW/hypo_slab
HEAD=$RES/ab_royer/SW/hypo/loc_headline.csv
mkdir -p "$OUT" locks logs
exec 8>locks/ab_hypo_slab_SW.lock
if ! flock -n 8; then echo "another SW slab follow-up holds the lock; exiting"; exit 0; fi
echo $$ > "$OUT/PID"
rm -f "$OUT/DONE" "$OUT/FAILED"
hb_loop() { while :; do date -u +%FT%TZ > "$OUT/HEARTBEAT"; sleep 120; done; }
hb_loop & HB=$!
STAGE=wait
finish() {
  rc=$?
  kill -TERM $HB 2>/dev/null
  if [ ! -f "$OUT/DONE" ]; then
    { echo "stage=$STAGE"; echo "rc=$rc"; date -u +%FT%TZ; } > "$OUT/FAILED"
  fi
}
trap finish EXIT
trap 'STAGE="$STAGE(signalled)"; exit 143' TERM INT HUP
QPID=$(cat "$RES/AB_ROYER_PID" 2>/dev/null)
echo "[slabSW] $(date -u +%FT%TZ) waiting for $HEAD (queue pid $QPID)"
while [ ! -s "$HEAD" ]; do
  if [ -f "$RES/AB_ROYER_FAILED" ]; then echo "[slabSW] queue FAILED marker; giving up"; exit 2; fi
  if [ -n "$QPID" ] && ! kill -0 "$QPID" 2>/dev/null && [ ! -s "$HEAD" ]; then
    sleep 60; [ -s "$HEAD" ] || { echo "[slabSW] queue pid gone, no headline; giving up"; exit 3; }
  fi
  sleep 300
done
# the headline csv is written in one to_csv call; let it settle
sleep 60
STAGE=locate
echo "[slabSW] $(date -u +%FT%TZ) headline present; running ab_hypo_slab.py --sub SW"
python scripts/ab_hypo_slab.py --sub SW || exit $?
date -u +%FT%TZ > "$OUT/DONE"
echo "[slabSW] $(date -u +%FT%TZ) DONE"
