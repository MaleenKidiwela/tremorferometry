#!/usr/bin/env bash
# GPU-GATE WATCHDOG -- a dead Step 1b must never strand Step 1c.
#
# `w_step1c_queue.sh` Phase B waits for results/STEP1B_DONE or results/STEP1B_FAILED before it
# touches the GPU.  The Step-1b queue now traps TERM/INT/HUP and writes STEP1B_FAILED itself, but
# a SIGKILL (or a vanished container) still leaves no marker.  This watchdog implements the
# release rule the marker cannot:
#
#   RELEASE RULE -- if, for RELEASE_AFTER_S (600 s) of continuous polling,
#       (i)  no `w_step1b_queue.sh` process exists, AND
#       (ii) neither STEP1B_DONE nor STEP1B_FAILED exists, AND
#       (iii) the GPU is idle (utilisation < 5 % and < 2000 MiB in use),
#   then write results/STEP1B_FAILED ("released by the watchdog: no Step-1b process, GPU idle")
#   so Step 1c proceeds.  The watchdog never kills anything and never touches results it did not
#   create.
#
# It is deliberately a SEPARATE script: w_step1c_queue.sh was already running when this rule was
# added, and editing a bash script that bash is currently reading corrupts its execution.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
POLL_S=${POLL_S:-60}
RELEASE_AFTER_S=${RELEASE_AFTER_S:-600}
idle_for=0
echo "[watchdog] started pid $$ (poll ${POLL_S}s, release after ${RELEASE_AFTER_S}s idle)"
while true; do
  if [ -f results/STEP1B_DONE ] || [ -f results/STEP1B_FAILED ]; then
    echo "[watchdog] $(date -u +%FT%TZ) Step-1b marker present -- nothing to do, exiting"
    exit 0
  fi
  if ! pgrep -f "[w]_step1b_queue.sh" >/dev/null 2>&1; then
    used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1)
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null | head -1)
    used=${used:-0}; util=${util:-0}
    if [ "$used" -lt 2000 ] && [ "$util" -lt 5 ]; then
      idle_for=$((idle_for + POLL_S))
      echo "[watchdog] $(date -u +%FT%TZ) no Step-1b process, GPU idle (${used} MiB, ${util} %) -- ${idle_for}/${RELEASE_AFTER_S}s"
    else
      idle_for=0
    fi
  else
    idle_for=0
  fi
  if [ "$idle_for" -ge "$RELEASE_AFTER_S" ]; then
    { echo "step: <watchdog release>";
      echo "reason: no w_step1b_queue.sh process and an idle GPU for ${RELEASE_AFTER_S}s, and no";
      echo "        STEP1B_DONE/STEP1B_FAILED marker -- Step 1b died without running its trap";
      echo "        (SIGKILL or a vanished container).  Released so Step 1c is not stranded.";
      echo "time: $(date -u +%FT%TZ)"; echo "watchdog_pid: $$";
      echo "--- last 60 lines of logs/step1b_queue.log ---";
      tail -60 logs/step1b_queue.log 2>/dev/null; } > results/STEP1B_FAILED
    echo "[watchdog] $(date -u +%FT%TZ) RELEASED Step 1c (wrote results/STEP1B_FAILED)"
    exit 0
  fi
  sleep "$POLL_S"
done
