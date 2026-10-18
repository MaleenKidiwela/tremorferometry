#!/bin/bash
# Detached queue for scripts/wa_dl2.py (tiers C, D, D2, T). Launch:
#   nohup setsid -f bash scripts/wa_dl2_queue.sh > /dev/null 2>&1
# Markers in results/wa_dl2/: WA_DL2_PID, WA_DL2_HEARTBEAT (python progress), WA_DL2_ALIVE (shell, 5 min),
# WA_DL2_DONE | WA_DL2_FAILED (+ STOPPED_STORAGE if the 2.6 TB cap was hit). Log logs/wa_dl2.log.
cd /home/jovyan/tremorferometry/lfe_extension_wa || exit 1
O=results/wa_dl2; L=logs/wa_dl2.log; X=logs/wa_dl2.out
echo $$ > $O/WA_DL2_PID
rm -f $O/WA_DL2_DONE $O/WA_DL2_FAILED
( while kill -0 $$ 2>/dev/null; do echo "$(date -u +%FT%TZ) queue alive pid $$" > $O/WA_DL2_ALIVE; sleep 300; done ) &
echo "$(date -u +%FT%TZ) [queue] start pid $$" >> $L
# stage 1: tiers C + D (planned and smoke-tested first)
python scripts/wa_dl2.py run --tiers C,D --workers 4 --rate 2 --cap-tb 2.6 >> $X 2>&1; rc=$?
# stage 2: everything in the plan as it is THEN (tier T rows are added to plan.csv by `wa_dl2.py plan` once the
# temporary-network inventory selection exists), then D2; the manifest skips everything already done
if [ $rc -eq 0 ]; then
  echo "$(date -u +%FT%TZ) [queue] stage 2 (re-read plan: T, D2)" >> $L
  python scripts/wa_dl2.py run --workers 4 --rate 2 --cap-tb 2.6 >> $X 2>&1; rc=$?
fi
if [ $rc -eq 0 ]; then
  echo "$(date -u +%FT%TZ) [queue] retry pass for errors" >> $L
  python scripts/wa_dl2.py run --retry-errors-only --workers 3 --rate 1 --cap-tb 2.6 >> $X 2>&1; rc=$?
fi
python scripts/wa_dl2.py summary >> $X 2>&1
if [ $rc -eq 3 ]; then          # the 2.6 TB cap stopped the run: not a failure, the rest waits for approval
  echo "$(date -u +%FT%TZ) [queue] DONE (stopped at the storage cap, see STOPPED_STORAGE)" >> $L
  echo "$(date -u +%FT%TZ) stopped at storage cap: $(cat $O/STOPPED_STORAGE)" > $O/WA_DL2_DONE
elif [ $rc -eq 0 ]; then
  echo "$(date -u +%FT%TZ) [queue] DONE" >> $L; date -u +%FT%TZ > $O/WA_DL2_DONE
else
  echo "$(date -u +%FT%TZ) [queue] FAILED rc=$rc" >> $L; echo "$(date -u +%FT%TZ) rc=$rc" > $O/WA_DL2_FAILED
fi
