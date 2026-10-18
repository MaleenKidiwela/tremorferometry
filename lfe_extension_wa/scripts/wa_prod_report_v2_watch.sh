#!/bin/bash
# Wait (read-only, every 10 min, at most 21 days) for the production queue's WA_PROD_DONE marker, which the queue writes
# right after its final report; then append the post-processing v2 link section (scripts/wa_prod_report_v2.py).
# Touches no queue marker. Log: results/wa_prod/post_v2/report_v2_watch.log
cd /home/jovyan/tremorferometry/lfe_extension_wa || exit 1
O=results/wa_prod
LOGW=$O/post_v2/report_v2_watch.log
echo "$(date -u +%FT%TZ) watcher start pid $$" >> $LOGW
for i in $(seq 1 3024); do
  if [ -f $O/WA_PROD_DONE ]; then
    sleep 60   # let the queue finish writing
    python scripts/wa_prod_report_v2.py >> $LOGW 2>&1
    echo "$(date -u +%FT%TZ) watcher done" >> $LOGW
    exit 0
  fi
  sleep 600
done
echo "$(date -u +%FT%TZ) watcher gave up after 21 days" >> $LOGW
