#!/bin/bash
# runs the post-hoc kernel-error model test (s11_wa_kerr.py) after the WA-NW queue ends. Marker wa_nw/KERR_DONE.
cd /home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/s11_depth/wa_nw
while pgrep -f s11_wa_queue.sh >/dev/null; do sleep 60; done
[ -f $O/KERR_DONE ] && exit 0
S11_REGION=WA_NW python s11_wa_kerr.py --workers 8 >> $O/logs/kerr_stage.out 2>&1 && touch $O/KERR_DONE
echo "$(TZ=America/Los_Angeles date '+%F %-I:%M %p PT') kerr chain end rc=$?" >> $O/logs/queue.log
