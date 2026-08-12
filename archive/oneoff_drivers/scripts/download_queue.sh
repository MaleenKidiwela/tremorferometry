#!/bin/bash
cd /home/jovyan/tremorferometry
for STA in B928 B004 B001 B003 B005 B013 B014; do
  echo "=== [queue] downloading $STA $(date) ==="
  /home/jovyan/envs/tremorferometry/bin/python scripts/download_borehole_3comp.py \
    --net PB --sta "$STA" --start 2007-01-01 --end 2026-05-31 --workers 10
  echo "=== [queue] $STA done $(date) ==="
done
echo "=== [queue] ALL DOWNLOADS COMPLETE ==="
