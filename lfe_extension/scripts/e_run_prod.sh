#!/bin/bash
# Track E — the per-site pass over the WHOLE production record, 2 concurrent GPU passes.
# Resumable: a month whose persite_<ym>_meta.json exists is skipped.
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=3
N=${1:-2}
ls results/trackC/detections_*.parquet | sed 's/.*detections_//; s/\.parquet//' | sort > /tmp/e_months.txt
: > logs/e_prod_queue.txt
while read ym; do
  [ -f results/trackE/persite_${ym}_meta.json ] || echo $ym >> logs/e_prod_queue.txt
done < /tmp/e_months.txt
echo "[$(date -u +%FT%TZ)] $(wc -l < logs/e_prod_queue.txt) months to run, $N workers"
cat logs/e_prod_queue.txt | xargs -P $N -I{} sh -c \
  'python scripts/e_persite.py --mode prod --ym {} --tag {} > logs/e_persite_{}.log 2>&1 || echo "FAIL {}" >> logs/e_prod_fail.txt'
echo "[$(date -u +%FT%TZ)] DONE $(ls results/trackE/persite_*_meta.json | wc -l) meta files"
