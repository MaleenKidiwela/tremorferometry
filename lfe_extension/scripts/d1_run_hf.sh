#!/bin/bash
# Track D, D1(a) -- the frozen HF-ratio veto applied to the four epoch-ablation passes, so that the
# per-family per-epoch recall is measured at EXACTLY the production operating point (tier AND hf_pass).
# Same script, same band (PGC 10-19 Hz / 2-8 Hz) and same cut (3.6384, the BH epochs) as Track A step 1(b).
cd /home/jovyan/tremorferometry/lfe_extension
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=4
O=/home/jovyan/tremorferometry/lfe_extension/results
run () { python scripts/a1_hf_veto.py --det $O/e5r_detections_$1.parquet --tag $1 --L 12 --year $2 \
           > logs/d1_hf_$1.log 2>&1; echo "$1 exit=$? $(date -u +%FT%TZ)" >> logs/d1_hf.log; }
run d1_abl2012_E2 2012 &
run d1_abl2012_E3 2012 &
wait
run d1_abl2010_E2 2010 &
run d1_abl2010_E3 2010 &
wait
echo "[$(date -u +%FT%TZ)] D1(a) HF passes DONE" >> logs/d1_hf.log
