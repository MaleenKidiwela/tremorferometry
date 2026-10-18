#!/bin/bash
# Amendment R2b: debug the low overlap (self-detection, stacks at Table-S1 times, overlap diagnostics, processing grid).
# Waits for R2_MF_DONE (noise templates of the chosen setting are needed for the selection-matched control).
# 7 workers (+1 for r2_queue2's stage-2 analysis) = 8. Markers R2B_SELF_DONE R2B_OVERLAP_DONE R2B_GRID_DONE R2B_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while [ ! -f R2_MF_DONE ]; do
  if ! pgrep -f "bash scripts/r2_queue2.sh" > /dev/null; then echo "queue2 ended without R2_MF_DONE"; exit 1; fi
  sleep 30
done
( python scripts/hg_r2b.py self > logs/r2b_self.log 2>&1 && touch R2B_SELF_DONE
  python scripts/hg_r2b.py overlap > logs/r2b_overlap.log 2>&1 && touch R2B_OVERLAP_DONE ) &
( python scripts/hg_r2b.py grid > logs/r2b_grid.log 2>&1 && touch R2B_GRID_DONE ) &
for K in 1 2 3 4 5; do python scripts/hg_r2b.py stack $K > logs/r2b_stack_T$K.log 2>&1 & done
wait
python scripts/hg_r2b.py stacksum > logs/r2b_stacksum.log 2>&1
n=$(ls R2B_SELF_DONE R2B_OVERLAP_DONE R2B_GRID_DONE data/r2b/stack_T?.parquet 2>/dev/null | wc -l)
[ "$n" = "8" ] && touch R2B_DONE
echo "$(TZ=America/Los_Angeles date) r2b done ($n/8 parts)"
