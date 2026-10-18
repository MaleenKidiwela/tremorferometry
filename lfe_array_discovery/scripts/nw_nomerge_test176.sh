#!/usr/bin/env bash
# Step 1 of Merlin 09-28 plan: no-merge growth of the 176 NW_mute it0 seeds with >=2 impulsive sites, 3 iterations.
set -u -o pipefail
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p locks logs
exec 9>locks/nw_nomerge_test.lock; flock -n 9 || exit 0
T=results/ab_royer/NW_nomerge/test176
export AB_MUTE_DIR=$PWD/results/ab_royer/NW_mute/mute
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
rm -f $T/TEST_DONE $T/TEST_FAILED; echo $$ > $T/PID
python scripts/ab_grow_fast.py --sub NW --tpl $T/seeds176.npz --out $T/grow --no-merge \
  --max-iters 3 --tol 0.0 --prep-workers 6 --restack-workers 8 || { echo grow > $T/TEST_FAILED; exit 1; }
python scripts/ab_restack28.py --sub NW --det $T/grow/det_it3.parquet --out $T/stacks28_it3.npz || { echo r28 > $T/TEST_FAILED; exit 1; }
date -u +%FT%TZ > $T/TEST_DONE
