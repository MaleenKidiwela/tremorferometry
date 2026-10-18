#!/bin/bash
# Track D — rebuild the whole product from whatever Track C months are on disk.
# Idempotent: every step overwrites its own outputs under results/trackD/ and figures/d*; nothing else
# is touched.  The GPU epoch ablation (scripts/d1_run_ablation.sh) and its HF passes are NOT repeated
# here — they depend only on the referee years and are already on disk.
set -u
cd /home/jovyan/tremorferometry/lfe_extension
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=4
P=/opt/conda/bin/python
L=logs/d_run_all.log
echo "[$(date -u +%FT%TZ)] D chain START" >> $L
$P scripts/d1_score_epochs.py                >> $L 2>&1; echo "  d1_score_epochs=$?"   >> $L
$P scripts/d1_mc_model.py                    >> $L 2>&1; echo "  d1_mc_model=$?"       >> $L
$P scripts/d1_mc_apply.py                    >> $L 2>&1; echo "  d1_mc_apply=$?"       >> $L
$P scripts/d1_floor.py                       >> $L 2>&1; echo "  d1_floor=$?"          >> $L
$P scripts/d2_assemble.py                    >> $L 2>&1; echo "  d2_assemble=$?"       >> $L
$P scripts/d3_sse_windows.py --out results/trackD/sse_windows.csv                >> $L 2>&1
$P scripts/d3_sse_windows.py --max-gap-d 5 --out results/trackD/sse_windows_sub.csv >> $L 2>&1
$P - <<'PY' >> $L 2>&1
import pandas as pd
for f in ("results/trackD/sse_windows.csv", "results/trackD/sse_windows_sub.csv"):
    W = pd.read_csv(f); W["major"] = W.n_active_days >= 10; W.to_csv(f, index=False)
PY
$P scripts/d3_run_sse.py --workers "${DW:-8}" --null-draws "${DN:-20}" >> $L 2>&1; echo "  d3_run_sse=$?" >> $L
# The thinned reference is a 2012-REFEREE comparison against Bostock, so it must use the referee
# operating point (8.75 / 8.25) -- the same one `d3_run_sse.py` builds `ref2012_ours_*` at.  The
# per-epoch production thresholds do not enter here.
$P scripts/d3_thinned_reference.py --T 8.75 --tier T_core --epoch E0 --seeds 20 --workers 8 >> $L 2>&1
$P scripts/d3_thinned_reference.py --T 8.25 --tier T_ext  --epoch E0 --seeds 20 --workers 8 >> $L 2>&1
$P scripts/d1d_season_diagnostic.py          >> $L 2>&1; echo "  d1d_diagnostic=$?"    >> $L
$P scripts/d1_figures.py                     >> $L 2>&1; echo "  d1_figures=$?"        >> $L
$P scripts/d3_figures.py --tier core         >> $L 2>&1; echo "  d3_figures=$?"        >> $L
$P scripts/d3_figures.py --tier core_hp      >> $L 2>&1
$P scripts/d3_figures.py --tier ext          >> $L 2>&1
echo "[$(date -u +%FT%TZ)] D chain DONE" >> $L
