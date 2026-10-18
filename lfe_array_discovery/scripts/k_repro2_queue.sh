#!/usr/bin/env bash
# K_REPRO2: rebuild the 7 numerically-broken hours, regroup at full budget, grow EVERY group
# through 3 iterations of the literal chain, dedupe, score with Bostock's criterion.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/repro
exec 9>locks/k_repro2.lock
flock -n 9 || { echo "already running"; exit 0; }
rm -f results/K_REPRO2_DONE results/K_REPRO2_FAILED
STAGE=start
trap 'echo "[q2] TRAP exit $?"; [ -f results/K_REPRO2_DONE ] || echo "stage=$STAGE" > results/K_REPRO2_FAILED' EXIT
echo "[q2] start $(date -u +%FT%TZ) pid=$$"

STAGE=redo_bad_hours
python scripts/k_auto.py --dates 2005-08-30 --hours 11 --out results/repro/auto_fix || exit 1
python scripts/k_auto.py --dates 2005-09-02 --hours 2  --out results/repro/auto_fix || exit 1
python scripts/k_auto.py --dates 2005-09-07 --hours 16 --out results/repro/auto_fix || exit 1
python scripts/k_auto.py --dates 2005-09-14 --hours 20 --out results/repro/auto_fix || exit 1
python scripts/k_auto.py --dates 2005-09-19 --hours 16,17 --out results/repro/auto_fix || exit 1
python scripts/k_auto.py --dates 2005-09-21 --hours 22 --out results/repro/auto_fix || exit 1
python - <<'PY' || exit 1
import glob, pandas as pd, numpy as np
A=[pd.read_parquet('results/repro/auto_full_clean.parquet')]
for f in sorted(glob.glob('results/repro/auto_fix/*.parquet')):
    x=pd.read_parquet(f)
    if len(x): A.append(x)
R=pd.concat(A,ignore_index=True)
R=R[np.isfinite(R.nmad)].reset_index(drop=True)
R.to_parquet('results/repro/auto_all.parquet',index=False)
print("merged detections", len(R), "hours", R.groupby(['date','hour']).ngroups)
PY
echo "[q2] hours fixed $(date -u +%FT%TZ)"

STAGE=group
python scripts/k_group.py --detfile results/repro/auto_all.parquet --budget 1000000 \
  --min-members 3 --out results/repro/gb || exit 1

STAGE=grow
python scripts/k_grow.py --nodes results/repro/gb_nodes.csv --summary results/repro/gb_summary.csv \
  --iters 3 --out results/repro/grow || exit 1

STAGE=score
for it in 0 1 2 3; do
  [ -s results/repro/grow/det_it$it.parquet ] || continue
  python scripts/k_dedupe.py --det results/repro/grow/det_it$it.parquet \
     --out results/repro/grow/det_it${it}_dd.parquet
  echo "===== iteration $it (deduped)"
  python scripts/k_score_fam.py --det results/repro/grow/det_it${it}_dd.parquet \
     --out results/repro/grow/score_it${it}.csv
done 2>&1 | tee logs/k_score_iters.log

STAGE=loss
python - <<'PY'
import pandas as pd, glob, numpy as np, sys
sys.path.insert(0,'scripts'); import k_common as K, k_score_pairs as KS
R=pd.read_parquet('results/repro/auto_all.parquet')
C=pd.read_csv('results/repro/s1all_times.csv'); C['t_corr']=C.t_origin+C.dt
TT=K.fam_tt(K.CLEAN11)
TTM={int(b):float(np.median([TT[(int(b),s)] for s in K.CLEAN11])) for b in C.bkey.unique()}
M=KS.score(R,C,TTM,0.6,8.0); M.to_csv('results/repro/auto_all_matched.csv',index=False)
print("pair-stage matched", len(M), "of", len(R), "over", M.bkey.nunique(), "kept families")
PY
LAST=$(ls results/repro/grow/score_it*.csv | tail -1)
python scripts/k_loss.py --matched results/repro/auto_all_matched.csv \
  --purity results/repro/gb_purity.csv --score "$LAST" --out results/repro/loss.csv 2>&1 | tee logs/k_loss.log

echo "OK $(date -u +%FT%TZ)" > results/K_REPRO2_DONE
echo "[q2] DONE $(date -u +%FT%TZ)"
