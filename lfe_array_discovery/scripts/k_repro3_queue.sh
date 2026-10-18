#!/usr/bin/env bash
# K_REPRO3: 15-station autocorrelation over the whole episode, regroup, grow every group,
# score PRE-MERGE (the two published uniqueness tests were measured not to separate families),
# plus the catalogue-built ceiling and a looser-first-pass variant.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/repro
exec 9>locks/k_repro3.lock
flock -n 9 || { echo "already running"; exit 0; }
rm -f results/K_REPRO3_DONE results/K_REPRO3_FAILED
STAGE=start
trap 'echo "[q3] TRAP exit $?"; [ -f results/K_REPRO3_DONE ] || echo "stage=$STAGE" > results/K_REPRO3_FAILED' EXIT
S15="GLBC,TWKB,MGCB,TSJB,TWBB,PGC,SILB,SSIB,LCBC,B003,B005,JRBC,LZB,VGZ,SNB"
DATES=$(python - <<'PY'
import pandas as pd
print(",".join(d.strftime("%Y-%m-%d") for d in pd.date_range("2005-08-29","2005-09-21")))
PY
)
echo "[q3] start $(date -u +%FT%TZ) pid=$$"

STAGE=autocorr15
python scripts/k_auto.py --dates "$DATES" --hours all --stas "$S15" \
  --out results/repro/auto15 > logs/k_auto15_full.log 2>&1 || exit 1
python - <<'PY' || exit 1
import glob, pandas as pd, numpy as np
A=[pd.read_parquet(f) for f in sorted(glob.glob('results/repro/auto15/*.parquet'))]
R=pd.concat([x for x in A if len(x)],ignore_index=True)
S=pd.read_csv('results/repro/auto15/_stats.csv')
bad={(r.date,int(r.hour)) for r in S.itertuples() if not (r.nmad_in_sd>3)}
R=R[~pd.Series(list(zip(R.date,R.hour)),index=R.index).isin(bad)]
R=R[np.isfinite(R.nmad)].reset_index(drop=True)
R.to_parquet('results/repro/auto15_all.parquet',index=False)
print("15-station detections",len(R),"bad hours",len(bad))
PY
echo "[q3] autocorr15 done $(date -u +%FT%TZ)"

STAGE=score_pairs15
python - <<'PY' > logs/k_score_pairs15.log 2>&1
import sys; sys.path.insert(0,'scripts')
import numpy as np, pandas as pd, k_common as K, k_score_pairs as KS
R=pd.read_parquet('results/repro/auto15_all.parquet')
C=pd.read_csv('results/repro/s1all_times.csv'); C['t_corr']=C.t_origin+C.dt
TT=K.fam_tt(K.CLEAN11)
TTM={int(b):float(np.median([TT[(int(b),s)] for s in K.CLEAN11])) for b in C.bkey.unique()}
M=KS.score(R,C,TTM,0.6,8.0); N=KS.score(R,C,TTM,0.6,8.0,shift=300.)
M.to_csv('results/repro/auto15_matched.csv',index=False)
print(f"{len(R)} det; matched {len(M)} ({len(M)/len(R):.2%}) fams {M.bkey.nunique()}; chance {len(N)/len(R):.2%}")
for q in (100,500,2000,10000,50000):
    if len(R)>=q:
        Rt=R.nlargest(q,'nmad'); Mt=KS.score(Rt,C,TTM,0.6,8.0)
        print(f"top {q}: {len(Mt)/q:.1%} fams {Mt.bkey.nunique()}")
PY
cat logs/k_score_pairs15.log

STAGE=group15
python scripts/k_group.py --detfile results/repro/auto15_all.parquet --budget 1000000 \
  --min-members 3 --out results/repro/gb15 > logs/k_group15.log 2>&1 || exit 1
grep -E "nodes,|purity>=0.5" logs/k_group15.log

STAGE=grow15
python scripts/k_grow.py --nodes results/repro/gb15_nodes.csv --summary results/repro/gb15_summary.csv \
  --stas "$S15" --iters 1 --out results/repro/grow15 > logs/k_grow15.log 2>&1 || exit 1
for it in 0 1; do
  python scripts/k_score_fam.py --det results/repro/grow15/det_it$it.parquet \
    --out results/repro/grow15/score_pre_it$it.csv > logs/k_score15_it$it.log 2>&1
  echo "== grow15 iter $it"; grep "grown groups" logs/k_score15_it$it.log
done

STAGE=grow15_loose
python scripts/k_grow.py --nodes results/repro/gb15_nodes.csv --summary results/repro/gb15_summary.csv \
  --stas "$S15" --iters 1 --nmad-first 6.0 --out results/repro/grow15L > logs/k_grow15L.log 2>&1 || true
for it in 0 1; do
  [ -s results/repro/grow15L/det_it$it.parquet ] || continue
  python scripts/k_score_fam.py --det results/repro/grow15L/det_it$it.parquet \
    --out results/repro/grow15L/score_pre_it$it.csv > logs/k_score15L_it$it.log 2>&1
  echo "== grow15 LOOSE-first iter $it"; grep "grown groups" logs/k_score15L_it$it.log
done

STAGE=ceiling
python scripts/k_ceiling.py --stas "$S15" --iters 1 --out results/repro/ceiling15 \
  > logs/k_ceiling15.log 2>&1 || true
for it in 0 1; do
  [ -s results/repro/ceiling15/det_it$it.parquet ] || continue
  python scripts/k_score_fam.py --det results/repro/ceiling15/det_it$it.parquet \
    --out results/repro/ceiling15/score_it$it.csv > logs/k_ceiling_score_it$it.log 2>&1
  echo "== CEILING iter $it"; grep -A3 "grown groups" logs/k_ceiling_score_it$it.log | head -4
  grep "per-family medians" logs/k_ceiling_score_it$it.log
done

STAGE=loss
python scripts/k_loss.py --matched results/repro/auto15_matched.csv \
  --purity results/repro/gb15_purity.csv --score results/repro/grow15/score_pre_it0.csv \
  --out results/repro/loss15.csv > logs/k_loss15.log 2>&1
cat logs/k_loss15.log

echo "OK $(date -u +%FT%TZ)" > results/K_REPRO3_DONE
echo "[q3] DONE $(date -u +%FT%TZ)"
