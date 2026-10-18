#!/usr/bin/env bash
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs; exec 9>locks/k_repro6.lock; flock -n 9 || exit 0
rm -f results/K_REPRO6_DONE results/K_REPRO6_FAILED
STAGE=start
trap 'echo "[q6] TRAP $?"; [ -f results/K_REPRO6_DONE ] || echo "stage=$STAGE" > results/K_REPRO6_FAILED' EXIT
D0=$(python -c "
import glob,os,pandas as pd
ds=sorted({int(os.path.basename(p)[:3]) for p in glob.glob('/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz/*/2010/*.npz')})
print((pd.Timestamp('2010-01-01')+pd.Timedelta(days=ds[0]-1)).strftime('%Y-%m-%d'))")
D1=$(python -c "
import glob,os,pandas as pd
ds=sorted({int(os.path.basename(p)[:3]) for p in glob.glob('/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz/*/2010/*.npz')})
print((pd.Timestamp('2010-01-01')+pd.Timedelta(days=ds[-1]-1)).strftime('%Y-%m-%d'))")
S="PGC,VGZ,SNB,NLLB,LRIV,LZB,B011,B926"
E="K_YEAR=2010 K_DAY0=$D0 K_DAY1=$D1 K_OUT=results/repro2010"
DATES=$(python - "$D0" "$D1" <<'PY'
import sys, pandas as pd
print(",".join(d.strftime("%Y-%m-%d") for d in pd.date_range(sys.argv[1], sys.argv[2])))
PY
)
echo "[q6] 2010 bed $D0..$D1 stations $S"
STAGE=auto; env $E python scripts/k_auto.py --dates "$DATES" --hours all --stas "$S" \
  --out results/repro2010/auto2 > logs/k2_2010_auto.log 2>&1 || exit 1
STAGE=merge; env $E python - <<'PY' > logs/k2_2010_merge.log 2>&1 || exit 1
import glob, pandas as pd, numpy as np
A=[pd.read_parquet(f) for f in sorted(glob.glob("results/repro2010/auto2/*.parquet"))]
R=pd.concat([x for x in A if len(x)],ignore_index=True)
S=pd.read_csv("results/repro2010/auto2/_stats.csv")
bad={(r.date,int(r.hour)) for r in S.itertuples() if not (r.nmad_in_sd>3)}
R=R[~pd.Series(list(zip(R.date,R.hour)),index=R.index).isin(bad)]
R=R[np.isfinite(R.nmad)].reset_index(drop=True)
R.to_parquet("results/repro2010/auto_all.parquet",index=False)
print("det",len(R),"hours",R.groupby(['date','hour']).ngroups,"bad",len(bad))
PY
cat logs/k2_2010_merge.log
STAGE=pairs; env $E python - "$S" <<'PY' > logs/k2_2010_pairs.log 2>&1
import sys; sys.path.insert(0,'scripts')
import numpy as np, pandas as pd, k_common as K, k_score_pairs as KS
ST=sys.argv[1].split(",")
R=pd.read_parquet("results/repro2010/auto_all.parquet")
C=pd.read_csv("results/repro2010/s1_times.csv"); C['t_corr']=C.t_origin+C.dt
TT=K.fam_tt(ST); TTM={int(b):float(np.median([TT[(int(b),s)] for s in ST])) for b in C.bkey.unique()}
M=KS.score(R,C,TTM,0.6,8.0); N=KS.score(R,C,TTM,0.6,8.0,shift=300.)
M.to_csv("results/repro2010/auto_matched.csv",index=False)
print(f"{len(R)} det; matched {len(M)} ({len(M)/max(len(R),1):.2%}) fams {M.bkey.nunique()}; chance {len(N)/max(len(R),1):.2%}")
PY
cat logs/k2_2010_pairs.log
STAGE=group; env $E python scripts/k_group.py --detfile results/repro2010/auto_all.parquet \
  --budget 1000000 --min-members 3 --stas "$S" --times results/repro2010/s1_times.csv \
  --out results/repro2010/gb > logs/k2_2010_group.log 2>&1 || exit 1
grep -E "nodes,|purity>=0.5" logs/k2_2010_group.log
STAGE=grow; env $E python scripts/k_grow.py --nodes results/repro2010/gb_nodes.csv \
  --summary results/repro2010/gb_summary.csv --stas "$S" --iters 0 \
  --out results/repro2010/grow2 > logs/k2_2010_grow.log 2>&1 || exit 1
env $E python scripts/k_score_fam.py --det results/repro2010/grow2/det_it0.parquet \
  --times results/repro2010/s1_times.csv --loud results/repro2010/s1_loud.csv \
  --out results/repro2010/score_it0.csv > logs/k2_2010_score.log 2>&1
echo "== 2010 grown:"; grep "grown groups" logs/k2_2010_score.log
STAGE=ceiling; env $E python scripts/k_ceiling.py --times results/repro2010/s1_times.csv \
  --stas "$S" --iters 0 --out results/repro2010/ceiling2 > logs/k2_2010_ceiling.log 2>&1
env $E python scripts/k_score_fam.py --det results/repro2010/ceiling2/det_it0.parquet \
  --times results/repro2010/s1_times.csv --loud results/repro2010/s1_loud.csv \
  --out results/repro2010/ceiling_score.csv > logs/k2_2010_ceiling_score.log 2>&1
echo "== 2010 ceiling:"; grep -E "grown groups|per-family" logs/k2_2010_ceiling_score.log
STAGE=loss; env $E python scripts/k_loss.py --matched results/repro2010/auto_matched.csv \
  --purity results/repro2010/gb_purity.csv --score results/repro2010/score_it0.csv \
  --times results/repro2010/s1_times.csv --out results/repro2010/loss.csv > logs/k2_2010_loss.log 2>&1
cat logs/k2_2010_loss.log
echo "OK $(date -u +%FT%TZ)" > results/K_REPRO6_DONE
