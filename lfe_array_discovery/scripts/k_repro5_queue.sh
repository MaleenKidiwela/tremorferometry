#!/usr/bin/env bash
# K_REPRO5: finish the three tasks properly.
#  - hidden corrected times + hidden CEILING (2005 bed, env explicitly pinned in every call)
#  - kept CEILING re-run so the template waveforms are saved (needed by the classifier)
#  - classification of every 2005 template against the CORRECT hidden score
#  - 2010 permanent bed with the stations that actually have coordinates
# Every stage pins its own environment; nothing is exported across stages (that leak voided the
# first blind run).
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/repro
exec 9>locks/k_repro5.lock
flock -n 9 || { echo "already running"; exit 0; }
rm -f results/K_REPRO5_DONE results/K_REPRO5_FAILED
STAGE=start
trap 'echo "[q5] TRAP exit $?"; [ -f results/K_REPRO5_DONE ] || echo "stage=$STAGE" > results/K_REPRO5_FAILED' EXIT
S15="GLBC,TWKB,MGCB,TSJB,TWBB,PGC,SILB,SSIB,LCBC,B003,B005,JRBC,LZB,VGZ,SNB"
E2005="K_YEAR=2005 K_DAY0=2005-08-29 K_DAY1=2005-09-21"
echo "[q5] start $(date -u +%FT%TZ) pid=$$"

STAGE=hidden_corrected_times
env $E2005 python scripts/k_s1_loud.py --nfam 30 --hidden --stas "$S15" \
  --out results/repro/s1hid2005 > logs/k_s1hid2005.log 2>&1 || exit 1
tail -12 logs/k_s1hid2005.log
env $E2005 python scripts/k_score_fam.py --det results/repro/grow15/det_it0.parquet \
  --times results/repro/s1hid2005_times.csv --loud results/repro/s1hid2005_loud.csv \
  --out results/repro/score_hidden_corr.csv > logs/k_blind_corr.log 2>&1 || true
echo "== blind, CORRECTED times (confirmatory; the reported number stays the raw-time one):"
grep -E "grown groups|per-family" logs/k_blind_corr.log

STAGE=kept_ceiling_with_templates
env $E2005 python scripts/k_ceiling.py --times results/repro/s1all_times.csv --stas "$S15" \
  --iters 0 --out results/repro/ceiling15b > logs/k_ceiling15b.log 2>&1 || exit 1

STAGE=hidden_ceiling
env $E2005 python scripts/k_ceiling.py --times results/repro/s1hid2005_times.csv --stas "$S15" \
  --iters 0 --out results/repro/ceiling_hid2005 > logs/k_ceiling_hid2005.log 2>&1 || true

STAGE=classify
env $E2005 python scripts/k_classify.py --loc results/repro/loc15_it0.csv \
  --score-kept results/repro/grow15/score_pre_it0.csv \
  --score-hidden results/repro/score_hidden2005_raw.csv \
  --tpl results/repro/grow15/tpl_it0.npz \
  --tpl-kept-ceiling results/repro/ceiling15b/tpl_it0.npz \
  --tpl-hidden-ceiling results/repro/ceiling_hid2005/tpl_it0.npz \
  --out results/repro/classify15.csv > logs/k_classify2.log 2>&1 || exit 1
tail -45 logs/k_classify2.log

STAGE=bed2010
D0=$(python - <<'PY'
import glob, os, pandas as pd
ds=sorted({int(os.path.basename(p)[:3]) for p in glob.glob('/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz/*/2010/*.npz')})
print((pd.Timestamp('2010-01-01')+pd.Timedelta(days=ds[0]-1)).strftime('%Y-%m-%d'))
PY
)
D1=$(python - <<'PY'
import glob, os, pandas as pd
ds=sorted({int(os.path.basename(p)[:3]) for p in glob.glob('/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz/*/2010/*.npz')})
print((pd.Timestamp('2010-01-01')+pd.Timedelta(days=ds[-1]-1)).strftime('%Y-%m-%d'))
PY
)
S2010=$(python - <<'PY'
import sys; sys.path.insert(0,'scripts')
import k_common as K
want="PGC,VGZ,SNB,NLLB,LRIV,LZB,B011,B926,PFB,YOUB,OCP,STW".split(",")
ST=K.station_table()
ok=[s for s in want if s in ST.index and float(ST.lat.get(s, float('nan')))==float(ST.lat.get(s, float('nan')))]
print(",".join(ok))
PY
)
echo "[q5] 2010 bed $D0..$D1 stations: $S2010"
E2010="K_YEAR=2010 K_DAY0=$D0 K_DAY1=$D1 K_OUT=results/repro2010"
mkdir -p results/repro2010
DATES=$(python - "$D0" "$D1" <<'PY'
import sys, pandas as pd
print(",".join(d.strftime("%Y-%m-%d") for d in pd.date_range(sys.argv[1], sys.argv[2])))
PY
)
env $E2010 python scripts/k_s1_loud.py --nfam 100 --stas "$S2010" --out results/repro2010/s1 \
  > logs/k_2010_s1.log 2>&1 || { echo "[q5] 2010 s1 failed"; tail -5 logs/k_2010_s1.log; }
if [ -s results/repro2010/s1_times.csv ]; then
  env $E2010 python scripts/k_auto.py --dates "$DATES" --hours all --stas "$S2010" \
    --out results/repro2010/auto > logs/k_2010_auto.log 2>&1
  env $E2010 python - <<'PY' > logs/k_2010_merge.log 2>&1
import glob, pandas as pd, numpy as np
A=[pd.read_parquet(f) for f in sorted(glob.glob("results/repro2010/auto/*.parquet"))]
R=pd.concat([x for x in A if len(x)],ignore_index=True)
S=pd.read_csv("results/repro2010/auto/_stats.csv")
bad={(r.date,int(r.hour)) for r in S.itertuples() if not (r.nmad_in_sd>3)}
R=R[~pd.Series(list(zip(R.date,R.hour)),index=R.index).isin(bad)]
R=R[np.isfinite(R.nmad)].reset_index(drop=True)
R.to_parquet("results/repro2010/auto_all.parquet",index=False); print("det",len(R))
PY
  env $E2010 python - "$S2010" <<'PY' > logs/k_2010_pairs.log 2>&1
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
  cat logs/k_2010_pairs.log
  env $E2010 python scripts/k_group.py --detfile results/repro2010/auto_all.parquet \
    --budget 1000000 --min-members 3 --stas "$S2010" --times results/repro2010/s1_times.csv \
    --out results/repro2010/gb > logs/k_2010_group.log 2>&1
  grep -E "nodes,|purity>=0.5" logs/k_2010_group.log
  env $E2010 python scripts/k_grow.py --nodes results/repro2010/gb_nodes.csv \
    --summary results/repro2010/gb_summary.csv --stas "$S2010" --iters 0 \
    --out results/repro2010/grow > logs/k_2010_grow.log 2>&1
  env $E2010 python scripts/k_score_fam.py --det results/repro2010/grow/det_it0.parquet \
    --times results/repro2010/s1_times.csv --loud results/repro2010/s1_loud.csv \
    --out results/repro2010/score_it0.csv > logs/k_2010_score.log 2>&1
  echo "== BED 2010 grown:"; grep "grown groups" logs/k_2010_score.log
  env $E2010 python scripts/k_ceiling.py --times results/repro2010/s1_times.csv --stas "$S2010" \
    --iters 0 --out results/repro2010/ceiling > logs/k_2010_ceiling.log 2>&1
  env $E2010 python scripts/k_score_fam.py --det results/repro2010/ceiling/det_it0.parquet \
    --times results/repro2010/s1_times.csv --loud results/repro2010/s1_loud.csv \
    --out results/repro2010/ceiling_score.csv > logs/k_2010_ceiling_score.log 2>&1
  echo "== BED 2010 ceiling:"; grep -E "grown groups|per-family" logs/k_2010_ceiling_score.log
  env $E2010 python scripts/k_loss.py --matched results/repro2010/auto_matched.csv \
    --purity results/repro2010/gb_purity.csv --score results/repro2010/score_it0.csv \
    --times results/repro2010/s1_times.csv --out results/repro2010/loss.csv \
    > logs/k_2010_loss.log 2>&1
  cat logs/k_2010_loss.log
fi

echo "OK $(date -u +%FT%TZ)" > results/K_REPRO5_DONE
echo "[q5] DONE $(date -u +%FT%TZ)"
