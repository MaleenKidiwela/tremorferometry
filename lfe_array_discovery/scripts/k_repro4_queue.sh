#!/usr/bin/env bash
# K_REPRO4: (1) the SAME literal recipe on the PERMANENT-station beds (2012, 2010);
#           (2) the ONE blind check of the 30 sealed families on the unchanged 2005 15-station run;
#           (3) locations and the new-family classification.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/repro
exec 9>locks/k_repro4.lock
flock -n 9 || { echo "already running"; exit 0; }
rm -f results/K_REPRO4_DONE results/K_REPRO4_FAILED
STAGE=start
trap 'echo "[q4] TRAP exit $?"; [ -f results/K_REPRO4_DONE ] || echo "stage=$STAGE" > results/K_REPRO4_FAILED' EXIT
echo "[q4] start $(date -u +%FT%TZ) pid=$$"
S15="GLBC,TWKB,MGCB,TSJB,TWBB,PGC,SILB,SSIB,LCBC,B003,B005,JRBC,LZB,VGZ,SNB"

run_bed () {           # $1 year  $2 day0  $3 day1  $4 stations  $5 outdir
  local Y=$1 D0=$2 D1=$3 ST=$4 OD=$5
  export K_YEAR=$Y K_DAY0=$D0 K_DAY1=$D1 K_OUT=$OD
  mkdir -p "$OD"
  local DATES; DATES=$(python - "$D0" "$D1" <<'PY'
import sys, pandas as pd
print(",".join(d.strftime("%Y-%m-%d") for d in pd.date_range(sys.argv[1], sys.argv[2])))
PY
)
  [ -s "$OD/s1_times.csv" ] || python scripts/k_s1_loud.py --nfam 100 --stas "$ST" --out "$OD/s1" \
      > logs/k_${Y}_s1.log 2>&1 || return 1
  python scripts/k_auto.py --dates "$DATES" --hours all --stas "$ST" --out "$OD/auto" \
      > logs/k_${Y}_auto.log 2>&1 || return 1
  python - "$OD" <<'PY' || return 1
import sys, glob, pandas as pd, numpy as np
OD=sys.argv[1]
A=[pd.read_parquet(f) for f in sorted(glob.glob(f"{OD}/auto/*.parquet"))]
R=pd.concat([x for x in A if len(x)],ignore_index=True)
S=pd.read_csv(f"{OD}/auto/_stats.csv")
bad={(r.date,int(r.hour)) for r in S.itertuples() if not (r.nmad_in_sd>3)}
R=R[~pd.Series(list(zip(R.date,R.hour)),index=R.index).isin(bad)]
R=R[np.isfinite(R.nmad)].reset_index(drop=True)
R.to_parquet(f"{OD}/auto_all.parquet",index=False)
print("detections",len(R),"bad hours",len(bad),"hours",R.groupby(['date','hour']).ngroups)
PY
  python scripts/k_score_pairs.py --det "$OD/auto" --times "$OD/s1_times.csv" --stas "$ST" \
      > logs/k_${Y}_pairs.log 2>&1 || true
  python - "$OD" "$ST" <<'PY' > logs/k_${Y}_pairs2.log 2>&1
import sys; sys.path.insert(0,'scripts')
import numpy as np, pandas as pd, k_common as K, k_score_pairs as KS
OD,ST=sys.argv[1],sys.argv[2].split(",")
R=pd.read_parquet(f"{OD}/auto_all.parquet")
C=pd.read_csv(f"{OD}/s1_times.csv"); C['t_corr']=C.t_origin+C.dt
TT=K.fam_tt(ST); TTM={int(b):float(np.median([TT[(int(b),s)] for s in ST])) for b in C.bkey.unique()}
M=KS.score(R,C,TTM,0.6,8.0); N=KS.score(R,C,TTM,0.6,8.0,shift=300.)
M.to_csv(f"{OD}/auto_matched.csv",index=False)
print(f"{len(R)} det; matched {len(M)} ({len(M)/max(len(R),1):.2%}) fams {M.bkey.nunique()}; chance {len(N)/max(len(R),1):.2%}")
PY
  cat logs/k_${Y}_pairs2.log
  python scripts/k_group.py --detfile "$OD/auto_all.parquet" --budget 1000000 --min-members 3 \
      --stas "$ST" --times "$OD/s1_times.csv" --out "$OD/gb" > logs/k_${Y}_group.log 2>&1 || return 1
  grep -E "nodes,|purity>=0.5" logs/k_${Y}_group.log
  python scripts/k_grow.py --nodes "$OD/gb_nodes.csv" --summary "$OD/gb_summary.csv" \
      --stas "$ST" --iters 0 --out "$OD/grow" > logs/k_${Y}_grow.log 2>&1 || return 1
  python scripts/k_score_fam.py --det "$OD/grow/det_it0.parquet" --times "$OD/s1_times.csv" \
      --loud "$OD/s1_loud.csv" --out "$OD/score_it0.csv" > logs/k_${Y}_score.log 2>&1 || true
  echo "== BED $Y grown:"; grep "grown groups" logs/k_${Y}_score.log
  python scripts/k_ceiling.py --times "$OD/s1_times.csv" --stas "$ST" --iters 0 \
      --out "$OD/ceiling" > logs/k_${Y}_ceiling.log 2>&1 || true
  python scripts/k_score_fam.py --det "$OD/ceiling/det_it0.parquet" --times "$OD/s1_times.csv" \
      --loud "$OD/s1_loud.csv" --out "$OD/ceiling_score.csv" > logs/k_${Y}_ceiling_score.log 2>&1 || true
  echo "== BED $Y ceiling:"; grep -E "grown groups|per-family medians" logs/k_${Y}_ceiling_score.log
  python scripts/k_loss.py --matched "$OD/auto_matched.csv" --purity "$OD/gb_purity.csv" \
      --score "$OD/score_it0.csv" --times "$OD/s1_times.csv" --out "$OD/loss.csv" \
      > logs/k_${Y}_loss.log 2>&1 || true
  cat logs/k_${Y}_loss.log
  unset K_YEAR K_DAY0 K_DAY1 K_OUT
}

STAGE=bed2012_permanent
run_bed 2012 2012-08-15 2012-10-05 \
  "PGC,VGZ,SNB,NLLB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926" \
  results/repro2012 || echo "[q4] 2012 bed FAILED (continuing)"

STAGE=bed2010_permanent
D2010=$(python - <<'PY'
import glob, os, pandas as pd
ds=set()
for p in glob.glob('/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz/*/2010/*.npz'):
    ds.add(int(os.path.basename(p)[:3]))
ds=sorted(ds)
a=pd.Timestamp('2010-01-01')+pd.Timedelta(days=ds[0]-1)
b=pd.Timestamp('2010-01-01')+pd.Timedelta(days=ds[-1]-1)
print(a.strftime('%Y-%m-%d'), b.strftime('%Y-%m-%d'))
PY
)
run_bed 2010 $(echo $D2010) "PGC,VGZ,SNB,NLLB,LRIV,LZB,PFB,YOUB,B011,B926,OCP,STW" \
  results/repro2010 || echo "[q4] 2010 bed FAILED (continuing)"

STAGE=blind_hidden
python scripts/k_s1_loud.py --nfam 30 --hidden --stas "$S15" --out results/repro/s1hidden \
  > logs/k_s1hidden.log 2>&1 || exit 1
tail -12 logs/k_s1hidden.log
python scripts/k_score_fam.py --det results/repro/grow15/det_it0.parquet \
  --times results/repro/s1hidden_times.csv --loud results/repro/s1hidden_loud.csv \
  --out results/repro/score_hidden_it0.csv > logs/k_blind_hidden.log 2>&1 || exit 1
echo "===== BLIND CHECK (30 sealed families, unchanged 2005 15-station run) ====="
grep -E "grown groups|per-family medians" logs/k_blind_hidden.log
python - <<'PY' | tee logs/k_blind_summary.log
import pandas as pd, sys
sys.path.insert(0,'/home/jovyan/tremorferometry/lfe_picker/scripts')
ISO=[8,10,12,15,29,62,87,113,115,119]
S=pd.read_csv('results/repro/score_hidden_it0.csv')
ok=S[S.frac>=0.30]
print("HIDDEN rediscovered:", ok.bkey.nunique(), "of 30 ; templates passing:", len(ok))
print("ISOLATED stratum:", len(set(ok.bkey)&set(ISO)), "of 10 ->", sorted(set(ok.bkey)&set(ISO)))
per=ok.sort_values('n_match',ascending=False).drop_duplicates('bkey')
if len(per):
    print("per-family medians: recall_all %.3f recall_60 %.3f recall_loud %.3f purity %.3f"%(
      per.recall_all.median(),per.recall_60.median(),per.recall_loud.median(),per.frac.median()))
PY
echo "BLIND DONE $(date -u +%FT%TZ)" > results/K_BLIND_DONE

STAGE=hidden_ceiling
python scripts/k_ceiling.py --times results/repro/s1hidden_times.csv --stas "$S15" --iters 0 \
  --out results/repro/ceiling_hidden > logs/k_ceiling_hidden.log 2>&1 || true

STAGE=locate
[ -s results/repro/loc15_it0.csv ] || python scripts/k_locate.py --tpl results/repro/grow15/tpl_it0.npz \
  --det results/repro/grow15/det_it0.parquet --score results/repro/grow15/score_pre_it0.csv \
  --out results/repro/loc15_it0.csv > logs/k_locate15.log 2>&1
python scripts/k_classify.py --loc results/repro/loc15_it0.csv \
  --score-kept results/repro/grow15/score_pre_it0.csv \
  --score-hidden results/repro/score_hidden_it0.csv \
  --tpl results/repro/grow15/tpl_it0.npz \
  --tpl-kept-ceiling results/repro/ceiling15/tpl_it0.npz \
  --tpl-hidden-ceiling results/repro/ceiling_hidden/tpl_it0.npz \
  --out results/repro/classify15.csv > logs/k_classify.log 2>&1 || true
tail -45 logs/k_classify.log

echo "OK $(date -u +%FT%TZ)" > results/K_REPRO4_DONE
echo "[q4] DONE $(date -u +%FT%TZ)"
