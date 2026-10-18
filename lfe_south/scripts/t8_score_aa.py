#!/usr/bin/env python
"""One-to-one scoring (a copy of t8_score.py) with the DENOMINATOR FAULT FIXED.

t8_score.py takes the published catalogue over the whole --t0..--t1 span while our detections exist
only on the held-out days the matched filter was actually run on.  For Shelly that means 61,441
published events in the denominator of which only 5,021 (per-family split) / a subset (source union)
could ever be matched: the recall CAP is 0.082, so the reported "pooled recall 0.026" was 0.32 of
the achievable maximum, not 0.026 of it.  Here the truth is restricted to exactly the days scored.

  --days-mode source : the union of the source's held-out days (every family was run on all of
                       them, so this is the honest denominator)  [default]
  --days-mode fam    : each family's own held-out days only
  --days-mode all    : t8_score.py behaviour, kept so the old number is reproducible

Also sweeps the detection statistic: with --sweep the score is reported at every net_mad threshold
in the detections file, so the threshold is a ranking column and not a gate.
"""
import os,sys,argparse,json
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
from t8_score import one_to_one,score

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ours",required=True); ap.add_argument("--truth",default="shelly",choices=["ducellier","shelly"])
    ap.add_argument("--t0",required=True); ap.add_argument("--t1",required=True)
    ap.add_argument("--tol",type=float,default=2.0); ap.add_argument("--out",required=True)
    ap.add_argument("--days-file",default=""); ap.add_argument("--days-mode",default="source",choices=["source","fam","all"])
    ap.add_argument("--source",default="Shelly2025")
    ap.add_argument("--sweep",default="",help="comma list of net_mad thresholds to score at")
    a=ap.parse_args()
    O=pd.read_parquet(a.ours); O["fam"]=O.fam.astype(str)
    if a.truth=="ducellier":
        T=pd.read_parquet(f"{RES}/ducellier_detections.parquet"); T=T[T.cat=="net"].copy()
    else:
        S=pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv")
        S.columns=[c.strip('%') for c in S.columns]
        S["t"]=(pd.to_datetime(dict(year=S.year,month=S.month,day=S.day,hour=S.hour,minute=S.minute,
                second=S.second.astype(int))).astype("datetime64[ns]")+pd.to_timedelta(S.second%1,unit="s"))
        T=S.rename(columns={"idnum":"fam","meancc":"cc"})[["fam","t","cc"]]
    T["fam"]=T.fam.astype(str)
    t0,t1=pd.Timestamp(a.t0),pd.Timestamp(a.t1)
    T=T[(T.t>=t0)&(T.t<t1)].copy(); O=O[(O.t>=t0)&(O.t<t1)].copy()
    DF={}
    if a.days_file and a.days_mode!="all":
        dd=json.load(open(a.days_file))
        for k,v in dd.items():
            s_,f_=k.split("|")
            if s_!=a.source: continue
            DF[f_]={pd.Timestamp(x) for x in v}
        if a.days_mode=="source":
            U=set().union(*DF.values()) if DF else set()
            DF={f_:U for f_ in DF}
        T["day"]=T.t.dt.floor("D")
        T=T[[ (f in DF and d in DF[f]) for f,d in zip(T.fam,T.day)]].copy()
    thr=[None]+[float(x) for x in a.sweep.split(",") if x.strip()] if a.sweep else [None]
    OUT=[]
    for tm in thr:
        Oc=O if tm is None else O[O.net_mad>=tm]
        rows=[]
        for f in sorted(set(T.fam)|set(Oc.fam)):
            oo=Oc[Oc.fam==f].t.values.astype("datetime64[ns]").astype(np.int64)/1e9
            tt=T[T.fam==f]
            r=score(oo,tt.t.values.astype("datetime64[ns]").astype(np.int64)/1e9,
                    tt.cc.values if "cc" in tt else None,a.tol)
            r.update(fam=f,tmin=(np.nan if tm is None else tm)); rows.append(r)
        D=pd.DataFrame(rows); D["bostock30"]=(D.precision>=0.30).astype(int)
        OUT.append(D)
        tot=dict(th=int(D.n_theirs.sum()),ou=int(D.n_ours.sum()),ma=int(D.n_match.sum()))
        log(f"[score-aa {a.truth} days={a.days_mode} tmin={tm}] fam {len(D)} | theirs {tot['th']} "
            f"ours {tot['ou']} matched {tot['ma']} | pooled recall {tot['ma']/max(tot['th'],1):.3f} "
            f"precision {tot['ma']/max(tot['ou'],1):.4f} | median recall {D.recall_all.median():.3f} "
            f"| fam with >=1 match {int((D.n_match>0).sum())}/{len(D)}")
    A=pd.concat(OUT,ignore_index=True)
    A=A[["tmin","fam","n_theirs","n_ours","n_match","recall_all","n_indep","recall_indep","n_loud",
         "recall_loud","precision","bostock30"]].sort_values(["tmin","recall_all"],ascending=[True,False])
    A.to_csv(a.out,index=False)
if __name__=="__main__": main()
