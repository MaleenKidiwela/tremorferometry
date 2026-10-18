#!/usr/bin/env python
"""One-to-one scoring against a published catalogue (Bostock's criterion, REPRO_BOSTOCK §11).

A detection matches a published event when |dt| <= TOL and the assignment is ONE-TO-ONE
(greedy by |dt| within the family).  Reported against three denominators:
  all       every published event of that family in the window
  indep     published events >= 60 s from the family's previous published event
  loud      the top decile of the family's published events by the catalogue's own cc
Plus PRECISION = fraction of OUR detections that match a published event, and Bostock's
">= 30 % of the grown family's detections" as a REPORTED column, never a gate.
Every family keeps its row, scored or not.
"""
import os,sys,argparse
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
def one_to_one(ours,theirs,tol):
    """greedy one-to-one match by |dt|; returns (n_match, idx_ours_matched, idx_theirs_matched)."""
    if len(ours)==0 or len(theirs)==0: return 0,np.array([],int),np.array([],int)
    o=np.sort(ours); t=np.sort(theirs)
    j=np.searchsorted(t,o)
    cand=[]
    for i,x in enumerate(o):
        for k in (j[i]-1,j[i]):
            if 0<=k<len(t) and abs(t[k]-x)<=tol: cand.append((abs(t[k]-x),i,k))
    cand.sort()
    uo=set();ut=set()
    for dt,i,k in cand:
        if i in uo or k in ut: continue
        uo.add(i);ut.add(k)
    return len(uo),np.array(sorted(uo),int),np.array(sorted(ut),int)
def score(ours,theirs_t,theirs_cc,tol=2.0):
    n,io,it=one_to_one(ours,theirs_t,tol)
    th=np.sort(theirs_t)
    indep=np.ones(len(th),bool)
    if len(th)>1: indep[1:]=np.diff(th)>=60.0
    loud=np.zeros(len(th),bool)
    if theirs_cc is not None and len(th)>=10:
        order=np.argsort(np.argsort(theirs_cc[np.argsort(theirs_t)]))
        loud=order>=int(0.9*len(th))
    m=np.zeros(len(th),bool); m[it]=True
    return dict(n_ours=len(ours),n_theirs=len(th),n_match=n,
        recall_all=n/max(len(th),1),
        n_indep=int(indep.sum()),recall_indep=float(m[indep].sum())/max(int(indep.sum()),1),
        n_loud=int(loud.sum()),recall_loud=(float(m[loud].sum())/max(int(loud.sum()),1) if loud.any() else np.nan),
        precision=n/max(len(ours),1))
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ours",required=True); ap.add_argument("--truth",required=True,
        choices=["ducellier","shelly"])
    ap.add_argument("--t0",required=True); ap.add_argument("--t1",required=True)
    ap.add_argument("--tol",type=float,default=2.0); ap.add_argument("--out",required=True)
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
    T=T[(T.t>=t0)&(T.t<t1)]; O=O[(O.t>=t0)&(O.t<t1)]
    rows=[]
    fams=sorted(set(T.fam)|set(O.fam))
    for f in fams:
        oo=O[O.fam==f].t.values.astype("datetime64[ns]").astype(np.int64)/1e9
        tt=T[T.fam==f]
        r=score(oo,tt.t.values.astype("datetime64[ns]").astype(np.int64)/1e9,
                tt.cc.values if "cc" in tt else None,a.tol)
        r.update(fam=f); rows.append(r)
    D=pd.DataFrame(rows)
    D["bostock30"]=(D.precision>=0.30).astype(int)      # REPORTED reference line, not a gate
    D=D[["fam","n_theirs","n_ours","n_match","recall_all","n_indep","recall_indep","n_loud",
         "recall_loud","precision","bostock30"]].sort_values("recall_all",ascending=False)
    D.to_csv(a.out,index=False)
    tot=dict(n_theirs=int(D.n_theirs.sum()),n_ours=int(D.n_ours.sum()),n_match=int(D.n_match.sum()))
    log(f"[score {a.truth}] families {len(D)} | events theirs {tot['n_theirs']} ours {tot['n_ours']} "
        f"matched {tot['n_match']} | pooled recall {tot['n_match']/max(tot['n_theirs'],1):.3f} "
        f"precision {tot['n_match']/max(tot['n_ours'],1):.3f} | median per-family recall_all "
        f"{D.recall_all.median():.3f}, recall_loud {D.recall_loud.median():.3f} | "
        f"families with >=1 match {int((D.n_match>0).sum())}/{len(D)} (reference line: "
        f"precision>=0.30 -> {int(D.bostock30.sum())})")
if __name__=="__main__": main()
