#!/usr/bin/env python
"""Score blind southern discovery against Ducellier's 66 families (one-to-one +-2 s, Bostock's
>=30 % criterion as a REPORTED line), plus the CEILING: his own templates through the identical
filter (that number is the bed A reproduction of stage R1).  Nothing is filtered."""
import os,sys,argparse
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
from t8_score import one_to_one
ap=argparse.ArgumentParser(); ap.add_argument("--grown",required=True); ap.add_argument("--out",required=True)
ap.add_argument("--t0",default="2008-04-15"); ap.add_argument("--t1",default="2008-05-06")
a=ap.parse_args()
G=pd.read_parquet(a.grown)
T=pd.read_parquet(f"{RES}/ducellier_detections.parquet"); T=T[T.cat=="net"]
T=T[(T.t>=a.t0)&(T.t<a.t1)]
tcol=[c for c in ("t","time","t_epoch") if c in G.columns][0]
gt=(G[tcol].values.astype("datetime64[ns]").astype(np.int64)/1e9 if np.issubdtype(G[tcol].dtype,np.datetime64)
    else G[tcol].values.astype(float))
gcol=[c for c in ("grp","group","gid","fam") if c in G.columns][0]
rows=[]
for g,sub in pd.DataFrame(dict(g=G[gcol].values,t=gt)).groupby("g"):
    best=None
    for f,tf in T.groupby("fam"):
        tt=tf.t.values.astype("datetime64[ns]").astype(np.int64)/1e9
        n,io,it=one_to_one(sub.t.values,tt,2.0)
        frac=n/max(len(sub),1)
        if best is None or frac>best[1]: best=(f,frac,n,len(tt))
    rows.append(dict(grp=g,n_det=len(sub),best_family=best[0],frac_matched=best[1],
                     n_match=best[2],n_family_events=best[3]))
D=pd.DataFrame(rows).sort_values("frac_matched",ascending=False)
D.to_csv(a.out.replace(".csv","_groups.csv"),index=False)
hit=D[D.frac_matched>=0.30]
txt=[f"* grown templates: **{len(D)}**; distinct Ducellier families matched at Bostock's >=30 % "
     f"criterion: **{hit.best_family.nunique()} of {T.fam.nunique()}** (reference line, not a filter).",
     f"* median best-family matched fraction over all grown templates: **{D.frac_matched.median():.3f}**; "
     f"templates with >=1 matched event: **{int((D.n_match>0).sum())}/{len(D)}**.",
     f"* the CEILING on this bed (his own templates through the identical filter) is the bed A "
     f"reproduction in section 1.\n"]
open(a.out,"w").write("\n".join(txt)+"\n")
print("\n".join(txt))
