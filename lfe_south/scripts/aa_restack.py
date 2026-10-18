#!/usr/bin/env python
"""Shelly's template ITERATION (agent AA): restack each station's template at OUR OWN strongest
detection times, the step Shelly et al. 2025 sec.2 describes as "stacking the detections with the
highest correlation sum across the network ... The stacked waveforms are then used to construct new
templates with higher signal-to-noise ratios, with the process repeating (2-3 iterations)".

Input  : a detections parquet produced by t11_mf_aa.py on the BUILD days only (so the held-out test
         days stay unseen), plus the template table that produced it.
Output : a new templates npz/csv in the same format as t9_templates.py, same keys, same tt_s/best_tau
         columns (the per-station offsets are unchanged -- only the waveform is restacked).

  python scripts/aa_restack.py --det results/aa_shelly/mfA_build.parquet \
      --templates results/templates_south_build --nstack 200 \
      --out results/aa_shelly/templatesD_build.npz
"""
import os,sys,argparse,warnings
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--det",required=True); ap.add_argument("--templates",required=True)
    ap.add_argument("--out",required=True); ap.add_argument("--nstack",type=int,default=200)
    ap.add_argument("--source",default="Shelly2025")
    a=ap.parse_args(); warnings.filterwarnings("ignore")
    M=pd.read_csv(a.templates+".csv"); M=M[M.source==a.source].copy(); M["fam"]=M.fam.astype(str)
    Z=np.load(a.templates+".npz")
    D=pd.read_parquet(a.det); D["fam"]=D.fam.astype(str)
    D=D.sort_values("net_mad",ascending=False).groupby("fam").head(a.nstack).copy()
    D["day"]=D.t.dt.floor("D")
    log(f"[restack] {len(D)} detections over {D.fam.nunique()} families, <= {a.nstack}/family")
    PRE=float(M.pre_s.iloc[0]); LEN=float(M.len_s.iloc[0]); L=int(round(LEN*FS))
    TPL={}; rows=[]
    for (net,sta),Ms in M.groupby(["net","sta"]):
        acc={}; n={}
        fams=set(Ms.fam)
        sub=D[D.fam.isin(fams)]
        off={r.fam:(float(r.tt_s)+float(r.best_tau)) for r in Ms.itertuples()}
        C={r.fam:int(r.ncomp) for r in Ms.itertuples()}
        for day,Td in sub.groupby("day"):
            p=npz_path(net,sta,day)
            if not os.path.exists(p): continue
            try: z=np.load(p); X=z["x"].astype(np.float32); t0=float(z["t0"])
            except Exception: continue
            for f,Tf in Td.groupby("fam"):
                if f not in off or X.shape[0]!=C.get(f): continue
                st0=(Tf.t.values.astype("datetime64[ns]").astype(np.int64)/1e9)+off[f]-PRE-t0
                idx=np.round(st0*FS).astype(int); idx=idx[(idx>=0)&(idx+L<=X.shape[1])]
                if len(idx)==0: continue
                w=X[:,idx[:,None]+np.arange(L)]
                r=np.sqrt((w.astype(np.float64)**2).mean(axis=(0,2)))+1e-20
                s=(w/r[None,:,None]).sum(1)
                if f in acc and acc[f].shape!=s.shape: continue
                acc[f]=s if f not in acc else acc[f]+s; n[f]=n.get(f,0)+len(idx)
        for f,s in acc.items():
            e=Ms[Ms.fam==f].iloc[0]
            TPL[e.key]=(s/max(n[f],1)).astype(np.float32)
            d=e._asdict() if hasattr(e,"_asdict") else dict(e)
            d["n_stack"]=n[f]; rows.append(d)
        log(f"  {net}.{sta}: {len(acc)} restacked templates")
    if not rows:
        log("[restack] nothing restacked"); sys.exit(3)
    T=pd.DataFrame(rows)
    keep=[c for c in ["source","fam","net","sta","netsta","key","ncomp","ncomp_station","n_days",
        "n_days_offcount","n_stack","pre_s","len_s","best_tau","excess_db","null_pctile","z_vs_null",
        "hypo_km","tt_s"] if c in T.columns]
    T[keep].to_csv(a.out.replace(".npz",".csv"),index=False)
    np.savez_compressed(a.out,**TPL)
    log(f"[restack] {len(TPL)} templates -> {a.out}")
if __name__=="__main__": main()
