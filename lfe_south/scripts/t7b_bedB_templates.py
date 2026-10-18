#!/usr/bin/env python
"""BED B -- 'his set plus every station installed since'.  Ducellier released templates only for
HIS stations, so a new station can only join the network if it is given a template.  Each new
station's template is built from ITS OWN stack at HIS published detection times in the window,
cut to the same length and the same time origin as that family's own templates (so the detection
time convention is unchanged).  A waveform is never transferred between stations.

  results/duc_bedB_templates.npz : "<family>|<STA>|<CHA>" -> (L,) float32
"""
import os,sys,pickle,argparse,warnings
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
DUC=f"{ROOT}/catalogs/ducellier_catalog"; CACHE_DUC=f"{WS}/cache_duc"
ap=argparse.ArgumentParser(); ap.add_argument("--t0",default="2008-04-01"); ap.add_argument("--t1",default="2008-05-10")
ap.add_argument("--min-events",type=int,default=30); ap.add_argument("--out",default=f"{RES}/duc_bedB_templates.npz")
a=ap.parse_args(); warnings.filterwarnings("ignore")
B=pd.read_csv(f"{RES}/duc_bedB_stations.csv"); NEW=B[B.in_bedA==0]
T=pd.read_parquet(f"{RES}/ducellier_detections.parquet"); T=T[(T.cat=="net")&(T.t>=a.t0)&(T.t<a.t1)]
T["day"]=T.t.dt.floor("D")
LEN={}
for f in sorted(os.listdir(f"{DUC}/src/templates_v1")):
    fs=[x for x in os.listdir(f"{DUC}/src/templates_v1/{f}") if x.endswith(".pkl")]
    if fs: LEN[f]=pickle.load(open(f"{DUC}/src/templates_v1/{f}/{fs[0]}","rb"))[0].stats.npts
OUT={};rows=[]
for r in NEW.itertuples():
    for day,Td in T.groupby("day"):
        pass
    p_any=f"{CACHE_DUC}/{r.net}.{r.sta}"
    if not os.path.isdir(p_any): continue
    for fam,Tf in T.groupby("fam"):
        L=LEN.get(fam)
        if L is None or len(Tf)<a.min_events: continue
        acc=None;n=0;chans=None
        for day,Td in Tf.groupby("day"):
            d=pd.Timestamp(day); p=f"{CACHE_DUC}/{r.net}.{r.sta}/{d.year}/{d.dayofyear:03d}.npz"
            if not os.path.exists(p): continue
            z=np.load(p); X=z["x"]; V=z["valid"]; t0=float(z["t0"]); chans=list(map(str,z["chans"]))
            idx=np.round(((Td.t.values.astype("datetime64[ns]").astype(np.int64)/1e9)-t0)*FS).astype(int)
            idx=idx[(idx>=0)&(idx+L<=X.shape[1])]
            if len(idx)==0: continue
            w=X[:,idx[:,None]+np.arange(L)]
            ok=V[:,idx[:,None]+np.arange(L)].all(axis=(0,2))
            w=w[:,ok]; 
            if w.shape[1]==0: continue
            rr=np.sqrt((w.astype(np.float64)**2).mean(axis=(0,2)))+1e-20
            s=(w/rr[None,:,None]).sum(1); acc=s if acc is None else acc+s; n+=w.shape[1]
        if acc is None or n<a.min_events: continue
        for i,c in enumerate(chans):
            OUT[f"{fam}|{r.sta}|{c}"]=(acc[i]/n).astype(np.float32)
        rows.append(dict(fam=fam,net=r.net,sta=r.sta,n_stack=n,nchan=len(chans),L=L,
                         min_hypo_km=r.min_hypo_km))
np.savez_compressed(a.out,**OUT)
pd.DataFrame(rows).to_csv(a.out.replace(".npz",".csv"),index=False)
log(f"[bedB templates] {len(OUT)} channel-templates for {len(set(k.split('|')[1] for k in OUT))} new "
    f"stations over {len(set(k.split('|')[0] for k in OUT))} families")
