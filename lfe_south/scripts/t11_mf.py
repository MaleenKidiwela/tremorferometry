#!/usr/bin/env python
"""STAGE 3a -- network matched filter with the stage-2 southern templates.

Per family: per-component normalised cc of each station's own template against that station's own
20 Hz / 1-8 Hz cache, SUMMED over all admitted channels of all admitted stations (the published
plain channel sum, as in REPRO_BOSTOCK §0 and Ducellier's find_all_LFEs), then
net_mad = (cc - median(cc)) / MAD(cc) over the day, detections at net_mad >= --tmin, declustered at
4 s (Bostock 2015 §2).  No family is dropped: a family with no detections keeps a zero row.

  python scripts/t11_mf.py --source Shelly2025 --t0 2018-01-01 --t1 2018-02-01
"""
import os,sys,time,argparse,warnings
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
try:
    import cupy as cp; XP=cp; GPU=True
except Exception: XP=np; GPU=False
def to_np(a): return cp.asnumpy(a) if GPU else a
def ncc(data,tpl):
    M=tpl.size; N=data.size
    L=int(2**np.ceil(np.log2(N+M)))
    dot=XP.fft.irfft(XP.fft.rfft(data,L)*XP.conj(XP.fft.rfft(tpl,L)),L)[:N-M+1]
    c1=XP.concatenate([XP.zeros(1),XP.cumsum(data)]); c2=XP.concatenate([XP.zeros(1),XP.cumsum(data*data)])
    S1=c1[M:N+1]-c1[:N-M+1]; S2=c2[M:N+1]-c2[:N-M+1]
    st=XP.sqrt(M*XP.sum(tpl*tpl)-XP.sum(tpl)**2)
    sd=XP.sqrt(XP.maximum(M*S2-S1*S1,1e-30))
    return (M*dot-XP.sum(tpl)*S1)/(st*sd+1e-30)
def decluster(t,v,gap):
    o=np.argsort(-v); keep=[]
    taken=np.zeros(len(t),bool)
    for i in o:
        if taken[i]: continue
        keep.append(i); taken |= np.abs(t-t[i])<gap
    return np.array(sorted(keep))
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--source",required=True)
    ap.add_argument("--t0",required=True); ap.add_argument("--t1",required=True)
    ap.add_argument("--tmin",type=float,default=8.0); ap.add_argument("--gap",type=float,default=4.0)
    ap.add_argument("--stations",default=""); ap.add_argument("--out",required=True)
    ap.add_argument("--days-file",default="",help="JSON {source|fam: [days]} - score only those days")
    ap.add_argument("--templates",default=f"{RES}/templates_south")
    a=ap.parse_args(); warnings.filterwarnings("ignore")
    M=pd.read_csv(a.templates+".csv"); M=M[M.source==a.source].copy()
    if a.stations: M=M[M.netsta.isin(a.stations.split(","))]
    Z=np.load(a.templates+".npz")
    TPL={r.key:Z[r.key] for r in M.itertuples() if r.key in Z}
    fams=sorted(M.fam.astype(str).unique()); stas=sorted(M.netsta.unique())
    log(f"[mf {a.source}] {len(fams)} families, {len(stas)} stations, {len(TPL)} templates, GPU={GPU}")
    DAYS=None
    if a.days_file:
        import json as _j
        dd=_j.load(open(a.days_file)); DAYS=set()
        for k,v in dd.items():
            if k.split("|")[0]==a.source: DAYS|={pd.Timestamp(x) for x in v}
        log(f"  scoring restricted to {len(DAYS)} held-out days")
    rows=[];cov=[];MISM={}
    for d in pd.date_range(a.t0,a.t1,freq="D",inclusive="left"):
        if DAYS is not None and d not in DAYS: continue
        DAT={}
        for ns in stas:
            net,sta=ns.split(".")
            p=npz_path(net,sta,d)
            if not os.path.exists(p): continue
            try: z=np.load(p); x=z["x"].astype(np.float64)
            except Exception: continue
            if x.shape[1]<0.5*86400*FS or not np.any(x): continue
            DAT[ns]=(XP.asarray(x),float(z["t0"]))
        cov.append(dict(day=str(d.date()),n_sta=len(DAT)))
        if not DAT: continue
        for f in fams:
            sub=M[(M.fam.astype(str)==f)&M.netsta.isin(DAT)]
            cc=None;nch=0;t0=None
            for r in sub.itertuples():
                X,tt=DAT[r.netsta]
                if t0 is None: t0=tt
                T=TPL.get(r.key)
                if T is None: continue
                if T.shape[0]!=X.shape[0]:
                    # a station whose channel count changed since the template was built (two
                    # co-located sensors in the same band, e.g. BK.YBH) contributes nothing that
                    # day -- counted and reported, never silently dropped
                    MISM[(r.netsta,T.shape[0],int(X.shape[0]))]=MISM.get((r.netsta,T.shape[0],int(X.shape[0])),0)+1
                    continue
                for c in range(X.shape[0]):
                    v=ncc(X[c],XP.asarray(T[c].astype(np.float64)))
                    cc=v if cc is None else (cc[:min(cc.size,v.size)]+v[:min(cc.size,v.size)])
                    nch+=1
            if cc is None or nch==0: continue
            c=to_np(cc)
            med=np.median(c); mad=np.median(np.abs(c-med))+1e-30
            nm=(c-med)/mad
            idx=np.where(nm>=a.tmin)[0]
            if idx.size==0: continue
            tt_=t0+idx/FS
            keep=decluster(tt_,nm[idx],a.gap)
            for i in keep:
                rows.append((f,tt_[i],float(nm[idx][i]),float(c[idx][i]/nch),nch,len(sub)))
    for (ns,ct,cd),k in sorted(MISM.items()):
        log(f"  [channels] {ns}: template has {ct} channels, {k} day(s) cached with {cd} -- not stacked")
    D=pd.DataFrame(rows,columns=["fam","t_epoch","net_mad","meancc","nchan","nsta"])
    if len(D): D["t"]=pd.to_datetime(D.t_epoch,unit="s").astype("datetime64[ns]")
    D.to_parquet(a.out)
    pd.DataFrame(cov).to_csv(a.out.replace(".parquet","_coverage.csv"),index=False)
    log(f"[mf {a.source}] {len(D)} detections {a.t0}..{a.t1} -> {a.out}")
if __name__=="__main__": main()
