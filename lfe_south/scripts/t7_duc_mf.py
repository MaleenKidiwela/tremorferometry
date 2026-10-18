#!/usr/bin/env python
"""STAGE R1c -- run DUCELLIER'S OWN TEMPLATES on our download of HIS OWN permanent stations with
HIS OWN parameters, and score one-to-one against his published detection times.

Published recipe (catalogs/ducellier_catalog/src/find_all_LFEs_parallel.py __main__, verbatim):
  templates  templates_v1/<family>/<STA>_<CHA>.pkl      (20 Hz, response-removed, 1.5-9 Hz)
  data       his stations_permanent.txt set, his channels, 20 Hz, same chain (t6_duc_cache.py)
  cc         per-CHANNEL normalised cross-correlation (correlate.optimized), then the MEAN over
             all channels that have a COMPLETE gapless record of the block  (his `nchannel` rule)
  threshold  8 x MAD, MAD = median(|meancc - mean(meancc)|) of the block
  decluster  clean_LFEs: runs of samples within floor(1/(dt*freq0)) = 20 samples (1 s) collapse to
             their maximum-cc sample
  blocks     12 h (the block length of his production __main__)
  time       block start + lag*dt  = the sample at which the template's first sample aligns --
             the same convention his catalogue is written in, so the comparison is like-for-like
Nothing is filtered: every family keeps its row whatever it scores.

  python scripts/t7_duc_mf.py --bed A --t0 2008-04-01 --t1 2008-05-10
"""
import os,sys,time,argparse,pickle,warnings,json
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
DUC=f"{ROOT}/catalogs/ducellier_catalog"
CACHE_DUC=f"{WS}/cache_duc"
BLOCK_H=12; THR=8.0; FREQ0=1.0
try:
    import cupy as cp; XP=cp; GPU=True
except Exception:
    XP=np; GPU=False
def to_np(a): return cp.asnumpy(a) if GPU else a
def ncc(data,tpl):
    """normalised cc of `tpl` (M,) over `data` (N,) -> (N-M+1,), Ducellier correlate.optimized."""
    M=tpl.size; N=data.size
    L=int(2**np.ceil(np.log2(N+M)))
    F=XP.fft.rfft(data,L)*XP.conj(XP.fft.rfft(tpl,L))
    dot=XP.fft.irfft(F,L)[:N-M+1]
    c1=XP.concatenate([XP.zeros(1),XP.cumsum(data)])
    c2=XP.concatenate([XP.zeros(1),XP.cumsum(data*data)])
    S1=c1[M:N+1]-c1[:N-M+1]; S2=c2[M:N+1]-c2[:N-M+1]
    st=XP.sqrt(M*XP.sum(tpl*tpl)-XP.sum(tpl)**2)
    sd=XP.sqrt(XP.maximum(M*S2-S1*S1,1e-30))
    return (M*dot-XP.sum(tpl)*S1)/(st*sd+1e-30)
def decluster(idx,cc,maxdiff):
    """clean_LFEs: chain consecutive indices within maxdiff, keep the max-cc sample of each run."""
    out=[];run=[idx[0]]
    for i in idx[1:]:
        if i-run[-1]<=maxdiff: run.append(i)
        else: out.append(run); run=[i]
    out.append(run)
    return np.array([r[int(np.argmax(cc[r]))] for r in out])
def load_templates(fam_stations):
    T={}
    for fam,stas in fam_stations.items():
        d=f"{DUC}/src/templates_v1/{fam}"
        if not os.path.isdir(d): continue
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".pkl"): continue
            sta,cha=fn[:-4].split("_",1)
            if sta not in stas: continue
            tr,ori=pickle.load(open(os.path.join(d,fn),"rb"))
            T.setdefault(fam,{})[(sta,cha)]=np.asarray(tr.data,dtype=np.float64)
    return T
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--bed",default="A"); ap.add_argument("--t0",default="2008-04-01")
    ap.add_argument("--t1",default="2008-05-10"); ap.add_argument("--out",default=None)
    ap.add_argument("--extra-templates",default="",help="bed B: npz of '<fam>|<STA>|<CHA>' templates")
    a=ap.parse_args(); warnings.filterwarnings("ignore")
    out=a.out or f"{RES}/duc_repro_bed{a.bed}_detections.parquet"
    fp=pd.read_csv(f"{DUC}/data/Ducellier/families_permanent.txt",sep=r"\s+",header=None,
                   names=["family","stations","duration"])
    fam_stations={r.family:set(r.stations.split(",")) for r in fp.itertuples()}
    if a.bed=="B":                                   # bed B = his set + every new station near a family
        B=pd.read_csv(f"{RES}/duc_bedB_stations.csv")
        for f in fam_stations: fam_stations[f]|=set(B[B.in_bedA==0].sta)
    T=load_templates(fam_stations)
    if a.extra_templates and os.path.exists(a.extra_templates):
        Z=np.load(a.extra_templates); nx=0
        for k in Z.files:
            fam,sta,cha=k.split("|")
            if fam in T: T[fam][(sta,cha)]=Z[k].astype(np.float64); nx+=1
        log(f"  bed B: +{nx} channel-templates from stations not in his set")
    log(f"[duc-mf bed {a.bed}] {len(T)} families with templates, GPU={GPU}")
    A=pd.read_csv(f"{RES}/duc_bedA_stations.csv")
    NETOF=dict(zip(A.sta,A.net))
    Bs=pd.read_csv(f"{RES}/duc_bedB_stations.csv"); NETOF.update(dict(zip(Bs.sta,Bs.net)))
    need={}
    for fam,dd in T.items():
        for (sta,cha) in dd: need.setdefault(sta,set()).add(cha)
    days=pd.date_range(a.t0,a.t1,freq="D")
    maxdiff=int(np.floor(1.0/((1/FS)*FREQ0)))
    rows=[];nblk=0;t00=time.time()
    for d in days:
        DAT={}
        for sta,chs in need.items():
            net=NETOF.get(sta)
            if net is None: continue
            p=f"{CACHE_DUC}/{net}.{sta}/{d.year}/{d.dayofyear:03d}.npz"
            if not os.path.exists(p): continue
            try: z=np.load(p)
            except Exception: continue
            ch=list(z["chans"]); X=z["x"]; V=z["valid"]
            for c in chs:
                if c in ch:
                    i=ch.index(c); DAT[(sta,c)]=(X[i],V[i])
        if not DAT: continue
        nspb=int(round(BLOCK_H*3600*FS))
        for b in range(int(24/BLOCK_H)):
            s0,s1=b*nspb,(b+1)*nspb
            bt0=d.timestamp()+b*BLOCK_H*3600
            FFTC={}
            for k,(x,v) in DAT.items():
                if not v[s0:s1].all(): continue          # his rule: complete gapless block only
                seg=x[s0:s1].astype(np.float64)
                if not np.any(seg): continue
                FFTC[k]=XP.asarray(seg)
            if not FFTC: continue
            nblk+=1
            for fam,dd in T.items():
                cs=[ncc(FFTC[k],XP.asarray(tpl)) for k,tpl in dd.items() if k in FFTC]
                nch=len(cs)
                if nch==0: continue
                n=min(int(c.size) for c in cs)          # templates in a family share a duration
                mcc=to_np(sum(c[:n] for c in cs)/nch)
                mad=np.median(np.abs(mcc-mcc.mean()))
                idx=np.where(mcc>=THR*mad)[0]
                if idx.size==0: continue
                keep=decluster(idx,mcc,maxdiff)
                for i in keep:
                    rows.append((fam,bt0+i/FS,float(mcc[i]),nch,float(mad)))
        if len(days)>1 and (d.dayofyear%5==0):
            log(f"  {d.date()} blocks={nblk} det={len(rows)} {time.time()-t00:.0f}s")
    D=pd.DataFrame(rows,columns=["fam","t_epoch","cc","nchannel","mad"])
    D["t"]=pd.to_datetime(D.t_epoch,unit="s").astype("datetime64[ns]")
    D.to_parquet(out)
    log(f"[duc-mf bed {a.bed}] {len(D)} detections over {nblk} blocks, {time.time()-t00:.0f}s -> {out}")
if __name__=="__main__": main()
