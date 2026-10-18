#!/usr/bin/env python
"""20 Hz / 1-8 Hz cache, frozen chain (lfe_picker/scripts/o_cache_month.cache_one), widened only
to accept 1-C and 2-C stations (the southern permanent network is largely vertical-only
short-period before ~2012 -- Ducellier's own catalogue used exactly those stations) and
non-integer fs (resample_poly with the exact rational fs/20, the documented twin in u_cache).

  cache20hz/<NET.STA>/<YYYY>/<DOY>.npz : x (C,N) float32, t0 float64, chans (C,), fs 20.0
  python scripts/t3_cache.py results/stage1_download_plan.csv [--workers 24]
"""
import os,sys,time,argparse,warnings
import numpy as np, pandas as pd
from fractions import Fraction
from scipy.signal import butter,sosfiltfilt,resample_poly
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
BANDPREF=["HH","BH","EH","SH","HN"]
def bp(x,fs,f1=1.0,f2=8.0):
    return sosfiltfilt(butter(4,[f1/(fs/2),f2/(fs/2)],btype="band",output="sos"),x,axis=-1)
def cache_one(job):
    net,sta,day,force=job
    netsta=f"{net}.{sta}"; d=pd.Timestamp(day)
    dst=npz_path(net,sta,day); src=day_path(TR,net,sta,day)
    if os.path.exists(dst) and not force: return dict(netsta=netsta,day=day,status="cached")
    if not os.path.exists(src): return dict(netsta=netsta,day=day,status="nosrc")
    try:
        from obspy import read,UTCDateTime
        st=read(src)
        for tr in st:                                   # snap near-integer rates (CN/NC fragments)
            fs_=float(tr.stats.sampling_rate); fr=round(fs_)
            if fr>0 and fs_!=fr and abs(fs_-fr)/fr<1e-4: tr.stats.sampling_rate=float(fr)
        st.merge(method=1,fill_value=0)
        t0=UTCDateTime(str(d)); st.trim(t0,t0+86400,pad=True,fill_value=0)
        bands={}
        for tr in st: bands.setdefault(tr.stats.channel[:2],[]).append(tr)
        band=None
        for b in BANDPREF:
            if b in bands and float(bands[b][0].stats.sampling_rate)>=20: band=b; break
        if band is None: return dict(netsta=netsta,day=day,status="noband")
        trs=sorted(bands[band],key=lambda t:t.stats.channel)
        fs0=float(trs[0].stats.sampling_rate)
        npts=min(t.stats.npts for t in trs)
        X=np.vstack([t.data[:npts].astype(np.float64) for t in trs])
        if not np.isfinite(X).all(): X=np.nan_to_num(X)
        X=bp(X,fs0)
        fr=Fraction(FS/fs0).limit_denominator(1000)
        if fr!=1: X=resample_poly(X,fr.numerator,fr.denominator,axis=-1)
        chans=[t.stats.channel for t in trs]
        os.makedirs(os.path.dirname(dst),exist_ok=True)
        tmp=dst+f".tmp{os.getpid()}.npz"
        np.savez(tmp,x=X.astype(np.float32),t0=np.float64(trs[0].stats.starttime.timestamp),
                 chans=np.array(chans),fs=np.float64(FS))
        os.replace(tmp,dst)
        return dict(netsta=netsta,day=day,status="built",band=band,nc=len(chans),
                    fs0=fs0,frac_nonzero=float((np.abs(X).sum(0)>0).mean()))
    except Exception as e:
        return dict(netsta=netsta,day=day,status="error",err=str(e)[:150])
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("plan"); ap.add_argument("--workers",type=int,default=24)
    ap.add_argument("--force",action="store_true"); ap.add_argument("--out",default=f"{RES}/cache_stage1.parquet")
    a=ap.parse_args(); warnings.filterwarnings("ignore")
    P=pd.read_csv(a.plan).drop_duplicates(["net","sta","day"])
    jobs=[(r.net,r.sta,str(r.day),a.force) for r in P.itertuples()]
    log(f"[cache] {len(jobs)} station-days, {a.workers} workers")
    import multiprocessing as mp
    t0=time.time()
    with mp.Pool(a.workers) as pool:
        rows=[]
        for i,r in enumerate(pool.imap_unordered(cache_one,jobs,chunksize=4),1):
            rows.append(r)
            if i%500==0: log(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s")
    D=pd.DataFrame(rows); D.to_parquet(a.out)
    log(f"[cache] {D.status.value_counts().to_dict()} in {time.time()-t0:.0f}s")
if __name__=="__main__": main()
