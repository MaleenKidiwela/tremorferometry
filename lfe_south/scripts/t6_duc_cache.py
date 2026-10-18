#!/usr/bin/env python
"""STAGE R1b -- Ducellier's OWN preprocessing chain, verbatim from
catalogs/ducellier_catalog/src/get_data.py (get_from_IRIS / get_from_NCEDC):

  detrend('linear') -> taper(hann, max_length=5) -> remove_response(VEL,
  pre_filt=(0.2,0.5,10,15), water_level=80) -> bandpass 1.5-9 Hz zerophase ->
  interpolate(ratio.denominator*fs, lanczos a=10) -> decimate(ratio.numerator, no_filter)
  with ratio = Fraction(int(fs), 20)            [chain=rational when denominator > 1]

Single-component vertical stations are kept: each channel is cached and correlated on its own,
which is what find_all_LFEs_parallel.py does (it loops over CHANNELS, not stations).
A per-sample VALID mask records where a continuous trace segment existed before any fill, so the
matched filter can reproduce his rule that a channel contributes to a block only when it has a
complete, gapless record of it.

  cache_duc/<NET>.<STA>/<YYYY>/<DOY>.npz : x (C,N) f32, valid (C,N) bool, chans, t0, fs, chain
  python scripts/t6_duc_cache.py results/duc_download_plan.csv [--workers 16]
"""
import os,sys,time,argparse,warnings
from fractions import Fraction
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
CACHE_DUC=f"{WS}/cache_duc"; RESP=f"{WS}/resp"
FILT=(1.5,9.0); BANDPREF=["BH","HH","EH","SH"]
_HIS=None
def his_channels(sta):
    """the channels HIS stations_permanent.txt names for this station (union over his epochs).
    Using his list is a fidelity requirement: his templates are keyed <STA>_<CHA>, and a station
    like BK.YBH carries both BH (40 Hz) and HH (80 Hz) -- only his choice reproduces his network."""
    global _HIS
    if _HIS is None:
        _HIS={}
        for ln in open(f"{ROOT}/catalogs/ducellier_catalog/data/Ducellier/stations_permanent.txt"):
            q=ln.split()
            if len(q)>=3: _HIS.setdefault(q[0],set()).update(q[2].split(","))
    return _HIS.get(sta)
def resp_path(net,sta): return f"{RESP}/{net}.{sta}.xml"
def get_resp(net,sta,prov):
    p=resp_path(net,sta)
    if os.path.exists(p): return p
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    os.makedirs(RESP,exist_ok=True)
    for pv in ([prov] if prov else [])+["NCEDC","IRIS"]:
        try:
            inv=Client(pv,timeout=180).get_stations(network=net,station=sta,level="response",
                    starttime=UTCDateTime("2000-01-01"),endtime=UTCDateTime("2026-09-23"))
            inv.write(p,format="STATIONXML"); return p
        except Exception: continue
    return None
def cache_one(job):
    net,sta,day,prov,force=job
    d=pd.Timestamp(day); dst=f"{CACHE_DUC}/{net}.{sta}/{d.year}/{d.dayofyear:03d}.npz"
    if os.path.exists(dst) and not force: return dict(net=net,sta=sta,day=day,status="cached")
    src=day_path(TR,net,sta,day)
    if not os.path.exists(src): return dict(net=net,sta=sta,day=day,status="nosrc")
    rp=resp_path(net,sta)
    if not os.path.exists(rp): return dict(net=net,sta=sta,day=day,status="noresp")
    try:
        from obspy import read,read_inventory,UTCDateTime
        inv=read_inventory(rp,format="STATIONXML")
        D=read(src).select(channel="[BHES]H?")
        want=his_channels(sta)
        if want:
            sel=[tr for tr in D if tr.stats.channel in want]
        else:                                    # bed-B station: his list has none -> best band
            have={tr.stats.channel[:2] for tr in D}
            b=next((x for x in BANDPREF if x in have),None)
            sel=[tr for tr in D if tr.stats.channel[:2]==b] if b else []
        from obspy.core.stream import Stream
        D=Stream(sel)
        if len(D)==0: return dict(net=net,sta=sta,day=day,status="nochan")
        t0=UTCDateTime(str(d))-60; t1=UTCDateTime(str(d))+86400+60
        D.trim(t0,t1)
        segs={}                                        # per-channel valid segments (pre-fill)
        for tr in D: segs.setdefault(tr.stats.channel,[]).append(
            (tr.stats.starttime.timestamp,tr.stats.endtime.timestamp))
        D.detrend(type="linear")
        D.taper(type="hann",max_percentage=None,max_length=5.0)
        D.remove_response(inventory=inv,output="VEL",pre_filt=(0.2,0.5,10.0,15.0),water_level=80.0)
        D.filter("bandpass",freqmin=FILT[0],freqmax=FILT[1],zerophase=True)
        from obspy.core.stream import Stream
        rates=sorted({float(tr.stats.sampling_rate) for tr in D})
        out=Stream(); chains=set(); fs0=rates[0]
        for fr in rates:                          # one ratio per native rate, his formula exactly
            G=Stream([tr for tr in D if float(tr.stats.sampling_rate)==fr])
            ratio=Fraction(int(fr),int(FS))
            chains.add("frozen" if ratio.denominator==1 else "rational")
            G.interpolate(ratio.denominator*fr,method="lanczos",a=10)
            G.decimate(ratio.numerator,no_filter=True)
            out+=G
        D=out; chain="+".join(sorted(chains))
        D.merge(method=1,fill_value=0)
        tA=UTCDateTime(str(d)); tB=tA+86400
        D.trim(tA,tB,pad=True,fill_value=0)
        chans=sorted({tr.stats.channel for tr in D})
        N=int(round(86400*FS))
        X=np.zeros((len(chans),N),np.float32); V=np.zeros((len(chans),N),bool)
        for i,c in enumerate(chans):
            tr=D.select(channel=c)[0]
            n=min(tr.stats.npts,N); X[i,:n]=tr.data[:n]
            for (a0,a1) in segs.get(c,[]):
                i0=max(0,int(round((a0-tA.timestamp)*FS))); i1=min(N,int(round((a1-tA.timestamp)*FS)))
                if i1>i0: V[i,i0:i1]=True
        os.makedirs(os.path.dirname(dst),exist_ok=True)
        tmp=dst+f".tmp{os.getpid()}.npz"
        np.savez(tmp,x=X,valid=V,chans=np.array(chans),t0=np.float64(tA.timestamp),
                 fs=np.float64(FS),chain=np.array(chain),fs0=np.float64(fs0))
        os.replace(tmp,dst)
        return dict(net=net,sta=sta,day=day,status="built",nc=len(chans),fs0=fs0,chain=chain,
                    valid_frac=float(V.mean()))
    except Exception as e:
        return dict(net=net,sta=sta,day=day,status="error",err=str(e)[:160])
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("plan"); ap.add_argument("--workers",type=int,default=16)
    ap.add_argument("--force",action="store_true"); ap.add_argument("--out",default=f"{RES}/cache_duc.parquet")
    a=ap.parse_args(); warnings.filterwarnings("ignore")
    P=pd.read_csv(a.plan).drop_duplicates(["net","sta","day"])
    for (net,sta),_ in P.groupby(["net","sta"]):
        if not os.path.exists(resp_path(net,sta)):
            r=get_resp(net,sta,PROV.get(net,"IRIS")); log(f"  resp {net}.{sta} {'ok' if r else 'FAIL'}")
    jobs=[(r.net,r.sta,str(r.day),PROV.get(r.net,"IRIS"),a.force) for r in P.itertuples()]
    log(f"[duc-cache] {len(jobs)} station-days, {a.workers} workers")
    import multiprocessing as mp
    t0=time.time(); rows=[]
    with mp.Pool(a.workers) as pool:
        for i,r in enumerate(pool.imap_unordered(cache_one,jobs,chunksize=2),1):
            rows.append(r)
            if i%200==0: log(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s")
    D=pd.DataFrame(rows); D.to_parquet(a.out)
    log(f"[duc-cache] {D.status.value_counts().to_dict()} in {time.time()-t0:.0f}s")
    if "chain" in D: log(f"  chain: {D.chain.value_counts(dropna=True).to_dict()}")
if __name__=="__main__": main()
