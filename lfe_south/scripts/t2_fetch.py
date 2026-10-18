#!/usr/bin/env python
"""Download the station-days of a plan CSV (net,sta,day) at <=8 concurrent FDSN requests.

Native sample rate, no filtering/merging: the bytes that come back are what is written, exactly
as lfe_picker/scripts/r_download_day.fetch_day does.  `.nd` marks a verified empty day.
Appends to results/download_manifest.csv so the run is resumable.
  python scripts/t2_fetch.py results/stage1_download_plan.csv [--workers 8]
"""
import os,sys,csv,time,argparse,threading
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
CH="BH?,HH?,EH?,SH?,HN?".replace(",HN?","")           # velocity channels only
_lock=threading.Lock(); _n={"data":0,"nd":0,"err":0,"bytes":0}
_clients={}
def client(prov):
    import threading as th
    from obspy.clients.fdsn import Client
    k=(prov,th.get_ident())
    if k not in _clients: _clients[k]=Client(prov,timeout=180)
    return _clients[k]
def fetch(job):
    net,sta,day=job
    from obspy import UTCDateTime
    p=day_path(TR,net,sta,day)
    if os.path.exists(p): return dict(net=net,sta=sta,day=day,status="data",bytes=os.path.getsize(p),note="resumed")
    if os.path.exists(p+".nd"): return dict(net=net,sta=sta,day=day,status="nodata",bytes=0,note="resumed")
    t0=UTCDateTime(str(day)); t1=t0+86400
    prov=PROV.get(net,"IRIS")
    for attempt in range(3):
        try:
            st=client(prov).get_waveforms(net,sta,"*",CH,t0-60,t1+60)
            if len(st)==0: raise ValueError("empty")
            os.makedirs(os.path.dirname(p),exist_ok=True)
            tmp=p+f".tmp{os.getpid()}"
            st.write(tmp,format="MSEED"); os.replace(tmp,p)
            with _lock: _n["data"]+=1; _n["bytes"]+=os.path.getsize(p)
            return dict(net=net,sta=sta,day=day,status="data",bytes=os.path.getsize(p),
                        note=",".join(sorted({tr.stats.channel for tr in st})))
        except Exception as e:
            msg=str(e)[:120]
            if "No data" in msg or "empty" in msg or "204" in msg:
                os.makedirs(os.path.dirname(p),exist_ok=True); open(p+".nd","w").close()
                with _lock: _n["nd"]+=1
                return dict(net=net,sta=sta,day=day,status="nodata",bytes=0,note=msg)
            time.sleep(2*(attempt+1))
    with _lock: _n["err"]+=1
    return dict(net=net,sta=sta,day=day,status="error",bytes=0,note=msg)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("plan"); ap.add_argument("--workers",type=int,default=8)
    ap.add_argument("--manifest",default=f"{RES}/download_manifest.csv")
    a=ap.parse_args()
    P=pd.read_csv(a.plan).drop_duplicates(["net","sta","day"])
    jobs=[(r.net,r.sta,str(r.day)) for r in P.itertuples()]
    log(f"[fetch] {len(jobs)} station-days, {a.workers} threads")
    new=not os.path.exists(a.manifest)
    f=open(a.manifest,"a",newline=""); w=csv.DictWriter(f,fieldnames=["net","sta","day","status","bytes","note"])
    if new: w.writeheader()
    t0=time.time()
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for i,r in enumerate(ex.map(fetch,jobs),1):
            w.writerow(r)
            if i%200==0:
                f.flush(); log(f"  {i}/{len(jobs)} data={_n['data']} nd={_n['nd']} err={_n['err']} "
                               f"{_n['bytes']/1e9:.1f} GB {time.time()-t0:.0f}s")
    f.close(); log(f"[fetch] done data={_n['data']} nodata={_n['nd']} err={_n['err']} {_n['bytes']/1e9:.2f} GB")
if __name__=="__main__": main()
