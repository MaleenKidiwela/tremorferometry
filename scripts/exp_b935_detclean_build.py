#!/usr/bin/env python
"""B935 detection-level cleaning — build RAW and CLEANED daily-stack npz for target families.
Cleaning = keep only detections whose single-window coda (2-4s) correlates with the family all-time
reference coda above THR (probe showed GOLD real-LFE detections sit +0.18 over random; FAIL +0.05).
RAW = all QC-passed detections (reproduces long_window_daily_B935.npz). Both saved in the
long_window_daily format (stacks,patches,dates,t,fs) so dvv_roll30cal.py runs on either unchanged.
FAIL families are the circularity control: cleaning must NOT manufacture a dv/v signal for them.
"""
import os, sys, glob
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ.setdefault(v,"1")
import numpy as np, pandas as pd
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

NET,STA="PB","B935"; FS=40.0; PRE,POST=3.0,10.0
PRE_N,POST_N=int(PRE*FS),int(POST*FS); WIN_N=PRE_N+POST_N
T=(np.arange(WIN_N)-PRE_N)/FS; CM=(T>=2)&(T<=4)
FMIN,FMAX=2.0,8.0; THR=0.20; MIN_DET=3
FAMS=['40.800_-123.450__c18','40.050_-122.600__c130',      # GOLD
      '40.550_-122.800__c406','40.050_-122.600__c193']     # FAIL (control)

d=np.load('data/long_window_daily_B935.npz',allow_pickle=True)
pat=d['patches'].astype(str); stk=d['stacks']
REFC={}
for fam in FAMS:
    r=stk[pat==fam].astype(float).mean(0); rc=r[CM]; rc=rc-rc.mean(); REFC[fam]=rc/(np.linalg.norm(rc)+1e-12)

def _day(arg):
    day_str, by_fam = arg
    from obspy import read
    day=datetime.fromisoformat(day_str)
    p=f"data/waveforms/{NET}.{STA}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"
    if not os.path.exists(p): return []
    try:
        st=read(p).select(component="Z")
        if not len(st): return []
        for tr in st:
            if abs(tr.stats.sampling_rate-FS)>1e-6: tr.resample(FS)
        st.merge(fill_value=0); st.detrend("demean")
        st.filter("bandpass",freqmin=FMIN,freqmax=FMAX,corners=4,zerophase=True)
    except Exception: return []
    tr=st[0]; data=tr.data.astype(np.float64); start=tr.stats.starttime.datetime
    refc=REFC; out=[]
    for fam,tlist in by_fam.items():
        sraw=np.zeros(WIN_N); nraw=0; scln=np.zeros(WIN_N); ncln=0
        for tus in tlist:
            dt=datetime.utcfromtimestamp(tus/1e6); off=(dt-start).total_seconds()
            i0=int(round(off*FS))-PRE_N; i1=i0+WIN_N
            if i0<0 or i1>data.size: continue
            seg=data[i0:i1]; nrm=np.linalg.norm(seg)
            if nrm==0 or not np.isfinite(nrm): continue
            seg=seg/nrm; sraw+=seg; nraw+=1
            c=seg[CM]; c=c-c.mean(); cn=np.linalg.norm(c)
            if cn==0: continue
            if float(np.dot(c/cn, refc[fam]))>THR: scln+=seg; ncln+=1
        rec=[]
        if nraw>=MIN_DET:
            w=sraw/nraw; nn=np.linalg.norm(w)
            if nn>0: rec.append(('RAW',fam,day_str,(w/nn).astype(np.float32),nraw))
        if ncln>=MIN_DET:
            w=scln/ncln; nn=np.linalg.norm(w)
            if nn>0: rec.append(('CLN',fam,day_str,(w/nn).astype(np.float32),ncln))
        out.extend(rec)
    return out

def main():
    # detection times for all target families
    rows=[]
    for f in sorted(glob.glob('data/mf_b935_20*.csv')):
        m=pd.read_csv(f,usecols=['template','time']); m=m[m.template.isin(FAMS)]
        if len(m): rows.append(m)
    M=pd.concat(rows,ignore_index=True); M['t']=pd.to_datetime(M.time)
    M['day']=M.t.dt.strftime('%Y-%m-%d'); M['tus']=(M.t.values.view('int64')//1000)
    print(f"{len(M):,} detections across {M.day.nunique()} days, {M.template.nunique()} fams",flush=True)
    jobs={}
    for (day,fam),g in M.groupby(['day','template']):
        jobs.setdefault(day,{}).setdefault(fam,[]).extend(int(x) for x in g.tus.values)
    items=list(jobs.items())
    print(f"processing {len(items)} day-files (workers=24)...",flush=True)
    raw={'stacks':[],'patches':[],'dates':[],'n':[]}; cln={'stacks':[],'patches':[],'dates':[],'n':[]}
    done=0
    with ProcessPoolExecutor(max_workers=24,mp_context=mp.get_context("spawn")) as ex:
        for r in ex.map(_day, items, chunksize=4):
            for kind,fam,ds,w,n in r:
                tgt=raw if kind=='RAW' else cln
                tgt['stacks'].append(w); tgt['patches'].append(fam); tgt['dates'].append(ds); tgt['n'].append(n)
            done+=1
            if done%500==0: print(f"  {done}/{len(items)} days",flush=True)
    for tag,obj,fn in [('RAW',raw,'data/long_window_daily_B935_RAWexp.npz'),
                       ('CLN',cln,'data/long_window_daily_B935_CLEANexp.npz')]:
        S=np.array(obj['stacks'],dtype=np.float32); P=np.array(obj['patches']); D=np.array(obj['dates'])
        N=np.array(obj['n']); o=np.lexsort((D,P));
        np.savez(fn,stacks=S[o],patches=P[o],dates=D[o],n_det=N[o],t=T.astype(np.float32),fs=FS)
        print(f"{tag}: {len(S)} daily stacks -> {fn}",flush=True)
        for fam in FAMS:
            mm=P==fam; print(f"   {fam}: {mm.sum()} days, median {int(np.median(N[mm])) if mm.any() else 0} det/day")
    print("BUILD_DONE")

if __name__=="__main__": main()
