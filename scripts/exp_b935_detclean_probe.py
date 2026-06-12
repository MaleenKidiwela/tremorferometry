#!/usr/bin/env python
"""B935 detection-level cleaning PROBE (autonomous test, 2026-06-11).
Question: can we clean a family's detection list at the DETECTION level by coda coherence?
Only feasible if real-LFE detections have separably higher single-detection coda-CC (2-4s) with the
family reference than noise/chance detections do. Single detections are low-SNR, so this is NOT obvious
-> measure it empirically before committing to a full cleaning+restack run.

For GOLD c18 and FAIL c406 (B935): harvest detection windows the SAME way build_long_window_daily did
(Z, 40Hz, demean, bandpass 2-8, L2-normalize), compute each detection's coda-CC vs the family all-time
reference coda (from long_window_daily_B935.npz), and compare to RANDOM-time windows.
"""
import os, sys, glob
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ.setdefault(v,"1")
import numpy as np, pandas as pd
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

NET,STA="PB","B935"; FS=40.0; PRE,POST=3.0,10.0
PRE_N,POST_N=int(PRE*FS),int(POST*FS); WIN_N=PRE_N+POST_N
T=(np.arange(WIN_N)-PRE_N)/FS
CM=(T>=2)&(T<=4)            # coda window 2-4 s
FMIN,FMAX=2.0,8.0
GOLD='40.800_-123.450__c18'; FAIL='40.550_-122.800__c406'
NDAYS=45; SEED=7

# references (all-time mean daily stack, coda portion, zero-mean unit-norm)
d=np.load('data/long_window_daily_B935.npz',allow_pickle=True)
pat=d['patches'].astype(str); stk=d['stacks']
def ref_coda(fam):
    S=stk[pat==fam].astype(float); r=S.mean(0); rc=r[CM]; rc=rc-rc.mean(); return rc/ (np.linalg.norm(rc)+1e-12)
REF={GOLD:ref_coda(GOLD), FAIL:ref_coda(FAIL)}

def _day(arg):
    day_str, jobs = arg          # jobs: list of (fam, t_us); plus random probes fam='RAND'
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
    tr=st[0]; data=tr.data.astype(np.float64); start=tr.stats.starttime.datetime; out=[]
    for fam,tus in jobs:
        dt=datetime.utcfromtimestamp(tus/1e6); off=(dt-start).total_seconds()
        i0=int(round(off*FS))-PRE_N; i1=i0+WIN_N
        if i0<0 or i1>data.size: continue
        seg=data[i0:i1]; nrm=np.linalg.norm(seg)
        if nrm==0 or not np.isfinite(nrm): continue
        seg=seg/nrm
        c=seg[CM]; c=c-c.mean(); cn=np.linalg.norm(c)
        if cn==0: continue
        c=c/cn
        rfam = GOLD if fam in ('RAND_G',) else FAIL if fam in ('RAND_F',) else fam
        cc=float(np.dot(c, REF[rfam]))
        out.append((fam, cc))
    return out

def load_times(fam, n):
    rng=np.random.RandomState(SEED); rows=[]
    for f in sorted(glob.glob('data/mf_b935_20*.csv')):
        m=pd.read_csv(f,usecols=['template','time']); m=m[m.template==fam]
        if len(m): rows.append(m)
    M=pd.concat(rows,ignore_index=True); M['t']=pd.to_datetime(M.time)
    M['day']=M.t.dt.strftime('%Y-%m-%d')
    days=rng.choice(M.day.unique(), min(NDAYS,M.day.nunique()), replace=False)
    return M[M.day.isin(days)]

def main():
    rng=np.random.RandomState(SEED)
    jobs={}   # day -> list[(fam,t_us)]
    for fam in (GOLD,FAIL):
        M=load_times(fam,NDAYS)
        for day,g in M.groupby('day'):
            jobs.setdefault(day,[]).extend((fam,int(t.value/1e3)) for t in g.t)
        # random-time probes on the SAME days (calibration): one random time per detection
        for day,g in M.groupby('day'):
            base=pd.Timestamp(day)
            rt=base+pd.to_timedelta(rng.uniform(0,86400,len(g)),unit='s')
            tag='RAND_G' if fam==GOLD else 'RAND_F'
            jobs.setdefault(day,[]).extend((tag,int(t.value/1e3)) for t in rt)
    items=list(jobs.items())
    print(f"harvesting {sum(len(v) for v in jobs.values())} windows over {len(items)} days",flush=True)
    res=[]
    with ProcessPoolExecutor(max_workers=20,mp_context=mp.get_context("spawn")) as ex:
        for r in ex.map(_day, items, chunksize=2): res.extend(r)
    R=pd.DataFrame(res,columns=['fam','codaCC'])
    print(f"\nharvested {len(R)} windows\n")
    for tag,label in [(GOLD,'GOLD c18 detections'),('RAND_G','  GOLD random-time'),
                      (FAIL,'FAIL c406 detections'),('RAND_F','  FAIL random-time')]:
        x=R[R.fam==tag].codaCC.values
        if not len(x): print(f"{label}: none"); continue
        print(f"{label:24s} n={len(x):5d}  codaCC mean {x.mean():+.3f} med {np.median(x):+.3f} "
              f"std {x.std():.3f} | frac>0.2 {np.mean(x>0.2):.2f} frac>0.4 {np.mean(x>0.4):.2f}")
    # separability: how much higher is real-detection codaCC than its random baseline?
    for fam,rnd,lab in [(GOLD,'RAND_G','GOLD'),(FAIL,'RAND_F','FAIL')]:
        a=R[R.fam==fam].codaCC.values; b=R[R.fam==rnd].codaCC.values
        if len(a) and len(b):
            from scipy.stats import mannwhitneyu
            try: u,pp=mannwhitneyu(a,b,alternative='greater')
            except Exception: pp=np.nan
            print(f"[{lab}] det-mean {a.mean():+.3f} vs random-mean {b.mean():+.3f}  "
                  f"shift {a.mean()-b.mean():+.3f}  (MWU p>{pp:.1e})")
    R.to_csv('data/exp_b935_detclean_probe.csv',index=False)
    print("\n-> data/exp_b935_detclean_probe.csv\nPROBE_DONE")

if __name__=="__main__": main()
