#!/usr/bin/env python
"""STAGE 3a-FIXED (agent AA, 2026-09-25) -- network matched filter at SHELLY (2025) parameters.

A COPY of t11_mf.py, never an edit: t11_mf.py is in use by the running t_south_queue.sh extension
stage (pid 1509294) and must not change under it.

TWO MEASURED FAULTS IN t11_mf.py THAT THIS COPY FIXES
-----------------------------------------------------
F1 -- cc normalisation blows up on zero-filled gaps.  `ncc()` divides by
sd = sqrt(max(M*S2 - S1*S1, 1e-30)) built from cumulative sums.  The 20 Hz cache is written with
obspy `st.merge(method=1, fill_value=0)` + `trim(pad=True, fill_value=0)` (`t3_cache.py:31-32`) and
stores NO validity mask, so every gap is an exactly-constant run.  In such a run the true window
variance is 0, the cumsum difference goes slightly NEGATIVE from round-off, the clamp gives
sd = 1e-15, and the "normalised" cc comes out at 1e10-1e15 instead of <=1.
MEASURED, Shelly family 23, 2022-05-11: BK.RVIT and NC.KIP are 33-50 % constant samples,
|cc|max = 2.4e15 and 2.8e10; the channel sum then has MAD 25.7 (a normalised cc sum over 20
channels cannot exceed 20) with 16.6 % of all samples above 8*MAD -> after 4 s declustering every
4 s bin of every day fires: 128,922 detections on that ONE day for 27 families.
Ducellier's path (`t7_duc_mf.py:113`) is immune only because his cache carries `valid` and he
requires a gapless block -- which is why the SAME chain reproduced his catalogue at recall 0.836.
  FIX: a window is live for a channel only if it contains at least one sample-to-sample change and
  its running std exceeds 1e-4 of the day's median window std.  Dead windows contribute nothing and
  are not counted as a live channel.

F2 -- the network stack was summed WITHOUT the per-station moveout.  Each station's template is cut
at that station's own offset tt_s + tau_s (`t9_templates.py:93`), so a cc peak at lag i means the
event time was t0 + i/FS - (tt_s + tau_s - pre), a DIFFERENT event time for every station.  Summing
the raw cc series at a common lag therefore stacks stations that are offset from one another by
tt_s + tau_s, which for Shelly family 1 spans 1.49 s (PB.B934) to 5.85 s (NC.KBN) -- many coda
cycles at 1.5-9 Hz -- so true events do not stack and the detection time carries a station-average
offset relative to the published time.  Ducellier's path is immune because his templates are all cut
at the SAME absolute catalogue time on every station (the moveout lives inside his 30-60 s
templates), which is the second reason his reproduction worked and this one did not.
  FIX: shift each channel's cc series by round((tt_s + tau_s - pre)*FS) before summing, so index j
  of the stack is the published-catalogue time axis t0 + j/FS, identical for every station.

PUBLISHED PARAMETERS NOW USED (Shelly et al. 2025, GRL 52, e2025GL116116)
  Fig. 4a  "a correlation sum detection threshold of 9 times daily median absolute deviation"
           -> --tmin 9, MAD over the DAY.
  Sec. 3   "we retain only the detection with highest correlation sum within each 4 s window"
           -> --decl bin, fixed 4 s bins (not a greedy exclusion radius).
  Fig. 2   "Waveforms are filtered 2-9 Hz" -> --hp2 adds a 2 Hz high-pass to the 1.5-9 Hz cache and
           to the template, so the effective band is Shelly's.
  Fig. 2a  "Original candidate LFE waveform (6 s)" -> a 6 s template is one of the rows we run.

NOTHING IS FILTERED.  Every family keeps its row; a threshold SWEEP (6,7,8,9,10,12,15) is written
per family-day, together with the statistic's value at every published detection time on that day,
so the threshold is reported as a separation curve rather than chosen.

  python scripts/t11_mf_aa.py --source Shelly2025 --templates results/templates_south_build \
      --days-file results/stage2_days_test.json --truth shelly --tmin 9 \
      --t0 2018-01-01 --t1 2024-12-01 --out results/aa_shelly/mf12_fixed.parquet
"""
import os,sys,time,argparse,warnings,json
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
try:
    import cupy as cp; XP=cp; GPU=True
except Exception: XP=np; GPU=False
def to_np(a): return cp.asnumpy(a) if GPU else a
SWEEP=[6.0,7.0,8.0,9.0,10.0,12.0,15.0]

def ncc_masked(data,tpl,relfloor=1e-4):
    """normalised cc (N-M+1,) and a live mask; flat / near-flat windows -> cc 0, live False."""
    M=tpl.size; N=data.size
    L=int(2**np.ceil(np.log2(N+M)))
    dot=XP.fft.irfft(XP.fft.rfft(data,L)*XP.conj(XP.fft.rfft(tpl,L)),L)[:N-M+1]
    c1=XP.concatenate([XP.zeros(1),XP.cumsum(data)]); c2=XP.concatenate([XP.zeros(1),XP.cumsum(data*data)])
    S1=c1[M:N+1]-c1[:N-M+1]; S2=c2[M:N+1]-c2[:N-M+1]
    st=XP.sqrt(XP.maximum(M*XP.sum(tpl*tpl)-XP.sum(tpl)**2,1e-30))
    sd=XP.sqrt(XP.maximum(M*S2-S1*S1,0.0))
    chg=(XP.abs(XP.diff(data))>0).astype(XP.float64)
    cch=XP.concatenate([XP.zeros(1),XP.cumsum(chg)])
    nchg=cch[M-1:N]-cch[:N-M+1]
    good=sd>0
    if bool(XP.any(good)):
        med=XP.median(sd[good]); live=(nchg>0)&(sd>relfloor*med)
    else:
        live=XP.zeros(sd.shape,dtype=bool)
    cc=XP.where(live,(M*dot-XP.sum(tpl)*S1)/(st*sd+1e-30),0.0)
    return cc,live

def shift_into(dst,dstlive,src,srclive,sh):
    """dst[j] += src[j+sh] (zero outside); dst is on the published-time axis."""
    n=dst.size
    if sh>=0:
        m=min(n,src.size-sh)
        if m<=0: return
        dst[:m]+=src[sh:sh+m]; dstlive[:m]+=srclive[sh:sh+m]
    else:
        s=-sh; m=min(n-s,src.size)
        if m<=0: return
        dst[s:s+m]+=src[:m]; dstlive[s:s+m]+=srclive[:m]

def decl_bin(nm,thr,gap_n):
    """Shelly sec.3: keep the highest value in each fixed 4 s window that clears thr."""
    n=(nm.size//gap_n)*gap_n
    if n==0: return np.array([],int)
    B=nm[:n].reshape(-1,gap_n)
    am=B.argmax(1); mx=B.max(1)
    ok=np.where(mx>=thr)[0]
    return ok*gap_n+am[ok]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--source",required=True)
    ap.add_argument("--t0",required=True); ap.add_argument("--t1",required=True)
    ap.add_argument("--tmin",type=float,default=9.0); ap.add_argument("--gap",type=float,default=4.0)
    ap.add_argument("--stations",default=""); ap.add_argument("--out",required=True)
    ap.add_argument("--days-file",default=""); ap.add_argument("--templates",default=f"{RES}/templates_south")
    ap.add_argument("--truth",default="",choices=["","shelly"],help="also record the statistic at "
                    "every published detection time (the separation curve)")
    ap.add_argument("--hp2",action="store_true",help="2 Hz high-pass -> Shelly's 2-9 Hz band")
    ap.add_argument("--align",default="moveout",choices=["moveout","none"],
                    help="moveout: shift each channel by tt+tau-pre (F2 fix); none: t11_mf.py behaviour")
    a=ap.parse_args(); warnings.filterwarnings("ignore")
    os.makedirs(os.path.dirname(a.out) or ".",exist_ok=True)
    SOS=None
    if a.hp2:
        from scipy.signal import butter; SOS=butter(4,2.0/(FS/2),btype="highpass",output="sos")
    def hp(x):
        if SOS is None: return x
        from scipy.signal import sosfiltfilt; return sosfiltfilt(SOS,x,axis=-1)
    M=pd.read_csv(a.templates+".csv"); M=M[M.source==a.source].copy()
    if a.stations: M=M[M.netsta.isin(a.stations.split(","))]
    Z=np.load(a.templates+".npz")
    TPL={r.key:Z[r.key] for r in M.itertuples() if r.key in Z}
    fams=sorted(M.fam.astype(str).unique()); stas=sorted(M.netsta.unique())
    PRE=float(M.pre_s.iloc[0]) if "pre_s" in M.columns else 2.0
    log(f"[mf-aa {a.source}] {len(fams)} fam, {len(stas)} sta, {len(TPL)} tpl, GPU={GPU}, "
        f"tmin={a.tmin} gap={a.gap} hp2={a.hp2} align={a.align} pre={PRE}")
    TRU=None
    if a.truth=="shelly":
        S=pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv")
        S.columns=[c.strip('%') for c in S.columns]
        S["t"]=(pd.to_datetime(dict(year=S.year,month=S.month,day=S.day,hour=S.hour,minute=S.minute,
                second=S.second.astype(int))).astype("datetime64[ns]")+pd.to_timedelta(S.second%1,unit="s"))
        TRU=S.rename(columns={"idnum":"fam","meancc":"cc"})[["fam","t","cc"]]
        TRU["fam"]=TRU.fam.astype(str); TRU["day"]=TRU.t.dt.floor("D")
        log(f"  truth: {len(TRU)} published detections, {TRU.fam.nunique()} families")
    DAYS=None
    if a.days_file:
        dd=json.load(open(a.days_file)); DAYS=set()
        for k,v in dd.items():
            if k.split("|")[0]==a.source: DAYS|={pd.Timestamp(x) for x in v}
        log(f"  restricted to {len(DAYS)} held-out days")
    gap_n=int(round(a.gap*FS))
    rows=[];cov=[];diag=[];tru=[];MISM={};t00=time.time()
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
            DAT[ns]=(XP.asarray(hp(x)),float(z["t0"]))
        cov.append(dict(day=str(d.date()),n_sta=len(DAT)))
        if not DAT: continue
        for f in fams:
            sub=M[(M.fam.astype(str)==f)&M.netsta.isin(DAT)]
            NS=None; SUM=None; LIVE=None; nch=0; t0=None; flat=[]
            for r in sub.itertuples():
                X,tt=DAT[r.netsta]
                if t0 is None:
                    t0=tt; NS=int(X.shape[1])
                    SUM=np.zeros(NS); LIVE=np.zeros(NS)
                T=TPL.get(r.key)
                if T is None: continue
                if T.shape[0]!=X.shape[0]:
                    MISM[(r.netsta,T.shape[0],int(X.shape[0]))]=MISM.get((r.netsta,T.shape[0],int(X.shape[0])),0)+1
                    continue
                sh=int(round((float(r.tt_s)+float(r.best_tau)-PRE)*FS)) if a.align=="moveout" else 0
                for c in range(X.shape[0]):
                    tp=XP.asarray(hp(T[c].astype(np.float64)) if SOS is not None else T[c].astype(np.float64))
                    v,lv=ncc_masked(X[c],tp)
                    shift_into(SUM,LIVE,to_np(v),to_np(lv).astype(np.float64),sh)
                    flat.append(float(to_np(XP.mean((~lv).astype(XP.float64))))); nch+=1
            if SUM is None or nch==0: continue
            ok=LIVE>0
            if not ok.any(): continue
            stat=np.zeros(NS); stat[ok]=SUM[ok]/LIVE[ok]
            med=np.median(stat[ok]); mad=np.median(np.abs(stat[ok]-med))+1e-30
            nm=np.where(ok,(stat-med)/mad,-np.inf)
            keep=decl_bin(nm,a.tmin,gap_n)
            for i in keep:
                rows.append((f,t0+i/FS,float(nm[i]),float(stat[i]),int(LIVE[i]),nch,len(sub)))
            row=dict(fam=f,day=str(d.date()),nchan=nch,nsta=len(sub),med_live=float(np.median(LIVE)),
                flat_frac=float(np.mean(flat)) if flat else np.nan,mad=float(mad),
                nm_p50=float(np.percentile(nm[ok],50)),nm_p99=float(np.percentile(nm[ok],99)),
                nm_p999=float(np.percentile(nm[ok],99.9)),nm_max=float(nm[ok].max()),
                stat_max=float(stat[ok].max()),frac_over=float((nm[ok]>=a.tmin).mean()),
                ndet=int(len(keep)))
            for s_ in SWEEP: row[f"n{int(s_)}"]=int(len(decl_bin(nm,s_,gap_n)))
            if TRU is not None:
                th=TRU[(TRU.fam==f)&(TRU.day==d)]
                row["n_theirs"]=int(len(th))
                for _,e in th.iterrows():
                    j=int(round((e.t.value/1e9-t0)*FS))
                    if 0<=j<NS:
                        lo=max(0,j-int(2*FS)); hi=min(NS,j+int(2*FS)+1)
                        tru.append(dict(fam=f,day=str(d.date()),t=e.t,cc_theirs=e.cc,
                            nm_at=float(nm[j]) if ok[j] else np.nan,
                            nm_pm2s=float(np.nanmax(np.where(ok[lo:hi],nm[lo:hi],np.nan))),
                            dt_best=float((np.nanargmax(np.where(ok[lo:hi],nm[lo:hi],np.nan))+lo-j)/FS),
                            nm_p999=row["nm_p999"],mad=row["mad"]))
            diag.append(row)
        if len(cov)%25==0: log(f"  {d.date()} days={len(cov)} det={len(rows)} {time.time()-t00:.0f}s")
    for (ns,ct,cd),k in sorted(MISM.items()):
        log(f"  [channels] {ns}: template {ct} ch, {k} day(s) cached {cd} ch -- not stacked")
    D=pd.DataFrame(rows,columns=["fam","t_epoch","net_mad","meancc","nlive","nchan","nsta"])
    if len(D): D["t"]=pd.to_datetime(D.t_epoch,unit="s").astype("datetime64[ns]")
    D.to_parquet(a.out)
    pd.DataFrame(cov).to_csv(a.out.replace(".parquet","_coverage.csv"),index=False)
    pd.DataFrame(diag).to_csv(a.out.replace(".parquet","_diag.csv"),index=False)
    if tru: pd.DataFrame(tru).to_csv(a.out.replace(".parquet","_truthstat.csv"),index=False)
    log(f"[mf-aa {a.source}] {len(D)} detections {a.t0}..{a.t1} in {time.time()-t00:.0f}s -> {a.out}")
if __name__=="__main__": main()
