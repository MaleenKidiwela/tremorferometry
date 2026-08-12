#!/usr/bin/env python
"""Direct single-station Z-channel AUTOCORRELATION dv/v — NO templates, NO catalog.
Per day: bandpass 2-8 Hz Z, resample_poly to 40 Hz, cut into hourly windows, FFT-autocorrelate each,
normalize+stack -> daily ACF (positive lags 0-10 s). Then 30-cal-day rolling stack of daily ACFs,
stretch each on the 2-4 s LAG window vs the all-time-mean ACF -> dv/v(date).
This is the textbook passive single-station autocorrelation monitor. If the LFE 'families' are accidentally
doing this, this independent dv/v should reproduce their seasonal.
"""
import os, sys, glob, argparse
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ.setdefault(v,"1")
import numpy as np, pandas as pd
from scipy.signal import resample_poly, butter, sosfiltfilt
sys.path.insert(0,"src")
from tremorferometry.dvv import stretch_dvv

FS=40.0; LAGS=int(10*FS)  # 10 s of positive lag
WHITEN=True
def daily_acf(args):
    f,=args
    try:
        import obspy
        st=obspy.read(f); st.merge(fill_value=0); tr=st[0]
        x=tr.data.astype(float)
        if tr.stats.sampling_rate<=0 or len(x)<tr.stats.sampling_rate*3600: return None
        # resample to 40 Hz via poly
        sr=int(round(tr.stats.sampling_rate))
        if sr!=40:
            from math import gcd
            g=gcd(sr,40); x=resample_poly(x,40//g,sr//g)
        x=x-x.mean()
        sos=butter(4,[2/(FS/2),8/(FS/2)],btype='band',output='sos'); x=sosfiltfilt(sos,x)
        win=int(3600*FS); acf=np.zeros(LAGS); n=0
        from scipy.ndimage import uniform_filter1d
        for i in range(0,len(x)-win,win):
            w=x[i:i+win]; w=w-w.mean();
            if w.std()==0: continue
            w=w/w.std()
            X=np.fft.rfft(w,2*win)
            if WHITEN:
                # spectral whitening: flatten amplitude over a smoothed spectrum (removes seasonal
                # noise-SOURCE spectral variation that otherwise reads as fake dv/v; Fable flag, 06-10),
                # keep phase / Green's-function structure. smoothing ~0.1 Hz; 1% waterlevel.
                amp=np.abs(X); nb=max(int(0.1*2*win/FS),11)
                sm=uniform_filter1d(amp,size=nb); wl=sm.mean()*0.01
                X=X/(sm+wl)
            ac=np.fft.irfft(X*np.conj(X))[:LAGS]
            if ac[0]!=0: acf+=ac/ac[0]; n+=1
        if n<6: return None
        date=os.path.basename(os.path.dirname(os.path.dirname(f)))  # not reliable; use jday
        return (f,acf/n)
    except Exception:
        return None

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--net",required=True); ap.add_argument("--sta",required=True)
    ap.add_argument("--y0",type=int,default=2013); ap.add_argument("--y1",type=int,default=2022)
    ap.add_argument("--workers",type=int,default=16); ap.add_argument("--out",required=True)
    ap.add_argument("--no-whiten",action="store_true",help="disable spectral whitening (default ON)")
    a=ap.parse_args()
    globals()['WHITEN']=not a.no_whiten
    base=f"data/waveforms/{a.net}.{a.sta}"
    files=[]
    for y in range(a.y0,a.y1+1):
        for f in sorted(glob.glob(f"{base}/{y}/*.mseed")):
            jday=int(os.path.basename(f).split('.')[0]); files.append((f,y,jday))
    print(f"{a.sta}: {len(files)} day-files {a.y0}-{a.y1}",flush=True)
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    dates=[]; acfs=[]
    with ProcessPoolExecutor(max_workers=a.workers,mp_context=mp.get_context("spawn")) as ex:
        for (f,y,jd),r in zip(files,ex.map(daily_acf,[((f,)) for f,_,_ in files],chunksize=8)):
            if r is None: continue
            dates.append(pd.Timestamp(f"{y}-01-01")+pd.Timedelta(days=jd-1)); acfs.append(r[1])
    A=np.array(acfs); D=pd.to_datetime(dates)
    o=np.argsort(D); A=A[o]; D=D[o]
    print(f"  {len(A)} daily ACFs",flush=True)
    # 30-cal-day rolling stack + stretch on 2-4 s LAG window vs all-time mean
    Dn=D.values.astype('datetime64[D]'); idx=np.arange(len(A))
    los=np.searchsorted(Dn,Dn-np.timedelta64(30,'D'),side='right'); cnt=idx+1-los
    cs=np.vstack([np.zeros((1,A.shape[1])),np.cumsum(A,0)])
    R=(cs[idx+1]-cs[los])/cnt[:,None]; keep=cnt>=5; R=R[keep]; Dk=D[keep]
    ref=R.mean(0); t=np.arange(LAGS)/FS
    rows=[]
    for i in range(len(R)):
        try:
            r=stretch_dvv(ref,R[i],fs=FS,t_min=2.0,t_max=4.0,eps_max=0.02,n_eps=201)
            rows.append((Dk[i],r.dvv,r.cc_max))
        except Exception: pass
    out=pd.DataFrame(rows,columns=['date','dvv','cc_max']); out.to_csv(a.out,index=False)
    print(f"  wrote {a.out}: {len(out)} days, mean cc {out.cc_max.mean():.3f}",flush=True)

if __name__=="__main__": main()
