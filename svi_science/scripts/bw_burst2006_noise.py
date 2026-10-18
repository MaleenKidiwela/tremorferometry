import os, sys
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import numpy as np, pandas as pd
from scipy.signal import butter, sosfiltfilt
import bw_p4env_detect as E
FS=20.0; SOS=butter(4,[2/10,8/10],"band",output="sos")
stas=["PGC","B010","B004","VGZ","NLLB","SNB","TWKB","SILB","KLNB","SSIB","TSJB","LZB","MGCB","TWBB"]
days=[str(x.date()) for x in pd.date_range("2006-03-18","2006-04-02")]+["2006-01-16","2006-01-10","2006-05-26","2006-05-20"]
rows=[]
for d in days:
    for s in stas:
        p=E.cache_path(s,d)
        if p is None: continue
        z=np.load(p); x=z["x"].astype(float); t0=float(z["t0"]); day0=pd.Timestamp(d).value/1e9
        n=int(86400*FS); X=np.zeros((3,n)); i0=int(round((t0-day0)*FS)); a,b=max(i0,0),min(i0+x.shape[1],n)
        X[:,a:b]=x[:,a-i0:b-i0]; X[~np.isfinite(X)]=0; gap=~np.any(X!=0,axis=0)
        Y=sosfiltfilt(SOS,X,axis=1); e=np.sqrt((Y**2).sum(0)); e[gap]=np.nan
        k=int(10*FS); r=np.sqrt(np.nanmean(e.reshape(-1,k)**2,1))   # 10-s rms
        h=r.reshape(24,360)
        for hh in range(24):
            seg=h[hh]; ok=np.isfinite(seg).mean()
            rows.append(dict(date=d,sta=s,hour=hh,p10=np.nanpercentile(seg,10) if ok>0.2 else np.nan,med=np.nanmedian(seg) if ok>0.2 else np.nan,cov=ok))
pd.DataFrame(rows).to_csv("/home/jovyan/tremorferometry/svi_science/results/backward/first_look/burst2006/station_rms_hourly_abs.csv",index=False,float_format="%.4g")
