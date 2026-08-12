#!/usr/bin/env python
"""Complete Merlin's decision table for B926 H2: (a) noise-floor PSD H1/H2/Z, (b) azimuth spread,
plus a skeptical add: H2-vs-Z coda correlation (independent horizontal info vs Z-bleed)."""
import os
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"): os.environ.setdefault(v,"1")
import numpy as np, pandas as pd, glob, obspy
from scipy.signal import welch, butter, sosfiltfilt

STA_LAT, STA_LON, EPI_BOX, EPI_MIN = 48.82, -124.131, 0.4, 3
cat = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time","lat","lon"])
cat = cat[cat.lat.between(STA_LAT-EPI_BOX,STA_LAT+EPI_BOX) & cat.lon.between(STA_LON-EPI_BOX,STA_LON+EPI_BOX)]
cat["d"] = pd.to_datetime(cat.time).dt.strftime("%Y-%m-%d")
g = cat.groupby("d").size(); epi_days = sorted(g[g>=EPI_MIN].index)
cert = list(pd.read_csv("data/b926_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)

# ---- (b) azimuth spread of certified families ----
def azi(f):
    p=f.split("_"); la,lo=float(p[0]),float(p[1])
    import math
    dlon=math.radians(lo-STA_LON)
    y=math.sin(dlon)*math.cos(math.radians(la))
    x=math.cos(math.radians(STA_LAT))*math.sin(math.radians(la))-math.sin(math.radians(STA_LAT))*math.cos(math.radians(la))*math.cos(dlon)
    return (math.degrees(math.atan2(y,x))+360)%360
az = np.array([azi(f) for f in cert])
print(f"(b) certified-family azimuth: span {az.min():.0f}-{az.max():.0f} deg, "
      f"range {az.max()-az.min():.0f}, std {az.std():.0f}")

# ---- H2-vs-Z coda correlation per certified family (bleed check) ----
dZ=np.load("data/long_window_daily_B926p90f40_Z.npz",allow_pickle=True)
dH=np.load("data/long_window_daily_B926p90f40_H2.npz",allow_pickle=True)
t=dZ["t"]; coda=(t>=2)&(t<=4)
epi=set(epi_days)
def grand(d,fam):
    p,dt,nd,S=d["patches"],d["dates"],d["n_det"].astype(float),d["stacks"]
    m=p==fam; dm=np.array([x in epi for x in dt[m]]); w=nd[m][dm]
    return (S[m][dm]*w[:,None]).sum(0)/w.sum() if w.sum()>200 else None
bleed=[]
for f in cert:
    gz,gh=grand(dZ,f),grand(dH,f)
    if gz is None or gh is None: continue
    a=gz[coda]-gz[coda].mean(); b=gh[coda]-gh[coda].mean()
    bleed.append(abs(a@b/((np.linalg.norm(a)*np.linalg.norm(b))+1e-30)))
bleed=np.array(bleed)
print(f"(add) H2-vs-Z coda |cc| per family: median {np.median(bleed):.2f} 90th {np.quantile(bleed,.9):.2f}")
print("   -> HIGH (>0.7)=H2 is Z-bleed (not independent); LOW/MID=independent horizontal coda")

# ---- (a) raw PSD floor H1/H2/Z on sampled episode days ----
sos=butter(4,[2/50.,8/50.],btype='band',output='sos')
psd={'EH1':[], 'EH2':[], 'EHZ':[]}
picked=0
for d in epi_days[::max(1,len(epi_days)//25)]:
    yr,doy = d[:4], pd.Timestamp(d).dayofyear
    f=f"data/waveforms/PB.B926/{yr}/{doy:03d}.mseed"
    if not os.path.exists(f): continue
    try:
        st=obspy.read(f)
        for tr in st:
            ch=tr.stats.channel
            if ch not in psd: continue
            x=tr.data.astype(float); x=x-x.mean()
            xf=sosfiltfilt(sos,x)
            psd[ch].append(np.sqrt(np.mean(xf**2)))   # 2-8Hz RMS as floor proxy
        picked+=1
    except Exception: pass
    if picked>=20: break
med={k:np.median(v) for k,v in psd.items() if v}
print(f"(a) raw 2-8Hz RMS floor (n={picked} days): "
      f"EHZ {med.get('EHZ',0):.1f}  EH1 {med.get('EH1',0):.1f}  EH2 {med.get('EH2',0):.1f}")
if med.get('EH2') and med.get('EH1'):
    print(f"   H1/H2 floor ratio = {med['EH1']/med['EH2']:.2f}  (>~3 => H1 just a noisier axis, H2 survives)")
