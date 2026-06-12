#!/usr/bin/env python
"""T2a shallow-share slope for B935 (the NON-CIRCULAR depth discriminator; FAMILY_TRUST_TEST.md §2).
Regress each family's 2-4s dv/v on the template-free station Z-autocorrelation dv/v (the pure shallow
monitor). SLOPE not correlation (the stable deep fault makes everyone correlate). Theil-Sen robust slope.
  slope ~ 1/3 -> genuine two-lobed family (deep source + shallow receiver; deep stable) = TRUSTED-eligible
  slope ~ 1   -> the family IS the shallow monitor = SITE-CARRIER (real but shallow)
Report on RAW (with seasonal) and DESEASONED (annual sinusoid removed) so a seasonal-amplitude-ratio
coincidence can't masquerade as coupling. Tag families GOLD/FAIL/other from the trust table.
"""
import os
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ.setdefault(v,"1")
import numpy as np, pandas as pd
from scipy.stats import theilslopes

ACF=pd.read_csv('data/autocorr_dvv_B935.csv')           # date,dvv,cc_max  (shallow monitor)
FAM=pd.read_csv('data/daily_dvv_B935_2to4_calT.csv')    # patch,date,dvv,cc_max,n_stack
TR =pd.read_csv('data/family_trust_tier1_B935.csv'); TR=TR[TR.kind=='F']
fmax=pd.read_csv('data/family_trust_tier1_B935.csv').query("kind=='R'").coda_sigma.max()
gold=set(TR[TR.coda_sigma>fmax].fam); fail=set(TR[TR.verdict=='FAIL'].fam)
sig={r.fam:r.coda_sigma for r in TR.itertuples()}

ACF['date']=pd.to_datetime(ACF.date); ACF=ACF.dropna(subset=['dvv']).sort_values('date')
def deseason(d,y):
    doy=pd.to_datetime(d).dt.dayofyear.values
    X=np.c_[np.ones_like(doy,float),np.cos(2*np.pi*doy/365.25),np.sin(2*np.pi*doy/365.25),
            np.cos(4*np.pi*doy/365.25),np.sin(4*np.pi*doy/365.25)]
    b,_,_,_=np.linalg.lstsq(X,y,rcond=None); return y-X@b

rows=[]
for fam,g in FAM.groupby('patch'):
    g=g.copy(); g['date']=pd.to_datetime(g.date)
    m=pd.merge(g[['date','dvv']], ACF[['date','dvv']], on='date', suffixes=('_f','_a')).dropna()
    if len(m)<80: continue
    # robust slope family ~ acf
    sl_raw,_,lo,hi=theilslopes(m.dvv_f.values, m.dvv_a.values)
    fr=deseason(m.date, m.dvv_f.values); ar=deseason(m.date, m.dvv_a.values)
    sl_des,_,_,_=theilslopes(fr, ar)
    # variance share: how much of family var is explained by acf (R2 of OLS)
    A=np.c_[np.ones(len(m)), m.dvv_a.values]; bb,_,_,_=np.linalg.lstsq(A,m.dvv_f.values,rcond=None)
    pred=A@bb; r2=1-np.sum((m.dvv_f.values-pred)**2)/np.sum((m.dvv_f.values-m.dvv_f.mean())**2)
    cls='GOLD' if fam in gold else ('FAIL' if fam in fail else 'other')
    rows.append(dict(fam=fam,cls=cls,coda_sigma=round(sig.get(fam,np.nan),1),n=len(m),
                     slope_raw=round(sl_raw,3),slope_des=round(sl_des,3),
                     slope_lo=round(lo,3),slope_hi=round(hi,3),r2=round(r2,2)))
D=pd.DataFrame(rows)
D.to_csv('data/exp_b935_t2a_slope.csv',index=False)
print(f"ACF dv/v: {len(ACF)} days. Families with >=80 shared days: {len(D)}\n")
print("== slope summary by class (slope~1/3 deep two-lobed; ~1 shallow site-carrier) ==")
for c in ['GOLD','FAIL','other']:
    s=D[D.cls==c]
    if not len(s): continue
    print(f"{c:6s} n={len(s):2d} | slope_raw med {s.slope_raw.median():+.2f} (IQR {s.slope_raw.quantile(.25):+.2f}..{s.slope_raw.quantile(.75):+.2f}) "
          f"| slope_deseason med {s.slope_des.median():+.2f} | r2 med {s.r2.median():.2f}")
print("\n== GOLD families (the 6 that passed T1 tail rule) ==")
print(D[D.cls=='GOLD'].sort_values('coda_sigma',ascending=False).to_string(index=False))
print("\n-> data/exp_b935_t2a_slope.csv\nT2A_DONE")
