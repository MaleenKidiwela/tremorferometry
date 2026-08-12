#!/usr/bin/env python
"""T2a as a TRUSTWORTHY SITE FILTER (Fable-corrected 2026-06-11). Per family: regress DESEASONED
family 2-4s dv/v on the DESEASONED station Z-autocorrelation dv/v; OLS slope + R2; CI by 30-DAY BLOCK
BOOTSTRAP (honest for the 30-day-rolling autocorrelated series, ~120 indep pts not ~3600). Classify by
whether the slope CI excludes 0 (sign-aware):
  SHALLOW-COUPLED (ci_lo>0)  -> family dv/v tracks the shallow signal = candidate SITE-CARRIER
  DECOUPLED       (CI spans 0)-> NOT shallow (deep-candidate OR noise; T2a cannot separate)
  ANOMALOUS       (ci_hi<0)   -> robust NEGATIVE coupling, impossible for a mixing share = model violation
This is a contamination/site filter, NOT a depth certificate (slope magnitude geometry is NOT interpreted).
Usage: exp_t2a_fixed.py <STA>
"""
import os, sys
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"): os.environ.setdefault(v,"1")
import numpy as np, pandas as pd

STA=sys.argv[1] if len(sys.argv)>1 else 'B018'
NB=500; BLK=30; SEED=11
rng=np.random.RandomState(SEED)
ACF=pd.read_csv(f'data/autocorr_dvv_{STA}.csv'); ACF['date']=pd.to_datetime(ACF.date)
ACF=ACF.dropna(subset=['dvv']).sort_values('date')
FAM=pd.read_csv(f'data/daily_dvv_{STA}_2to4_calT.csv')
TR=pd.read_csv(f'data/family_trust_tier1_{STA}.csv'); fmax=TR[TR.kind=='R'].coda_sigma.max(); TR=TR[TR.kind=='F']
gold=set(TR[TR.coda_sigma>fmax].fam); fail=set(TR[TR.verdict=='FAIL'].fam)
sig={r.fam:r.coda_sigma for r in TR.itertuples()}

def deseason(dates,y):
    doy=pd.to_datetime(dates).dt.dayofyear.values
    X=np.c_[np.ones_like(doy,float),np.cos(2*np.pi*doy/365.25),np.sin(2*np.pi*doy/365.25),
            np.cos(4*np.pi*doy/365.25),np.sin(4*np.pi*doy/365.25)]
    b,_,_,_=np.linalg.lstsq(X,y,rcond=None); return y-X@b

def olsslope(x,y):
    x=x-x.mean(); vx=np.dot(x,x)
    return np.dot(x,y-y.mean())/vx if vx>0 else np.nan

def block_boot_ci(x,y):
    n=len(x); nbl=max(n//BLK,2); starts=np.arange(0,n-BLK+1)
    sl=np.empty(NB)
    for k in range(NB):
        picks=rng.choice(starts, nbl, replace=True)
        idx=np.concatenate([np.arange(s,s+BLK) for s in picks])
        sl[k]=olsslope(x[idx],y[idx])
    return np.nanpercentile(sl,2.5), np.nanpercentile(sl,97.5)

rows=[]
for fam,g in FAM.groupby('patch'):
    g=g.copy(); g['date']=pd.to_datetime(g.date)
    m=pd.merge(g[['date','dvv']], ACF[['date','dvv']], on='date', suffixes=('_f','_a')).dropna().sort_values('date')
    if len(m)<120: continue
    fr=deseason(m.date,m.dvv_f.values); ar=deseason(m.date,m.dvv_a.values)
    sl=olsslope(ar,fr);
    pred=fr.mean()+sl*(ar-ar.mean()); r2=1-np.sum((fr-pred)**2)/np.sum((fr-fr.mean())**2)
    lo,hi=block_boot_ci(ar,fr)
    cls=('SHALLOW-COUPLED' if lo>0 else 'ANOMALOUS' if hi<0 else 'DECOUPLED')
    tier='GOLD' if fam in gold else ('FAIL' if fam in fail else 'other')
    rows.append(dict(fam=fam,tier=tier,coda_sigma=round(sig.get(fam,np.nan),1),n=len(m),
                     slope=round(sl,3),ci_lo=round(lo,3),ci_hi=round(hi,3),r2=round(r2,2),cls=cls))
D=pd.DataFrame(rows); D.to_csv(f'data/exp_t2a_{STA}_fixed.csv',index=False)
print(f"[{STA}] T2a site-filter (block-bootstrap CI, 30-day blocks): {len(D)} families\n")
print("== class x tier ==")
print(pd.crosstab(D.tier,D.cls).to_string())
print("\n== GOLD families that are SHALLOW-COUPLED (candidate site-carriers — exclude/flag for deep use) ==")
sc=D[(D.tier=='GOLD')&(D.cls=='SHALLOW-COUPLED')].sort_values('slope',ascending=False)
print(sc[['fam','coda_sigma','slope','ci_lo','ci_hi','r2']].to_string(index=False) if len(sc) else "  none")
print(f"\n== GOLD summary: SHALLOW-COUPLED {sum((D.tier=='GOLD')&(D.cls=='SHALLOW-COUPLED'))}, "
      f"DECOUPLED {sum((D.tier=='GOLD')&(D.cls=='DECOUPLED'))}, ANOMALOUS {sum((D.tier=='GOLD')&(D.cls=='ANOMALOUS'))} "
      f"of {sum(D.tier=='GOLD')} GOLD ==")
print(f"-> data/exp_t2a_{STA}_fixed.csv\nT2A_FIXED_DONE")
