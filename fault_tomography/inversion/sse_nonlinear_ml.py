#!/usr/bin/env python
"""Nonlinear multivariate test: per deep patch, gradient-boosted regression of the patch dv/v 30-d change (t) on ALL model
variables (12) at lags 0, +-30, +-60, +-90, +-180 d (patch means, trailing frame) + calendar month; 5-fold contiguous
time-block CV; out-of-sample R^2 vs a null where the target series is rigidly shifted (30 draws). Products: interface,
family 2-4 s, family 3-5 s. Output: res_catalog_sse/nonlinear_ml.csv"""
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
SD="fault_tomography/inversion/res_catalog_sse"; RD="fault_tomography/inversion/res_catalog"; rng=np.random.RandomState(9)
zc=np.load(f"{SD}/components_at_cells_trailing30.npz",allow_pickle=True); days=pd.to_datetime(zc["days"]); G=len(days); cells=zc["cell"]; dep=zc["iface_depth_km"]
zt=np.load(f"{RD}/transmission_T.npz",allow_pickle=True); PAT=np.asarray(zt["patterns"],float); plat=zt["patch_lat"]; band=zt["patch_band"]
VARS={"cum":zc["slip_thrust_cum"].astype(float)*1e3,"inc":zc["slip_thrust_inc30"].astype(float)*1e3,"dil":zc["dilatation"].astype(float)*1e8,"coul":zc["coulomb"].astype(float)/1e3,"shear":zc["shear_stress"].astype(float)/1e3,"normal":zc["normal_stress"].astype(float)/1e3}
for c in ("exx","eyy","ezz","exy","exz","eyz"): VARS[c]=zc[f"strain_{c}"].astype(float)*1e8
zp=np.load(f"{RD}/inversion_4d_daily.npz",allow_pickle=True); md=pd.to_datetime(zp["days"]); Y=pd.DataFrame(np.asarray(zp["MODEL"],float),index=md).reindex(days).values
F24=np.load(f"{SD}/product_family24.npz",allow_pickle=True)["MODEL"].astype(float); F35=np.load(f"{SD}/product_family35.npz",allow_pickle=True)["MODEL"].astype(float)
win=np.asarray((days>="2011-01-01")&(days<="2026-07-20")); month=days.month.values
def pm(X,q): w=PAT[q]; fin=np.isfinite(X); num=np.nansum(np.where(fin,X,0)*w[None],1); den=(fin*w[None]).sum(1); o=num/np.where(den>0,den,np.nan); o[(fin*w[None]).sum(1)/w.sum()<0.5]=np.nan; return o
def deseason(y): s=pd.Series(y,index=days); return (s-s.groupby(s.index.month).transform("mean")).values
def d30(y): o=np.full_like(y,np.nan); o[30:]=y[30:]-y[:-30]; return o
LAGS=(0,30,-30,60,-60,90,-90,180,-180); rows=[]
for q in [k for k in np.argsort(plat) if band[k]=="deep" and 40.5<plat[k]<49.5]:
    feats={f"{v}_{L:+d}":np.roll(pm(X,q),L) for v,X in VARS.items() for L in LAGS}; feats["month"]=month.astype(float)
    Xf=np.column_stack(list(feats.values()))
    for pname,P in (("interface",Y),("family 2-4 s",F24),("family 3-5 s",F35)):
        y=d30(deseason(pm(P,q))); m=win&np.isfinite(y)&np.isfinite(Xf).all(1); Xm,ym=Xf[m],y[m]; n=len(ym)
        if n<2000: continue
        def cv_r2(target):
            folds=np.array_split(np.arange(n),5); pred=np.full(n,np.nan)
            for f in folds:
                tr=np.setdiff1d(np.arange(n),np.concatenate([f,np.arange(max(f[0]-60,0),f[0]),np.arange(f[-1]+1,min(f[-1]+61,n))]))   # 60-d buffers
                mdl=HistGradientBoostingRegressor(max_iter=150,learning_rate=0.05,max_depth=4,min_samples_leaf=50,random_state=0).fit(Xm[tr],target[tr]); pred[f]=mdl.predict(Xm[f])
            return 1-np.nanmean((target-pred)**2)/np.nanvar(target)
        r2=cv_r2(ym); nul=[cv_r2(np.roll(ym,rng.randint(400,n-400))) for _ in range(20)]
        rows.append(dict(patch=round(plat[q],1),product=pname,n=n,cv_r2=r2,null_mean=np.mean(nul),null_sd=np.std(nul),p=float(np.mean(np.array(nul)>=r2)))); print(f"patch {plat[q]:.1f}N {pname:13s}: CV R^2 {r2:+.4f} | shifted-target null {np.mean(nul):+.4f} +- {np.std(nul):.4f} -> p {rows[-1]['p']:.2f}",flush=True)
T=pd.DataFrame(rows); T.to_csv(f"{SD}/nonlinear_ml.csv",index=False); print(f"cells with p<0.05: {(T.p<0.05).sum()}/{len(T)}; max CV R^2 {T.cv_r2.max():+.4f}"); print("DONE_ML")
