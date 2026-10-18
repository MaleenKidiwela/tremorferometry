#!/usr/bin/env python
"""(1) CROSS-SEGMENT response matrix: regional dv/v 30-d change in segment X vs slip-in-window in segment Y, lags -365..+365
(15-d), max-statistic over lags, shift null (300). (2) LONG-LAG level scan: patch-mean dv/v vs cumulative slip and slip-in-
window at lags -730..+730 d. (3) VARIANCE test: per region, the day-to-day scatter (robust SD of daily residuals about a
30-d running median) of dv/v during slip windows (segment slip-in-window >= 3 mm) vs quiet, ratio vs shift null.
Products: interface (deseasoned), family 2-4 s, family 3-5 s.  Output: res_catalog_sse/crosseg_longlag_var.csv"""
import numpy as np, pandas as pd
SD="fault_tomography/inversion/res_catalog_sse"; RD="fault_tomography/inversion/res_catalog"; rng=np.random.RandomState(7)
zc=np.load(f"{SD}/components_at_cells_trailing30.npz",allow_pickle=True); days=pd.to_datetime(zc["days"]); G=len(days); cells=zc["cell"]; clat,dep,pdist=zc["cell_lat"],zc["iface_depth_km"],zc["patch_dist_km"]
inc=zc["slip_thrust_inc30"].astype(float)*1e3; cum=zc["slip_thrust_cum"].astype(float)*1e3; win=np.asarray((days>="2011-01-01")&(days<="2026-07-20")); month=days.month.values
zp=np.load(f"{RD}/inversion_4d_daily.npz",allow_pickle=True); md=pd.to_datetime(zp["days"]); Y=pd.DataFrame(np.asarray(zp["MODEL"],float),index=md).reindex(days).values
F24=np.load(f"{SD}/product_family24.npz",allow_pickle=True)["MODEL"].astype(float); F35=np.load(f"{SD}/product_family35.npz",allow_pickle=True)["MODEL"].astype(float)
ep=pd.read_csv(f"{SD}/episodes.csv"); mask=np.zeros(G,bool)
for i in ep[ep.ok].ipeak: mask[max(int(i)-45,0):int(i)+91]=True
deepc=np.isfinite(dep)&(dep>=30)&(dep<=45)&(pdist<=15); SEG={"S":(40,42.5),"SC":(42.5,45),"C":(45,47.3),"N":(47.3,50)}; REG={s:np.where(deepc&(clat>=a)&(clat<b))[0] for s,(a,b) in SEG.items()}
def rmean(X,cs): v=X[:,cs]; o=np.nanmean(v,1); o[np.isfinite(v).mean(1)<0.5]=np.nan; return o
def deseason(y):
    o=y.copy()
    for m in range(1,13): s=(month==m)&~mask&np.isfinite(y); o[month==m]-=np.nanmean(y[s]) if s.any() else 0
    return o
def d30(y): o=np.full_like(y,np.nan); o[30:]=y[30:]-y[:-30]; return o
def corr(a,b): m=np.isfinite(a)&np.isfinite(b)&win; return np.corrcoef(a[m],b[m])[0,1] if m.sum()>500 else np.nan
PROD={"interface":Y,"family 2-4 s":F24,"family 3-5 s":F35}; DV={p:{s:deseason(rmean(X,cs)) for s,cs in REG.items()} for p,X in PROD.items()}; SL={s:rmean(inc,cs) for s,cs in REG.items()}; CU={s:rmean(cum,cs) for s,cs in REG.items()}
rows=[]
LAG1=np.arange(-365,366,15)
print("(1) CROSS-SEGMENT: dv/v 30-d change in X vs slip-in-window in Y (best lag within +-365 d; p = max-stat over lags)")
for p in PROD:
    for X in SEG:
        for Yseg in SEG:
            y=d30(DV[p][X]); x=SL[Yseg]; xl=np.array([corr(y,np.roll(x,L)) for L in LAG1]); k=int(np.nanargmax(np.abs(xl))) if np.isfinite(xl).any() else 0
            nul=[np.nanmax(np.abs([corr(y,np.roll(np.roll(x,rng.randint(365,G-365)),L)) for L in LAG1[::3]])) for _ in range(150)]
            rows.append(dict(test="cross-segment",product=p,region=X,driver=Yseg,r=xl[k],lag=int(LAG1[k]),p=float(np.mean(np.array(nul)>=abs(xl[k])))))
    tb=pd.DataFrame([r for r in rows if r["test"]=="cross-segment" and r["product"]==p]).pivot(index="region",columns="driver",values="r").round(2); tp=pd.DataFrame([r for r in rows if r["test"]=="cross-segment" and r["product"]==p]).pivot(index="region",columns="driver",values="p").round(2)
    print(f"  {p}: r matrix (rows = dv/v region, cols = slip driver)\n{tb.to_string()}\n  p matrix\n{tp.to_string()}",flush=True)
print("\n(2) LONG LAGS (-730..+730 d, 30-d steps): patch/segment dv/v LEVEL vs cumulative slip and slip-in-window")
LAG2=np.arange(-730,731,30)
for p in PROD:
    for s in SEG:
        y=DV[p][s]
        for vn,x in (("cum slip",CU[s]),("slip in window",SL[s])):
            xl=np.array([corr(y,np.roll(x,L)) for L in LAG2]); k=int(np.nanargmax(np.abs(xl))) if np.isfinite(xl).any() else 0
            nul=[np.nanmax(np.abs([corr(y,np.roll(np.roll(x,rng.randint(800,G-800)),L)) for L in LAG2[::3]])) for _ in range(150)]
            rows.append(dict(test="long-lag",product=p,region=s,driver=vn,r=xl[k],lag=int(LAG2[k]),p=float(np.mean(np.array(nul)>=abs(xl[k])))))
            print(f"  {p:13s} {s:2s} vs {vn:15s}: best r {xl[k]:+.2f} at {LAG2[k]:+d} d, p {rows[-1]['p']:.2f}",flush=True)
print("\n(3) VARIANCE during slip vs quiet (robust SD of daily residual about a 30-d running median), ratio vs shift null")
for p in PROD:
    for s in SEG:
        y=DV[p][s]; res=y-pd.Series(y).rolling(30,center=True,min_periods=10).median().values; act=win&(SL[s]>=3); qui=win&(SL[s]<1)
        mad=lambda v: 1.4826*np.nanmedian(np.abs(v-np.nanmedian(v)))
        ratio=mad(res[act])/mad(res[qui]); nul=[]
        for _ in range(300):
            sl=np.roll(SL[s],rng.randint(365,G-365)); nul.append(mad(res[win&(sl>=3)])/mad(res[win&(sl<1)]))
        nul=np.array(nul); rows.append(dict(test="variance",product=p,region=s,driver="slip>=3mm vs <1mm",r=ratio,lag=0,p=float(np.mean(np.abs(nul-1)>=abs(ratio-1)))))
        print(f"  {p:13s} {s:2s}: scatter ratio active/quiet = {ratio:.3f} (null {nul.mean():.3f} +- {nul.std():.3f}) p {rows[-1]['p']:.2f}",flush=True)
pd.DataFrame(rows).to_csv(f"{SD}/crosseg_longlag_var.csv",index=False); print("DONE_CLV")
