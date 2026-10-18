#!/usr/bin/env python
"""Moving-patch inside-vs-outside test for any dv/v product npz (MODEL key on model days): per day, deep cells currently
active vs quiet in the same segment; active = model slip-in-window >= 3 mm (mode 'slip') or within 35 km of tremor in the
trailing 30 d (mode 'tremor'); statistic = mean dv/v 30-d change (lag L after the map date) inside minus outside, averaged
over days with >= 8 active and >= 8 quiet cells; rigid-shift null (100) and max-stat over lags 0..120 step 10.
Usage: python sse_moving_patch.py <product.npz> <label> [slip|tremor]"""
import numpy as np, pandas as pd, sys
SD="fault_tomography/inversion/res_catalog_sse"; RD="fault_tomography/inversion/res_catalog"; rng=np.random.RandomState(31)
prod, label, mode = sys.argv[1], sys.argv[2], (sys.argv[3] if len(sys.argv)>3 else "slip")
zc=np.load(f"{SD}/components_at_cells_trailing30.npz",allow_pickle=True); days=pd.to_datetime(zc["days"]); G=len(days); clat,clon,dep,pdist=zc["cell_lat"],zc["cell_lon"],zc["iface_depth_km"],zc["patch_dist_km"]; gi=pd.Series(np.arange(G),index=days)
km=lambda la,lo: np.column_stack([(lo+124)*111*np.cos(np.radians(47)),(la-47)*111]); XY=km(clat,clon); win=np.asarray((days>="2011-01-01")&(days<="2026-07-20")); month=days.month.values
if mode=="slip": ACT=zc["slip_thrust_inc30"].astype(float)*1e3; hi,lo=3.0,0.5
else:
    tr=pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv",usecols=["time","lat","lon"]); tr["t"]=pd.to_datetime(tr.time,errors="coerce"); tr=tr.dropna(subset=["t"]); tr["di"]=gi.reindex(tr.t.dt.floor("D")).values; tr=tr.dropna(subset=["di"]); tr["di"]=tr.di.astype(int)
    TXY=km(tr.lat.values,tr.lon.values); ACT=np.zeros((G,len(clat)),np.float32); di=tr.di.values
    from scipy.spatial import cKDTree
    tree=cKDTree(XY); near=tree.query_ball_point(TXY,35.0)
    for k,cs in enumerate(near):
        if cs: ACT[di[k],cs]+=1
    ACT=pd.DataFrame(ACT).rolling(30,min_periods=1).sum().values; hi,lo=20.0,0.0    # >=20 tremor detections within 35 km in 30 d = active; 0 = quiet
zp=np.load(prod,allow_pickle=True); Y=pd.DataFrame(np.asarray(zp["MODEL"],float),index=pd.to_datetime(zp["days"])).reindex(days).values.copy()
ep=pd.read_csv(f"{SD}/episodes.csv"); mask=np.zeros(G,bool)
for i in ep[ep.ok].ipeak: mask[max(int(i)-45,0):int(i)+91]=True
for m in range(1,13): sel=(month==m)&~mask; Y[month==m]-=np.nanmean(Y[sel],0)[None]
D=np.full_like(Y,np.nan); D[30:]=Y[30:]-Y[:-30]
deepc=np.isfinite(dep)&(dep>=30)&(dep<=45)&(pdist<=15); SEG={"N":(47.3,50),"C":(45,47.3),"S":(40,42.5)}; LAGS=np.arange(0,121,10)
print(f"MOVING PATCH [{mode}] on {label}: inside-minus-outside mean 30-d change; max-stat p over 13 lags")
for s,(a,b) in SEG.items():
    cs=deepc&(clat>=a)&(clat<b)
    def series(A,L):
        Dl=np.roll(D,-L,axis=0); out=[]
        for t in np.where(win)[0][::2]:
            sl=A[t,cs]; y=Dl[t,cs]; ins=(sl>=hi)&np.isfinite(y); outm=(sl<=lo)&np.isfinite(y)
            if ins.sum()>=8 and outm.sum()>=8: out.append(np.nanmean(y[ins])-np.nanmean(y[outm]))
        return np.nanmean(out) if len(out)>30 else np.nan
    obs=np.array([series(ACT,L) for L in LAGS]); nul=np.array([[series(np.roll(ACT,sh,axis=0),L) for L in LAGS] for sh in rng.randint(365,G-365,60)])
    z=(obs-np.nanmean(nul,0))/np.nanstd(nul,0); zmax=np.nanmax(np.abs((nul-np.nanmean(nul,0))/np.nanstd(nul,0)),1); k=int(np.nanargmax(np.abs(z))) if np.isfinite(z).any() else 0
    print(f"  {s}: extreme {obs[k]:+.3f}% at lag +{LAGS[k]} d (z {z[k]:+.1f}, p_lags {np.mean(zmax>=abs(z[k])):.3f}) | lags 0..120: "+" ".join(f"{v:+.2f}" for v in obs))
