#!/usr/bin/env python
"""Coarse 3D VOLUME monthly inversion (multi-window 1-3/2-4/3-5 s).
Model = 3D crustal cells (lon,lat,depth) at the RESOLVABLE coarse scale (~0.4 deg horizontal, 3 depth
zones: shallow/mid/deep) -- validated by volume3d_checkerboard (coarse corr 0.86, depth separates).
No lumped site terms: the shallow layer IS the near-surface, now resolved spatially. Each datum's row =
the single-scatter ellipsoid kernel over the 3D cells (the 3 windows give the depth leverage). Monthly LSQR
on the deseasoned calendar data -> 3D delta-beta/beta volume per month -> depth-slice movie with coastlines.
RESOLUTION IS COARSE (~150-200 km horizontal, 3 depth zones); fine detail is NOT resolved.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import cartopy.crs as ccrs, cartopy.feature as cfeature
sys.path.insert(0, 'fault_tomography/kernels')
from kernel import kernel_singlescatter

BETA=3.5; LSTAR=40.0; CCMIN=0.7; HGRID=0.4
DEPTHS=np.array([10.,26.,42.]); DLAB=['shallow ~10 km','mid ~26 km','deep ~42 km']
WINS=[('13',1.0,3.0,'1to3_cal_des'),('24',2.0,4.0,'2to4_cal_des'),('35',3.0,5.0,'3to5_cal_des')]
STA={'B927':(49.2188,-124.8113),'NLLB':(49.2271,-123.9882),'B928':(48.834,-125.134),'PGC':(48.6498,-123.4521),
 'B011':(48.65,-123.448),'B004':(48.202,-124.427),'B013':(47.813,-122.9108),'HDW':(47.6490,-123.0530),
 'GNW':(47.5641,-122.8250),'B014':(47.5133,-123.8125),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B020':(46.3827,-123.8445),'B201':(46.3033,-122.2648),'B204':(46.136,-122.169),'B023':(46.1112,-123.0787),
 'B022':(45.9546,-123.931),'B026':(45.3094,-123.8231),'COLT':(45.17044,-122.438152),'COR':(44.5855,-123.3046),
 'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),
 'B036':(42.5058,-123.3817),'B040':(41.8308,-122.4205),'B039':(41.4667,-122.4847),'B935':(40.4787,-123.5732)}

# ---- monthly tensor (source cell, station, window, month) ----
parts=[]
for s,(sla,slo) in STA.items():
    for wk,w1,w2,suf in WINS:
        f=f'data/daily_dvv_{s}_{suf}.csv'
        if not os.path.exists(f): continue
        d=pd.read_csv(f); d=d[d.cc_max>CCMIN]
        if not len(d): continue
        pre=d.patch.astype(str).str.split('__').str[0]
        d=d.assign(plat=(np.round(pre.str.split('_').str[0].astype(float)/0.1)*0.1).round(3),
                   plon=(np.round(pre.str.split('_',n=1).str[1].astype(float)/0.1)*0.1).round(3))
        d['src']=d.plat.astype(str)+'_'+d.plon.astype(str); d['ym']=pd.to_datetime(d.date).dt.to_period('M').astype(str)
        d['station']=s; d['sla']=sla; d['slo']=slo; d['win']=wk
        parts.append(d.groupby(['src','station','win','ym']).agg(dvv=('dvv','mean'),plat=('plat','first'),plon=('plon','first'),
                     sla=('sla','first'),slo=('slo','first')).reset_index())
ten=pd.concat(parts,ignore_index=True)
cat=ten[['src','station','win','plat','plon','sla','slo']].drop_duplicates(['src','station','win']).reset_index(drop=True)
slab=pd.read_csv('data/cas_slab2_input_04-18.csv',low_memory=False)[['lat','lon','depth','etype']].dropna(subset=['depth'])
slab=slab[(~slab.etype.isin(['BA','TO']))&(slab.depth.between(10,70))]
sd=griddata(slab[['lat','lon']].values,slab.depth.values,cat[['plat','plon']].values,'linear')
sd=np.where(np.isnan(sd),griddata(slab[['lat','lon']].values,slab.depth.values,cat[['plat','plon']].values,'nearest'),sd)
cat['sdep']=np.clip(sd,12,55)
print(f'{len(ten):,} (src,sta,win,month) rows; {len(cat)} unique (src,sta,win) pairs')

# ---- 3D coarse model + G ----
lat0,lon0=cat.plat.mean(),cat.plon.mean()
def proj(lat,lon): return (np.asarray(lon)-lon0)*111.0*np.cos(np.radians(lat0)),(np.asarray(lat)-lat0)*111.0
cat['ax'],cat['ay']=proj(cat.plat.values,cat.plon.values); cat['bx'],cat['by']=proj(cat.sla.values,cat.slo.values)
lons=np.arange(cat.plon.min()-0.4,cat.plon.max()+0.6,HGRID); lats=np.arange(cat.plat.min()-0.4,cat.plat.max()+0.4,HGRID)
GL,GA,GD=np.meshgrid(lons,lats,DEPTHS,indexing='ij'); cx,cy=proj(GA.ravel(),GL.ravel()); cz=GD.ravel()
M0=len(cz); WB={wk:(w1,w2) for wk,w1,w2,_ in WINS}; Xk=np.column_stack([cx,cy])
G=np.zeros((len(cat),M0))
for i,r in cat.iterrows():
    w1,w2=WB[r.win]; G[i]=kernel_singlescatter(Xk,cz,np.array([r.ax,r.ay]),r.sdep,np.array([r.bx,r.by]),BETA,w1,w2,ell=LSTAR)
keep=(G>0).sum(0)>=3; G=-G[:,keep]
GLk,GAk,GDk,cxk,cyk,czk=GL.ravel()[keep],GA.ravel()[keep],GD.ravel()[keep],cx[keep],cy[keep],cz[keep]
M=int(keep.sum()); print(f'3D model: {M0}->{M} illuminated cells, {len(DEPTHS)} depth layers')
P=np.column_stack([cxk,cyk,czk*3.0]); tree=cKDTree(P); rr,cc,vv=[],[],[]
for i in range(M):
    _,nn=tree.query(P[i],7)
    for j in nn[1:]: rr+=[i,i]; cc+=[i,j]; vv+=[1.,-1.]
L=sparse.coo_matrix((vv,(rr,cc)),shape=(M,M)).tocsr()
key=cat[['src','station','win']].copy(); key['row']=np.arange(len(cat))

# ---- monthly inversion ----
lam=0.8; reg=lam*L
months=sorted(ten.ym.unique()); MV=np.full((M,len(months)),np.nan); ND=np.zeros(len(months),int)
for j,E in enumerate(months):
    te=ten[ten.ym==E][['src','station','win','dvv']]; mm=key.merge(te,on=['src','station','win'])
    if len(mm)<60 or mm.station.nunique()<10: continue
    rows=mm.row.values; dd=mm.dvv.values*100
    A=sparse.vstack([sparse.csr_matrix(G[rows]),reg]).tocsr(); b=np.concatenate([dd,np.zeros(M)])
    MV[:,j]=lsqr(A,b,atol=1e-7,btol=1e-7,iter_lim=2000)[0]; ND[j]=len(mm)
ok=ND>0; print(f'inverted {ok.sum()}/{len(months)} months')
np.savez('fault_tomography/inversion/fault_volume3d.npz',months=np.array(months),lon=GLk,lat=GAk,depth=GDk,MV=MV,ok=ok,DEPTHS=DEPTHS)

# ---- depth-slice movie (3-month smoothed) ----
df=pd.DataFrame(MV.T,index=pd.PeriodIndex(months,freq='M')); MVs=df.rolling(3,center=True,min_periods=2).mean().values.T
idx=np.where(ok)[0]; slat=[v[0] for v in STA.values()]; slon=[v[1] for v in STA.values()]
fig=plt.figure(figsize=(15,9)); axs=[fig.add_subplot(1,3,k+1,projection=ccrs.PlateCarree()) for k in range(3)]
sm={}
for k,(ax,dep,lab) in enumerate(zip(axs,DEPTHS,DLAB)):
    ax.set_extent([-126.2,-121,40,49.8]); ax.add_feature(cfeature.OCEAN.with_scale('50m'),facecolor='#eaf1f6'); ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='#f6f4ef')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'),lw=.6); ax.add_feature(cfeature.BORDERS.with_scale('50m'),lw=.4)
    m=np.abs(GDk-dep)<1
    sm[k]=ax.scatter(GLk[m],GAk[m],c=MVs[m,idx[0]],cmap='RdBu_r',vmin=-0.2,vmax=0.2,s=70,marker='s',transform=ccrs.PlateCarree())
    ax.plot(slon,slat,'k^',ms=4,transform=ccrs.PlateCarree()); ax.set_title(lab,fontsize=12)
cb=fig.colorbar(sm[0],ax=axs,shrink=0.5,label=u'δβ/β (%)')
sup=fig.suptitle('',y=0.86,fontsize=14)
def upd(t):
    j=idx[t]
    for k,dep in enumerate(DEPTHS):
        m=np.abs(GDk-dep)<1; sm[k].set_array(MVs[m,j])
    sup.set_text('Coarse 3D interface velocity change (3 depth zones, ~150 km resolution)  '+months[j]); return list(sm.values())+[sup]
anim=FuncAnimation(fig,upd,frames=len(idx),interval=130,blit=False)
anim.save('figures/fault_movie_volume3d.mp4',writer=FFMpegWriter(fps=7,bitrate=2400))
print('saved figures/fault_movie_volume3d.mp4 (%d frames)'%len(idx))
