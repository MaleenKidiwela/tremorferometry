#!/usr/bin/env python
"""Synthetic recovery (resolution) test for the multi-window geometry.
Plant isolated delta-beta/beta spikes up and down the margin, forward-model the dv/v that every real
(cell,station,window) datum would see, add realistic noise, run the SAME joint fault+site inversion,
and compare recovered vs planted -> where can we actually resolve a fault change, and at what amplitude.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
sys.path.insert(0, 'fault_tomography/kernels')
from kernel import kernel_singlescatter

GRID=0.10; BETA=3.5; LSTAR=40.0; CCMIN=0.7
WINS=[('13',1.0,3.0,'1to3_cal_des'),('24',2.0,4.0,'2to4_cal_des'),('35',3.0,5.0,'3to5_cal_des')]
STA={'B927':(49.2188,-124.8113),'NLLB':(49.2271,-123.9882),'B928':(48.834,-125.134),'PGC':(48.6498,-123.4521),
 'B011':(48.65,-123.448),'B004':(48.202,-124.427),'B013':(47.813,-122.9108),'HDW':(47.6490,-123.0530),
 'GNW':(47.5641,-122.8250),'B014':(47.5133,-123.8125),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B020':(46.3827,-123.8445),'B201':(46.3033,-122.2648),'B204':(46.136,-122.169),'B023':(46.1112,-123.0787),
 'B022':(45.9546,-123.931),'B026':(45.3094,-123.8231),'COLT':(45.17044,-122.438152),'COR':(44.5855,-123.3046),
 'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),
 'B036':(42.5058,-123.3817),'B040':(41.8308,-122.4205),'B039':(41.4667,-122.4847),'B935':(40.4787,-123.5732)}

# ---- geometry (cells + multi-window G_f + per-(station,window) site) from the deseasoned calendar data ----
parts=[]
for s,(sla,slo) in STA.items():
    for wk,w1,w2,suf in WINS:
        f=f'data/daily_dvv_{s}_{suf}.csv'
        if not os.path.exists(f): continue
        d=pd.read_csv(f); d=d[d.cc_max>CCMIN]
        if not len(d): continue
        pre=d.patch.astype(str).str.split('__').str[0]
        plat=pre.str.split('_').str[0].astype(float); plon=pre.str.split('_',n=1).str[1].astype(float)
        cell=(np.round(plat/GRID)*GRID).round(3).astype(str)+'_'+(np.round(plon/GRID)*GRID).round(3).astype(str)
        parts.append(pd.DataFrame({'cell':cell,'clat':(np.round(plat/GRID)*GRID).round(3),'clon':(np.round(plon/GRID)*GRID).round(3),
                                   'station':s,'sla':sla,'slo':slo,'win':wk}).drop_duplicates())
cat=pd.concat(parts,ignore_index=True).drop_duplicates(['cell','station','win']).reset_index(drop=True)
cells=cat[['cell','clat','clon']].drop_duplicates('cell').reset_index(drop=True)
slab=pd.read_csv('data/cas_slab2_input_04-18.csv',low_memory=False)[['lat','lon','depth','etype']].dropna(subset=['depth'])
slab=slab[(~slab.etype.isin(['BA','TO']))&(slab.depth.between(10,70))]
q=cells[['clat','clon']].values
dl=griddata(slab[['lat','lon']].values,slab.depth.values,q,'linear'); dn=griddata(slab[['lat','lon']].values,slab.depth.values,q,'nearest')
cells['depth_km']=np.clip(np.where(np.isnan(dl),dn,dl),12,55)
cells['n_stations']=cat.groupby('cell')['station'].nunique().reindex(cells.cell).values
lat0,lon0=cells.clat.mean(),cells.clon.mean()
def proj(lat,lon): return (np.asarray(lon)-lon0)*111.0*np.cos(np.radians(lat0)),(np.asarray(lat)-lat0)*111.0
cells['x'],cells['y']=proj(cells.clat.values,cells.clon.values); M=len(cells); cidx={c:i for i,c in enumerate(cells.cell)}
Xk=cells[['x','y']].values; Zk=cells.depth_km.values
cat['depth']=cat.cell.map(cells.set_index('cell').depth_km)
cat['ax'],cat['ay']=proj(cat.clat.values,cat.clon.values); cat['bx'],cat['by']=proj(cat.sla.values,cat.slo.values)
WB={wk:(w1,w2) for wk,w1,w2,_ in WINS}
Gf=np.zeros((len(cat),M))
for i,r in cat.iterrows():
    w1,w2=WB[r.win]
    Gf[i]=kernel_singlescatter(Xk,Zk,np.array([r.ax,r.ay]),r.depth,np.array([r.bx,r.by]),BETA,w1,w2,ell=LSTAR)
Gf=-Gf
sw=(cat.station+'|'+cat.win); swu=sorted(sw.unique()); swidx={x:i for i,x in enumerate(swu)}; B=len(swu)
Gs=sparse.csr_matrix((np.ones(len(cat)),(np.arange(len(cat)),sw.map(swidx).values)),shape=(len(cat),B))
def graph_laplacian(P,k=6):
    tree=cKDTree(P); m=len(P); rr,cc,vv=[],[],[]
    for i in range(m):
        _,nn=tree.query(P[i],k+1)
        for j in nn[1:]: rr+=[i,i]; cc+=[i,j]; vv+=[1.0,-1.0]
    return sparse.coo_matrix((vv,(rr,cc)),shape=(m,m)).tocsr()
L=graph_laplacian(cells[['x','y']].values)
print(f'geometry: {M} cells, {len(cat)} data rows, {B} site terms')

# ---- plant spikes up the margin, forward + noise, invert ----
rng=np.random.RandomState(0)
targets=[(48.5,'N dense'),(47.5,'N-cent'),(46.2,'cent'),(44.6,'Oregon'),(43.0,'S Oregon'),(41.0,'S sparse')]
AMP=-0.3; NOISE=0.05
m_true=np.zeros(M); planted=[]
for tlat,name in targets:
    # nearest resolved-ish cell to (tlat, mid-band lon)
    band=cells[(cells.clat.between(tlat-0.06,tlat+0.06))]
    if not len(band): continue
    ci=band.iloc[(band.clat-tlat).abs().argmin()]; k=cidx[ci.cell]
    # 3-cell patch around it
    dd=np.hypot(cells.x-cells.x[k],cells.y-cells.y[k]); patch=np.where(dd<16)[0]
    m_true[patch]=AMP; planted.append((k,name,ci.clat,ci.clon,cells.n_stations[k]))
lam_f,lam_s=0.4,0.05
d=Gf@m_true + NOISE*rng.randn(len(cat))
A=sparse.vstack([sparse.hstack([sparse.csr_matrix(Gf),Gs]),
                 sparse.hstack([lam_f*L,sparse.csr_matrix((M,B))]),
                 sparse.hstack([sparse.csr_matrix((B,M)),lam_s*sparse.identity(B)])]).tocsr()
b=np.concatenate([d,np.zeros(M),np.zeros(B)])
sol=lsqr(A,b,atol=1e-7,btol=1e-7,iter_lim=3000)[0]; m_rec=sol[:M]
print('\\nplanted -0.30%% patches  ->  recovered (peak in patch), and amplitude recovery:')
for k,name,cla,clo,nst in planted:
    dd=np.hypot(cells.x-cells.x[k],cells.y-cells.y[k]); patch=dd<16
    rec=m_rec[patch][np.argmax(np.abs(m_rec[patch]))]
    print(f'  {name:9s} {cla:.1f}N  ({int(nst)} sta)  recovered {rec:+.3f}%%  ({100*rec/AMP:.0f}%% of planted)')

# ---- figure: planted vs recovered on the margin ----
fig=plt.figure(figsize=(13,11))
for p,(arr,ttl) in enumerate([(m_true,'PLANTED (truth)'),(m_rec,'RECOVERED (inversion)')]):
    ax=fig.add_subplot(1,2,p+1,projection=ccrs.PlateCarree()); ax.set_extent([-126.2,-121.0,40.0,49.8])
    ax.add_feature(cfeature.OCEAN.with_scale('50m'),facecolor='#eaf1f6'); ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='#f6f4ef')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'),lw=0.6); ax.add_feature(cfeature.BORDERS.with_scale('50m'),lw=0.4)
    sc=ax.scatter(cells.clon,cells.clat,c=arr,cmap='RdBu_r',vmin=-0.3,vmax=0.3,s=30,marker='s',transform=ccrs.PlateCarree())
    for k,name,cla,clo,nst in planted: ax.text(clo+0.15,cla,'%dst'%int(nst),fontsize=7,transform=ccrs.PlateCarree())
    ax.set_title(ttl,fontsize=12)
fig.colorbar(sc,ax=fig.axes,shrink=0.4,label=u'δβ/β (%)')
fig.suptitle('Synthetic spike-recovery test (multi-window): can we resolve a -0.3%% fault change?',y=0.95,fontsize=13)
plt.savefig('figures/smoke_synth_resolution.png',dpi=120,bbox_inches='tight'); print('\\nsaved figures/smoke_synth_resolution.png')
