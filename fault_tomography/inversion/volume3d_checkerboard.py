#!/usr/bin/env python
"""3D VOLUME resolution test for the multi-window geometry (1-3/2-4/3-5 s).
Model = a 3D grid of crustal cells (lon,lat,depth), NOT just the interface. For each datum
(source family on interface, station at surface, window) the row of G is the single-scatter ellipsoid
kernel evaluated over the 3D volume. Plant a 3D checkerboard, forward + noise, invert with 3D smoothness,
recover -> report recovery corr BY DEPTH LAYER (which depths the geometry can actually resolve).
This decides whether a 3D volume movie is meaningful before we build it.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, 'fault_tomography/kernels')
from kernel import kernel_singlescatter

BETA=3.5; LSTAR=40.0; CCMIN=0.7
WINS=[('13',1.0,3.0,'1to3_cal_des'),('24',2.0,4.0,'2to4_cal_des'),('35',3.0,5.0,'3to5_cal_des')]
STA={'B927':(49.2188,-124.8113),'NLLB':(49.2271,-123.9882),'B928':(48.834,-125.134),'PGC':(48.6498,-123.4521),
 'B011':(48.65,-123.448),'B004':(48.202,-124.427),'B013':(47.813,-122.9108),'HDW':(47.6490,-123.0530),
 'GNW':(47.5641,-122.8250),'B014':(47.5133,-123.8125),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B020':(46.3827,-123.8445),'B201':(46.3033,-122.2648),'B204':(46.136,-122.169),'B023':(46.1112,-123.0787),
 'B022':(45.9546,-123.931),'B026':(45.3094,-123.8231),'COLT':(45.17044,-122.438152),'COR':(44.5855,-123.3046),
 'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),
 'B036':(42.5058,-123.3817),'B040':(41.8308,-122.4205),'B039':(41.4667,-122.4847),'B935':(40.4787,-123.5732)}

# ---- data pairs (source family-cell on interface, station, window) ----
parts=[]
for s,(sla,slo) in STA.items():
    for wk,w1,w2,suf in WINS:
        f=f'data/daily_dvv_{s}_{suf}.csv'
        if not os.path.exists(f): continue
        d=pd.read_csv(f,usecols=['patch','cc_max']); d=d[d.cc_max>CCMIN]
        if not len(d): continue
        pre=d.patch.astype(str).str.split('__').str[0]
        plat=pre.str.split('_').str[0].astype(float); plon=pre.str.split('_',n=1).str[1].astype(float)
        sub=pd.DataFrame({'plat':(np.round(plat/0.1)*0.1).round(3),'plon':(np.round(plon/0.1)*0.1).round(3)}).drop_duplicates()
        sub['station']=s; sub['sla']=sla; sub['slo']=slo; sub['win']=wk; parts.append(sub)
cat=pd.concat(parts,ignore_index=True).drop_duplicates(['plat','plon','station','win']).reset_index(drop=True)
slab=pd.read_csv('data/cas_slab2_input_04-18.csv',low_memory=False)[['lat','lon','depth','etype']].dropna(subset=['depth'])
slab=slab[(~slab.etype.isin(['BA','TO']))&(slab.depth.between(10,70))]
cat['sdep']=np.clip(np.where(np.isnan(griddata(slab[['lat','lon']].values,slab.depth.values,cat[['plat','plon']].values,'linear')),
    griddata(slab[['lat','lon']].values,slab.depth.values,cat[['plat','plon']].values,'nearest'),
    griddata(slab[['lat','lon']].values,slab.depth.values,cat[['plat','plon']].values,'linear')),12,55)
lat0,lon0=cat.plat.mean(),cat.plon.mean()
def proj(lat,lon): return (np.asarray(lon)-lon0)*111.0*np.cos(np.radians(lat0)),(np.asarray(lat)-lat0)*111.0
cat['ax'],cat['ay']=proj(cat.plat.values,cat.plon.values); cat['bx'],cat['by']=proj(cat.sla.values,cat.slo.values)
print(f'{len(cat)} data rows (source,station,window)')

# ---- 3D model grid (lon,lat,depth) over the data footprint ----
HGRID=0.2; DEPTHS=np.array([3.,12.,22.,32.,40.,48.])           # km, layer centers (spans both lobes + the trough)
lons=np.arange(cat.plon.min()-0.3,cat.plon.max()+0.5,HGRID); lats=np.arange(cat.plat.min()-0.3,cat.plat.max()+0.3,HGRID)
GL,GA,GD=np.meshgrid(lons,lats,DEPTHS,indexing='ij')
cx,cy=proj(GA.ravel(),GL.ravel()); cz=GD.ravel()
M0=len(cz); WB={wk:(w1,w2) for wk,w1,w2,_ in WINS}
# build G over all cells, keep illuminated ones
G=np.zeros((len(cat),M0))
Xk=np.column_stack([cx,cy])
for i,r in cat.iterrows():
    w1,w2=WB[r.win]
    G[i]=kernel_singlescatter(Xk,cz,np.array([r.ax,r.ay]),r.sdep,np.array([r.bx,r.by]),BETA,w1,w2,ell=LSTAR)
hit=(G>0).sum(0)
keep=hit>=3                                                    # cells illuminated by >=3 data
G=-G[:,keep]; cz_k=cz[keep]; GLk=GL.ravel()[keep]; GAk=GA.ravel()[keep]; GDk=GD.ravel()[keep]; cxk=cx[keep]; cyk=cy[keep]
M=int(keep.sum())
print(f'3D model: {M0} cells -> {M} illuminated (>=3 data); depth layers {DEPTHS.tolist()} km')
# 3D Laplacian
P=np.column_stack([cxk,cyk,cz_k*3.0]); tree=cKDTree(P); rr,cc,vv=[],[],[]
for i in range(M):
    _,nn=tree.query(P[i],7)
    for j in nn[1:]: rr+=[i,i]; cc+=[i,j]; vv+=[1.,-1.]
L=sparse.coo_matrix((vv,(rr,cc)),shape=(M,M)).tocsr()
# checkerboard at SEVERAL scales (reuse G,L) to find where it starts to resolve
rng=np.random.RandomState(0)
SCALES=[(0.5,18,'fine 0.5deg/18km'),(1.0,30,'med 1.0deg/30km'),(1.8,60,'coarse 1.8deg/60km'),(3.0,200,'shallow-vs-deep zones')]
print('\\n=== 3D checkerboard recovery vs SCALE (corr by depth) ===')
results={}
for hdeg,dkm,nm in SCALES:
    mt=np.sign(np.sin(np.pi*GAk/hdeg)*np.sin(np.pi*GLk/hdeg)*np.sin(np.pi*GDk/dkm))*0.01
    d=G@mt + 0.05*np.std(G@mt)*rng.randn(len(cat))
    A=sparse.vstack([sparse.csr_matrix(G),0.5*L]); b=np.concatenate([d,np.zeros(M)])
    mr=lsqr(A,b,atol=1e-6,btol=1e-6,iter_lim=3000)[0]
    overall=np.corrcoef(mt,mr)[0,1]
    bydep=[(np.corrcoef(mt[np.abs(cz_k-dep)<1],mr[np.abs(cz_k-dep)<1])[0,1] if (np.abs(cz_k-dep)<1).sum()>3 else np.nan) for dep in DEPTHS]
    results[nm]=(mt,mr); print('  %-24s overall %.2f | depths %s'%(nm,overall,['%.2f'%c for c in bydep]))
# --- DEPTH-SEPARATION test (the real 3D question): shallow + / deep -, COARSE horizontal ---
env=np.sign(np.sin(np.pi*GAk/2.5)*np.sin(np.pi*GLk/2.5))
mt=env*np.where(cz_k<18,0.01,np.where(cz_k>34,-0.01,0.0))
d=G@mt + 0.05*np.std(G@mt)*rng.randn(len(cat))
mr=lsqr(sparse.vstack([sparse.csr_matrix(G),0.5*L]).tocsr(),np.concatenate([d,np.zeros(M)]),atol=1e-6,btol=1e-6,iter_lim=3000)[0]
print('\\n=== DEPTH-SEPARATION test (planted: shallow +0.01 / deep -0.01, coarse horizontal) ===')
print('  if recovered flips sign shallow->deep, depth IS resolved; if same sign throughout, it smears:')
for dep in DEPTHS:
    sl=np.abs(cz_k-dep)<1
    print('    %2.0f km: planted %+.4f  recovered %+.4f'%(dep,(mt[sl]*env[sl]).mean(),(mr[sl]*env[sl]).mean()))
# plot the coarsest (zones) recovery by depth
mt,mr=results['shallow-vs-deep zones']
fig,axs=plt.subplots(2,len(DEPTHS),figsize=(3.1*len(DEPTHS),7),sharex=True,sharey=True)
for di,dep in enumerate(DEPTHS):
    sl=np.abs(cz_k-dep)<1; corr=np.corrcoef(mt[sl],mr[sl])[0,1] if sl.sum()>3 else np.nan
    for row,(arr,nm2) in enumerate([(mt,'true'),(mr,'rec')]):
        ax=axs[row,di]; ax.scatter(GLk[sl],GAk[sl],c=arr[sl],cmap='RdBu_r',vmin=-0.01,vmax=0.01,s=14,marker='s')
        if row==0: ax.set_title('%.0f km  corr %.2f'%(dep,corr),fontsize=10)
        if di==0: ax.set_ylabel(nm2,fontsize=11)
fig.suptitle('3D resolution, COARSEST test (shallow-vs-deep zones): top=planted, bottom=recovered',y=0.98,fontsize=13)
plt.tight_layout(); plt.savefig('figures/smoke_volume3d_checkerboard.png',dpi=110,bbox_inches='tight'); print('\\nsaved (coarsest-scale figure)')
