#!/usr/bin/env python
"""A2 (2-4 s DESEASONED): same tensor build as a2_all_tensor.py but ingests the corrected
deseasoned 2-4 s rolling-stack dv/v (data/daily_dvv_{S}_2to4_des.csv) for all 28 stations.
Outputs prefix tomo24des_* (does NOT touch the 1-4 s tomo/borehole baselines)."""
import os
import numpy as np, pandas as pd
from scipy.interpolate import griddata

GRID = 0.10; OUT = 'fault_tomography/data_assembly'

STA = {  # station: (lat, lon, tier)
 'B011':(48.65,-123.448,'A'),'B013':(47.813,-122.9108,'A'),'B941':(46.9868,-122.219,'A'),
 'B018':(46.9795,-123.0203,'A'),'B023':(46.1112,-123.0787,'A'),'B928':(48.834,-125.134,'A'),
 'B004':(48.202,-124.427,'A'),'B026':(45.3094,-123.8231,'A'),'B014':(47.5133,-123.8125,'A'),
 'B204':(46.136,-122.169,'A'),'B028':(44.4937,-122.9638,'A'),'B030':(43.9713,-122.7717,'A'),
 'B032':(43.668,-123.3923,'A'),'B033':(43.2917,-123.1245,'A'),'B036':(42.5058,-123.3817,'A'),
 'B039':(41.4667,-122.4847,'A'),'B040':(41.8308,-122.4205,'A'),'B927':(49.2188,-124.8113,'A'),
 'B020':(46.3827,-123.8445,'A'),'B022':(45.9546,-123.931,'A'),'B935':(40.4787,-123.5732,'A'),
 'B201':(46.3033,-122.2648,'A'),
 'GNW':(47.5641,-122.8250,'A'),'HDW':(47.6490,-123.0530,'A'),'NLLB':(49.2271,-123.9882,'A'),
 'COR':(44.5855,-123.3046,'A'),'PGC':(48.6498,-123.4521,'B'),'COLT':(45.17044,-122.438152,'A'),
}

rows = []
for s,(sla,slo,tier) in STA.items():
    f = f'data/daily_dvv_{s}_2to4_des.csv'
    if not os.path.exists(f): print('  MISSING', f); continue
    d = pd.read_csv(f); pre = d['patch'].astype(str).str.split('__').str[0]
    d['plat']=pre.str.split('_').str[0].astype(float); d['plon']=pre.str.split('_',n=1).str[1].astype(float)
    d['station']=s; d['sta_lat']=sla; d['sta_lon']=slo; d['tier']=tier
    if 'n_det' not in d: d['n_det']=1
    rows.append(d[['station','sta_lat','sta_lon','tier','plat','plon','patch','date','dvv','cc_max','n_det']])
df = pd.concat(rows, ignore_index=True); df['date']=pd.to_datetime(df['date'])
df['cell_lat']=(np.round(df.plat/GRID)*GRID).round(3); df['cell_lon']=(np.round(df.plon/GRID)*GRID).round(3)
df['cell']=df.cell_lat.astype(str)+'_'+df.cell_lon.astype(str)
print(f'loaded {len(df):,} rows across {df.station.nunique()} stations (2-4s deseasoned)')

slab=pd.read_csv('data/cas_slab2_input_04-18.csv',low_memory=False)[['lat','lon','depth','etype']].dropna(subset=['depth'])
slab=slab[(~slab.etype.isin(['BA','TO']))&(slab.depth.between(10,70))]
cells=df[['cell','cell_lat','cell_lon']].drop_duplicates('cell').reset_index(drop=True)
q=cells[['cell_lat','cell_lon']].values
dl=griddata(slab[['lat','lon']].values,slab.depth.values,q,'linear'); dn=griddata(slab[['lat','lon']].values,slab.depth.values,q,'nearest')
cells['depth_km']=np.clip(np.where(np.isnan(dl),dn,dl),12,55)

cat=(df.groupby(['cell','station']).agg(tier=('tier','first'),n_families=('patch','nunique'),n_obs=('dvv','size'),
        sta_lat=('sta_lat','first'),sta_lon=('sta_lon','first'),date_min=('date','min'),date_max=('date','max')).reset_index())
cat=cat.merge(cells,on='cell')
cells=cells.merge(cat.groupby('cell')['station'].nunique().rename('n_stations'),on='cell')

df['ym']=df.date.dt.to_period('M').astype(str); df['_wd']=df.dvv*df.n_det.clip(lower=1)
g=df.groupby(['cell','station','ym'])
ten=g.agg(sta_lat=('sta_lat','first'),sta_lon=('sta_lon','first'),cell_lat=('cell_lat','first'),
          cell_lon=('cell_lon','first'),tier=('tier','first'),n=('dvv','size')).reset_index()
ten['dvv']=(g['_wd'].sum()/g.apply(lambda x:x.n_det.clip(lower=1).sum())).values
ten=ten.merge(cells[['cell','depth_km','n_stations']],on='cell')

cat.to_csv(f'{OUT}/tomo24des_patch_catalog.csv',index=False); cells.to_csv(f'{OUT}/tomo24des_cells.csv',index=False)
ten.to_parquet(f'{OUT}/tomo24des_dvv_monthly.parquet',index=False)
n3=(cells.n_stations>=3).sum(); n4=(cells.n_stations>=4).sum()
print(f'fault cells: {len(cells)} | >=3 stations: {n3} | >=4: {n4}')
print(f'monthly tensor: {len(ten):,} (cell,station,month); {ten.ym.nunique()} months {ten.ym.min()}..{ten.ym.max()}')
print(f'wrote {OUT}/tomo24des_*.csv/parquet')
