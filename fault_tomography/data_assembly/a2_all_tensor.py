#!/usr/bin/env python
"""A2 (ALL tomography stations): same as a2_borehole_tensor.py but adds the 6 non-borehole stations,
each with its FINALIZED / clean-subset dv/v file (GNW & PGC restricted to their clean modern eras).
Outputs tomo_cells.csv, tomo_patch_catalog.csv, tomo_dvv_monthly.parquet (the full-roster tensor).
CPW excluded (co-located with B018). PGC/COLT are lower-tier (noisy / accelerometer)."""
import os
import numpy as np, pandas as pd
from scipy.interpolate import griddata

GRID = 0.10; OUT = 'fault_tomography/data_assembly'; os.makedirs(OUT, exist_ok=True)

# station -> (lat, lon, dvv csv).  Boreholes: all-time 1-4s.  Non-boreholes: finalized / clean subset.
BORE = {
 'B011':(48.65,-123.448),'B013':(47.813,-122.9108),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B023':(46.1112,-123.0787),'B928':(48.834,-125.134),'B004':(48.202,-124.427),'B026':(45.3094,-123.8231),
 'B014':(47.5133,-123.8125),'B204':(46.136,-122.169),'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),
 'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),'B036':(42.5058,-123.3817),'B039':(41.4667,-122.4847),
 'B040':(41.8308,-122.4205),'B927':(49.2188,-124.8113),'B020':(46.3827,-123.8445),'B022':(45.9546,-123.931),
 'B935':(40.4787,-123.5732),'B201':(46.3033,-122.2648),
}
NONBORE = {  # station: (lat, lon, dvv file, tier)
 'GNW': (47.5641,-122.8250,'data/daily_dvv_GNW_tomo_e23.csv','A'),         # era 2+3 only
 'HDW': (47.6490,-123.0530,'data/daily_dvv_HDW_coda_1to4_full.csv','A'),
 'NLLB':(49.2271,-123.9882,'data/daily_dvv_NLLB_perera.csv','A'),
 'COR': (44.5855,-123.3046,'data/daily_dvv_COR_coda_1to4.csv','A'),        # cleanest, fills Oregon
 'PGC': (48.6498,-123.4521,'data/daily_dvv_PGC_tomo.csv','B'),             # >=2017-08, genuinely noisy (cc 0.87) -> lower tier
 'COLT':(45.17044,-122.438152,'data/daily_dvv_COLT_coda_1to4.csv','A'),    # accelerometer but cc 0.979 (unit-agnostic dv/v) -> full tier
}
FILES = {s:(la,lo,f'data/daily_dvv_{s}_coda_1to4.csv','A') for s,(la,lo) in BORE.items()}
FILES.update(NONBORE)

rows = []
for s,(sla,slo,f,tier) in FILES.items():
    if not os.path.exists(f): print('  MISSING', f); continue
    d = pd.read_csv(f); pre = d['patch'].astype(str).str.split('__').str[0]
    d['plat']=pre.str.split('_').str[0].astype(float); d['plon']=pre.str.split('_',n=1).str[1].astype(float)
    d['station']=s; d['sta_lat']=sla; d['sta_lon']=slo; d['tier']=tier
    rows.append(d[['station','sta_lat','sta_lon','tier','plat','plon','patch','date','dvv','cc_max','n_det']])
df = pd.concat(rows, ignore_index=True); df['date']=pd.to_datetime(df['date'])
df['cell_lat']=(np.round(df.plat/GRID)*GRID).round(3); df['cell_lon']=(np.round(df.plon/GRID)*GRID).round(3)
df['cell']=df.cell_lat.astype(str)+'_'+df.cell_lon.astype(str)
nlow = df[df.tier == 'B'].station.nunique()
print(f'loaded {len(df):,} rows across {df.station.nunique()} stations ({nlow} lower-tier)')

# Slab2 interface depth (interface-tracking constraints only)
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

cat.to_csv(f'{OUT}/tomo_patch_catalog.csv',index=False); cells.to_csv(f'{OUT}/tomo_cells.csv',index=False)
ten.to_parquet(f'{OUT}/tomo_dvv_monthly.parquet',index=False)
n3=(cells.n_stations>=3).sum(); n4=(cells.n_stations>=4).sum()
print(f'fault cells: {len(cells)} | >=3 stations: {n3} (was 123 borehole-only) | >=4: {n4}')
print(f'monthly tensor: {len(ten):,} (cell,station,month); {ten.ym.nunique()} months {ten.ym.min()}..{ten.ym.max()}')
print(f'wrote {OUT}/tomo_cells.csv, tomo_patch_catalog.csv, tomo_dvv_monthly.parquet')
