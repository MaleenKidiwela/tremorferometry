#!/usr/bin/env python
"""A2 (borehole-only first pass): assemble the (fault-cell x station x epoch) dv/v data tensor
from the 22 PB borehole stations -- all processed identically (single-era EHZ, all-time 1-4s ref,
--no-deconv, --despike-mad 8), so no per-era / co-location complications.

Outputs (in fault_tomography/data_assembly/):
  borehole_patch_catalog.csv   -- one row per (fault-cell, station): coverage backbone for G's rows
  borehole_dvv_monthly.parquet -- (cell, station, year-month) -> mean dv/v, std, n  (the 4-D tensor)
  borehole_cells.csv           -- one row per fault cell: lat, lon, Slab2 depth, #stations (the mesh seed)
Fault cell = 0.10 deg (~11 km), the resolution A1 settled. Depth interpolated from the Slab2 input DB.
"""
import os
import numpy as np, pandas as pd
from scipy.interpolate import griddata

GRID = 0.10  # fault-cell size (deg)
OUT = 'fault_tomography/data_assembly'
os.makedirs(OUT, exist_ok=True)

# 22 PB borehole stations (lat, lon) -- the uniform first-pass set
STA = {
 'B011':(48.65,-123.448),'B013':(47.813,-122.9108),'B941':(46.9868,-122.219),
 'B018':(46.9795,-123.0203),'B023':(46.1112,-123.0787),'B928':(48.834,-125.134),
 'B004':(48.202,-124.427),'B026':(45.3094,-123.8231),'B014':(47.5133,-123.8125),
 'B204':(46.136,-122.169),'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),
 'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),'B036':(42.5058,-123.3817),
 'B039':(41.4667,-122.4847),'B040':(41.8308,-122.4205),'B927':(49.2188,-124.8113),
 'B020':(46.3827,-123.8445),'B022':(45.9546,-123.931),'B935':(40.4787,-123.5732),
 'B201':(46.3033,-122.2648),
}

# --- load all borehole dv/v, parse patch fault coords, snap to cells ---
rows = []
for s,(sla,slo) in STA.items():
    f = f'data/daily_dvv_{s}_coda_1to4.csv'
    if not os.path.exists(f):
        print('  MISSING', f); continue
    d = pd.read_csv(f)
    pre = d['patch'].astype(str).str.split('__').str[0]
    d['plat'] = pre.str.split('_').str[0].astype(float)
    d['plon'] = pre.str.split('_', n=1).str[1].astype(float)
    d['station'] = s; d['sta_lat'] = sla; d['sta_lon'] = slo
    rows.append(d[['station','sta_lat','sta_lon','plat','plon','patch','date','dvv','cc_max','n_det']])
df = pd.concat(rows, ignore_index=True)
df['date'] = pd.to_datetime(df['date'])
df['cell_lat'] = (np.round(df.plat/GRID)*GRID).round(3)
df['cell_lon'] = (np.round(df.plon/GRID)*GRID).round(3)
df['cell'] = df.cell_lat.astype(str)+'_'+df.cell_lon.astype(str)
print(f'loaded {len(df):,} (station,patch,day) dv/v rows across {df.station.nunique()} boreholes')

# --- Slab2 interface depth at each fault cell ---
# NB: the smooth Slab2 OUTPUT grid is firewall-blocked; we approximate from the INPUT constraint DB.
# That DB mixes constraint types with very different depths, so use only interface-TRACKING ones:
# drop BA (seafloor bathymetry, 1-5 km) and TO (deep tomography, 60-440 km), and keep constraints
# in the interface band (10-70 km) to reject shallow crustal seismicity. Clip to the physical
# Cascadia tremor-zone interface range [12,55] km. (Approximate -- replace with Slab2 output if obtainable.)
slab = pd.read_csv('data/cas_slab2_input_04-18.csv', low_memory=False)[['lat','lon','depth','etype']].dropna(subset=['depth'])
slab = slab[(~slab.etype.isin(['BA','TO'])) & (slab.depth.between(10,70))]
cells = df[['cell','cell_lat','cell_lon']].drop_duplicates('cell').reset_index(drop=True)
pts = slab[['lat','lon']].values; dep = slab['depth'].values
q = cells[['cell_lat','cell_lon']].values
d_lin = griddata(pts, dep, q, method='linear')
d_near = griddata(pts, dep, q, method='nearest')
cells['depth_km'] = np.clip(np.where(np.isnan(d_lin), d_near, d_lin), 12, 55)

# --- per-(cell, station) coverage catalog ---
cat = (df.groupby(['cell','station'])
         .agg(n_families=('patch','nunique'), n_obs=('dvv','size'),
              sta_lat=('sta_lat','first'), sta_lon=('sta_lon','first'),
              date_min=('date','min'), date_max=('date','max')).reset_index())
cat = cat.merge(cells, on='cell')
nsta_per_cell = cat.groupby('cell')['station'].nunique().rename('n_stations')
cells = cells.merge(nsta_per_cell, on='cell')

# --- monthly (cell, station, epoch) tensor: family-aggregated dv/v ---
df['ym'] = df.date.dt.to_period('M').astype(str)
w = df.n_det.clip(lower=1).astype(float)            # weight by detections
df['_wd'] = df.dvv * w
g = df.groupby(['cell','station','ym'])
ten = g.agg(sta_lat=('sta_lat','first'), sta_lon=('sta_lon','first'),
            cell_lat=('cell_lat','first'), cell_lon=('cell_lon','first'),
            dvv_std=('dvv','std'), n=('dvv','size')).reset_index()
ten['dvv'] = (g['_wd'].sum()/g.apply(lambda x: x.n_det.clip(lower=1).sum())).values
ten = ten.merge(cells[['cell','depth_km','n_stations']], on='cell')

cat.to_csv(f'{OUT}/borehole_patch_catalog.csv', index=False)
cells.to_csv(f'{OUT}/borehole_cells.csv', index=False)
ten.to_parquet(f'{OUT}/borehole_dvv_monthly.parquet', index=False)

# --- report ---
n3 = (cells.n_stations>=3).sum(); n4=(cells.n_stations>=4).sum()
print(f'fault cells: {len(cells)}  | >=3 stations: {n3}  >=4: {n4}')
print(f'depth range of >=3-station cells: {cells[cells.n_stations>=3].depth_km.min():.0f}-{cells[cells.n_stations>=3].depth_km.max():.0f} km')
print(f'monthly tensor: {len(ten):,} (cell,station,month) cells; {ten.ym.nunique()} months {ten.ym.min()}..{ten.ym.max()}')
print(f'wrote -> {OUT}/borehole_patch_catalog.csv, borehole_cells.csv, borehole_dvv_monthly.parquet')
