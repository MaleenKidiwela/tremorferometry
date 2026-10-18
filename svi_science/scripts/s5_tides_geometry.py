"""S5 setup, step 1: fault geometry at each SVI LFE family.

Strike/dip from the Slab2 depth grid (cas_slab2_dep_02.24.18.grd) at the family epicentre
(plane fit to slab2 depths within 10 km). Local ENU frame (x=east, y=north, z=up).
  n      = unit normal pointing UP into the hanging wall (North America)
  s_conv = hanging-wall slip direction: horizontal azimuth 236 deg (NA relative to JdF,
           the "fixed 236 deg" direction used in S1) projected onto the fault plane
  s_up   = pure up-dip (rake 90) direction
  s_str  = Aki-Richards strike direction (rake 0), strike = downdip azimuth - 90
Rake of s_conv is reported (Aki-Richards). Also copies the slip-model rake from
lfe_extension/results/slip_compare/slip_svi_family_meta.csv for comparison.
Output: results/s5_setup/family_fault_geometry.csv
"""
import numpy as np, pandas as pd, xarray as xr
from pathlib import Path

ROOT = Path('/home/jovyan/tremorferometry')
OUT = ROOT / 'svi_science/results/s5_setup'
OUT.mkdir(parents=True, exist_ok=True)
CONV_AZ = 236.0

fam = pd.read_csv(ROOT / 'svi_science/results/s0/families.csv')
slab = xr.open_dataset(ROOT / 'data/slab_models/slab2/cas_slab2_dep_02.24.18.grd')
X, Y = np.meshgrid(slab.x.values - 360.0, slab.y.values)
Z = -slab.z.values.astype(float)  # depth km positive down

rows = []
for _, f in fam.iterrows():
    lat0, lon0 = f.lat, f.lon
    dx = (X - lon0) * 111.195 * np.cos(np.radians(lat0))
    dy = (Y - lat0) * 111.195
    m = (np.hypot(dx, dy) <= 10.0) & np.isfinite(Z)
    A = np.c_[np.ones(m.sum()), dx[m], dy[m]]
    c, *_ = np.linalg.lstsq(A, Z[m], rcond=None)
    gx, gy = c[1], c[2]                       # d(depth)/dx, d(depth)/dy  (km/km)
    dip = np.degrees(np.arctan(np.hypot(gx, gy)))
    az_dd = np.degrees(np.arctan2(gx, gy)) % 360   # downdip azimuth (depth increases)
    strike = (az_dd - 90) % 360
    d, a = np.radians(dip), np.radians(az_dd)
    n = np.array([np.sin(d) * np.sin(a), np.sin(d) * np.cos(a), np.cos(d)])
    v_dd = np.array([np.cos(d) * np.sin(a), np.cos(d) * np.cos(a), -np.sin(d)])
    s_up = -v_dd
    s_str = np.array([np.sin(np.radians(strike)), np.cos(np.radians(strike)), 0.0])
    h = np.array([np.sin(np.radians(CONV_AZ)), np.cos(np.radians(CONV_AZ)), 0.0])
    s_conv = h - h.dot(n) * n
    s_conv /= np.linalg.norm(s_conv)
    rake = np.degrees(np.arctan2(s_conv.dot(s_up), s_conv.dot(s_str)))
    rows.append(dict(fid=int(f.fid), lat=lat0, lon=lon0, dep_cat_km=f.dep,
                     slab2_depth_km=f.slab2_depth_km, strike=strike, dip=dip, az_downdip=az_dd,
                     rake_conv236=rake, n_e=n[0], n_n=n[1], n_u=n[2],
                     s_e=s_conv[0], s_n=s_conv[1], s_u=s_conv[2],
                     sup_e=s_up[0], sup_n=s_up[1], sup_u=s_up[2],
                     sstr_e=s_str[0], sstr_n=s_str[1], sstr_u=s_str[2]))
g = pd.DataFrame(rows)
try:
    meta = pd.read_csv(ROOT / 'lfe_extension/results/slip_compare/slip_svi_family_meta.csv')
    meta = meta[['lat', 'lon', 'rake_deg']].rename(columns={'rake_deg': 'rake_slipmodel'})
    g = g.merge(meta.drop_duplicates(['lat', 'lon']), on=['lat', 'lon'], how='left')
except Exception as e:
    print('slip-model rake not merged:', e)
g.to_csv(OUT / 'family_fault_geometry.csv', index=False)
print(g[['strike', 'dip', 'az_downdip', 'rake_conv236', 'rake_slipmodel', 'dep_cat_km', 'slab2_depth_km']].describe().round(2))
