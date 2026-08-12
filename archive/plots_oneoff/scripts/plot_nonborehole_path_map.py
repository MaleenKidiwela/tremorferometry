#!/usr/bin/env python
"""Non-borehole (UW short-period + CN broadband) dv/v PATH MAP: each original station
connected to the LFE families that actually produced its FINALIZED dv/v curve, one color
per station. Companion to scripts/plot_pb_path_map.py (which does the PB boreholes).

Patch ids embed coords ('<lat>_<lon>__c<fam>'), same as the PB stations, so we parse the
family location straight from the dv/v csv -- but each original station has its own
finalized dv/v filename (per-era / clean / despike variants), listed below.
Re-run: PYTHONPATH=src python scripts/plot_nonborehole_path_map.py
"""
import os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.cm as cm, matplotlib.lines as ml
import cartopy.crs as ccrs, cartopy.feature as cfeature

# station -> (lat, lon, finalized dv/v csv, net/type label)
STATIONS = {
    'GNW':  (47.5641, -122.8250, 'data/daily_dvv_GNW_perera.csv', 'UW short-period'),
    'HDW':  (47.6490, -123.0530, 'data/daily_dvv_HDW_coda_1to4_perera.csv', 'UW short-period'),
    'CPW':  (46.9733, -123.1382, 'data/daily_dvv_CPW_clean_1to4.csv',       'UW short-period'),
    'PGC':  (48.6498, -123.4521, 'data/daily_dvv_PGC_perera_1to3.csv',      'CN broadband'),
    'NLLB': (49.2271, -123.9882, 'data/daily_dvv_NLLB_perera.csv',          'CN broadband'),
}


def patch_coords(csv):
    d = pd.read_csv(csv, usecols=['patch'])
    pts = []
    for p in d['patch'].drop_duplicates():
        pre = str(p).split('__')[0]
        pts.append((float(pre.split('_')[0]), float(pre.split('_', 1)[1])))
    return pd.DataFrame(pts, columns=['lat', 'lon'])


done = [s for s in STATIONS if os.path.exists(STATIONS[s][2])]
done.sort(key=lambda s: -STATIONS[s][0])   # north to south
patches = {s: patch_coords(STATIONS[s][2]) for s in done}
npatch = {s: len(patches[s]) for s in done}
print('non-borehole stations on path map:', done)

cmap = cm.get_cmap('Dark2')
colors = {s: cmap(i % 8) for i, s in enumerate(done)}

lats = [STATIONS[s][0] for s in done] + [v for s in done for v in patches[s].lat]
lons = [STATIONS[s][1] for s in done] + [v for s in done for v in patches[s].lon]
pad = 0.4
ext = [min(lons) - pad, max(lons) + pad, min(lats) - pad, max(lats) + pad]

fig = plt.figure(figsize=(10, 12)); ax = plt.axes(projection=ccrs.PlateCarree()); ax.set_extent(ext)
ax.add_feature(cfeature.LAND, facecolor='0.95'); ax.add_feature(cfeature.OCEAN, facecolor='#eaf3fb')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5, edgecolor='0.5')
ax.add_feature(cfeature.STATES.with_scale('50m'), lw=0.3, edgecolor='0.75')
gl = ax.gridlines(draw_labels=True, lw=0.2, color='0.85'); gl.top_labels = False; gl.right_labels = False
tk = dict(transform=ccrs.PlateCarree())

for s in done:
    la, lo, _, _ = STATIONS[s]; c = colors[s]; f = patches[s]
    for _, r in f.iterrows():
        ax.plot([lo, r.lon], [la, r.lat], color=c, lw=0.4, alpha=0.33, **tk)
    ax.scatter(f.lon, f.lat, s=12, color=c, alpha=0.7, edgecolor='none', **tk)
for s in done:
    la, lo, _, _ = STATIONS[s]
    ax.scatter([lo], [la], marker='*', s=320, color=colors[s], edgecolor='k', lw=1.0, zorder=10, **tk)
    ax.annotate(s, (lo, la), xytext=(5, 5), textcoords='offset points', fontsize=10, fontweight='bold')

hand = [ml.Line2D([], [], color=colors[s], marker='*', lw=1.2, markeredgecolor='k',
                  label=f'{s} ({npatch[s]} patches, {STATIONS[s][3]})') for s in done]
ax.legend(handles=hand, loc='lower left', fontsize=9, framealpha=0.92, title='station (dv/v patches)')
ax.set_title(f'Non-borehole (UW/CN) LFE-coda dv/v path map ({len(done)} stations):\n'
             'each station -> the families that produced its finalized dv/v')
plt.tight_layout(); plt.savefig('figures/smoke_nonborehole_path_map.png', dpi=130)
print('wrote figures/smoke_nonborehole_path_map.png; patches per station:', npatch)
