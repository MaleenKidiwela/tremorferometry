#!/usr/bin/env python
"""PB borehole dv/v PATH MAP: each done PB station connected to its coverage-selected LFE families
(the families used to generate that station's dv/v), one color per station. Auto-discovers every PB
station that has BOTH a <sta>_coverage_selection.summary.csv and a daily_dvv_<STA>_coda_1to4.csv.
Re-run after each new station completes: PYTHONPATH=src python scripts/plot_pb_path_map.py
"""
import os, glob
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.cm as cm, matplotlib.lines as ml
import cartopy.crs as ccrs, cartopy.feature as cfeature

# known PB borehole coords (lat, lon)
COORDS = {
    'B011': (48.65, -123.448), 'B013': (47.813, -122.9108), 'B941': (46.9868, -122.219),
    'B018': (46.9795, -123.0203), 'B023': (46.1112, -123.0787),
    'B928': (48.834, -125.134), 'B004': (48.202, -124.427),
    'B026': (45.3094, -123.8231), 'B014': (47.5133, -123.8125),
    'B204': (46.136, -122.169),
    'B028': (44.4937, -122.9638), 'B030': (43.9713, -122.7717),
    'B032': (43.668, -123.3923), 'B033': (43.2917, -123.1245), 'B036': (42.5058, -123.3817),
    'B039': (41.4667, -122.4847), 'B040': (41.8308, -122.4205), 'B927': (49.2188, -124.8113),
}

# discover done PB stations: have coverage selection + dv/v
done = []
for s in COORDS:
    cov = f'data/{s.lower()}_coverage_selection.summary.csv'
    dvv = f'data/daily_dvv_{s}_coda_1to4.csv'
    if os.path.exists(cov) and os.path.exists(dvv):
        done.append(s)
done.sort(key=lambda s: -COORDS[s][0])   # north to south
print('PB stations on path map:', done)

cmap = cm.get_cmap('tab10')
colors = {s: cmap(i % 10) for i, s in enumerate(done)}

lats = [COORDS[s][0] for s in done] + [r for s in done for r in pd.read_csv(f'data/{s.lower()}_coverage_selection.summary.csv').lat]
lons = [COORDS[s][1] for s in done] + [r for s in done for r in pd.read_csv(f'data/{s.lower()}_coverage_selection.summary.csv').lon]
pad = 0.4
ext = [min(lons)-pad, max(lons)+pad, min(lats)-pad, max(lats)+pad]

fig = plt.figure(figsize=(10, 12)); ax = plt.axes(projection=ccrs.PlateCarree()); ax.set_extent(ext)
ax.add_feature(cfeature.LAND, facecolor='0.95'); ax.add_feature(cfeature.OCEAN, facecolor='#eaf3fb')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), lw=0.5, edgecolor='0.5')
ax.add_feature(cfeature.STATES.with_scale('50m'), lw=0.3, edgecolor='0.75')
gl = ax.gridlines(draw_labels=True, lw=0.2, color='0.85'); gl.top_labels = False; gl.right_labels = False
tk = dict(transform=ccrs.PlateCarree())

nfam = {}
for s in done:
    la, lo = COORDS[s]; c = colors[s]
    f = pd.read_csv(f'data/{s.lower()}_coverage_selection.summary.csv'); nfam[s] = len(f)
    for _, r in f.iterrows():
        ax.plot([lo, r.lon], [la, r.lat], color=c, lw=0.4, alpha=0.33, **tk)
    ax.scatter(f.lon, f.lat, s=10, color=c, alpha=0.7, edgecolor='none', **tk)
for s in done:
    la, lo = COORDS[s]
    ax.scatter([lo], [la], marker='^', s=160, color=colors[s], edgecolor='k', lw=1.0, zorder=10, **tk)
    ax.annotate(s, (lo, la), xytext=(4, 4), textcoords='offset points', fontsize=9, fontweight='bold')

hand = [ml.Line2D([], [], color=colors[s], marker='^', lw=1.2, markeredgecolor='k', label=f'{s} ({nfam[s]} fam)') for s in done]
ax.legend(handles=hand, loc='upper right', fontsize=9, framealpha=0.92, title='station (coverage families)')
ax.set_title(f'PB borehole dv/v path map ({len(done)} stations): each station -> its coverage LFE families')
plt.tight_layout(); plt.savefig('figures/smoke_pb_path_map.png', dpi=130)
print('wrote figures/smoke_pb_path_map.png; families per station:', nfam)
