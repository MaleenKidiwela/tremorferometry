#!/usr/bin/env python
"""A1 FEASIBILITY CHECK (fault_tomography Phase A): do our LFE families form the
(fault-patch x station) bipartite redundancy the inversion needs?

The deep sensitivity lobe is COMMON-MODE across every station that records the same
family (theory/01_kernel.md sec 6). A single station cannot separate a fault-zone change
from a change anywhere along its up-going tube; you need an AZIMUTHAL RING of stations
sharing the same family/patch (theory/02 sec 6). So the key feasibility metric is:

    how many fault grid-cells are imaged by >= 3 stations, and how well-spread in azimuth?

Each dv/v patch id is '<lat>_<lon>__c<n>'; the '<lat>_<lon>' 0.05-deg grid cell IS the
fault-patch location, shared across stations (discovery is PNSN-catalog-driven). So we
match families across stations purely on that prefix. Read-only on data/daily_dvv_*.
"""
import glob, os
from collections import defaultdict
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

# station coords (only needed for the azimuth-spread refinement); PB + originals we know
STA = {
 'B011':(48.65,-123.448),'B013':(47.813,-122.9108),'B941':(46.9868,-122.219),
 'B018':(46.9795,-123.0203),'B023':(46.1112,-123.0787),'B928':(48.834,-125.134),
 'B004':(48.202,-124.427),'B026':(45.3094,-123.8231),'B014':(47.5133,-123.8125),
 'B204':(46.136,-122.169),'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),
 'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),'B036':(42.5058,-123.3817),
 'B039':(41.4667,-122.4847),'B040':(41.8308,-122.4205),'B927':(49.2188,-124.8113),
 'B020':(46.3827,-123.8445),'B022':(45.9546,-123.931),
 'GNW':(47.5641,-122.825),'PGC':(48.6498,-123.4521),'NLLB':(49.2271,-123.9882),
}

def cell_lalo(prefix):
    la = float(prefix.split('_')[0]); lo = float(prefix.split('_',1)[1]); return la, lo

files = sorted(glob.glob('data/daily_dvv_*_coda_1to4.csv'))
cell_st = defaultdict(set); sta_cells = {}
for f in files:
    sta = os.path.basename(f).replace('daily_dvv_','').replace('_coda_1to4.csv','')
    cells = set(str(p).split('__')[0] for p in pd.read_csv(f, usecols=['patch'])['patch'].unique())
    sta_cells[sta] = cells
    for c in cells: cell_st[c].add(sta)

counts = np.array([len(v) for v in cell_st.values()])
print(f'stations with coda_1to4 dv/v : {len(sta_cells)}  ({", ".join(sorted(sta_cells))})')
print(f'distinct fault grid-cells     : {len(cell_st)}')
for k in (1,2,3,4,5):
    print(f'  cells seen by >= {k} stations : {(counts>=k).sum():4d}  ({100*(counts>=k).mean():3.0f}%)')
print(f'  max stations on one cell     : {counts.max()}')

# azimuth spread for the well-covered (>=3) cells
def azspread(cell):
    la, lo = cell_lalo(cell)
    azs = []
    for s in cell_st[cell]:
        if s not in STA: continue
        sla, slo = STA[s]
        azs.append(np.degrees(np.arctan2((slo-lo)*np.cos(np.radians(la)), sla-la)) % 360)
    if len(azs) < 2: return 0.0
    azs = np.sort(azs); gaps = np.diff(np.r_[azs, azs[0]+360])
    return 360 - gaps.max()   # angular coverage = 360 minus largest empty wedge

well = [c for c in cell_st if len(cell_st[c]) >= 3]
spreads = np.array([azspread(c) for c in well]) if well else np.array([])
if len(spreads):
    print(f'  of >=3-station cells, azimuth coverage: median {np.median(spreads):.0f} deg, '
          f'{(spreads>=180).sum()} cells span >=180 deg')

# map: cells colored by number of imaging stations
os.makedirs('fault_tomography/figures', exist_ok=True)
las = np.array([cell_lalo(c)[0] for c in cell_st]); los = np.array([cell_lalo(c)[1] for c in cell_st])
fig, ax = plt.subplots(figsize=(8, 11))
sc = ax.scatter(los, las, c=counts, cmap='turbo', s=14, vmin=1, vmax=max(4, counts.max()))
for s,(la,lo) in STA.items():
    if s in sta_cells: ax.plot(lo, la, marker='^', ms=8, color='k', mfc='white', mew=1.2)
plt.colorbar(sc, label='# stations imaging this fault cell')
ax.set_xlabel('lon'); ax.set_ylabel('lat'); ax.set_aspect(1/np.cos(np.radians(46)))
ax.set_title('A1: fault-patch x station overlap\n(color = # stations sharing each LFE patch; >=3 = deep-lobe resolvable)')
plt.tight_layout(); plt.savefig('fault_tomography/figures/a1_patch_overlap.png', dpi=130)
print('wrote fault_tomography/figures/a1_patch_overlap.png')
