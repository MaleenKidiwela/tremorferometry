#!/usr/bin/env python
"""Multi-window joint inversion for the interface delta-beta/beta field, using 1-3, 2-4, 3-5 s.
One SHARED fault field m_f; each window contributes its own single-scatter kernel (different ellipsoid
size -> depth leverage) and its own per-station SITE term (later coda = more near-receiver, so the surface
contribution is window-specific).  The fault is the across-window common mode; the site absorbs the
window-specific surface signal.  Monthly joint LSQR -> 4-D m_f + spatial temporal-RMS.  Deseasoned
calendar (true 30-cal-day) data.  Also a checkerboard with the combined 3-window geometry.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
sys.path.insert(0, 'fault_tomography/kernels')
from kernel import kernel_singlescatter

GRID = 0.10; BETA = 3.5; LSTAR = 40.0; CCMIN = 0.7
DES = os.environ.get('DESEASON', '1') == '1'          # DESEASON=0 -> raw (seasonal kept; should land in site terms)
_sfx = '_cal_des' if DES else '_cal'
WINS = [('13', 1.0, 3.0, '1to3'+_sfx), ('24', 2.0, 4.0, '2to4'+_sfx), ('35', 3.0, 5.0, '3to5'+_sfx)]
STA = {
 'B927':(49.2188,-124.8113),'NLLB':(49.2271,-123.9882),'B928':(48.834,-125.134),'PGC':(48.6498,-123.4521),
 'B011':(48.65,-123.448),'B004':(48.202,-124.427),'B013':(47.813,-122.9108),'HDW':(47.6490,-123.0530),
 'GNW':(47.5641,-122.8250),'B014':(47.5133,-123.8125),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B020':(46.3827,-123.8445),'B201':(46.3033,-122.2648),'B204':(46.136,-122.169),'B023':(46.1112,-123.0787),
 'B022':(45.9546,-123.931),'B026':(45.3094,-123.8231),'COLT':(45.17044,-122.438152),'COR':(44.5855,-123.3046),
 'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),
 'B036':(42.5058,-123.3817),'B040':(41.8308,-122.4205),'B039':(41.4667,-122.4847),'B935':(40.4787,-123.5732),
}

# ---- 1. assemble combined monthly tensor across the 3 windows ----
parts = []
for s, (sla, slo) in STA.items():
    for wk, w1, w2, suf in WINS:
        f = f'data/daily_dvv_{s}_{suf}.csv'
        if not os.path.exists(f): continue
        d = pd.read_csv(f); d = d[d.cc_max > CCMIN]
        if not len(d): continue
        pre = d.patch.astype(str).str.split('__').str[0]
        plat = pre.str.split('_').str[0].astype(float); plon = pre.str.split('_', n=1).str[1].astype(float)
        d = d.assign(clat=(np.round(plat/GRID)*GRID).round(3), clon=(np.round(plon/GRID)*GRID).round(3))
        d['cell'] = d.clat.astype(str) + '_' + d.clon.astype(str)
        d['ym'] = pd.to_datetime(d.date).dt.to_period('M').astype(str)
        d['station'] = s; d['sta_lat'] = sla; d['sta_lon'] = slo; d['win'] = wk
        g = d.groupby(['cell', 'station', 'win', 'ym']).agg(
            dvv=('dvv', 'mean'), clat=('clat', 'first'), clon=('clon', 'first'),
            sta_lat=('sta_lat', 'first'), sta_lon=('sta_lon', 'first')).reset_index()
        parts.append(g)
ten = pd.concat(parts, ignore_index=True)
print(f'combined tensor: {len(ten):,} (cell,station,win,month) rows | {ten.cell.nunique()} cells | '
      f'{ten.station.nunique()} stations | windows {sorted(ten.win.unique())}')

cells = ten[['cell', 'clat', 'clon']].drop_duplicates('cell').reset_index(drop=True)
slab = pd.read_csv('data/cas_slab2_input_04-18.csv', low_memory=False)[['lat', 'lon', 'depth', 'etype']].dropna(subset=['depth'])
slab = slab[(~slab.etype.isin(['BA', 'TO'])) & (slab.depth.between(10, 70))]
q = cells[['clat', 'clon']].values
dl = griddata(slab[['lat', 'lon']].values, slab.depth.values, q, 'linear')
dn = griddata(slab[['lat', 'lon']].values, slab.depth.values, q, 'nearest')
cells['depth_km'] = np.clip(np.where(np.isnan(dl), dn, dl), 12, 55)
nst = ten.groupby('cell')['station'].nunique().rename('n_stations')
cells = cells.merge(nst, on='cell')
lat0, lon0 = cells.clat.mean(), cells.clon.mean()
def proj(lat, lon): return (np.asarray(lon)-lon0)*111.0*np.cos(np.radians(lat0)), (np.asarray(lat)-lat0)*111.0
cells['x'], cells['y'] = proj(cells.clat.values, cells.clon.values)
M = len(cells); cidx = {c: i for i, c in enumerate(cells.cell)}
Xk = cells[['x', 'y']].values; Zk = cells.depth_km.values

# ---- 2. data pairs (cell,station,win) and per-window G_f + per-(station,win) site ----
cat = ten[['cell', 'station', 'win', 'clat', 'clon', 'sta_lat', 'sta_lon']].drop_duplicates(['cell', 'station', 'win']).reset_index(drop=True)
cat['a'] = cat.cell.map(cidx); cat['depth'] = cat.cell.map(cells.set_index('cell').depth_km)
cat['ax'], cat['ay'] = proj(cat.clat.values, cat.clon.values)
cat['bx'], cat['by'] = proj(cat.sta_lat.values, cat.sta_lon.values)
WB = {wk: (w1, w2) for wk, w1, w2, _ in WINS}
Gf = np.zeros((len(cat), M))
for i, r in cat.iterrows():
    w1, w2 = WB[r.win]
    K = kernel_singlescatter(Xk, Zk, np.array([r.ax, r.ay]), r.depth, np.array([r.bx, r.by]), BETA, w1, w2, ell=LSTAR)
    Gf[i] = K                                            # +K (sums to 1); d = -Gf m  -> flip below
Gf = -Gf
sw = (cat.station + '|' + cat.win); swu = sorted(sw.unique()); swidx = {x: i for i, x in enumerate(swu)}; B = len(swu)
Gs = sparse.csr_matrix((np.ones(len(cat)), (np.arange(len(cat)), sw.map(swidx).values)), shape=(len(cat), B))
key = cat[['cell', 'station', 'win']].copy(); key['row'] = np.arange(len(cat))
print(f'model: {M} fault cells | {B} site terms (station x window) | {len(cat)} data rows/full-coverage')

def graph_laplacian(P, k=6):
    tree = cKDTree(P); m = len(P); rr, cc, vv = [], [], []
    for i in range(m):
        _, nn = tree.query(P[i], k+1)
        for j in nn[1:]:
            rr += [i, i]; cc += [i, j]; vv += [1.0, -1.0]
    return sparse.coo_matrix((vv, (rr, cc)), shape=(m, m)).tocsr()
L = graph_laplacian(cells[['x', 'y']].values)

# ---- 3. monthly joint inversion ----
lam_f, lam_s = 0.4, 0.05
regf = sparse.hstack([lam_f*L, sparse.csr_matrix((M, B))])
regs = sparse.hstack([sparse.csr_matrix((B, M)), lam_s*sparse.identity(B)])
months = sorted(ten.ym.unique())
MF = np.full((M, len(months)), np.nan); SITE = np.full((B, len(months)), np.nan)
VR = np.zeros(len(months)); ND = np.zeros(len(months), int)
for j, E in enumerate(months):
    te = ten[ten.ym == E][['cell', 'station', 'win', 'dvv']]
    mm = key.merge(te, on=['cell', 'station', 'win'])
    if len(mm) < 60 or mm.station.nunique() < 10: continue
    rows = mm.row.values; dd = mm.dvv.values*100.0
    top = sparse.hstack([sparse.csr_matrix(Gf[rows]), Gs[rows]])
    Aaug = sparse.vstack([top, regf, regs]).tocsr()
    baug = np.concatenate([dd, np.zeros(M), np.zeros(B)])
    sol = lsqr(Aaug, baug, atol=1e-7, btol=1e-7, iter_lim=2000)[0]
    mf = sol[:M]; pred = Gf[rows] @ mf + Gs[rows] @ sol[M:]
    MF[:, j] = mf; SITE[:, j] = sol[M:]; VR[j] = 1 - np.var(dd - pred)/np.var(dd); ND[j] = len(mm)

ok = ND > 0
cell_var = np.nanstd(MF[:, ok], axis=1)
well = (cells.n_stations >= 3).values
print(f'inverted {ok.sum()}/{len(months)} months; mean VR {100*np.nanmean(VR[ok]):.0f}%')
print(f'resolved cells (>=3 sta): {well.sum()} | temporal RMS among resolved: median {np.median(cell_var[well]):.3f}% max {cell_var[well].max():.3f}%')
fidx = np.nanmean(MF[well][:, ok], axis=0)
print(f'network fault index range [{np.nanmin(fidx):.3f},{np.nanmax(fidx):.3f}] %')
OUT = 'fault_tomography/inversion/fault_4d_multiwin%s.npz' % ('' if DES else '_raw')
np.savez(OUT, months=np.array(months), lat=cells.clat.values, lon=cells.clon.values, depth=cells.depth_km.values,
         nsta=cells.n_stations.values, MF=MF, VR=VR, ND=ND, ok=ok, well=well, cell_var=cell_var, fault_idx=fidx,
         SITE=SITE, site_labels=np.array(swu), deseasoned=DES)
print('wrote', OUT, '(deseasoned=%s)' % DES)

# ---- checkerboard with the combined 3-window geometry ----
m_true = np.sign(np.sin(np.pi*cells.clat.values/0.4)*np.sin(np.pi*cells.clon.values/0.4))*0.01
d = (-Gf) @ m_true; d += 0.05*np.std(d)*np.sin(np.arange(len(d)))
A = sparse.vstack([sparse.csr_matrix(-Gf), 0.3*L]); b = np.concatenate([d, np.zeros(M)])
m_rec = lsqr(A, b, atol=1e-6, btol=1e-6, iter_lim=2000)[0]
corr = np.corrcoef(m_true[well], m_rec[well])[0, 1]
print(f'CHECKERBOARD (3-window combined geometry, >=3-sta cells): corr {corr:.2f} (n={well.sum()})')
