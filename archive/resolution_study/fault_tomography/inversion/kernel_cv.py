#!/usr/bin/env python
"""Kernel selection by CROSS-VALIDATION on the fixed 2-4 s deseasoned data.
Hold out 20% of (cell,station,month) rows, joint fault+site invert on 80%, predict the 20%.
The kernel with the best HELD-OUT variance reduction matches the physics (overfitting can't help).
Candidates: single-scatter (matched 2-4 / mismatched 1-4 / different mean-free-paths) + diffusion.
"""
import sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree
sys.path.insert(0, 'fault_tomography/kernels')
from kernel import kernel_singlescatter, kernel_diffusion

PFX = 'fault_tomography/data_assembly/tomo24des'
BETA, D = 3.5, 47.0
cells = pd.read_csv(f'{PFX}_cells.csv'); cat = pd.read_csv(f'{PFX}_patch_catalog.csv')
lat0, lon0 = cells.cell_lat.mean(), cells.cell_lon.mean()
def proj(lat, lon): return (lon-lon0)*111.0*np.cos(np.radians(lat0)), (lat-lat0)*111.0
cells['x'], cells['y'] = proj(cells.cell_lat.values, cells.cell_lon.values); cells = cells.reset_index(drop=True)
cidx = {c: i for i, c in enumerate(cells.cell)}; M = len(cells)
Xk = cells[['x', 'y']].values; Zk = cells.depth_km.values
cat = cat[cat.cell.isin(cidx)].copy()
cat['ax'], cat['ay'] = proj(cat.cell_lat.values, cat.cell_lon.values)
cat['bx'], cat['by'] = proj(cat.sta_lat.values, cat.sta_lon.values)

def buildGf(kind, w1, w2, ell):
    rows = []
    for _, r in cat.iterrows():
        s_h = np.array([r.ax, r.ay]); d = r.depth_km; r_h = np.array([r.bx, r.by])
        if kind == 'ss':
            K = kernel_singlescatter(Xk, Zk, s_h, d, r_h, BETA, w1, w2, ell=ell)
        else:
            K = kernel_diffusion(Xk, Zk, s_h, d, r_h, D, (w1+w2)/2 + 2.5)
        rows.append(-K)
    return np.array(rows)

# graph Laplacian
P = cells[['x', 'y']].values; tree = cKDTree(P); rr, cc, vv = [], [], []
for i in range(M):
    _, nn = tree.query(P[i], 7)
    for j in nn[1:]:
        rr += [i, i]; cc += [i, j]; vv += [1.0, -1.0]
L = sparse.coo_matrix((vv, (rr, cc)), shape=(M, M)).tocsr()

stations = sorted(cat.station.unique()); sidx = {s: i for i, s in enumerate(stations)}; B = len(stations)
Gs = sparse.csr_matrix((np.ones(len(cat)), (np.arange(len(cat)), cat.station.map(sidx).values)), shape=(len(cat), B))
key = cat[['cell', 'station']].copy(); key['row'] = np.arange(len(cat))
ten = pd.read_parquet(f'{PFX}_dvv_monthly.parquet'); months = sorted(ten.ym.unique())
lam_f, lam_s = 0.4, 0.05
regf = sparse.hstack([lam_f*L, sparse.csr_matrix((M, B))])
regs = sparse.hstack([sparse.csr_matrix((B, M)), lam_s*sparse.identity(B)])

def cv(Gf):
    Gf = -Gf  # +K
    obs_all, prd_all = [], []
    for E in months:
        te = ten[ten.ym == E][['cell', 'station', 'dvv']]
        mm = key.merge(te, on=['cell', 'station'])
        if len(mm) < 60 or mm.station.nunique() < 10: continue
        rows = mm.row.values; d = mm.dvv.values*100.0
        test = (np.arange(len(mm)) % 5 == 0)           # deterministic 20% hold-out
        if test.sum() < 5 or (~test).sum() < 40: continue
        tr_rows = rows[~test]; dtr = d[~test]
        top = sparse.hstack([sparse.csr_matrix(Gf[tr_rows]), Gs[tr_rows]])
        A = sparse.vstack([top, regf, regs]).tocsr()
        b = np.concatenate([dtr, np.zeros(M), np.zeros(B)])
        sol = lsqr(A, b, atol=1e-7, btol=1e-7, iter_lim=1500)[0]
        pred = Gf[rows[test]] @ sol[:M] + Gs[rows[test]] @ sol[M:]
        obs_all.append(d[test]); prd_all.append(pred)
    o = np.concatenate(obs_all); p = np.concatenate(prd_all)
    return 1 - np.var(o - p)/np.var(o), len(o)

cands = [('single-scatter 2-4 ell40 (matched)', 'ss', 2.0, 4.0, 40),
         ('single-scatter 1-4 ell40 (mismatch)', 'ss', 1.0, 4.0, 40),
         ('single-scatter 2-4 ell20', 'ss', 2.0, 4.0, 20),
         ('single-scatter 2-4 ell80', 'ss', 2.0, 4.0, 80),
         ('diffusion D47', 'diff', 2.0, 4.0, 0)]
print(f'CROSS-VALIDATION (held-out VR), 2-4 s deseasoned data, {len(months)} months:')
for name, kind, w1, w2, ell in cands:
    Gf = buildGf(kind, w1, w2, ell)
    vr, n = cv(Gf)
    print(f'  {name:38s}  held-out VR = {100*vr:5.1f}%   (n_test={n})')
