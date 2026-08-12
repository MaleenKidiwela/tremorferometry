#!/usr/bin/env python
"""Phase E (borehole-only first pass): 4-D loop. Joint fault+site inversion of EVERY month ->
time series of the interface delta-beta/beta map. Then surface where/when a coherent fault-zone
change appears (the science: does any resolved patch move together across its stations, esp. during ETS).
Builds G_f once (expensive); each month is a cheap LSQR over the rows present that month.
"""
import sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
sys.path.insert(0, 'fault_tomography/inversion')
from forward import build, graph_laplacian

cells, cat, Gf, eta, M = build()
Gf = -Gf                                          # +K
cat = cat.reset_index(drop=True)
stations = sorted(cat.station.unique()); sidx = {s: i for i, s in enumerate(stations)}; B = len(stations)
Gs = sparse.csr_matrix((np.ones(len(cat)), (np.arange(len(cat)), cat.station.map(sidx).values)), shape=(len(cat), B))
L = graph_laplacian(cells)
key = cat[['cell', 'station']].copy(); key['row'] = np.arange(len(cat))

import os
PFX = os.environ.get('TOMO_PFX', 'fault_tomography/data_assembly/tomo')
SUF = os.path.basename(PFX)   # 'tomo' / 'borehole' / 'tomo24des' -> fault_4d_{SUF}.npz (no clobber)
ten = pd.read_parquet(f'{PFX}_dvv_monthly.parquet')
months = sorted(ten.ym.unique())
lam_f, lam_s = 0.4, 0.05
regf = sparse.hstack([lam_f * L, sparse.csr_matrix((M, B))])
regs = sparse.hstack([sparse.csr_matrix((B, M)), lam_s * sparse.identity(B)])

MF = np.full((M, len(months)), np.nan); VR = np.zeros(len(months)); ND = np.zeros(len(months), int)
for j, E in enumerate(months):
    te = ten[ten.ym == E][['cell', 'station', 'dvv']]
    mm = key.merge(te, on=['cell', 'station'])
    if len(mm) < 50 or mm.station.nunique() < 10:
        continue
    rows = mm.row.values; d = mm.dvv.values * 100.0
    top = sparse.hstack([sparse.csr_matrix(Gf[rows]), Gs[rows]])
    Aaug = sparse.vstack([top, regf, regs]).tocsr()
    baug = np.concatenate([d, np.zeros(M), np.zeros(B)])
    sol = lsqr(Aaug, baug, atol=1e-7, btol=1e-7, iter_lim=2000)[0]
    mf = sol[:M]; pred = Gf[rows] @ mf + Gs[rows] @ sol[M:]
    MF[:, j] = mf; VR[j] = 1 - np.var(d - pred) / np.var(d); ND[j] = len(mm)

t = pd.to_datetime([m + '-15' for m in months])
well = (cells.n_stations >= 3).values
ok = ND > 0
# per-cell temporal std among resolved cells -> where does the fault map vary most over time?
cell_var = np.nanstd(MF[:, ok], axis=1)
print(f'inverted {ok.sum()}/{len(months)} months (>=10 sta, >=50 data); mean variance reduction {100*np.nanmean(VR[ok]):.0f}%')
print(f'resolved cells: {well.sum()} | temporal RMS of fault dv/v among resolved cells: '
      f'median {np.median(cell_var[well]):.3f}%  max {np.max(cell_var[well]):.3f}%')
# network-coherent fault index: mean over resolved cells, per month
fault_idx = np.nanmean(MF[well][:, ok], axis=0)
print(f'network fault index (mean resolved dv/v) range over time: [{np.nanmin(fault_idx):.3f}, {np.nanmax(fault_idx):.3f}] %')
np.savez(f'fault_tomography/inversion/fault_4d_{SUF}.npz', months=np.array(months), t=t.astype(str),
         lat=cells.cell_lat, lon=cells.cell_lon, nsta=cells.n_stations, depth=cells.depth_km,
         MF=MF, VR=VR, ND=ND, well=well, cell_var=cell_var, fault_idx=fault_idx, ok=ok)
print('wrote fault_tomography/inversion/fault_4d.npz')
