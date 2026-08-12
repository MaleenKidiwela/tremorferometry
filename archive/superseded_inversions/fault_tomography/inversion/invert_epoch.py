#!/usr/bin/env python
"""Phase D (real data): joint fault+site inversion of one epoch of borehole dv/v.
   d = G_f m_f + G_s m_s + eps   (theory 02 eq.6, fault-vs-site split)
m_f = delta-beta/beta on the interface (target); m_s = per-station shallow site terms (absorb the
common-mode-over-families seasonal/site signal so it does not contaminate the fault).
Sign: our 'dvv' column IS delta-v/v, and dvv_ab = integral K (dbeta/beta) => G_f = +K (= -Gf from forward.py).
"""
import sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
sys.path.insert(0, 'fault_tomography/inversion')
from forward import build, graph_laplacian

cells, cat, Gf, eta, M = build()
Gf = -Gf                                         # +K : d(dvv) = Gf @ (dbeta/beta)
cat = cat.reset_index(drop=True)
stations = sorted(cat.station.unique()); sidx = {s: i for i, s in enumerate(stations)}; B = len(stations)
Gs_full = sparse.csr_matrix((np.ones(len(cat)), (np.arange(len(cat)), cat.station.map(sidx).values)),
                            shape=(len(cat), B))
L = graph_laplacian(cells)

ten = pd.read_parquet('fault_tomography/data_assembly/tomo_dvv_monthly.parquet')
# choose the best-covered month that also spans many stations
cov = ten.groupby('ym').agg(n=('dvv','size'), nsta=('station','nunique'))
E = cov[(cov.nsta>=15)].sort_values('n').index[-1]
te = ten[ten.ym==E][['cell','station','dvv']]

# align epoch data to the (cell,station) rows of cat
key = cat[['cell','station']].copy(); key['row'] = np.arange(len(cat))
m = key.merge(te, on=['cell','station'])
rows = m.row.values
d = m.dvv.values * 100.0                          # percent
A_f = Gf[rows]; A_s = Gs_full[rows]
print(f'epoch {E}: {len(d)} (cell,station) data over {m.station.nunique()} stations, {m.cell.nunique()} cells')

# joint inverse: [A_f A_s; lam_f L 0; 0 lam_s I] [m_f; m_s] = [d; 0; 0]
lam_f, lam_s = 0.4, 0.05
top = sparse.hstack([sparse.csr_matrix(A_f), A_s])
regf = sparse.hstack([lam_f * L, sparse.csr_matrix((M, B))])
regs = sparse.hstack([sparse.csr_matrix((B, M)), lam_s * sparse.identity(B)])
Aaug = sparse.vstack([top, regf, regs]).tocsr()
baug = np.concatenate([d, np.zeros(M), np.zeros(B)])
sol = lsqr(Aaug, baug, atol=1e-7, btol=1e-7, iter_lim=4000)[0]
m_f = sol[:M]; m_s = sol[M:]

# variance reduction
pred = A_f @ m_f + (A_s @ m_s)
vr = 1 - np.var(d - pred)/np.var(d)
print(f'variance reduction {100*vr:.0f}%  | fault m_f rms {np.std(m_f):.3f}%  site m_s rms {np.std(m_s):.3f}%')
print(f'site terms range [{m_s.min():.3f}, {m_s.max():.3f}] %  (these absorb the common-mode/site signal)')

np.savez('fault_tomography/inversion/epoch_map.npz', epoch=E,
         lat=cells.cell_lat, lon=cells.cell_lon, nsta=cells.n_stations, depth=cells.depth_km,
         m_f=m_f, stations=np.array(stations), m_s=m_s, vr=vr)
print('wrote fault_tomography/inversion/epoch_map.npz')
