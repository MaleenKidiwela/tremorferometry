#!/usr/bin/env python
"""Phase C + A5: fault mesh, medium params, and the forward operator G_f for the borehole set,
then a Phase-D CHECKERBOARD resolution test (the right first milestone: validates the whole
kernel->forward->inverse chain and maps where the borehole geometry resolves the interface).

Model space = the fault cells from A2 (each a node at its lat/lon/Slab2-depth). For each datum
(family-cell a, station b) the row of G_f is the deep-source/surface-receiver kernel sourced at a,
received at b, evaluated at every fault node k and restricted to the interface (theory 02 eq.1,3):
    G_f[(a,b), k] = -K(X_k; X_a, R_b, t_ab) * h * A_k          (deep lobe sits at k=a: self-imaging)
Lapse time t_ab = hypocentral_dist(a,b)/beta + coda_offset.  D, beta from A4 (below).
"""
import os, sys
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
sys.path.insert(0, 'fault_tomography/kernels')
from kernel import kernel_diffusion, kernel_singlescatter

# ---- A4: medium parameters (literature-based first pass; refine from coda envelopes later) ----
BETA = 3.5          # S speed [km/s], crustal average
LSTAR = 40.0        # transport mean free path [km], typical crustal scattering
D = BETA * LSTAR / 3.0   # diffusivity [km^2/s]  (~47)
CODA_OFF = 2.5      # lapse offset into the coda window [s]
H_THICK = 3.0       # fault-zone thickness h [km] (scales amplitude, not pattern)
CODA_W1 = float(os.environ.get('CODA_W1', '1.0'))   # single-scatter kernel coda window (env-overridable;
CODA_W2 = float(os.environ.get('CODA_W2', '4.0'))   # set 2.0/4.0 to match the corrected 2-4 s data)

def project(lat, lon, lat0, lon0):
    x = (lon - lon0) * 111.0 * np.cos(np.radians(lat0))
    y = (lat - lat0) * 111.0
    return x, y

def build():
    PFX = os.environ.get('TOMO_PFX', 'fault_tomography/data_assembly/tomo')  # set TOMO_PFX=.../borehole for borehole-only
    cells = pd.read_csv(f'{PFX}_cells.csv')
    cat   = pd.read_csv(f'{PFX}_patch_catalog.csv')
    lat0, lon0 = cells.cell_lat.mean(), cells.cell_lon.mean()
    cells['x'], cells['y'] = project(cells.cell_lat.values, cells.cell_lon.values, lat0, lon0)
    cells = cells.reset_index(drop=True)
    cidx = {c: i for i, c in enumerate(cells.cell)}
    M = len(cells)
    Xk = cells[['x', 'y']].values; Zk = cells.depth_km.values
    Ak = (0.10 * 111.0) * (0.10 * 111.0 * np.cos(np.radians(lat0)))   # cell area [km^2]

    # data pairs (cell a, station b)
    cat = cat[cat.cell.isin(cidx)].copy()
    cat['ax'], cat['ay'] = project(cat.cell_lat.values, cat.cell_lon.values, lat0, lon0)
    cat['bx'], cat['by'] = project(cat.sta_lat.values, cat.sta_lon.values, lat0, lon0)
    cat['a'] = cat.cell.map(cidx)   # cat already carries depth_km (cell depth) from A2

    # assemble G_f (n_data x M)
    rows = []
    for _, r in cat.iterrows():
        s_h = np.array([r.ax, r.ay]); d = r.depth_km; r_h = np.array([r.bx, r.by])
        # 1-4 s coda window = EARLY coda -> single-scattering ellipsoid kernel (theory 01 eq.15),
        # NOT the diffusion kernel (which over-smears early coda by sqrt(2Dt)~40 km).
        K = kernel_singlescatter(Xk, Zk, s_h, d, r_h, BETA, CODA_W1, CODA_W2, ell=LSTAR)   # (M,), sums to 1
        rows.append(-K)
    Gf = np.array(rows)                                            # (Ndata, M)
    # row-normalise the captured on-fault sensitivity -> eta (theory 02 eq.2); keep eta as a weight proxy
    eta = np.abs(Gf).sum(1)
    return cells, cat, Gf, eta, M

def graph_laplacian(cells, k=6):
    from scipy.spatial import cKDTree
    P = cells[['x', 'y']].values
    tree = cKDTree(P); M = len(P)
    rows, cols, vals = [], [], []
    for i in range(M):
        dd, nn = tree.query(P[i], k + 1)
        for j in nn[1:]:
            rows += [i, i]; cols += [i, j]; vals += [1.0, -1.0]
    return sparse.coo_matrix((vals, (rows, cols)), shape=(M, M)).tocsr()

def checkerboard(cells, Gf, L, lam=0.3, wl_deg=0.4, noise=0.05):
    # plant alternating +/- 1% on a ~wl_deg checkerboard
    m_true = np.sign(np.sin(np.pi*cells.cell_lat/wl_deg) * np.sin(np.pi*cells.cell_lon/wl_deg)) * 0.01
    d = Gf @ m_true
    rng_scale = noise * np.std(d)
    d_noisy = d + rng_scale * np.sin(np.arange(len(d)))      # deterministic pseudo-noise (no RNG in env)
    # solve  [Gf; lam L] m = [d; 0]
    A = sparse.vstack([sparse.csr_matrix(Gf), lam * L])
    b = np.concatenate([d_noisy, np.zeros(L.shape[0])])
    m_rec = lsqr(A, b, atol=1e-6, btol=1e-6, iter_lim=2000)[0]
    return m_true, m_rec

if __name__ == '__main__':
    cells, cat, Gf, eta, M = build()
    print(f'model nodes (fault cells): {M} | data pairs (cell,station): {len(cat)} | beta={BETA} D={D:.0f} km^2/s')
    print(f'G_f shape {Gf.shape} | mean captured-sensitivity eta {eta.mean():.2f}')
    L = graph_laplacian(cells)
    m_true, m_rec = checkerboard(cells, Gf, L)
    well = cells.n_stations >= 3
    corr = np.corrcoef(m_true[well], m_rec[well])[0, 1]
    print(f'CHECKERBOARD recovery (>=3-station cells): corr(true,rec) = {corr:.2f}  '
          f'(n={well.sum()})')
    np.savez('fault_tomography/inversion/checkerboard.npz',
             lat=cells.cell_lat, lon=cells.cell_lon, nsta=cells.n_stations,
             depth=cells.depth_km, m_true=m_true, m_rec=m_rec)
    print('wrote fault_tomography/inversion/checkerboard.npz')
