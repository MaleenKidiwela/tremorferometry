#!/usr/bin/env python
"""CORRECTED input-vs-recovered checkerboard grid (oracle lambda, replaces the buggy frozen-lambda grid).
Rows = scales; cols = INPUT | recovered 2011 | recovered 2025. Measurement-noise (resolution) tier sigma=0.13%.
Per (epoch, scale): sweep lambda, pick max-recovery, show the mean recovered map at that lambda. Deep cells >30 km.
"""
import os
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = "fault_tomography/inversion/res_catalog"
SCALES = [44, 70, 100, 140, 250]; SIG = 0.13; LAMS = np.logspace(-2, 1, 10); NREAL = 6
VIZ = ["2011-06", "2025-06"]
d = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = d["G"].astype(np.float64); f = d["captured"].astype(np.float64); cxy = d["cxy"]; depth = d["depth_km"]
clat = d["cell_lat"]; clon = d["cell_lon"]; ptag = d["pair_tag"]; pcell = d["pair_cell"]; Mf = G.shape[1]
keep = f > 0; Gc = (G*f[:, None])[keep]; ptag, pcell = ptag[keep], pcell[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
cells = pd.read_csv(f"{OUT}/cells.csv").sort_values("cell").reset_index(drop=True)
mask = (depth > 30) & (cells.n_sites.values >= 3); deep = depth > 30
_, nbr = cKDTree(cxy).query(cxy, k=7); L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T @ L
pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]

def checker(s): return np.sign(np.sin(np.pi*cxy[:, 0]/s)*np.sin(np.pi*cxy[:, 1]/s))
grids = {}
for ym in VIZ:
    g = pm[pm.ym == ym]; rws = g.row.values; tags = g.tag.values
    ut = pd.unique(tags); tc = {t: k for k, t in enumerate(ut)}
    S = np.zeros((len(rws), len(ut)))
    for r, t in enumerate(tags):
        S[r, tc[t]] = 1.0
    A = np.hstack([Gc[rws], S]); w = 1.0/SIG; Aw = A*w; AtA = Aw.T@Aw; AtW2 = A.T*w**2
    for s in SCALES:
        mt = checker(s); dclean = Gc[rws] @ mt; best = (-1, None)
        for lam in LAMS:
            reg = np.zeros_like(AtA); reg[:Mf, :Mf] = lam**2*LtL; reg[np.arange(Mf), np.arange(Mf)] += 1e-6
            si = np.arange(Mf, Mf+len(ut)); reg[si, si] += (lam*0.125)**2
            Minv = inv(AtA + reg); acc = np.zeros(Mf); cs = []
            for r in range(NREAL):
                rng = np.random.RandomState(r+s); dn = dclean + SIG*rng.randn(len(rws))
                m = (Minv @ (AtW2 @ dn))[:Mf]; acc += m; cs.append(np.corrcoef(mt[mask], m[mask])[0, 1])
            c = float(np.nanmean(cs))
            if c > best[0]: best = (c, acc/NREAL)
        grids[(ym, s)] = (mt, best[1], best[0])
    print(f"{ym}: " + "  ".join(f"{s}km:{grids[(ym,s)][2]:.2f}" for s in SCALES))

fig, ax = plt.subplots(len(SCALES), 3, figsize=(10.2, 3.0*len(SCALES)), squeeze=False)
lo, la = clon[deep], clat[deep]
def panel(a, vals, title):
    sc = a.scatter(lo, la, c=vals[deep], cmap="RdBu_r", vmin=-1, vmax=1, s=26, marker="s"); a.set_title(title, fontsize=10); a.set_xticks([]); a.set_yticks([]); return sc
for ri, s in enumerate(SCALES):
    sc = panel(ax[ri][0], grids[(VIZ[0], s)][0], (f"INPUT — {s} km" if ri == 0 else f"{s} km"))
    for ci, ym in enumerate(VIZ):
        mt, mrec, c = grids[(ym, s)]
        panel(ax[ri][ci+1], mrec, (f"RECOVERED {ym[:4]}\nr={c:.2f}" if ri == 0 else f"{ym[:4]}  r={c:.2f}"))
fig.suptitle("Deep-interface checkerboard: input vs recovered (CORRECTED, oracle λ, σ=0.13% resolution tier)", fontsize=12, y=1.004)
cax = fig.add_axes([1.01, 0.3, 0.015, 0.4]); plt.colorbar(sc, cax=cax, label="δβ/β (%)")
fig.tight_layout(); fig.savefig(f"{OUT}/checkerboard_grid.png", dpi=130, bbox_inches="tight")
print("wrote (overwrote) checkerboard_grid.png with oracle-lambda recovery")
