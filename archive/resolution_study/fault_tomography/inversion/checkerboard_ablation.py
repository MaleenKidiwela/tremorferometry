#!/usr/bin/env python
"""Apples-to-apples: run the checkerboard on the SAME (new, full) network under the OLD regime vs the HONEST
regime, to show the recovery drop is the METHOD corrections, not the network. All-time coverage (like the old
static test), 2-D interface, deep well-covered cells (>30 km, >=3 sites).
  OLD   = unit-sum kernel, NO site terms, inverse-crime noise (5% of signal), lam=0.3 (the forward.py recipe).
  +noise= OLD but with REAL absolute noise (0.25% dv/v).
  +cap  = capture-weighted rows (calibrated G), real noise.
  HONEST= capture-weighted + per-station site terms + real noise + PSF-floor lam.
"""
import os
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree

OUT = "fault_tomography/inversion/res_catalog"
d = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = d["G"].astype(np.float64); f = d["captured"].astype(np.float64); cxy = d["cxy"]; depth = d["depth_km"]
ptag = d["pair_tag"]; psite = d["pair_site"]; Mf = G.shape[1]
keep = f > 0; G = G[keep]; f = f[keep]; ptag = ptag[keep]
cells = pd.read_csv(f"{OUT}/cells.csv").sort_values("cell").reset_index(drop=True)
deepwell = (depth > 30) & (cells.n_sites.values >= 3)
_, nbr = cKDTree(cxy).query(cxy, k=7)
L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T @ L
ut = pd.unique(ptag); tc = {t: k for k, t in enumerate(ut)}
S = np.zeros((len(ptag), len(ut)))
for r, t in enumerate(ptag):
    S[r, tc[t]] = 1.0
SIG_REAL = 0.25   # % dv/v, the time-varying residual floor (Gate B)


def checker(size):
    return np.sign(np.sin(np.pi*cxy[:, 0]/size)*np.sin(np.pi*cxy[:, 1]/size))*1.0


def run(size, regime, nreal=8):
    mt = checker(size); corrs = []
    if regime == "OLD":
        Gop, lam, site, invcrime = G, 0.3, False, True
    elif regime == "+noise":
        Gop, lam, site, invcrime = G, 0.3, False, False
    elif regime == "+cap":
        Gop, lam, site, invcrime = G*f[:, None], 2.069, False, False
    else:  # HONEST
        Gop, lam, site, invcrime = G*f[:, None], 2.069, True, False
    dclean = Gop @ mt
    if site:
        A = np.hstack([Gop, S]); reg = np.zeros((A.shape[1],)*2); reg[:Mf, :Mf] = lam**2*LtL
        reg[np.arange(Mf), np.arange(Mf)] += 1e-4
        si = np.arange(Mf, Mf+len(ut)); reg[si, si] += (lam*0.125)**2
        w = np.full(len(dclean), 1.0/SIG_REAL)
        M = (A*w[:, None]).T @ (A*w[:, None]) + reg; Minv = inv(M)
    else:
        reg = lam**2*LtL + 1e-4*np.eye(Mf); Minv = inv(Gop.T@Gop + reg)
    for r in range(nreal):
        rng = np.random.RandomState(r + int(size))
        nstd = 0.05*np.std(dclean) if invcrime else SIG_REAL
        dn = dclean + nstd*rng.randn(len(dclean))
        if site:
            m = (Minv @ ((A*w[:, None]).T @ (w*dn)))[:Mf]
        else:
            m = Minv @ (Gop.T @ dn)
        corrs.append(np.corrcoef(mt[deepwell], m[deepwell])[0, 1])
    return float(np.nanmean(corrs))


print(f"deep well-covered cells: {int(deepwell.sum())} | pairs {len(ptag)} | sites {len(ut)}")
print(f"{'scale':>8} | {'OLD(unit-sum,no-site,invcrime)':>30} | {'+real noise':>12} | {'+capture-wt':>12} | {'HONEST(full)':>12}")
rows = []
for size in [44, 250]:
    r = {k: run(size, k) for k in ["OLD", "+noise", "+cap", "HONEST"]}
    rows.append(dict(scale_km=size, **r))
    print(f"{size:>6}km | {r['OLD']:>30.2f} | {r['+noise']:>12.2f} | {r['+cap']:>12.2f} | {r['HONEST']:>12.2f}")
pd.DataFrame(rows).to_csv(f"{OUT}/checkerboard_ablation.csv", index=False)
print("\n-> the NEW network reproduces ~old recovery UNDER OLD ASSUMPTIONS; the drop is the 3 honest corrections.")
