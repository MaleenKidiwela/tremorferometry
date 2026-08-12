#!/usr/bin/env python
"""CORRECTED interface checkerboard (Merlin audit): the earlier PSF-floor lam=2.069 was transplanted from the
atlas's per-pair-whitened monthly system into the checkerboards -> over-smoothed. Fixes:
 - ORACLE lam per (epoch, tier, scale): sweep lam, report the MAX recovery (a resolution test reports the best
   achievable under the estimator class; production lam is re-tunable and lives only in sensitivity_atlas.py).
 - THREE noise tiers: noise-free (pure geometry ceiling), sigma=0.13% (measurement floor), sigma=0.25% (total
   time-varying residual = detectability). Consistent whitening w=1/sigma.
 - Scales {44,70,100,140,180,250} km; report the CROSSING scale (corr first >0.7) per epoch per tier + the
   noise-free ceiling. Also print the data-space signal std per checker (anchors +/-1% to the observed 0.043% ETS).
Deep well-covered cells (>30 km, all-time >=3 sites) so epochs are compared on one cell set.
Outputs res_catalog/checkerboard_corrected.{csv,png}."""
import os
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = "fault_tomography/inversion/res_catalog"
SCALES = [44, 70, 100, 140, 180, 250]
TIERS = [("noise-free", 0.0), ("meas 0.13%", 0.13), ("real 0.25%", 0.25)]
EPOCHS = [("2011", "2011-06"), ("2016", "2016-06"), ("2021", "2021-06"), ("2025", "2025-06"), ("all", None)]
LAMS = np.logspace(-2, 1.0, 12)
NREAL = 8

d = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = d["G"].astype(np.float64); f = d["captured"].astype(np.float64); cxy = d["cxy"]; depth = d["depth_km"]
ptag = d["pair_tag"]; pcell = d["pair_cell"]; Mf = G.shape[1]; keep = f > 0
Gc = (G*f[:, None])[keep]; ptag = ptag[keep]; pcell = pcell[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
cells = pd.read_csv(f"{OUT}/cells.csv").sort_values("cell").reset_index(drop=True)
mask = (depth > 30) & (cells.n_sites.values >= 3)
_, nbr = cKDTree(cxy).query(cxy, k=7); L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T @ L
allut = pd.unique(ptag); tcall = {t: k for k, t in enumerate(allut)}
pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]


def checker(size):
    return np.sign(np.sin(np.pi*cxy[:, 0]/size)*np.sin(np.pi*cxy[:, 1]/size))*1.0


def epoch_rows(ym):
    if ym is None:
        return np.arange(len(ptag)), ptag
    g = pm[pm.ym == ym]
    return g.row.values, ptag[g.row.values]


sig_data = {s: float(np.std(Gc @ checker(s))) for s in SCALES}   # checker is +/-1 (% dbeta/beta); Gc dimensionless -> d in %
print("data-space signal std per +/-1% checker: " + " ".join(f"{s}km:{sig_data[s]:.3f}%" for s in SCALES))
print("(anchor: observed mirror-corrected ETS transient ~0.043% -> +/-1% interface checker is the observed class)\n")

rows_out = []
for ename, ym in EPOCHS:
    rws, tags = epoch_rows(ym)
    if len(rws) < 50:
        continue
    ut = pd.unique(tags); tc = {t: k for k, t in enumerate(ut)}
    S = np.zeros((len(rws), len(ut)))
    for r, t in enumerate(tags):
        S[r, tc[t]] = 1.0
    A = np.hstack([Gc[rws], S]); nsite = len(ut)
    for tname, sig in TIERS:
        w = 1.0/max(sig, 0.13)                      # noise-free uses measurement-scale weighting
        Aw = A*w; AtA = Aw.T@Aw; AtW2 = A.T*w**2
        best = {s: -1.0 for s in SCALES}
        for lam in LAMS:
            reg = np.zeros_like(AtA); reg[:Mf, :Mf] = lam**2*LtL
            reg[np.arange(Mf), np.arange(Mf)] += 1e-6
            si = np.arange(Mf, Mf+nsite); reg[si, si] += (lam*0.125)**2
            Minv = inv(AtA + reg)
            for s in SCALES:
                mt = checker(s); dclean = Gc[rws] @ mt; cs = []
                for r in range(NREAL):
                    rng = np.random.RandomState(r + s + int(1000*sig))
                    dn = dclean + (sig*rng.randn(len(rws)) if sig > 0 else 0.0)
                    m = (Minv @ (AtW2 @ dn))[:Mf]
                    cs.append(np.corrcoef(mt[mask], m[mask])[0, 1])
                    if sig == 0:
                        break                       # noise-free is deterministic
                best[s] = max(best[s], float(np.nanmean(cs)))
        cross = next((s for s in SCALES if best[s] >= 0.7), None)
        rows_out.append(dict(epoch=ename, nsite=nsite, tier=tname, sigma=sig,
                             cross_km=cross, **{f"r{s}": round(best[s], 2) for s in SCALES}))
        print(f"{ename:>4} {tname:>11} sites={nsite:3d}  " + " ".join(f"{s}:{best[s]:.2f}" for s in SCALES) +
              f"  cross>=0.7: {cross}km")

R = pd.DataFrame(rows_out); R.to_csv(f"{OUT}/checkerboard_corrected.csv", index=False)

# figure: crossing scale vs year (2 noise tiers) + recovery-vs-scale for the recent epoch
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
yrs = [e for e, _ in EPOCHS if e != "all"]
for tname, sig in TIERS[1:]:
    cr = [R[(R.epoch == e) & (R.tier == tname)].cross_km.iloc[0] for e in yrs]
    cr = [c if c is not None else 300 for c in cr]
    ax[0].plot(yrs, cr, "o-", label=f"{tname}")
ax[0].axhline(44, ls=":", color="green", lw=1.2, label="geometry limit (noise-free = 44 km cell)")
ax[0].set_ylabel("finest resolved deep scale (km, corr>0.7)"); ax[0].set_xlabel("epoch"); ax[0].invert_yaxis()
ax[0].set_title("Resolution length through time (corrected)\nlower=finer; noise, not geometry, is the limit")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
rec = R[R.epoch == "2025"]
for tname, sig in TIERS:
    s = rec[rec.tier == tname]
    if len(s):
        ax[1].plot(SCALES, [s[f"r{k}"].iloc[0] for k in SCALES], "o-", label=tname)
ax[1].axhline(0.7, ls=":", color="k", lw=1); ax[1].set_xlabel("checker size (km)"); ax[1].set_ylabel("deep recovery corr")
ax[1].set_title("2025 interface recovery vs scale (oracle lam)\nnoise-free=geometry ceiling"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3); ax[1].set_ylim(-.05, 1.05)
fig.tight_layout(); fig.savefig(f"{OUT}/checkerboard_corrected.png", dpi=130)
print(f"\nwrote {OUT}/checkerboard_corrected.{{csv,png}}")
