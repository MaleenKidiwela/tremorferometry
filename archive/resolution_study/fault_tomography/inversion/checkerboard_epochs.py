#!/usr/bin/env python
"""Multi-scale checkerboard on the HONEST calibrated operator, per epoch (Merlin ruling 4c communication product;
grounds the resolution result in the earlier large-checkerboard finding). For representative epochs, plant a
km-projected +/-1% checker at several sizes, forward through the calibrated fault+site operator, add REALISTIC
per-datum noise (time-varying residual sigma), invert at the frozen PSF-floor lambda (seeded, 10 realizations),
and score recovery on deep well-covered cells. Shows: large checkers recover, fine ones don't, and the
recoverable scale sharpens as the network grows. Outputs res_catalog/checkerboard_epochs.{csv,png}."""
import os
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = "fault_tomography/inversion/res_catalog"
LAM_RATIO, NN, LAM_F, DEEP = 0.05/0.40, 6, 2.069, 30.0   # LAM_F frozen by PSF-floor (sensitivity_atlas)
SIZES_KM = [35, 70, 140, 250]                             # fine -> coarse (250 ~ the old 3deg/200km scale)
EPOCHS = ["2011-06", "2016-06", "2021-06", "2025-06"]
NREAL = 10

d = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = d["G"].astype(np.float64); f = d["captured"].astype(np.float64); cxy = d["cxy"]; depth = d["depth_km"]
clat = d["cell_lat"]; clon = d["cell_lon"]; Mf = G.shape[1]; keep = f > 0
Gc = (G*f[:, None])[keep]; ptag = d["pair_tag"][keep]; pcell = d["pair_cell"][keep]; psig = d["pair_sigma"][keep].astype(float)
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
deep = depth > DEEP
_, nbr = cKDTree(cxy).query(cxy, k=NN+1)
L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T @ L

pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["site_mean"] = pm.groupby(["tag", "ym"]).dvv_month.transform("mean")
pm["resid"] = pm.dvv_month - pm.site_mean; pm["pair"] = pm.tag + "|" + pm.cell
tv = pm.groupby("pair").resid.apply(lambda r: (r-r.mean()).std() if len(r) >= 6 else np.nan)
TVM = float(np.nanmedian(tv)); sig_real = {p: (v if np.isfinite(v) and v > 0 else TVM) for p, v in tv.items()}
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]
pm["site"] = d["pair_site"][keep][pm.row.values]


def checker(size):
    m = np.sign(np.sin(np.pi*cxy[:, 0]/size) * np.sin(np.pi*cxy[:, 1]/size)) * 1.0   # +/-1%
    return m


def epoch_op(ym):
    g = pm[pm.ym == ym]; rows = g.row.values; nd = g.n_days.values.astype(float)
    sig = np.array([sig_real[p] for p in g.pair.values]) / np.sqrt(np.maximum(1.0, nd/30.0))
    tags = g.tag.values; ut = pd.unique(tags); tc = {t: k for k, t in enumerate(ut)}
    S = np.zeros((len(rows), len(ut)))
    for r, t in enumerate(tags):
        S[r, tc[t]] = 1.0
    A = np.hstack([Gc[rows], S]); w = 1.0/sig
    AtA = (A*w[:, None]).T @ (A*w[:, None])
    reg = np.zeros_like(AtA); reg[:Mf, :Mf] = LAM_F**2*LtL
    reg[np.arange(Mf), np.arange(Mf)] += 1e-4
    idx = np.arange(Mf, Mf+len(ut)); reg[idx, idx] += (LAM_F*LAM_RATIO)**2
    Minv = inv(AtA + reg)
    cells_hit = {c: n for c, n in g.groupby("cell").site.nunique().items()}   # sites per illuminated cell
    return A, w, S.shape[1], Minv, rows, cells_hit


# map cell-id string -> index (to mark illuminated)
cid = {f"{clat[i]:.2f}_{clon[i]:.2f}": i for i in range(Mf)}
res = []
grids = {}
for ym in EPOCHS:
    if not len(pm[pm.ym == ym]):
        continue
    A, w, nsite, Minv, rows, cells_hit = epoch_op(ym)
    ill = np.zeros(Mf, bool)
    for c, n in cells_hit.items():
        if c in cid and n >= 3:
            ill[cid[c]] = True
    mask = ill & deep
    AtW2 = (A.T * w**2)
    for size in SIZES_KM:
        mt = checker(size); dclean = Gc[rows] @ mt
        corrs = []; mrec_acc = np.zeros(Mf)
        for r in range(NREAL):
            rng = np.random.RandomState(1000*SIZES_KM.index(size)+r+hash(ym) % 100)
            dn = dclean + (1.0/w) * rng.randn(len(rows))
            x = Minv @ (AtW2 @ dn); mrec = x[:Mf]; mrec_acc += mrec
            if mask.sum() > 3:
                corrs.append(np.corrcoef(mt[mask], mrec[mask])[0, 1])
        c = float(np.nanmean(corrs)) if corrs else np.nan
        res.append((ym, nsite, size, mask.sum(), c))
        grids[(ym, size)] = (mt.copy(), (mrec_acc/NREAL).copy(), mask.copy())   # mean recovered over realizations
    print(f"{ym}: sites~{nsite} deep-illum {mask.sum()}  " + "  ".join(f"{s}km:{[r for r in res if r[0]==ym and r[2]==s][0][4]:.2f}" for s in SIZES_KM))

R = pd.DataFrame(res, columns=["ym", "nsite", "size_km", "n_deep_cells", "recovery_corr"])
R.to_csv(f"{OUT}/checkerboard_epochs.csv", index=False)

# --- summary curve: recoverable scale vs epoch ---
figc, axc = plt.subplots(figsize=(6.2, 4.6))
for ym in EPOCHS:
    s = R[R.ym == ym]
    if len(s):
        axc.plot(s.size_km, s.recovery_corr, "o-", label=f"{ym} ({int(s.nsite.iloc[0])} sites)")
axc.axhline(0.7, ls=":", color="k", lw=1); axc.set_xlabel("checker size (km)"); axc.set_ylabel("deep recovery corr")
axc.set_title("Recoverable deep scale vs epoch\n(large recovers, fine doesn't; sharpens with coverage)")
axc.legend(fontsize=8); axc.grid(alpha=.3); axc.set_ylim(-0.1, 1.05)
figc.tight_layout(); figc.savefig(f"{OUT}/checkerboard_epochs.png", dpi=130)

# --- CLASSIC checkerboard grid: INPUT vs RECOVERED at each scale, early vs late epoch (deep cells) ---
VIZ = [e for e in ["2011-06", "2025-06"] if any(k[0] == e for k in grids)]
ncol = 1 + len(VIZ)
fig, ax = plt.subplots(len(SIZES_KM), ncol, figsize=(3.4*ncol, 3.0*len(SIZES_KM)), squeeze=False)
lo, la = clon[deep], clat[deep]
def panel(a, vals, title):
    sc = a.scatter(lo, la, c=vals[deep], cmap="RdBu_r", vmin=-1, vmax=1, s=26, marker="s", edgecolors="none")
    a.set_title(title, fontsize=10); a.set_xticks([]); a.set_yticks([]); a.set_aspect("auto")
    return sc
for ri, size in enumerate(SIZES_KM):
    mt = grids[(VIZ[0] if VIZ else EPOCHS[0], size)][0]
    sc = panel(ax[ri][0], mt, (f"INPUT — {size} km checker" if ri == 0 else f"{size} km"))
    ax[ri][0].set_ylabel(f"{size} km", fontsize=11, rotation=90)
    for ci, ym in enumerate(VIZ):
        _, mrec, _ = grids[(ym, size)]
        c = R[(R.ym == ym) & (R.size_km == size)].recovery_corr.iloc[0]
        ns = int(R[R.ym == ym].nsite.iloc[0])
        panel(ax[ri][ci+1], mrec, (f"RECOVERED {ym[:4]} ({ns} sites)\nr={c:.2f}" if ri == 0 else f"{ym[:4]}  r={c:.2f}"))
fig.suptitle("Deep-interface checkerboard: input vs recovered at each scale (>30 km cells)", fontsize=12, y=1.005)
cax = fig.add_axes([1.01, 0.25, 0.015, 0.5]); plt.colorbar(sc, cax=cax, label="dbeta/beta (%)")
fig.tight_layout(); fig.savefig(f"{OUT}/checkerboard_grid.png", dpi=130, bbox_inches="tight")
print(f"\nwrote {OUT}/checkerboard_epochs.csv, checkerboard_epochs.png, checkerboard_grid.png")
print(R.pivot(index="ym", columns="size_km", values="recovery_corr").round(2).to_string())
