#!/usr/bin/env python
"""STEP 3-4 (Merlin recipe): per-MONTH model resolution matrix of the PRODUCTION fault+site estimator on the
CALIBRATED operator (Gc = unit-sum single-scatter rows x captured-fraction f), then aggregate to a
resolution-through-time curve. NO velocity inversion is run.

Estimator mirrored: m_est = M^-1 A^T W^2 d,  A=[Gc | S_site],  M = A^T W^2 A + blkdiag(lam_f^2 L^T L, lam_s^2 I),
W=diag(1/sigma_i), sigma_i = sigma_daily / sqrt(max(1, n_days/30))  (30-day overlapping stacks -> ~1 indep/mo).
Model resolution R = M^-1 (A^T W^2 A); posterior cov C = M^-1 (A^T W^2 A) M^-1.  Report fault block only.
lam re-selected ONCE by discrepancy principle (chi^2 ~ N) on the densest recent months, then FROZEN.

Outputs res_catalog/resolution_epochs.npz (months x cells: R_kk, psf_km, sigma_m, leak) + prints annual summary.
"""
import os
import numpy as np, pandas as pd
from numpy.linalg import inv

OUT = "fault_tomography/inversion/res_catalog"
LAM_RATIO = 0.05 / 0.40           # production lam_s/lam_f
MIN_DAYS, MIN_ROWS, MIN_SITES = 5, 50, 10
NN = 6                            # graph-Laplacian neighbours (matches forward.py)

d = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = d["G"].astype(np.float64); f = d["captured"].astype(np.float64)
cxy = d["cxy"]; depth = d["depth_km"]; cell_ids = d["cell"]
ptag = d["pair_tag"]; pcell = d["pair_cell"]; psite = d["pair_site"]; psig = d["pair_sigma"].astype(float)
Mf = G.shape[1]
keep = f > 0
Gc = (G * f[:, None])[keep]        # CALIBRATED rows (Merlin: production fix, not test-only)
ptag, pcell, psite, psig = ptag[keep], pcell[keep], psite[keep], psig[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
print(f"cells Mf={Mf} | calibrated pairs {Gc.shape[0]} (dropped {int((~keep).sum())} f=0)")

# graph Laplacian on fault cells (6-NN)
from scipy.spatial import cKDTree
_, nbr = cKDTree(cxy).query(cxy, k=NN + 1)
L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T @ L

pm = pd.read_parquet(f"{OUT}/pair_months.parquet")
pm = pm[pm.n_days >= MIN_DAYS].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]
pm = pm[pm.row >= 0]
months = sorted(pm.ym.unique())


def build_month(g):
    """rows, per-datum sigma, site design S, data d for one month-group g."""
    rows = g.row.values
    nd = g.n_days.values.astype(float)
    sig = psig[rows] / np.sqrt(np.maximum(1.0, nd / 30.0))
    tags = g.tag.values
    utag = pd.unique(tags); tcol = {t: k for k, t in enumerate(utag)}
    S = np.zeros((len(rows), len(utag)))
    for r, t in enumerate(tags):
        S[r, tcol[t]] = 1.0
    A = np.hstack([Gc[rows], S])
    return A, sig, g.dvv_month.values, len(utag)


def solve_month(g, lam_f, want_res=False):
    A, sig, dat, nsite = build_month(g)
    w2 = 1.0 / sig ** 2
    AtW2A = A.T * w2 @ A
    reg = np.zeros_like(AtW2A)
    reg[:Mf, :Mf] = lam_f ** 2 * LtL
    reg[np.arange(Mf), np.arange(Mf)] += 1e-4      # tiny ridge: stabilizes the Laplacian constant null-mode
    idx = np.arange(Mf, Mf + nsite)
    reg[idx, idx] += (lam_f * LAM_RATIO) ** 2
    Minv = inv(AtW2A + reg)
    m_est = Minv @ (A.T * w2 @ dat)
    resid = dat - A @ m_est
    chi2 = float(np.sum((resid / sig) ** 2)); N = len(dat)
    if not want_res:
        return chi2, N
    AtW2A_f = AtW2A                       # full normal matrix
    R = Minv @ AtW2A_f                    # model resolution
    C = Minv @ AtW2A_f @ Minv            # posterior covariance
    Rff = R[:Mf, :Mf]; Rfs = R[:Mf, Mf:]
    return dict(Rkk=np.diag(Rff).copy(), Rff=Rff, sigma_m=np.sqrt(np.clip(np.diag(C)[:Mf], 0, None)),
                leak=np.abs(Rfs).sum(1), nrow=N, nsite=nsite,
                nsite_distinct=g.site.nunique() if "site" in g else nsite)


# ---- lambda selection: discrepancy principle on the 3 densest months ----
pm = pm.merge(pd.DataFrame({"row": np.arange(len(psite)), "site": psite}), on="row", how="left")
dens = pm.groupby("ym").size().sort_values(ascending=False)
dense_months = dens.index[:3].tolist()
print(f"lambda selection on densest months: {dense_months} ({dens.iloc[:3].tolist()} rows)")
grid = np.logspace(-3, 0, 25)
best, curve = None, []
for lf in grid:
    rr = [solve_month(pm[pm.ym == m], lf) for m in dense_months]
    ratio = np.median([c / N for c, N in rr])
    curve.append((lf, ratio))
    if best is None or abs(np.log(ratio)) < abs(np.log(best[1])):
        best = (lf, ratio)
LAM_F = best[0]
print("  chi2/N vs lam_f:", "  ".join(f"{lf:.3f}:{r:.1f}" for lf, r in curve if lf in grid[::4]))
print(f"  -> chosen lam_f={LAM_F:.4f} (chi2/N={best[1]:.2f}), lam_s={LAM_F*LAM_RATIO:.4f} [FROZEN]")

# ---- per-month resolution matrix ----
Rkk = np.full((len(months), Mf), np.nan, np.float32)
PSF = np.full((len(months), Mf), np.nan, np.float32)
SIG = np.full((len(months), Mf), np.nan, np.float32)
LEAK = np.full((len(months), Mf), np.nan, np.float32)
meta = []
dist = np.sqrt(((cxy[:, None, :] - cxy[None, :, :]) ** 2).sum(-1))   # Mf x Mf cell-cell distance (km)
order = np.argsort(dist, axis=1)
for mi, m in enumerate(months):
    g = pm[pm.ym == m]
    if len(g) < MIN_ROWS or g.site.nunique() < MIN_SITES:
        meta.append((m, len(g), g.site.nunique(), 0)); continue
    R = solve_month(g, LAM_F, want_res=True)
    Rkk[mi] = R["Rkk"]; SIG[mi] = R["sigma_m"]; LEAK[mi] = R["leak"]
    ill = R["Rkk"] > 0.05
    for k in np.where(ill)[0]:
        w = np.abs(R["Rff"][k, order[k]]); cw = np.cumsum(w)
        if cw[-1] > 0:
            half = np.searchsorted(cw, 0.5 * cw[-1])
            PSF[mi, k] = dist[k, order[k, min(half, Mf - 1)]]
    meta.append((m, len(g), g.site.nunique(), int(ill.sum())))

meta = pd.DataFrame(meta, columns=["ym", "nrow", "nsite", "n_illum"])
np.savez_compressed(f"{OUT}/resolution_epochs.npz", months=np.array(months), Rkk=Rkk, PSF=PSF, SIG=SIG,
                    LEAK=LEAK, depth=depth, cell_lat=d["cell_lat"], cell_lon=d["cell_lon"], cxy=cxy,
                    lam_f=LAM_F, lam_s=LAM_F * LAM_RATIO, cell=cell_ids)
meta.to_csv(f"{OUT}/resolution_months_meta.csv", index=False)

# ---- annual aggregation (DEEP >30 km headline) ----
deep = depth > 30
SIGBAR = 0.5   # data-space detectability: interface dbeta/beta resolvable if sigma_m <= ~0.5% (Merlin ruling 2)
def resolved(mi):
    ok = (Rkk[mi] >= 0.3) & (PSF[mi] <= 50) & (SIG[mi] <= SIGBAR) & deep
    return int(np.nansum(ok)), np.nanmedian(PSF[mi][ok]) if np.nansum(ok) else np.nan
meta["ym_dt"] = pd.PeriodIndex(meta.ym, freq="M")
meta["year"] = meta.ym_dt.dt.year
rr = [resolved(mi) for mi in range(len(months))]
meta["deep_resolved"] = [a for a, _ in rr]; meta["deep_psf_med"] = [b for _, b in rr]
print("\n=== DEEP (>30km) RESOLUTION THROUGH TIME (annual) ===")
ann = meta.groupby("year").agg(months=("ym", "size"), med_sites=("nsite", "median"),
                               deep_resolved_med=("deep_resolved", "median"),
                               deep_resolved_max=("deep_resolved", "max"),
                               psf_km_med=("deep_psf_med", "median")).reset_index()
print(ann.to_string(index=False))
ann.to_csv(f"{OUT}/resolution_annual.csv", index=False)
print(f"\nlam_f={LAM_F:.4f} lam_s={LAM_F*LAM_RATIO:.4f} | sigma_m median (illum,deep,recent): "
      f"{np.nanmedian(SIG[-12:][:, deep]):.3f}% | wrote resolution_epochs.npz + annual csv")
