#!/usr/bin/env python
"""FINAL resolution-through-time deliverable (Merlin ruling 4, Gate outcome B).
The deep dv/v is a regional INDEX, not a cell map; report a SENSITIVITY/INFORMATION ATLAS THROUGH TIME:
 (b) HEADLINE  resolved-mode count vs year   = # singular values >1 of the whitened, site-projected deep operator
 (a) COMPANION index precision vs year       = posterior sigma of the deep-cell area-weighted mean (both sigma tiers)
 (c) SUPPLEMENT per-cell 2-sigma minimum-detectable-amplitude (MDA) maps per epoch
Operator = calibrated (Gc = unit-sum single-scatter x captured f) + per-station site terms; lambda by PSF-floor
(smallest lam_f s.t. median deep-cell PSF half-width >= 30 km = kernel footprint); realistic per-pair sigma =
TIME-VARYING residual std (Gate B), measurement sigma = detrended-MAD as the optimistic tier.
Outputs res_catalog/atlas_*.{npz,csv,png}."""
import os
import numpy as np, pandas as pd
from numpy.linalg import inv, svd
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "fault_tomography/inversion/res_catalog"
LAM_RATIO = 0.05 / 0.40
MIN_DAYS, MIN_ROWS, MIN_SITES = 5, 50, 10
NN, PSF_FLOOR, DEEP = 6, 30.0, 30.0

d = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = d["G"].astype(np.float64); f = d["captured"].astype(np.float64)
cxy = d["cxy"]; depth = d["depth_km"]; clat = d["cell_lat"]; clon = d["cell_lon"]
ptag = d["pair_tag"]; pcell = d["pair_cell"]; psig = d["pair_sigma"].astype(float)
Mf = G.shape[1]; keep = f > 0
Gc = (G * f[:, None])[keep]; ptag, pcell, psig = ptag[keep], pcell[keep], psig[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
deep = depth > DEEP
warea = deep.astype(float) / deep.sum()           # area-weight vector for the deep index (equal cells)

# graph Laplacian (6-NN)
_, nbr = cKDTree(cxy).query(cxy, k=NN + 1)
L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T @ L

pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= MIN_DAYS].copy()
# realistic per-pair sigma = time-varying residual std (Gate B); fallback = global median
pm["site_mean"] = pm.groupby(["tag", "ym"]).dvv_month.transform("mean")
pm["resid"] = pm.dvv_month - pm.site_mean
pm["pair"] = pm.tag + "|" + pm.cell
tv = pm.groupby("pair").resid.apply(lambda r: (r - r.mean()).std() if len(r) >= 6 else np.nan)
TV_MED = float(np.nanmedian(tv))
sig_real_pair = {p: (v if np.isfinite(v) and v > 0 else TV_MED) for p, v in tv.items()}
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]
pm = pm[pm.row >= 0]
pm["site"] = psig[pm.row.values] * 0 + d["pair_site"][keep][pm.row.values]  # distinct-site id
months = sorted(pm.ym.unique())
print(f"cells Mf={Mf} deep={int(deep.sum())} | pairs {Gc.shape[0]} | realistic sigma median {TV_MED:.3f}% | months {len(months)}")


def assemble(g, tier):
    rows = g.row.values; nd = g.n_days.values.astype(float)
    if tier == "meas":
        sig = psig[rows] / np.sqrt(np.maximum(1.0, nd / 30.0))
    else:                                            # realistic: time-varying residual per pair
        sig = np.array([sig_real_pair[p] for p in g.pair.values]) / np.sqrt(np.maximum(1.0, nd / 30.0))
    tags = g.tag.values; ut = pd.unique(tags); tc = {t: k for k, t in enumerate(ut)}
    S = np.zeros((len(rows), len(ut)))
    for r, t in enumerate(tags):
        S[r, tc[t]] = 1.0
    return rows, sig, S, len(ut)


def month_solve(g, lam_f, tier="real", want=("cov",)):
    rows, sig, S, nsite = assemble(g, tier)
    w = 1.0 / sig; A = np.hstack([Gc[rows], S]); Aw = A * w[:, None]
    AtA = Aw.T @ Aw
    reg = np.zeros_like(AtA); reg[:Mf, :Mf] = lam_f ** 2 * LtL
    reg[np.arange(Mf), np.arange(Mf)] += 1e-4
    idx = np.arange(Mf, Mf + nsite); reg[idx, idx] += (lam_f * LAM_RATIO) ** 2
    Minv = inv(AtA + reg)
    out = {}
    if "res" in want or "psf" in want:
        R = Minv @ AtA; out["Rff"] = R[:Mf, :Mf]
    if "cov" in want:
        Cov = Minv @ AtA @ Minv
        out["sigm"] = np.sqrt(np.clip(np.diag(Cov)[:Mf], 0, None))
        out["idx_var"] = float(warea @ Cov[:Mf, :Mf] @ warea)
    if "modes" in want:
        # whiten deep fault operator, project out site columns, SVD, count sv>1
        Gd = Aw[:, :Mf][:, deep]; Sw = Aw[:, Mf:]
        P = Sw @ np.linalg.pinv(Sw)                  # projector onto whitened site space
        Gperp = Gd - P @ Gd
        sv = svd(Gperp, compute_uv=False)
        out["nmode"] = int((sv > 1.0).sum()); out["sv"] = sv[:6]
    return out


# ---- lambda by PSF-floor on the densest month ----
dens = pm.groupby("ym").size().sort_values(ascending=False)
dm = pm[pm.ym == dens.index[0]]
def med_deep_psf(lam_f):
    r = month_solve(dm, lam_f, want=("res",))["Rff"]
    dist = np.sqrt(((cxy[:, None] - cxy[None]) ** 2).sum(-1)); order = np.argsort(dist, 1)
    psf = []
    for k in np.where(deep)[0]:
        w = np.abs(r[k, order[k]]); cw = np.cumsum(w)
        if r[k, k] > 0.05 and cw[-1] > 0:
            psf.append(dist[k, order[k, min(np.searchsorted(cw, 0.5 * cw[-1]), Mf - 1)]])
    return np.nanmedian(psf) if psf else np.nan
LAM_F = None
for lf in np.logspace(-3, 0.5, 20):
    if med_deep_psf(lf) >= PSF_FLOOR:
        LAM_F = lf; break
LAM_F = LAM_F or 1.0
print(f"lambda (PSF-floor >= {PSF_FLOOR}km): lam_f={LAM_F:.3f}  median deep PSF={med_deep_psf(LAM_F):.0f}km [FROZEN]")

# ---- per-month atlas ----
rec = []; MDA = np.full((len(months), Mf), np.nan, np.float32); svtop = np.zeros((len(months), 6))
for mi, m in enumerate(months):
    g = pm[pm.ym == m]
    if len(g) < MIN_ROWS or g.site.nunique() < MIN_SITES:
        rec.append((m, len(g), g.site.nunique(), 0, np.nan, np.nan)); continue
    R = month_solve(g, LAM_F, "real", want=("cov", "modes"))
    Rm = month_solve(g, LAM_F, "meas", want=("cov",))
    MDA[mi] = 2 * R["sigm"]
    svtop[mi, :len(R["sv"])] = R["sv"]
    rec.append((m, len(g), g.site.nunique(), R["nmode"], np.sqrt(R["idx_var"]), np.sqrt(Rm["idx_var"])))

M = pd.DataFrame(rec, columns=["ym", "nrow", "nsite", "nmode", "idx_sig_real", "idx_sig_meas"])
M["year"] = pd.PeriodIndex(M.ym, freq="M").year
ann = M.groupby("year").agg(months=("ym", "size"), sites=("nsite", "median"), nmode=("nmode", "median"),
                            idx_real_pct=("idx_sig_real", "median"), idx_meas_pct=("idx_sig_meas", "median")).reset_index()
np.savez_compressed(f"{OUT}/atlas.npz", months=np.array(months), MDA=MDA, svtop=svtop, clat=clat, clon=clon,
                    depth=depth, deep=deep, lam_f=LAM_F, nmode=M.nmode.values, idx_real=M.idx_sig_real.values,
                    idx_meas=M.idx_sig_meas.values, ym=M.ym.values)
M.to_csv(f"{OUT}/atlas_monthly.csv", index=False); ann.to_csv(f"{OUT}/atlas_annual.csv", index=False)
print("\n=== SENSITIVITY ATLAS THROUGH TIME (annual) ===")
print(ann.to_string(index=False))
print(f"\ncross-check: recent-year deep index sigma (realistic) = {ann.idx_real_pct.iloc[-3:].median():.3f}% "
      f"(network bound from ETS-null ~0.02%); measurement-floor tier = {ann.idx_meas_pct.iloc[-3:].median():.3f}%")

# ---- figures ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
ax[0].plot(ann.year, ann.nmode, "o-", color="#1a1a2e"); ax[0].set_title("(b) resolved deep modes (sv>1) vs year")
ax[0].set_xlabel("year"); ax[0].set_ylabel("# independent deep patterns"); ax[0].grid(alpha=.3)
ax2 = ax[0].twinx(); ax2.plot(ann.year, ann.sites, "s--", color="#d97706", alpha=.6); ax2.set_ylabel("median sites", color="#d97706")
ax[1].plot(ann.year, ann.idx_real_pct, "o-", label="realistic sigma", color="#dc2626")
ax[1].plot(ann.year, ann.idx_meas_pct, "s--", label="measurement floor", color="#2563eb")
ax[1].axhline(0.02, ls=":", color="k", lw=1, label="ETS-null network bound"); ax[1].set_yscale("log")
ax[1].set_title("(a) deep-index precision vs year"); ax[1].set_xlabel("year"); ax[1].set_ylabel("posterior sigma (%)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
# (c) MDA map for the most recent full year
lastfull = [i for i, m in enumerate(months) if m.startswith("2025")]
mda = np.nanmedian(MDA[lastfull], 0) if lastfull else np.nanmedian(MDA, 0)
sc = ax[2].scatter(clon[deep], clat[deep], c=np.clip(mda[deep], 0, 3), s=14, cmap="viridis_r", vmin=0, vmax=3)
ax[2].set_title("(c) deep 2-sigma min. detectable dbeta/beta (%), 2025"); ax[2].set_xlabel("lon"); ax[2].set_ylabel("lat")
plt.colorbar(sc, ax=ax[2], label="MDA (%)")
fig.tight_layout(); fig.savefig(f"{OUT}/atlas.png", dpi=130)
print(f"wrote {OUT}/atlas.npz, atlas_monthly.csv, atlas_annual.csv, atlas.png")
