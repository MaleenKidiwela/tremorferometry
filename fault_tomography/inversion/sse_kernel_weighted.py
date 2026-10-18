#!/usr/bin/env python
"""KERNEL-WEIGHTED model variables vs pair-level coda dv/v (2026-09-04). For each production (station, cell) pair, the
S-pinned differential residence kernel K_eff (res_catalog_kernels/G_diff_w24.npz, cells x weights, signed) is applied to
the model fields on the interface cells (trailing frame): apparent_pair(t) = sum_c K_eff[pair,c] * field(t, c), for
cumulative thrust slip, dilatation, Coulomb, shear, normal (10 km) and the six strain components (0.5-km offset via
strain_gf_offsets E[0] @ A_trail). Response: the pair's family-mean coda dv/v (2-4 s, -eps, %) per day, deseasoned.
Tests: (a) LEVEL correlation per pair vs rigid-shift null (fraction of pairs beyond the 95th pct, mean signed r, sign
consistency); (b) pooled OLS of pair 30-d dv/v change on the kernel-weighted 30-d changes (6 strain + slip; per-pair
demeaned), R^2 vs 200 rigid shifts of the amplitude functions; per-station coefficient sign map.
Outputs: res_catalog_sse/kernel_weighted_pairs.csv, printed summary.
"""
import numpy as np, pandas as pd, sys
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(91)
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src"); from sse_strain import strain as st
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells); gi = pd.Series(np.arange(G), index=days)
A_trail = np.nan_to_num(zc["A_trail"].astype(np.float32)); win = np.asarray((days >= "2011-01-01") & (days <= "2026-07-20")); month = days.month.values
zk = np.load("fault_tomography/inversion/res_catalog_kernels/G_diff_w24.npz", allow_pickle=True); K = -zk["G"].astype(np.float32)   # K_eff (signed residence), rows = pairs
ptag, pcell = zk["pair_tag"], zk["pair_cell"]; assert list(zk["cell"]) == list(cells); NP = len(ptag)
zg = np.load(f"{SD}/strain_gf_offsets.npz", allow_pickle=True); E0 = zg["E"][0].astype(np.float32)                        # (M, 6, 13) at 0.5 km above the interface
# kernel-weighted 13-vectors per pair for the 6 strain components: (NP, 6, 13)
KE = np.einsum("pc,cbk->pbk", K, E0)
FIELDS = {"cum slip (mm)": zc["slip_thrust_cum"].astype(np.float32) * 1e3, "slip in window (mm)": zc["slip_thrust_inc30"].astype(np.float32) * 1e3,
          "dilatation (1e-8)": zc["dilatation"].astype(np.float32) * 1e8, "Coulomb (kPa)": zc["coulomb"].astype(np.float32) / 1e3, "shear (kPa)": zc["shear_stress"].astype(np.float32) / 1e3, "normal (kPa)": zc["normal_stress"].astype(np.float32) / 1e3}
KW = {k: np.nan_to_num(v) @ K.T for k, v in FIELDS.items()}                                                                # (G, NP)
STR = {c: np.einsum("tk,pk->tp", A_trail, KE[:, b, :]) * 1e8 for b, c in enumerate(("exx", "eyy", "ezz", "exy", "exz", "eyz"))}
print(f"pairs {NP} | kernel row sums: median {np.median(K.sum(1)):+.3f}, negative rows {100*np.mean(K.sum(1)<0):.0f}%", flush=True)
# pair-level response: family-mean 2-4 s coda dv/v per day
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "cell", "date", "dvv"]); fd["v"] = -fd.dvv
pidx = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}; fd["pi"] = [pidx.get((t, c), -1) for t, c in zip(fd.tag, fd.cell)]; fd = fd[fd.pi >= 0]
g = fd.groupby(["pi", "date"]).v.mean().reset_index(); g["di"] = gi.reindex(g.date).values; g = g.dropna(subset=["di"]); del fd
R = np.full((G, NP), np.nan, np.float32); R[g.di.astype(int).values, g.pi.values] = g.v.values; del g
for m in range(1, 13): R[month == m] -= np.nanmean(R[month == m], 0)[None]
valid = np.isfinite(R[win]).mean(0) >= 0.3; print(f"pairs with >=30% coverage: {valid.sum()}", flush=True)
def corr_cols(Y, X):
    m = np.isfinite(Y) & np.isfinite(X) & win[:, None]; n = m.sum(0); Ym = np.where(m, Y, 0); Xm = np.where(m, X, 0)
    my = Ym.sum(0) / np.maximum(n, 1); mx = Xm.sum(0) / np.maximum(n, 1); yc = np.where(m, Y - my, 0); xc = np.where(m, X - mx, 0)
    r = (yc * xc).sum(0) / np.sqrt((yc ** 2).sum(0) * (xc ** 2).sum(0) + 1e-20); r[n < 500] = np.nan; return r
rows = []
print("\n(a) LEVEL correlation per pair (deseasoned pair dv/v vs kernel-weighted model field), 200 rigid shifts:")
for name, X in list(KW.items()) + [(f"strain {c} (1e-8)", v) for c, v in STR.items()]:
    r = corr_cols(R, X); nul = np.array([corr_cols(R, np.roll(X, rng.randint(365, G - 365), axis=0)) for _ in range(200)])
    hi = np.nanpercentile(np.abs(nul), 95, axis=0); beyond = np.nanmean(np.abs(r) > hi); pos = np.nanmean(r > 0)
    z = (r - np.nanmean(nul, 0)) / np.nanstd(nul, 0); zm = np.nanmean(z[valid]); zsd = np.nanstd(np.nanmean(nul, 1))
    print(f"  {name:22s}: mean r {np.nanmean(r[valid]):+.4f} | pairs beyond own 95% null {100*beyond:.1f}% (5% expected) | positive {100*pos:.0f}% | mean z {zm:+.2f}", flush=True)
    rows.append(dict(variable=name, mean_r=float(np.nanmean(r[valid])), frac_beyond95=float(beyond), frac_pos=float(pos), mean_z=float(zm)))
# (b) pooled OLS on 30-d changes: dy ~ 6 strain changes + kernel-weighted slip-in-window, per-pair demeaned
def d30(X): o = np.full_like(X, np.nan); o[30:] = X[30:] - X[:-30]; return o
dY = d30(R); preds = {c: d30(v) for c, v in STR.items()}; preds["slip"] = KW["slip in window (mm)"]
sel = valid; m = np.isfinite(dY[:, sel]) & win[:, None]
for v in preds.values(): m &= np.isfinite(v[:, sel])
def demean_cols(X): Xm = np.where(m, X[:, sel], 0); mu = Xm.sum(0) / np.maximum(m.sum(0), 1); return np.where(m, X[:, sel] - mu, 0)
y = demean_cols(dY).ravel()[m.ravel()]
def fit(pr):
    Xd = np.column_stack([demean_cols(v).ravel()[m.ravel()] for v in pr.values()]); b = np.linalg.lstsq(Xd, y, rcond=None)[0]; res = y - Xd @ b; return 1 - (res ** 2).sum() / (y ** 2).sum(), b
R2, beta = fit(preds); print(f"\n(b) pooled OLS (30-d changes, {m.sum():,} pair-days, {sel.sum()} pairs): R^2 = {R2:.5f}; coefficients {dict(zip(preds, np.round(beta, 4)))}", flush=True)
nul = []
for _ in range(60):
    sh = rng.randint(365, G - 365); pr = {c: np.roll(v, sh, axis=0) for c, v in preds.items()}; nul.append(fit(pr)[0])
print(f"    shift null R^2: mean {np.mean(nul):.5f} sd {np.std(nul):.5f} -> p = {np.mean(np.array(nul) >= R2):.3f}")
pd.DataFrame(rows).to_csv(f"{SD}/kernel_weighted_pairs.csv", index=False); print("wrote kernel_weighted_pairs.csv")
