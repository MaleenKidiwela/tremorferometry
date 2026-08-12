#!/usr/bin/env python
"""Illustrative: MONTHLY vs ANNUAL deep INDEX (spatial mean of deep delta-beta/beta). Shows why maps use annual
(monthly per-cell SNR is too low -> regularization art) while the spatial-mean INDEX survives monthly (averages
~440 deep cells). Same calibrated fault+site operator; simple per-pair deseason+common-window demean (illustrative,
not the full conditional version). Output res_catalog/monthly_vs_annual_index.png."""
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "fault_tomography/inversion/res_catalog"; LAM_F, LAM_RATIO = 0.231, 0.125
d = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = -d["G"].astype(float); f = d["captured"].astype(float); cxy = d["cxy"]; depth = d["depth_km"]
ptag = d["pair_tag"]; pcell = d["pair_cell"]; psite = d["pair_site"]; Mf = G.shape[1]
keep = f > 0; Gc = (G*f[:, None])[keep]; ptag, pcell, psite = ptag[keep], pcell[keep], psite[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}; deep = depth > 30
_, nbr = cKDTree(cxy).query(cxy, k=7); L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]: L[i, i] += 1; L[i, j] -= 1
LtL = L.T@L
pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]
pm["site"] = psite[pm.row.values]; pm["pair"] = pm.tag+"|"+pm.cell; pm["t"] = pd.PeriodIndex(pm.ym, freq="M"); pm["year"] = pm.t.dt.year
def harm(t): yr = 12.; return np.column_stack([np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr), np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
def dd(g):
    v = g.dvv_month.values.astype(float); t = (g.t-g.t.min()).apply(lambda x: x.n).values.astype(float)
    if len(v) >= 24: b, *_ = np.linalg.lstsq(np.column_stack([np.ones_like(t), harm(t)]), v, rcond=None); v = v-harm(t)@b[1:]
    inw = (g.year.values >= 2019) & (g.year.values <= 2024); base = v[inw].mean() if inw.sum() >= 6 else v.mean(); return v-base
pm = pm.sort_values(["pair", "t"]); pm["anom"] = np.concatenate([dd(g) for _, g in pm.groupby("pair")])
sig_pair = pm.groupby("pair").anom.std().to_dict(); TVM = float(np.nanmedian(list(sig_pair.values())))
rho = pm.groupby("pair").anom.apply(lambda a: float(np.clip(np.corrcoef(a.values[:-1], a.values[1:])[0, 1], 0, .95)) if len(a) > 8 and a.std() > 0 else 0.0).to_dict()
war = deep/deep.sum()
def solve_win(sub):
    a = sub.groupby(["row", "pair"]).agg(dvv=("anom", "mean"), nm=("anom", "size"), tag=("tag", "first"), site=("site", "first")).reset_index()
    if len(a) < 50 or a.site.nunique() < 10: return None
    nf = lambda p: (1-rho.get(p, 0))/(1+rho.get(p, 0))
    sig = np.array([sig_pair.get(p, TVM)/np.sqrt(max(1., nm*nf(p))) for p, nm in zip(a.pair, a.nm)])
    ut = pd.unique(a.tag.values); tc = {t: k for k, t in enumerate(ut)}; S = np.zeros((len(a), len(ut)))
    for r, t in enumerate(a.tag.values): S[r, tc[t]] = 1
    A = np.hstack([Gc[a.row.values], S]); w = 1/sig; Aw = A*w[:, None]; AtA = Aw.T@Aw; AtW2 = A.T*w**2
    reg = np.zeros_like(AtA); reg[:Mf, :Mf] = LAM_F**2*LtL; reg[np.arange(Mf), np.arange(Mf)] += 1e-4
    si = np.arange(Mf, Mf+len(ut)); reg[si, si] += (LAM_F*LAM_RATIO)**2; Minv = inv(AtA+reg)
    m = (Minv@(AtW2@a.dvv.values))[:Mf]; Cov = Minv@AtA@Minv
    return float(war@m), float(np.sqrt(max(0, war@Cov[:Mf, :Mf]@war)))
months = sorted(pm.ym.unique()); mi = []; msg = []
for mm in months:
    r = solve_win(pm[pm.ym == mm]); mi.append(r[0] if r else np.nan); msg.append(r[1] if r else np.nan)
yrs = list(range(2010, 2026)); yi = []; ysg = []
for y in yrs:
    r = solve_win(pm[pm.year == y]); yi.append(r[0] if r else np.nan); ysg.append(r[1] if r else np.nan)
mi, msg, yi, ysg = map(np.array, (mi, msg, yi, ysg))
print(f"index sigma: monthly median {np.nanmedian(msg):.3f}% | annual median {np.nanmedian(ysg):.3f}% -> annual {np.nanmedian(msg)/np.nanmedian(ysg):.1f}x tighter")
md = pd.PeriodIndex(months, freq="M").to_timestamp()
fig, ax = plt.subplots(figsize=(12, 4.6))
ax.fill_between(md, mi-2*msg, mi+2*msg, color="grey", alpha=.22, label="monthly ±2σ")
ax.plot(md, mi, color="#888", lw=.8, label="MONTHLY index")
ax.errorbar([pd.Timestamp(f"{y}-07-01") for y in yrs], yi, yerr=2*ysg, fmt="o-", color="#c0392b", capsize=3, lw=1.6, label="ANNUAL index ±2σ")
ax.axhline(0, color="k", lw=.5); ax.set_ylabel("deep network-mean δβ/β (%, model)"); ax.set_xlabel("year")
ax.set_title("Monthly vs annual deep INDEX — monthly σ is ~%.0f× larger (why MAPS use annual; index survives monthly)" % (np.nanmedian(msg)/np.nanmedian(ysg)))
ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(f"{OUT}/monthly_vs_annual_index.png", dpi=130)
print("wrote monthly_vs_annual_index.png")
