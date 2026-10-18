#!/usr/bin/env python
"""Shape program stage 2 (Merlin FE design): from data/shape_metrics_all.csv build
(1) per-metric fixed-effects solve: metric = station FE + b1*dist + b2*dep + medium(midpoint cell),
    ridge + Laplacian smoothing on a 0.2-deg midpoint-cell grid;
(2) mini-checkerboard identification test on the SAME design (inject known cell pattern + station FE,
    re-solve, report recovery corr overall and per margin third);
(3) A2-rescoped product: within-cell co-located source spread (centroid std above split-half floor).
Outputs: res_catalog/shape_maps.npz + figures/shape_maps.png. Maps from RESIDUAL design only
(anti-circularity: no raw station values are mapped)."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = pd.read_csv("data/shape_metrics_all.csv")
COORD = {r.sta: (r.lat, r.lon) for _, r in pd.read_csv("data/candidate_stations_post2020.csv").iterrows()}
COORD.update({"PGC": (48.6498, -123.4521), "SHB": (49.5990, -123.8770), "CLRS": (48.8203, -124.1309)})
D["stem"] = D.tag.str.replace(r"p\d+f\d+$", "", regex=True).str.upper()
D["sla"] = D.stem.map(lambda s: COORD.get(s, (np.nan, np.nan))[0])
D["slo"] = D.stem.map(lambda s: COORD.get(s, (np.nan, np.nan))[1])
D["mlat"] = (D.lat + D.sla) / 2; D["mlon"] = (D.lon + D.slo) / 2
D = D.dropna(subset=["mlat", "mlon"]).copy()          # stations without coords can't be mapped
D["cell"] = (D.mlat / 0.2).round().astype(int).astype(str) + "_" + (D.mlon / 0.2).round().astype(int).astype(str)

def solve(metric):
    ok = D.dropna(subset=[metric, "dist", "dep", "mlat", "mlon"]).copy()
    stations = sorted(ok.tag.unique()); cells = sorted(ok.cell.unique())
    si = {s: i for i, s in enumerate(stations)}; ci = {c: i for i, c in enumerate(cells)}
    nS, nC, nG = len(stations), len(cells), 2
    n = len(ok)
    rowsA, colsA, vals = [], [], []
    def add(r, c, v): rowsA.append(r); colsA.append(c); vals.append(v)
    for r in range(n):
        add(r, si[ok.tag.iloc[r]], 1.0)
        add(r, nS + 0, ok.dist.iloc[r] / 100)
        add(r, nS + 1, ok.dep.iloc[r] / 10)
        add(r, nS + nG + ci[ok.cell.iloc[r]], 1.0)
    from scipy.sparse import coo_matrix, identity, vstack as svstack
    from scipy.sparse.linalg import lsqr
    A = coo_matrix((vals, (rowsA, colsA)), shape=(n, nS + nG + nC)).tocsr()
    # Laplacian on cells (4-neighbor grid)
    cl = {c: tuple(map(int, c.split("_"))) for c in cells}
    Lr, Lc, Lv = [], [], []
    k = 0
    for c, (a, b) in cl.items():
        for d2 in ((a+1, b), (a, b+1)):
            c2 = f"{d2[0]}_{d2[1]}"
            if c2 in ci:
                Lr += [k, k]; Lc += [nS+nG+ci[c], nS+nG+ci[c2]]; Lv += [1.0, -1.0]; k += 1
    LAM_C, LAM_S = 1.0, 0.1
    L = coo_matrix((np.array(Lv)*LAM_C, (Lr, Lc)), shape=(k, nS+nG+nC)).tocsr()
    Rs = identity(nS + nG + nC, format="csr") * LAM_S
    y = ok[metric].values.astype(float)
    Afull = svstack([A, L, Rs])
    yfull = np.concatenate([y - y.mean(), np.zeros(k), np.zeros(nS + nG + nC)])
    x = lsqr(Afull, yfull, atol=1e-10, btol=1e-10, iter_lim=5000)[0]
    med = x[nS+nG:]
    pred = A @ x
    r2 = 1 - np.var(y - y.mean() - pred) / np.var(y)
    # mini-checkerboard: pattern on cells + random station FE, same design
    rng = np.random.RandomState(0)
    latc = np.array([cl[c][0]*0.2 for c in cells])
    patt = np.sign(np.sin(np.pi * latc / 1.0))            # ~1-deg checker in latitude
    x0 = np.zeros(nS+nG+nC); x0[:nS] = rng.randn(nS)*np.std(x[:nS]); x0[nS+nG:] = patt*np.std(med) if np.std(med)>0 else patt*0.1
    ysyn = A @ x0 + rng.randn(n)*np.std(y - y.mean() - pred)
    xs = lsqr(svstack([A, L, Rs]), np.concatenate([ysyn, np.zeros(k), np.zeros(nS+nG+nC)]),
              atol=1e-10, btol=1e-10, iter_lim=5000)[0]
    rec = np.corrcoef(x0[nS+nG:], xs[nS+nG:])[0, 1]
    recN = {}
    for nm, msk in (("south", latc < 44), ("central", (latc >= 44) & (latc < 47)), ("north", latc >= 47)):
        if msk.sum() > 10 and np.std(x0[nS+nG:][msk]) > 0:
            recN[nm] = float(np.corrcoef(x0[nS+nG:][msk], xs[nS+nG:][msk])[0, 1])
    return ok, cells, cl, med, r2, rec, recN

results = {}
for metric in ("decay", "cdr"):
    ok, cells, cl, med, r2, rec, recN = solve(metric)
    results[metric] = (ok, cells, cl, med, r2, rec, recN)
    print(f"{metric}: n={len(ok)} cells={len(cells)} R2={r2:.2f} | mini-checkerboard recovery "
          f"overall {rec:.2f} | by third {recN}", flush=True)

# A2r: within-cell source spread of centroid (co-located, same station)
grp = D.dropna(subset=["cent"]).groupby(["tag", "lat", "lon"])
sp = grp.cent.agg(["std", "count"])
sp = sp[sp["count"] >= 2].reset_index()
floor = np.nanmedian(np.abs(D.cent_o - D.cent_e)) / np.sqrt(2)
print(f"A2r: {len(sp)} co-located sets; median within-cell centroid std {sp['std'].median():.3f} Hz "
      f"(floor {floor:.3f})", flush=True)

fig, axs = plt.subplots(1, 3, figsize=(16, 6.5))
for ax, metric in zip(axs[:2], ("decay", "cdr")):
    ok, cells, cl, med, r2, rec, recN = results[metric]
    la = np.array([cl[c][0]*0.2 for c in cells]); lo = np.array([cl[c][1]*0.2 for c in cells])
    v = med - np.median(med); lim = np.percentile(np.abs(v), 95)
    sc = ax.scatter(lo, la, c=v, cmap="coolwarm", vmin=-lim, vmax=lim, s=22, marker="s")
    plt.colorbar(sc, ax=ax, shrink=0.7)
    ax.set_title(f"{metric} medium term (midpoint cells)\nR2 {r2:.2f}, checker recovery {rec:.2f}", fontsize=9)
sc = axs[2].scatter(sp.lon, sp.lat, c=sp["std"], cmap="magma_r", vmin=0, vmax=max(0.5, sp["std"].quantile(.95)), s=14)
plt.colorbar(sc, ax=axs[2], shrink=0.7, label="within-cell centroid std (Hz)")
axs[2].set_title(f"A2r: co-located SOURCE spread\n(floor {floor:.3f} Hz)", fontsize=9)
for ax in axs: ax.set_xlabel("lon"); ax.set_ylabel("lat")
fig.suptitle("Shape program stage 2 — fixed-effects medium maps + source-variation map "
             "(PROVISIONAL, pre-Merlin-audit)", fontsize=11)
fig.tight_layout()
fig.savefig("figures/shape_maps.png", dpi=130, bbox_inches="tight")
np.savez_compressed("fault_tomography/inversion/res_catalog/shape_maps.npz",
                    **{f"{m}_cells": np.array(results[m][1], object) for m in results},
                    **{f"{m}_med": results[m][3] for m in results},
                    a2r=sp.to_records(index=False), floor=floor)
print("wrote figures/shape_maps.png + res_catalog/shape_maps.npz")
