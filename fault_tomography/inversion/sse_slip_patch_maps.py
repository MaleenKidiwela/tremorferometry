#!/usr/bin/env python
"""Maps of the space-time connected slip patches (sse_slip_patches.py catalog) starting on/after T0 (default 2011-01-01).
Subset CSV: res_catalog_sse/slip_patches_{THR}mm_post{T0year}.csv (ALL components after T0, flags kept).
Figure: one panel per MAJOR patch (Mw >= MWMIN, not minor, not broad_pulse): mesh triangles in grey, footprint triangles
filled with the cumulative thrust slip over the active span (raw_cum[end] - raw_cum[start-30], mm), centroid marked;
title = start..end, Mw, peak rate.  -> figures/sse_slip_patch_maps_{THR}mm_post{T0year}.png.  Env: THR, T0, MWMIN.
"""
import os, sys, numpy as np, pandas as pd
REPO = "/home/jovyan/cascadia-sse-strain"; sys.path.insert(0, f"{REPO}/src")
from sse_strain import load_fault_ascii, load_solution
SD = "fault_tomography/inversion/res_catalog_sse"; THR = float(os.environ.get("THR", "5")); T0 = pd.Timestamp(os.environ.get("T0", "2011-01-01")); MWMIN = float(os.environ.get("MWMIN", "6.2")); WIN = 30
GAP = int(os.environ.get("GAP", "0")); BRIDGE_KM = float(os.environ.get("BRIDGE_KM", "0")); SUF = f"_gap{GAP}_r{BRIDGE_KM:g}" if (GAP or BRIDGE_KM) else ""   # catalog variant (sse_slip_patches.py bridging)
cat = pd.read_csv(f"{SD}/slip_patches_{THR:g}mm{SUF}.csv", parse_dates=["start", "end", "peak_date"]); post = cat[cat.start >= T0].copy()
post.to_csv(f"{SD}/slip_patches_{THR:g}mm{SUF}_post{T0.year}.csv", index=False)
maj = post[(post.Mw >= MWMIN) & ~post.minor & ~post.broad_pulse].sort_values("start")
print(f"after {T0.date()}: {len(post)} components, {(~post.minor).sum()} non-minor, {len(maj)} major (Mw>={MWMIN:g}, not minor/broad-pulse)")
zl = np.load(f"{SD}/slip_patches_{THR:g}mm{SUF}.npz", allow_pickle=True); L = zl["label"].T; days = pd.to_datetime(zl["days"])     # L: (742, T)
fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt"); sol = load_solution(f"{REPO}/data/2026-08-18", fault); sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns(); A_raw = np.nan_to_num(np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True)["A_raw"].astype(float))
raw_cum = -(Ld @ A_raw.T) * 1e3                                                     # (742, T) mm thrust-positive
tri = np.loadtxt(f"{REPO}/data/cascadia_mesh.txt", comments="#"); verts = tri.reshape(-1, 3, 3)[:, :, :2]   # (742, 3, [lon, lat])
gi = pd.Series(np.arange(len(days)), index=days)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.tri as mtri; from matplotlib.collections import PolyCollection
nc = 6; ROWS = 7; PER = nc * ROWS; npages = int(np.ceil(len(maj) / PER))
FP, DS = {}, {}
for _, r in maj.iterrows():
    ts, te = int(gi[r.start]), int(gi[r.end]); fp = np.unique(np.where(L[:, ts:te + 1] == r.id)[0]); FP[r.id] = fp; DS[r.id] = raw_cum[fp, te] - raw_cum[fp, max(ts - WIN, 0)]
vmax = float(np.ceil(np.percentile(np.concatenate(list(DS.values())), 98) / 5) * 5)
cen = fault.centroid_llh; trg = mtri.Triangulation(cen[:, 0], cen[:, 1])                          # slab-depth contours from triangle centroids
base = f"figures/sse_slip_patch_maps_{THR:g}mm{SUF}_post{T0.year}" + ("" if MWMIN == 6.2 else f"_Mw{MWMIN:g}"); outs = []
for pg in range(npages):
    sub = maj.iloc[pg * PER:(pg + 1) * PER]; n = len(sub); nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3.1 * nc, 4.6 * nr), constrained_layout=True); axs = np.atleast_2d(axs)
    for k, (_, r) in enumerate(sub.iterrows()):
        a = axs[k // nc, k % nc]; fp, ds = FP[r.id], DS[r.id]
        a.add_collection(PolyCollection(verts, facecolor="#e2e8f0", edgecolor="none"))
        a.tricontour(trg, -cen[:, 2], levels=[20, 30, 40, 50], colors="#94a3b8", linewidths=0.5)
        pc = PolyCollection(verts[fp], array=ds, cmap="YlOrRd", edgecolor="face", lw=0.3); pc.set_clim(0, vmax); a.add_collection(pc)
        a.plot(r.lon_centroid, r.lat_centroid, "k+", ms=8, mew=1.3)
        a.set_xlim(-127.3, -121.3); a.set_ylim(39.7, 50.3); a.set_aspect(1 / np.cos(np.radians(45)))
        a.set_xticks([-126, -124, -122]); a.set_yticks(range(40, 51, 2)); a.tick_params(labelsize=7, length=2); a.grid(alpha=.25, lw=.5)
        if k % nc: a.set_yticklabels([])
        if k // nc < nr - 1 and k + nc < n: a.set_xticklabels([])
        a.set_title(f"{r.start:%Y-%m-%d} → {r.end:%Y-%m-%d}\nMw {r.Mw:.2f} · peak {r.peak_mm_per30d:.0f} mm/30 d · {r.duration_d} d", fontsize=8, pad=4)
        a.text(0.03, 0.02, f"#{r.id}", transform=a.transAxes, fontsize=7, color="#475569")
    for k in range(n, nr * nc): axs[k // nc, k % nc].axis("off")
    cb = fig.colorbar(pc, ax=axs.ravel().tolist(), shrink=0.45, pad=0.01, aspect=35); cb.set_label(f"cumulative thrust slip over the active span (mm; colour scale capped at the 98th percentile, {vmax:g} mm)", fontsize=9)
    fig.suptitle(f"SSE model slip patches > {THR:g} mm / 30 d, space-time connected{f' (gaps <= {GAP} d / {BRIDGE_KM:g} km bridged)' if SUF else ''}, starting on/after {T0:%Y-%m-%d} — Mw ≥ {MWMIN:g}, not minor / broad pulse ({len(maj)}{f', page {pg + 1}/{npages}' if npages > 1 else ''}); grey contours = slab depth 20-50 km", fontsize=11)
    out = base + (f"_p{pg + 1}" if npages > 1 else "") + ".png"; fig.savefig(out, dpi=200); plt.close(fig); outs.append(out)
print(f"wrote {', '.join(outs)}, {SD}/slip_patches_{THR:g}mm{SUF}_post{T0.year}.csv")
print(maj[["id","start","end","duration_d","lat_min","lat_max","area_km2","peak_mm_per30d","Mw"]].assign(start=maj.start.dt.date, end=maj.end.dt.date).round(2).to_string(index=False))
