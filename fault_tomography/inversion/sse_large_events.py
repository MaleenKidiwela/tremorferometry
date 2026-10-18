#!/usr/bin/env python
"""LARGE-EVENT analysis (user 2026-09-04: large SSEs may matter even if every SSE doesn't).
For every model episode (all 32 ok, all segments): per-episode footprint contrast ON[+15,+60]-PRE[-90,-30] and its own
z-score against 400 rigid shifts of the SAME footprint (per-episode null), in three products (interface deseasoned,
family 2-4 s, family 3-5 s). Then: (i) response vs episode size (footprint slip sum, segment slip increment, tremor
count, footprint cells) — Spearman with a permutation p; (ii) composite of the LARGEST quartile (by footprint slip sum)
and of the top-5 northern ETS, with Null A/B; (iii) lag composites for the top quartile; (iv) delta-beta/beta 30-d
change maps at ON for the three largest northern events next to the slip map.
Outputs: res_catalog_sse/large_events.csv, figures/sse_large_events.png, figures/sse_large_event_maps.png
"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(71)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
clat, clon, dep = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"]; gi = pd.Series(np.arange(G), index=days); inc = zc["slip_thrust_inc30"].astype(float) * 1e3
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok].reset_index(drop=True); zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
E0, E1 = int(gi[pd.Timestamp("2011-04-01")]), int(gi[pd.Timestamp("2026-04-01")]); span = E1 - E0; wrap = lambda i: E0 + (i - E0) % span
EP = [dict(key=r.key, seg=r.seg, ion=int(gi[pd.Timestamp(r.fp_onset)]), ipeak=int(r.ipeak), fp=zf[f"fp_{r.key}"], size=r.fp_slip_sum_mm, seginc=r.seg_inc_mm, tremor=r.tremor30_max, nfp=r.n_fp) for _, r in okep.iterrows()]
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); md = pd.to_datetime(zp["days"]); Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=md).reindex(days).values
F24 = np.load(f"{SD}/product_family24.npz", allow_pickle=True)["MODEL"].astype(float); F35 = np.load(f"{SD}/product_family35.npz", allow_pickle=True)["MODEL"].astype(float)
month = days.month.values; mask = np.zeros(G, bool)
for e in EP: mask[max(e["ipeak"] - 45, 0):e["ipeak"] + 91] = True
def deseason(X):
    Xd = X.copy()
    for m in range(1, 13): sel = (month == m) & ~mask; Xd[month == m] -= np.nanmean(X[sel], 0)[None]
    return Xd
PROD = {"interface": deseason(Y), "family 2-4 s": deseason(F24), "family 3-5 s": deseason(F35)}
def contrast(X, e, shift=0):
    i = wrap(e["ion"] + shift) if shift else e["ion"]; a = np.nanmean(X[max(i + 15, 0):i + 61][:, e["fp"]], 0); b = np.nanmean(X[max(i - 90, 0):i - 29][:, e["fp"]], 0); m = np.isfinite(a) & np.isfinite(b)
    return np.nanmean(a[m] - b[m]) if m.sum() >= 3 else np.nan
rows = []
for e in EP:
    rec = dict(key=e["key"], seg=e["seg"], size_mm=e["size"], seginc_mm=e["seginc"], tremor=e["tremor"], n_fp=e["nfp"])
    for pname, X in PROD.items():
        c = contrast(X, e); nul = np.array([contrast(X, e, rng.randint(90, span - 90)) for _ in range(400)]); nul = nul[np.isfinite(nul)]
        rec[f"c_{pname}"] = c; rec[f"z_{pname}"] = (c - nul.mean()) / nul.std() if len(nul) > 50 and nul.std() > 0 else np.nan; rec[f"p_{pname}"] = float(np.mean(np.abs(nul - nul.mean()) >= abs(c - nul.mean()))) if len(nul) > 50 else np.nan
    rows.append(rec)
T = pd.DataFrame(rows); T.to_csv(f"{SD}/large_events.csv", index=False)
print("PER-EPISODE (sorted by footprint slip sum, largest first):")
print(T.sort_values("size_mm", ascending=False)[["key", "size_mm", "seginc_mm", "tremor", "n_fp", "c_interface", "z_interface", "c_family 2-4 s", "z_family 2-4 s", "c_family 3-5 s", "z_family 3-5 s"]].round(2).to_string(index=False))
for pname in PROD:
    z = T[f"z_{pname}"]; print(f"{pname}: episodes with |z|>2: {int((z.abs() > 2).sum())}/{int(z.notna().sum())} (expected ~{0.05*z.notna().sum():.1f}); with z<-2: {int((z < -2).sum())}, z>+2: {int((z > 2).sum())}")
    for sz in ("size_mm", "seginc_mm", "tremor", "n_fp"):
        ok = z.notna() & T[sz].notna(); r = spearmanr(T[sz][ok], z[ok]).correlation; ra = spearmanr(T[sz][ok], z[ok].abs()).correlation
        print(f"    Spearman(z, {sz}) {r:+.2f} | (|z|, {sz}) {ra:+.2f}")
# composites of the largest quartile and top-5 northern
def battery(name, eps):
    for pname, X in PROD.items():
        s = np.array([contrast(X, e) for e in eps]); K, N, A = int((s < 0).sum()), int(np.isfinite(s).sum()), np.nanmean(s)
        nA = [np.nanmean([contrast(X, e, rng.randint(90, span - 90)) for e in eps]) for _ in range(600)]; nB = [np.nanmean([contrast(X, e, int(rng.choice([-1, 1]) * rng.randint(1, 15) * 365 + rng.randint(-15, 16))) for e in eps]) for _ in range(600)]
        print(f"  {name} {pname:13s}: n={N} mean {A:+.3f}% (K={K}/{N} decrease) pA {np.mean(np.abs(nA) >= abs(A)):.3f} pB {np.mean(np.abs(nB) >= abs(A)):.3f} (sd {np.std(nA):.3f})")
big = sorted(EP, key=lambda e: -e["size"])[:len(EP) // 4]; print(f"\nLARGEST QUARTILE by footprint slip ({len(big)}): {[e['key'] for e in big]}"); battery("top-quartile", big)
topN = sorted([e for e in EP if e["seg"] == "N"], key=lambda e: -e["seginc"])[:5]; print(f"\nTOP-5 NORTHERN by segment slip: {[e['key'] for e in topN]}"); battery("top-5 N", topN)
topT = sorted(EP, key=lambda e: -e["tremor"])[:8]; print(f"\nTOP-8 by tremor count: {[e['key'] for e in topT]}"); battery("top-8 tremor", topT)
# figures
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
for k, pname in enumerate(PROD):
    ax[k].scatter(T.size_mm, T[f"z_{pname}"], c=[{"N": "#1d4ed8", "C": "#059669", "S": "#dc2626", "SC": "#94a3b8"}[s] for s in T.seg], s=28)
    for _, r in T.nlargest(5, "size_mm").iterrows(): ax[k].annotate(r.key[2:9], (r.size_mm, r[f"z_{pname}"]), fontsize=7)
    ax[k].axhline(0, color="k", lw=.5); ax[k].axhline(2, color="r", lw=.5, ls="--"); ax[k].axhline(-2, color="r", lw=.5, ls="--"); ax[k].set_xlabel("episode footprint slip sum (mm)"); ax[k].set_ylabel("per-episode z of footprint contrast"); ax[k].set_title(pname, fontsize=9)
fig.suptitle("Per-episode footprint response vs episode size (blue N, green C, red S, grey SC); |z|>2 = beyond its own shift null", fontsize=10)
fig.savefig("figures/sse_large_events.png", dpi=110, bbox_inches="tight")
# maps for the 3 largest northern events: interface 30-d change at ON minus PRE per cell vs slip in window at peak
top3 = sorted([e for e in EP if e["seg"] == "N"], key=lambda e: -e["size"])[:3]
fig, ax = plt.subplots(2, 3, figsize=(14, 10))
for k, e in enumerate(top3):
    i = e["ion"]; X = PROD["interface"]; a = np.nanmean(X[max(i + 15, 0):i + 61], 0); b = np.nanmean(X[max(i - 90, 0):i - 29], 0); dm = a - b
    sc = ax[0, k].scatter(clon, clat, c=inc[e["ipeak"]], cmap="magma", vmin=0, vmax=15, s=12); ax[0, k].set_title(f"{e['key']}: modeled slip in window at peak (mm)", fontsize=9); ax[0, k].set_aspect(1 / np.cos(np.radians(47)))
    sc2 = ax[1, k].scatter(clon, clat, c=dm, cmap="RdBu_r", vmin=-0.6, vmax=0.6, s=12); ax[1, k].scatter(clon[e["fp"]], clat[e["fp"]], facecolor="none", edgecolor="k", s=16, linewidth=.4); ax[1, k].set_title("interface δβ/β: ON−PRE per cell (%) (black = footprint)", fontsize=9); ax[1, k].set_aspect(1 / np.cos(np.radians(47)))
    for a_ in ax[:, k]: a_.set_xlim(-126.5, -121.5); a_.set_ylim(45.5, 50)
plt.colorbar(sc, ax=ax[0, :], shrink=.6, label="mm"); plt.colorbar(sc2, ax=ax[1, :], shrink=.6, label="%")
fig.suptitle("Three largest northern modeled episodes: slip map vs interface δβ/β change map", fontsize=10); fig.savefig("figures/sse_large_event_maps.png", dpi=110, bbox_inches="tight")
print("wrote large_events.csv, figures/sse_large_events.png, figures/sse_large_event_maps.png")
