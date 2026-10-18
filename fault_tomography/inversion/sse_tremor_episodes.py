#!/usr/bin/env python
"""Tremor-defined ETS episodes margin-wide (independent of the geodetic model's weak southern slip) and the footprint
response of coda dv/v to them, with the same pre-registered statistics as the model-based test.
Catalog: PNSN tremor, 0.5-deg latitude bands 40-50N; per band the trailing-30-d count; episodes = local maxima above
the band's 90th percentile, >= 90 d apart, 2011-04..2026-04; onset = first day the trailing count exceeds 25 % of the
peak (within 45 d before); footprint = deep cells (30-45 km) within 35 km of tremor epicentres of the episode window.
Model check: modeled thrust slip-in-window in the footprint at the peak (does the model see the tremor episode?).
Response products: interface (deseasoned), family 2-4 s, family 3-5 s, SITES of footprint stations (>=5 % rows).
Statistic: footprint mean ON[+15,+60]-PRE[-90,-30], per segment (S 40-42.5, SC 42.5-45, C 45-47.3, N 47.3-50) and pooled;
Null A = rigid circular shift per segment (1000), Null B = integer-year shift + jitter (1000); ring control.
Outputs: res_catalog_sse/tremor_episodes.csv, tremor_episode_test.csv, figures/sse_tremor_episode_test.png
"""
import numpy as np, pandas as pd, os
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(51)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
clat, clon, dep, pdist = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"], zc["patch_dist_km"]; gi = pd.Series(np.arange(G), index=days); inc = zc["slip_thrust_inc30"].astype(float) * 1e3
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
E0, E1 = pd.Timestamp("2011-04-01"), pd.Timestamp("2026-04-01"); i0, i1 = int(gi[E0]), int(gi[E1]); span = i1 - i0; wrap = lambda i: i0 + (i - i0) % span
SEG = {"S": (40.0, 42.5), "SC": (42.5, 45.0), "C": (45.0, 47.3), "N": (47.3, 50.0)}
tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"]); tr["t"] = pd.to_datetime(tr.time, errors="coerce"); tr = tr.dropna(subset=["t"]); tr["di"] = gi.reindex(tr.t.dt.floor("D")).values
tr = tr.dropna(subset=["di"]); tr["di"] = tr.di.astype(int); TXY = km(tr.lat.values, tr.lon.values)
deepc = np.isfinite(dep) & (dep >= 30) & (dep <= 45)
# ---- catalog ----
EP = []
for b0 in np.arange(40.0, 50.0, 0.5):
    sel = (tr.lat >= b0) & (tr.lat < b0 + 0.5); cnt = np.bincount(tr.di[sel], minlength=G).astype(float); c30 = pd.Series(cnt).rolling(30, min_periods=1).sum().values
    thr = np.percentile(c30[(days >= "2010-01-01")], 90)
    if thr < 50: continue
    cand = [i for i in range(1, G - 1) if c30[i] >= thr and c30[i] >= c30[i - 1] and c30[i] >= c30[i + 1] and E0 <= days[i] <= E1]
    cand = sorted(cand, key=lambda i: -c30[i]); keep = []
    for i in cand:
        if all(abs(i - k) >= 90 for k in keep): keep.append(i)
    for i in sorted(keep):
        on = i
        while on > i - 45 and c30[on - 1] > 0.25 * c30[i]: on -= 1
        w = sel & (tr.di >= on - 5) & (tr.di <= i + 5); pts = TXY[w.values]
        if len(pts) < 20: continue
        d = np.sqrt(((XY[:, None, :] - pts[None, :, :]) ** 2).sum(-1)).min(1); fp = np.where(deepc & (d <= 35))[0]
        if len(fp) < 5: continue
        seg = [s for s, (a, b) in SEG.items() if a <= b0 + 0.25 < b][0]
        EP.append(dict(key=f"T{b0:.1f}_{days[i].date()}", band=b0, seg=seg, ipeak=i, ion=on, peak=str(days[i].date()), onset=str(days[on].date()), count30=c30[i], n_fp=len(fp), fp=fp,
                       model_slip_mm=float(np.nanmean(inc[i, fp])), doy=int(days[i].dayofyear)))
cat = pd.DataFrame([{k: v for k, v in e.items() if k != "fp"} for e in EP]); cat.to_csv(f"{SD}/tremor_episodes.csv", index=False)
print(f"tremor episodes: {len(EP)} | per segment: {cat.groupby('seg').size().to_dict()} | model slip-in-window at peak (mm): median {cat.model_slip_mm.median():.1f}, per segment {cat.groupby('seg').model_slip_mm.median().round(1).to_dict()}")
# merge overlapping episodes across adjacent bands (same peak within 20 d and footprints overlapping) -> keep the larger count
EP.sort(key=lambda e: -e["count30"]); kept = []
for e in EP:
    if any(abs(e["ipeak"] - k["ipeak"]) <= 20 and len(np.intersect1d(e["fp"], k["fp"])) > 0.5 * min(len(e["fp"]), len(k["fp"])) for k in kept): continue
    kept.append(e)
EP = sorted(kept, key=lambda e: e["ipeak"]); print(f"after merging overlapping bands: {len(EP)} episodes; per segment {pd.Series([e['seg'] for e in EP]).value_counts().to_dict()}")
# ---- products ----
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); md = pd.to_datetime(zp["days"]); Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=md).reindex(days).values
SITES = pd.DataFrame(np.asarray(zp["SITES"], float), index=md, columns=list(zp["sta_tags"])).reindex(days)
F24 = np.load(f"{SD}/product_family24.npz", allow_pickle=True)["MODEL"].astype(float); F35 = np.load(f"{SD}/product_family35.npz", allow_pickle=True)["MODEL"].astype(float)
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "cell"]); cpos = {c: i for i, c in enumerate(cells)}; fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos); sta_rows = fd.groupby(["ci", "tag"]).size().reset_index(name="n"); del fd
month = days.month.values; mask = np.zeros(G, bool)
for e in EP: mask[max(e["ipeak"] - 45, 0):e["ipeak"] + 91] = True
def deseason(X):
    Xd = X.copy()
    for m in range(1, 13): sel = (month == m) & ~mask; Xd[month == m] -= np.nanmean(X[sel], 0)[None]
    return Xd
PROD = {"interface": deseason(Y), "family 2-4 s": deseason(F24), "family 3-5 s": deseason(F35)}
for e in EP:
    sr = sta_rows[sta_rows.ci.isin(e["fp"])].groupby("tag").n.sum(); sr = sr / sr.sum(); e["stas"] = [t for t in sr[sr >= 0.05].index if t in SITES.columns]
    d = np.sqrt(((XY[:, None, :] - XY[e["fp"]][None, :, :]) ** 2).sum(-1)).min(1); e["ring"] = np.where(deepc & (d >= 100) & (d <= 300))[0]
def contrast(X, e, shift=0, cellset=None):
    i = wrap(e["ion"] + shift) if shift else e["ion"]; cs = e["fp"] if cellset is None else cellset
    a = np.nanmean(X[max(i + 15, 0):i + 61][:, cs], 0); b = np.nanmean(X[max(i - 90, 0):i - 29][:, cs], 0); m = np.isfinite(a) & np.isfinite(b); return np.nanmean(a[m] - b[m]) if m.sum() >= 3 else np.nan
def sites_contrast(e, shift=0):
    if not e["stas"]: return np.nan
    i = wrap(e["ion"] + shift) if shift else e["ion"]; y = SITES[e["stas"]].mean(axis=1).values; a = y[max(i + 15, 0):i + 61]; b = y[max(i - 90, 0):i - 29]
    return np.nanmean(a) - np.nanmean(b) if np.isfinite(a).mean() > .6 and np.isfinite(b).mean() > .6 else np.nan
SITESd = SITES - SITES.groupby(SITES.index.month).transform("mean")
rows = []
for subset in ("ALL", "N", "C", "SC", "S", "N+C"):
    eps = [e for e in EP if subset == "ALL" or e["seg"] == subset or (subset == "N+C" and e["seg"] in ("N", "C"))]
    if len(eps) < 4: continue
    for pname, X in list(PROD.items()) + [("SITES", None)]:
        f = (lambda e, sh=0: sites_contrast(e, sh)) if X is None else (lambda e, sh=0: contrast(X, e, sh))
        s = np.array([f(e) for e in eps]); K, N, A = int((s < 0).sum()), int(np.isfinite(s).sum()), np.nanmean(s)
        nA = [np.nanmean([f(e, rng.randint(90, span - 90)) for e in eps]) for _ in range(600)]
        nB = [np.nanmean([f(e, int(rng.choice([-1, 1]) * rng.randint(1, 15) * 365 + rng.randint(-15, 16))) for e in eps]) for _ in range(600)]
        ring = np.nanmean([contrast(X, e, 0, e["ring"]) for e in eps if len(e["ring"]) >= 10]) if X is not None else np.nan
        rows.append(dict(subset=subset, product=pname, n=N, K_decrease=K, mean=A, pA=float(np.mean(np.abs(nA) >= abs(A))), pB=float(np.mean(np.abs(nB) >= abs(A))), sdA=float(np.std(nA)), ring=ring))
        print(f"{subset:4s} {pname:13s}: n={N:2d} mean {A:+.3f}% (K={K}/{N} decrease) pA {rows[-1]['pA']:.3f} pB {rows[-1]['pB']:.3f} (sd {rows[-1]['sdA']:.3f}) | ring {ring:+.3f}", flush=True)
T = pd.DataFrame(rows); T.to_csv(f"{SD}/tremor_episode_test.csv", index=False)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
LAGS = np.arange(-150, 241); L0 = 150
fig, ax = plt.subplots(1, 4, figsize=(18, 4), sharey=True)
for k, subset in enumerate(("N", "C", "SC", "S")):
    eps = [e for e in EP if e["seg"] == subset]
    for pname, col in zip(PROD, ("#1d4ed8", "#7c3aed", "#0891b2")):
        ys = []
        for e in eps:
            idx = e["ion"] + LAGS; bad = (idx < 0) | (idx >= G); v = np.nanmean(PROD[pname][np.clip(idx, 0, G - 1)][:, e["fp"]], 1); v[bad] = np.nan; ys.append(v)
        c = np.nanmean(ys, 0); ax[k].plot(LAGS, c - np.nanmean(c[L0 - 90:L0 - 30]), color=col, lw=1.3, label=pname)
    ax[k].axhline(0, color="k", lw=.5); ax[k].axvspan(-90, -30, color="#94a3b8", alpha=.15); ax[k].axvspan(15, 60, color="#f59e0b", alpha=.15); ax[k].set_title(f"{subset}: {len(eps)} tremor episodes", fontsize=9); ax[k].set_xlabel("days from tremor onset")
ax[0].set_ylabel("dv/v rel. PRE (%)"); ax[0].legend(fontsize=7); ax[0].set_ylim(-0.4, 0.4)
fig.suptitle("Coda dv/v composites around TREMOR-defined ETS episodes (footprint = deep cells within 35 km of tremor), by segment", fontsize=10)
fig.savefig("figures/sse_tremor_episode_test.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_tremor_episode_test.png, tremor_episodes.csv, tremor_episode_test.csv")
