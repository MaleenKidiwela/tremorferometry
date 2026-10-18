#!/usr/bin/env python
"""Verification of the 49.2N lead (positive dv/v vs slip-in-window correlation in all products) and the 47.2N SITES
anticorrelation: per-patch (a) episode composites (N-segment ok episodes whose footprint overlaps the patch, aligned on
footprint onset), (b) Null A (rigid shift) and Null B (integer-year shift + jitter) for the window contrast and for
the whole-record correlation, (c) controls: cc, rows/cell-day, day-of-year of the patch's episodes, (d) time split
2011-18 / 2019-26, (e) lag of maximum cross-correlation between patch dv/v and slip-in-window (-120..+120 d).
Products: interface (deseasoned), family TOTAL, family site-removed, SITES. Output: figures/sse_verify_492.png + print.
"""
import numpy as np, pandas as pd
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(9)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
gi = pd.Series(np.arange(G), index=days); inc = zc["slip_thrust_inc30"].astype(float) * 1e3
zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True); PAT = np.asarray(zt["patterns"], float); plat = zt["patch_lat"]; band = zt["patch_band"]
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); md = pd.to_datetime(zp["days"])
Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=md).reindex(days).values
SITES = pd.DataFrame(np.asarray(zp["SITES"], float), index=md, columns=list(zp["sta_tags"])).reindex(days)
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "cell", "date", "dvv", "cc"]); cpos = {c: i for i, c in enumerate(cells)}
fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos); fd["v"] = -fd.dvv; fd["dw"] = fd.v - fd.groupby(["tag", "date"]).v.transform("mean")
cd = fd.groupby(["ci", "date"]).agg(tot=("v", "mean"), dw=("dw", "mean"), cc=("cc", "mean"), n=("v", "size")).reset_index(); okd = cd.date.isin(days)
ARR = {}
for k in ("tot", "dw", "cc", "n"):
    A = np.full((G, M), np.nan); A[gi[cd.date[okd]].values, cd.ci[okd].values] = cd[k][okd].values; ARR[k] = A
sta_rows = fd.groupby(["ci", "tag"]).size().reset_index(name="n"); del fd, cd
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
E0, E1 = int(gi[pd.Timestamp("2011-04-01")]), int(gi[pd.Timestamp("2026-04-01")]); span = E1 - E0; wrap = lambda i: E0 + (i - E0) % span
win = (days >= "2011-01-01") & (days <= "2026-07-20")
def pm(X, q):
    w = PAT[q]; fin = np.isfinite(X); num = np.nansum(np.where(fin, X, 0) * w[None], 1); den = (fin * w[None]).sum(1); o = num / np.where(den > 0, den, np.nan)
    o[(fin * w[None]).sum(1) / w.sum() < 0.5] = np.nan; return o
def deseason(y, mask=None):
    s = pd.Series(y, index=days); free = np.ones(G, bool) if mask is None else ~mask
    clim = {m: np.nanmean(y[(days.month == m) & free]) for m in range(1, 13)}; return y - np.array([clim[m] for m in days.month])
def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b) & win; return np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 500 else np.nan
LAGS = np.arange(-150, 201); L0 = 150; PRE, ON = (-90, -30), (15, 60)
def series_at(y, i): idx = i + LAGS; bad = (idx < 0) | (idx >= G); v = y[np.clip(idx, 0, G - 1)].copy(); v[bad] = np.nan; return v
def contrast(v): a = v[L0 + ON[0]:L0 + ON[1] + 1]; b = v[L0 + PRE[0]:L0 + PRE[1] + 1]; return np.nanmean(a) - np.nanmean(b) if np.isfinite(a).mean() > .6 and np.isfinite(b).mean() > .6 else np.nan
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 3, figsize=(17, 8))
for row, (qlat, seg) in enumerate(((49.2, "N"), (47.2, "N"))):
    q = [k for k in range(len(plat)) if band[k] == "deep" and abs(plat[k] - qlat) < 0.05][0]; cs = np.where(PAT[q] > 0.3)[0]
    sr = sta_rows[sta_rows.ci.isin(cs)].groupby("tag").n.sum(); sr = sr / sr.sum(); stas = [t for t in sr[sr >= 0.05].index if t in SITES.columns]
    okep = ep[ep.ok & (ep.seg == seg)]
    eps = [(r.key, int(gi[pd.Timestamp(r.fp_onset)]), int(r.doy)) for _, r in okep.iterrows() if len(np.intersect1d(zf[f"fp_{r.key}"], cs)) >= 3]
    ions = np.array([e[1] for e in eps]); mask = np.zeros(G, bool)
    for i in ions: mask[max(i - 30, 0):i + 121] = True
    slip = pm(inc, q)
    prods = {"interface": deseason(pm(Y, q), mask), "family TOTAL": deseason(pm(ARR["tot"], q), mask), "family site-removed": deseason(pm(ARR["dw"], q), mask),
             "SITES": deseason(SITES[stas].mean(axis=1).values, mask)}
    ctrl = {"cc": pm(ARR["cc"], q), "rows/cell-day": pm(ARR["n"], q)}
    print(f"\n=== deep patch {qlat}N: {len(cs)} cells, {len(eps)} N episodes overlapping, {len(stas)} stations ===")
    R = float(np.abs(np.mean(np.exp(2j * np.pi * np.array([e[2] for e in eps]) / 365.25)))); print(f"  episode day-of-year concentration R = {R:.2f}")
    for pname, y in list(prods.items()) + list(ctrl.items()):
        s_obs = np.array([contrast(series_at(y, i)) for i in ions]); K, N, A = int((s_obs < 0).sum()), int(np.isfinite(s_obs).sum()), np.nanmean(s_obs)
        r_obs = corr(y, slip)
        nA, nB, rA, rB = [], [], [], []
        for _ in range(2000):
            dA = rng.randint(90, span - 90); nA.append(np.nanmean([contrast(series_at(y, wrap(i + dA))) for i in ions])); rA.append(corr(y, np.roll(slip, dA)))
            dB = rng.choice([-1, 1]) * rng.randint(1, 15) * 365 + rng.randint(-15, 16); nB.append(np.nanmean([contrast(series_at(y, wrap(i + dB))) for i in ions])); rB.append(corr(y, np.roll(slip, dB)))
        pA = np.mean(np.abs(np.nan_to_num(nA)) >= abs(A)); pB = np.mean(np.abs(np.nan_to_num(nB)) >= abs(A)); prA = np.nanmean(np.abs(rA) >= abs(r_obs)); prB = np.nanmean(np.abs(rB) >= abs(r_obs))
        # time split + lag of max cross-corr
        e1 = days <= "2018-12-31"; r1 = corr(np.where(e1, y, np.nan), slip); r2 = corr(np.where(~e1, y, np.nan), slip)
        lags = np.arange(-120, 121, 5); xc = [corr(y, np.roll(slip, -L)) for L in lags]; Lbest = lags[int(np.nanargmax(np.abs(xc)))]
        tag = "CONTROL" if pname in ctrl else "product"
        print(f"  {tag:7s} {pname:20s}: ON-PRE mean {A:+.3f} (K={K}/{N} decr) pA {pA:.3f} pB {pB:.3f} | r(slip-in-window) {r_obs:+.3f} pA {prA:.3f} pB {prB:.3f} | 2011-18 {r1:+.2f} / 2019-26 {r2:+.2f} | best lag {Lbest:+d} d (r {xc[int(np.nanargmax(np.abs(xc)))]:+.2f})")
    # figure: composite + full series overlay
    a = ax[row, 0]
    for pname, col in zip(("interface", "family TOTAL", "SITES"), ("#1d4ed8", "#7c3aed", "#f59e0b")):
        ys = np.array([series_at(prods[pname], i) for i in ions]); c = np.nanmean(ys, 0); a.plot(LAGS, c - np.nanmean(c[L0 - 90:L0 - 30]), color=col, lw=1.5, label=pname)
    a2 = a.twinx(); sc_ = np.nanmean([series_at(slip, i) for i in ions], 0); a2.fill_between(LAGS, 0, sc_, color="#f59e0b", alpha=.2); a2.set_ylabel("slip in window (mm)", fontsize=8)
    a.axhline(0, color="k", lw=.5); a.axvspan(*PRE, color="#94a3b8", alpha=.15); a.axvspan(*ON, color="#f59e0b", alpha=.12); a.legend(fontsize=7); a.set_title(f"{qlat}N episode composite (N={len(eps)}, rel. PRE)", fontsize=9); a.set_xlabel("days from footprint onset"); a.set_ylabel("dv/v (%)")
    a = ax[row, 1]; y = prods["interface"]; y90 = pd.Series(y, index=days).rolling("90D", min_periods=20).mean()
    a.plot(days, y, color="#cbd5e1", lw=.4); a.plot(y90.index, y90.values, color="#1d4ed8", lw=1.2); a.set_ylim(-1.2, 1.2); a.set_ylabel("interface δβ/β (%)")
    a2 = a.twinx(); a2.fill_between(days, 0, slip, color="#f59e0b", alpha=.4); a2.set_ylim(0, None); a2.set_ylabel("slip in window (mm)", color="#b45309", fontsize=8)
    for i in ions: a.axvline(days[i], color="#dc2626", lw=.5, ls=":")
    a.set_xlim(pd.Timestamp("2011-01-01"), pd.Timestamp("2026-08-01")); a.set_title(f"{qlat}N: interface δβ/β (deseasoned) vs modeled slip; dotted = admitted episodes", fontsize=9)
    a = ax[row, 2]; ys = prods["SITES"]; y90 = pd.Series(ys, index=days).rolling("90D", min_periods=20).mean()
    a.plot(days, ys, color="#fde68a", lw=.4); a.plot(y90.index, y90.values, color="#b45309", lw=1.2); a.set_ylim(-0.6, 0.6); a.set_ylabel("near-station dv/v (%)")
    a2 = a.twinx(); a2.fill_between(days, 0, slip, color="#94a3b8", alpha=.4); a2.set_ylim(0, None)
    a.set_xlim(pd.Timestamp("2011-01-01"), pd.Timestamp("2026-08-01")); a.set_title(f"{qlat}N: SITES (mean of {len(stas)} footprint stations) vs modeled slip", fontsize=9)
fig.suptitle("Verification of the whole-record leads: 49.2N (positive dv/v–slip correlation in every product) and 47.2N (near-station anticorrelation)", fontsize=10)
fig.savefig("figures/sse_verify_492.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_verify_492.png")
