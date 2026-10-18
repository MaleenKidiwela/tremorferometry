#!/usr/bin/env python
"""CUMULATIVE slow-slip model results vs regional coda dv/v (user 2026-09-04): how the rise and fall of the modeled
cumulative slip / stress state relates to a region's delta-beta/beta.
Regions: 4 segments (S 40-42.5, SC 42.5-45, C 45-47.3, N 47.3-50; deep cells 30-45 km) and the 10 deep 125-km patches.
Model regional series (trailing frame, gauge 2019-24, regional mean over the region's cells): cumulative thrust slip
(mm), dilatation, Coulomb, shear, normal stress (10 km, fixed receiver, kPa), plus slip-in-window.
dv/v regional series: interface (production), interface gateD0.3, family 2-4 s, family 3-5 s, SITES of region stations —
all deseasoned by episode-masked monthly climatology.
Tests per region x product x model variable: (a) LEVEL correlation at lags -180..+180 d (30-d steps) vs 400 rigid
shifts of the model series (max-statistic over lags); (b) PHASE test: mean 30-d change of dv/v during LOADING phases
(cumulative slip falling over 30 d) vs RELEASE phases (rising), difference vs the same shift null; (c) overlays.
Outputs: res_catalog_sse/cumulative_regional.csv, figures/sse_cumulative_regional.png (segments), _patches.png
"""
import numpy as np, pandas as pd
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(81)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
clat, clon, dep, pdist = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"], zc["patch_dist_km"]; gi = pd.Series(np.arange(G), index=days)
win = np.asarray((days >= "2011-01-01") & (days <= "2026-07-20"))
MODEL = {"cum slip (mm)": zc["slip_thrust_cum"].astype(float) * 1e3, "slip in window (mm)": zc["slip_thrust_inc30"].astype(float) * 1e3, "dilatation (1e-8)": zc["dilatation"].astype(float) * 1e8,
         "Coulomb (kPa)": zc["coulomb"].astype(float) / 1e3, "shear (kPa)": zc["shear_stress"].astype(float) / 1e3, "normal (kPa)": zc["normal_stress"].astype(float) / 1e3}
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); md = pd.to_datetime(zp["days"]); Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=md).reindex(days).values
SITES = pd.DataFrame(np.asarray(zp["SITES"], float), index=md, columns=list(zp["sta_tags"])).reindex(days)
zg = np.load("fault_tomography/inversion/splits/gateD0.3/inversion_4d_daily.npz", allow_pickle=True); Yg = pd.DataFrame(np.asarray(zg["MODEL"], float), index=pd.to_datetime(zg["days"])).reindex(days).values
F24 = np.load(f"{SD}/product_family24.npz", allow_pickle=True)["MODEL"].astype(float); F35 = np.load(f"{SD}/product_family35.npz", allow_pickle=True)["MODEL"].astype(float)
ep = pd.read_csv(f"{SD}/episodes.csv"); okep = ep[ep.ok]; month = days.month.values; mask = np.zeros(G, bool)
for i in okep.ipeak: mask[max(int(i) - 45, 0):int(i) + 91] = True
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "cell"]); cpos = {c: i for i, c in enumerate(cells)}; fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos); sta_rows = fd.groupby(["ci", "tag"]).size().reset_index(name="n"); del fd
zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True); PAT = np.asarray(zt["patterns"], float); plat = zt["patch_lat"]; band = zt["patch_band"]
deepc = np.isfinite(dep) & (dep >= 30) & (dep <= 45) & (pdist <= 15)
REG = {s: np.where(deepc & (clat >= a) & (clat < b))[0] for s, (a, b) in {"S": (40, 42.5), "SC": (42.5, 45), "C": (45, 47.3), "N": (47.3, 50)}.items()}
for q in np.argsort(plat):
    if band[q] == "deep" and 40.5 < plat[q] < 49.5: REG[f"patch {plat[q]:.1f}N"] = np.where((PAT[q] > 0.3) & np.isfinite(dep))[0]
def rmean(X, cs):
    v = X[:, cs]; fin = np.isfinite(v); out = np.nanmean(v, 1); out[fin.mean(1) < 0.5] = np.nan; return out
def deseason(y):
    out = y.copy()
    for m in range(1, 13): sel = (month == m) & ~mask & np.isfinite(y); out[month == m] -= np.nanmean(y[sel]) if sel.any() else 0
    return out
def d30(y): o = np.full_like(y, np.nan); o[30:] = y[30:] - y[:-30]; return o
def corr(a, b): m = np.isfinite(a) & np.isfinite(b) & win; return np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 500 else np.nan
LAGS = np.arange(-180, 181, 30)
rows = []; series = {}
for rname, cs in REG.items():
    sr = sta_rows[sta_rows.ci.isin(cs)].groupby("tag").n.sum(); sr = sr / sr.sum(); stas = [t for t in sr[sr >= 0.05].index if t in SITES.columns]
    dv = {"interface": deseason(rmean(Y, cs)), "interface gateD0.3": deseason(rmean(Yg, cs)), "family 2-4 s": deseason(rmean(F24, cs)), "family 3-5 s": deseason(rmean(F35, cs)),
          "SITES": deseason(SITES[stas].mean(axis=1).values) if stas else np.full(G, np.nan)}
    mod = {k: rmean(X, cs) for k, X in MODEL.items()}; series[rname] = (dv, mod)
    cum = mod["cum slip (mm)"]; dcum = d30(cum); loading = dcum < -0.5; release = dcum > 0.5          # falling = inter-episode loading (slip below trend), rising = release
    for pname, y in dv.items():
        dy = d30(y)
        if np.isfinite(y[win]).mean() < 0.3: continue
        # (b) phase test
        ph = np.nanmean(dy[win & release]) - np.nanmean(dy[win & loading]); nph = []
        for _ in range(400):
            sh = rng.randint(365, G - 365); dc = np.roll(dcum, sh); nph.append(np.nanmean(dy[win & (dc > 0.5)]) - np.nanmean(dy[win & (dc < -0.5)]))
        p_ph = float(np.mean(np.abs(np.nan_to_num(nph)) >= abs(ph)))
        rec = dict(region=rname, product=pname, n_stations=len(stas), phase_release_minus_loading=ph, p_phase=p_ph, mean_dy_release=np.nanmean(dy[win & release]), mean_dy_loading=np.nanmean(dy[win & loading]))
        # (a) level correlation vs model variables at lags, max-stat over lags
        for vname, x in mod.items():
            xl = np.array([corr(y, np.roll(x, L)) for L in LAGS]); k = int(np.nanargmax(np.abs(xl))) if np.isfinite(xl).any() else 0
            nulmax = []
            for _ in range(200):
                sh = rng.randint(365, G - 365); xs = np.roll(x, sh); nulmax.append(np.nanmax(np.abs([corr(y, np.roll(xs, L)) for L in LAGS[::2]])))
            rec[f"r_{vname}"] = xl[k]; rec[f"lag_{vname}"] = int(LAGS[k]); rec[f"r0_{vname}"] = xl[LAGS == 0][0]; rec[f"p_{vname}"] = float(np.mean(np.array(nulmax) >= abs(xl[k])))
        rows.append(rec)
    print(f"{rname} done", flush=True)
T = pd.DataFrame(rows); T.to_csv(f"{SD}/cumulative_regional.csv", index=False)
print("\nPHASE TEST (mean 30-d dv/v change during RELEASE minus during LOADING, %; shift null):")
print(T.pivot(index="region", columns="product", values="phase_release_minus_loading").round(3).to_string()); print(T.pivot(index="region", columns="product", values="p_phase").round(3).to_string())
print("\nLEVEL CORRELATION with cumulative slip (best lag within +-180 d; p = max-statistic over lags):")
print(T.pivot(index="region", columns="product", values="r_cum slip (mm)").round(2).to_string()); print(T.pivot(index="region", columns="product", values="p_cum slip (mm)").round(3).to_string())
for v in ("Coulomb (kPa)", "dilatation (1e-8)", "normal (kPa)", "shear (kPa)"):
    sig = T[T[f"p_{v}"] < 0.05]; print(f"{v}: region-product cells with max-stat p<0.05: {len(sig)}/{len(T)} -> " + "; ".join(f"{r.region}|{r['product']} r={r[f'r_{v}']:+.2f}@{r[f'lag_{v}']:+d}d" for _, r in sig.iterrows()))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
for fname, keys in (("figures/sse_cumulative_regional.png", ["N", "C", "SC", "S"]), ("figures/sse_cumulative_patches.png", [k for k in REG if k.startswith("patch")])):
    fig, ax = plt.subplots(len(keys), 1, figsize=(14, 2.3 * len(keys)), sharex=True)
    for a, rname in zip(np.atleast_1d(ax), keys):
        dv, mod = series[rname]
        for pname, col in (("interface", "#1d4ed8"), ("family 2-4 s", "#7c3aed"), ("family 3-5 s", "#0891b2")):
            y = pd.Series(dv[pname], index=days).rolling("90D", min_periods=20).mean(); a.plot(y.index, y.values, color=col, lw=1.1, label=pname)
        a.set_ylim(-1, 1); a.set_ylabel("dv/v (%)", fontsize=8); a.axhline(0, color="k", lw=.4)
        a2 = a.twinx(); a2.plot(days, mod["cum slip (mm)"], color="#f59e0b", lw=1.2, label="cum. thrust slip"); a2.set_ylabel("cum slip (mm)", color="#b45309", fontsize=8)
        a3 = a.twinx(); a3.spines["right"].set_position(("axes", 1.06)); a3.plot(days, mod["Coulomb (kPa)"], color="#dc2626", lw=.8, alpha=.7); a3.set_ylabel("Coulomb (kPa)", color="#dc2626", fontsize=8)
        a.set_title(f"{rname}: regional dv/v (90-d mean) vs modeled cumulative slip (orange) and Coulomb stress state (red)", fontsize=9)
    np.atleast_1d(ax)[0].legend(fontsize=7, loc="upper left"); np.atleast_1d(ax)[-1].set_xlim(pd.Timestamp("2011-01-01"), pd.Timestamp("2026-08-01"))
    fig.savefig(fname, dpi=110, bbox_inches="tight"); print("wrote", fname)
