#!/usr/bin/env python
"""Exploratory PATTERN SEARCH (Merlin stage 5/6, max-statistic FWER): where, when, in which product, and for which
episodes does coda dv/v respond to modeled slow slip?  One shared Null-A draw set for every hypothesis so the best hit
carries a family-wise p.  Products: interface dbeta/beta (deseasoned; raw), family-level coda dv/v at footprint cells
(TOTAL = -eps, no site removal; DW = station-common removed), near-station SITES.  Cell classes relative to each
episode: footprint (30-45 km), up-dip band (15-30 km, same along-strike span), down-dip (45-55 km), ring (>=100 km).
Windows: pulse ON[+15,+60], delayed [+60,+150], step [+90,+180], precursory [-30,+15], all vs PRE[-90,-30].
Alignments: local onset, footprint onset. Subsets: N+C, N, C, S, big/small (footprint slip), early/late years.
Plus: PCA of footprint responses in 30-d lag bins (eigenvalue + split-half nulls), per-cell sensitivity map (100-km
smoothed, z vs null), contrast-vs-episode-size, and the composite per class.
Outputs: res_catalog_sse/pattern_search.csv (+ .npz), figures/sse_pattern_search.png
"""
import numpy as np, pandas as pd, os, itertools
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"
rng = np.random.RandomState(23); NDRAW = int(os.environ.get("NDRAW", "1000"))
PRE = (-90, -30); LAGS = np.arange(-200, 251); L0 = -LAGS[0]
WINS = {"pulse[+15,+60]": (15, 60), "delayed[+60,+150]": (60, 150), "step[+90,+180]": (90, 180), "precursory[-30,+15]": (-30, 15)}
SEGLAT = {"S": (40.0, 42.5), "SC": (42.5, 45.0), "C": (45.0, 47.3), "N": (47.3, 49.5)}
E0, E1 = pd.Timestamp("2011-04-01"), pd.Timestamp("2026-04-01")
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True)
days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells); clat, clon = zc["cell_lat"], zc["cell_lon"]
dep = zc["iface_depth_km"]; pdist = zc["patch_dist_km"]; gi = pd.Series(np.arange(G), index=days)
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); mdays = pd.to_datetime(zp["days"])
Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=mdays).reindex(days).values
SITES = pd.DataFrame(np.asarray(zp["SITES"], float), index=mdays, columns=list(zp["sta_tags"])).reindex(days)
inc30 = zc["slip_thrust_inc30"].astype(float)
ep = pd.read_csv(f"{SD}/episodes.csv"); zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True); ep["key"] = ep.seg + "_" + ep.peak
okep = ep[ep.ok].reset_index(drop=True)
segcell = np.full(M, "", object)
for s, (a, b) in SEGLAT.items(): segcell[(clat >= a) & (clat < b)] = s
i0, i1 = int(gi[E0]), int(gi[E1]); span = i1 - i0
wrap = lambda i: i0 + (i - i0) % span
# family-level products at cells (TOTAL and station-common-removed), true dv/v sign
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "cell", "date", "dvv"]); cpos = {c: i for i, c in enumerate(cells)}
fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos); fd["v"] = -fd.dvv
fd["dw"] = fd.v - fd.groupby(["tag", "date"]).v.transform("mean")
cd = fd.groupby(["ci", "date"]).agg(tot=("v", "mean"), dw=("dw", "mean")).reset_index(); okd = cd.date.isin(days)
TOT = np.full((G, M), np.nan); DW = np.full((G, M), np.nan)
TOT[gi[cd.date[okd]].values, cd.ci[okd].values] = cd.tot[okd].values; DW[gi[cd.date[okd]].values, cd.ci[okd].values] = cd.dw[okd].values
sta_rows = fd.groupby(["ci", "tag"]).size().reset_index(name="n"); del fd, cd
# episode-masked deseason
mask_days = {s: np.zeros(G, bool) for s in SEGLAT}
for _, r in okep.iterrows(): mask_days[r.seg][max(int(r.ipeak) - 30, 0):int(r.ipeak) + 91] = True
month = days.month.values
def deseason(X):
    Xd = X.copy()
    for s in SEGLAT:
        cs = np.where(segcell == s)[0]
        if not len(cs): continue
        free = ~mask_days[s]
        for m in range(1, 13):
            sel = (month == m) & free; Xd[np.ix_(month == m, cs)] -= np.nanmean(X[sel][:, cs], 0)[None, :]
    return Xd
Yd, TOTd, DWd = deseason(Y), deseason(TOT), deseason(DW)
PROD = {"interface(deseason)": Yd, "interface(raw)": Y, "family TOTAL(deseason)": TOTd, "family site-removed(deseason)": DWd}
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
deep = np.isfinite(dep) & (pdist <= 15)
class Ep: pass
EPS = []
for _, r in okep.iterrows():
    e = Ep(); e.seg = r.seg; e.key = r.key; e.fp = zf[f"fp_{r.key}"]; e.lo = zf[f"lo_{r.key}"]; e.ion = int(gi[pd.Timestamp(r.fp_onset)])
    e.size = r.fp_slip_sum_mm; e.year = int(r.peak[:4]); e.doy = int(r.doy)
    d = np.sqrt(((XY[:, None, :] - XY[e.fp][None, :, :]) ** 2).sum(-1)).min(1)
    y0, y1 = XY[e.fp, 1].min() - 25, XY[e.fp, 1].max() + 25; strip = (XY[:, 1] >= y0) & (XY[:, 1] <= y1) & np.isfinite(dep)
    e.cls = {"footprint": e.fp, "updip(15-30km)": np.where(strip & (dep >= 15) & (dep < 30))[0], "downdip(45-55km)": np.where(strip & (dep > 45) & (dep <= 55))[0],
             "ring(>=100km)": np.where(deep & (dep >= 30) & (dep <= 45) & (d >= 100) & (d <= 300))[0]}
    sr = sta_rows[sta_rows.ci.isin(e.fp)].groupby("tag").n.sum(); sr = sr / sr.sum(); e.stas = [t for t in sr[sr >= 0.05].index if t in SITES.columns]
    EPS.append(e)
def aligned(X, e, cs, local, shift=0):
    on = (e.lo if local and len(cs) == len(e.fp) and np.array_equal(cs, e.fp) else np.full(len(cs), e.ion)) + shift
    idx = on[:, None] + LAGS[None, :]; bad = (idx < 0) | (idx >= G); v = X[np.clip(idx, 0, G - 1), cs[:, None]]; v[bad] = np.nan
    return np.nanmean(v, 0) if np.isfinite(v).any() else np.full(len(LAGS), np.nan)
def contrast(y, w):
    a = y[L0 + w[0]:L0 + w[1] + 1]; b = y[L0 + PRE[0]:L0 + PRE[1] + 1]
    return np.nanmean(a) - np.nanmean(b) if (np.isfinite(a).mean() >= 0.6 and np.isfinite(b).mean() >= 0.6) else np.nan
def sites_series(e, shift=0):
    if not e.stas: return np.full(len(LAGS), np.nan)
    y = SITES[e.stas].mean(axis=1).values; idx = e.ion + shift + LAGS; bad = (idx < 0) | (idx >= G); v = y[np.clip(idx, 0, G - 1)]; v[bad] = np.nan; return v
def shifted(e, sh):
    f = Ep(); f.__dict__.update(e.__dict__); f.lo = wrap(e.lo + sh); f.ion = wrap(e.ion + sh); return f
med_size = np.median([e.size for e in EPS if e.seg in ("N", "C")])
SUBSETS = {"N+C": lambda e: e.seg in ("N", "C"), "N": lambda e: e.seg == "N", "C": lambda e: e.seg == "C", "S": lambda e: e.seg == "S",
           "N+C big": lambda e: e.seg in ("N", "C") and e.size >= med_size, "N+C small": lambda e: e.seg in ("N", "C") and e.size < med_size,
           "N+C 2011-18": lambda e: e.seg in ("N", "C") and e.year <= 2018, "N+C 2019-26": lambda e: e.seg in ("N", "C") and e.year >= 2019}
# ---- hypothesis list: (product, class, window, alignment, subset) -> per-episode contrasts; statistic = mean ----
HYP = []
for p, c, w, al in itertools.product(PROD, ["footprint", "updip(15-30km)", "downdip(45-55km)", "ring(>=100km)"], WINS, ["local", "fp-onset"]):
    if al == "local" and c != "footprint": continue
    HYP.append((p, c, w, al))
for w in WINS: HYP.append(("SITES(near-station)", "footprint stations", w, "fp-onset"))
def all_contrasts(eps):
    out = np.full((len(HYP), len(eps)), np.nan)
    cache = {}
    for j, e in enumerate(eps):
        for h, (p, c, w, al) in enumerate(HYP):
            k = (p, c, al, j)
            if k not in cache:
                cache[k] = sites_series(e) if p.startswith("SITES") else aligned(PROD[p], e, e.cls[c], al == "local")
            out[h, j] = contrast(cache[k], WINS[w])
    return out
C_obs = all_contrasts(EPS)
subs = {nm: np.array([f(e) for e in EPS]) for nm, f in SUBSETS.items()}
def stat_table(C):
    T = np.full((len(HYP), len(SUBSETS)), np.nan)
    for k, (nm, m) in enumerate(subs.items()):
        v = C[:, m]; T[:, k] = np.nanmean(v, 1)
    return T
S_obs = stat_table(C_obs)
print(f"hypotheses {len(HYP)} x subsets {len(SUBSETS)} = {S_obs.size} statistics; episodes {len(EPS)}; null draws {NDRAW}", flush=True)
null = np.zeros((NDRAW,) + S_obs.shape)
for d in range(NDRAW):
    sh = {s: int(rng.randint(90, span - 90)) for s in SEGLAT}
    null[d] = stat_table(all_contrasts([shifted(e, sh[e.seg]) for e in EPS]))
    if d % 200 == 0: print(f"  draw {d}", flush=True)
sd = np.nanstd(null, 0); z = S_obs / np.where(sd > 0, sd, np.nan)
p_raw = (np.abs(null) >= np.abs(S_obs)[None]).mean(0)
zmax_null = np.nanmax(np.abs(null / np.where(sd > 0, sd, np.nan)[None]).reshape(NDRAW, -1), 1)
p_fwer = np.array([(zmax_null >= abs(zz)).mean() for zz in z.ravel()]).reshape(z.shape)
rows = []
for h, (p, c, w, al) in enumerate(HYP):
    for k, nm in enumerate(SUBSETS):
        n = int(np.isfinite(C_obs[h, subs[nm]]).sum())
        rows.append(dict(product=p, cells=c, window=w, align=al, subset=nm, n=n, mean=S_obs[h, k], null_sd=sd[h, k], z=z[h, k], p_raw=p_raw[h, k], p_fwer=p_fwer[h, k],
                         frac_neg=float(np.nanmean(C_obs[h, subs[nm]] < 0)) if n else np.nan))
tab = pd.DataFrame(rows).sort_values("p_raw"); tab.to_csv(f"{SD}/pattern_search.csv", index=False)
print("\nTOP 12 hypotheses (lowest raw p; p_fwer = family-wise over all %d):" % S_obs.size)
print(tab.head(12).round(3).to_string(index=False))
print(f"\nbest family-wise p = {tab.p_fwer.min():.3f}; number with p_raw<0.01: {(tab.p_raw < 0.01).sum()} (expected by chance ~{0.01 * S_obs.size:.1f})")
# ---- contrast vs episode size / year (N+C, interface footprint pulse, local) ----
from scipy.stats import spearmanr
h0 = HYP.index(("interface(deseason)", "footprint", "pulse[+15,+60]", "local")); m = subs["N+C"]
cs_, sz, yr = C_obs[h0, m], np.array([e.size for e in EPS])[m], np.array([e.year for e in EPS])[m]
ok = np.isfinite(cs_); rs = spearmanr(cs_[ok], sz[ok]).correlation; ry = spearmanr(cs_[ok], yr[ok]).correlation
print(f"\ncontrast vs footprint slip: Spearman {rs:+.2f}; vs year: {ry:+.2f}  (N+C, n={ok.sum()})")
# ---- PCA of footprint responses (Merlin Q4) ----
bins = np.arange(-90, 151, 30)
def resp_matrix(eps):
    R = []
    for e in eps:
        y = aligned(Yd, e, e.fp, True); pre = y[L0 - 90:L0 - 30]; base = np.nanmean(pre); s = np.nanstd(y[L0 - 180:L0 - 30])
        if not np.isfinite(base) or not s: continue
        R.append([(np.nanmean(y[L0 + b:L0 + b + 30]) - base) / s for b in bins[:-1]])
    return np.array(R)
def lead(R):
    R = R[np.isfinite(R).all(1)]
    if len(R) < 6: return np.nan, None, len(R)
    Rc = R - R.mean(0); w, v = np.linalg.eigh(Rc.T @ Rc / len(R)); return w[-1] / w.sum(), v[:, -1], len(R)
PR = [e for e in EPS if e.seg in ("N", "C")]
ev_obs, mode, nR = lead(resp_matrix(PR)); ev_null = []
for d in range(300):
    sh = {s: int(rng.randint(90, span - 90)) for s in SEGLAT}; ev_null.append(lead(resp_matrix([shifted(e, sh[e.seg]) for e in PR]))[0])
p_ev = float(np.nanmean(np.array(ev_null) >= ev_obs))
R = resp_matrix(PR); R = R[np.isfinite(R).all(1)]; odd, even = R[::2], R[1::2]
sh_corr = abs(np.corrcoef(lead(odd)[1], lead(even)[1])[0, 1]) if len(odd) >= 6 and len(even) >= 6 else np.nan
print(f"PCA (N+C footprint responses, {nR} episodes x {len(bins)-1} lag bins): leading eigenvalue fraction {ev_obs:.2f}, null p = {p_ev:.2f}; split-half mode corr {sh_corr:.2f}; mode shape {np.round(mode, 2) if mode is not None else None}")
# ---- per-cell sensitivity map (N+C, interface deseasoned, local onset, pulse window), 100-km smoothed z ----
def cellmap(eps):
    num = np.zeros(M); cnt = np.zeros(M)
    for e in eps:
        idx_on = e.lo[:, None] + np.arange(15, 61)[None, :]; idx_pre = e.lo[:, None] + np.arange(-90, -29)[None, :]
        a = np.nanmean(Yd[np.clip(idx_on, 0, G - 1), e.fp[:, None]], 1); b = np.nanmean(Yd[np.clip(idx_pre, 0, G - 1), e.fp[:, None]], 1)
        c = a - b; okc = np.isfinite(c); num[e.fp[okc]] += c[okc]; cnt[e.fp[okc]] += 1
    m = np.where(cnt > 0, num / np.maximum(cnt, 1), np.nan)
    W = np.exp(-0.5 * (np.sqrt(((XY[:, None, :] - XY[None, :, :]) ** 2).sum(-1)) / 50.0) ** 2); W[:, ~np.isfinite(m)] = 0
    return (W @ np.nan_to_num(m)) / np.maximum(W.sum(1), 1e-9), cnt
map_obs, cnt = cellmap(PR); maps_null = np.array([cellmap([shifted(e, sh[e.seg]) for e in PR])[0] for sh in [{s: int(rng.randint(90, span - 90)) for s in SEGLAT} for _ in range(200)]])
zmap = map_obs / np.where(maps_null.std(0) > 0, maps_null.std(0), np.nan); zmap[cnt < 3] = np.nan
zmax_map = np.nanmax(np.abs(maps_null / np.where(maps_null.std(0) > 0, maps_null.std(0), np.nan)), 1)
print(f"sensitivity map: max |z| {np.nanmax(np.abs(zmap)):.2f} at {cells[int(np.nanargmax(np.abs(zmap)))]}; FWER p = {(zmax_map >= np.nanmax(np.abs(zmap))).mean():.2f} (cells with >=3 episodes: {(cnt >= 3).sum()})")
np.savez_compressed(f"{SD}/pattern_search.npz", S_obs=S_obs, z=z, p_raw=p_raw, p_fwer=p_fwer, hyp=np.array(HYP, object), subsets=np.array(list(SUBSETS)),
                    zmap=zmap, cnt=cnt, ev_obs=ev_obs, p_ev=p_ev, sh_corr=sh_corr, mode=mode, C_obs=C_obs)
# ---- figure ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig = plt.figure(figsize=(17, 9)); gs = fig.add_gridspec(2, 3, width_ratios=[1.1, 1.3, 1])
ax = fig.add_subplot(gs[:, 0]); sc = ax.scatter(clon, clat, c=zmap, cmap="RdBu_r", vmin=-3, vmax=3, s=14); ax.set_aspect(1 / np.cos(np.radians(47)))
plt.colorbar(sc, ax=ax, shrink=.6, label="z of mean ON−PRE contrast (100-km smoothed)"); ax.set_title("per-cell sensitivity to modeled slip (N+C, local onset, pulse window)\nred = velocity increase during slip; blue = decrease", fontsize=9)
ax = fig.add_subplot(gs[0, 1])
for c_, col in zip(["footprint", "updip(15-30km)", "downdip(45-55km)", "ring(>=100km)"], ["#1d4ed8", "#dc2626", "#059669", "#94a3b8"]):
    ys = np.array([aligned(Yd, e, e.cls[c_], c_ == "footprint") for e in PR]); cc = np.nanmean(ys, 0); ax.plot(LAGS, cc - np.nanmean(cc[L0 - 90:L0 - 30]), color=col, lw=1.4, label=c_)
ax.axhline(0, color="k", lw=.5); ax.axvspan(*PRE, color="#94a3b8", alpha=.15); ax.legend(fontsize=8); ax.set_title("interface δβ/β composites by cell class (N+C, rel. PRE)", fontsize=9); ax.set_xlabel("days from onset"); ax.set_ylabel("%")
ax = fig.add_subplot(gs[1, 1])
for p_, col in zip(["interface(deseason)", "family TOTAL(deseason)", "family site-removed(deseason)"], ["#1d4ed8", "#7c3aed", "#0891b2"]):
    ys = np.array([aligned(PROD[p_], e, e.fp, True) for e in PR]); cc = np.nanmean(ys, 0); ax.plot(LAGS, cc - np.nanmean(cc[L0 - 90:L0 - 30]), color=col, lw=1.4, label=p_)
ss = np.nanmean([sites_series(e) for e in PR], 0); ax.plot(LAGS, ss - np.nanmean(ss[L0 - 90:L0 - 30]), color="#f59e0b", lw=1.4, label="SITES (near-station)")
ax.axhline(0, color="k", lw=.5); ax.legend(fontsize=7); ax.set_title("footprint composites by product (N+C, rel. PRE)", fontsize=9); ax.set_xlabel("days from onset"); ax.set_ylabel("%")
ax = fig.add_subplot(gs[0, 2]); ax.scatter(sz[ok], cs_[ok], c=["#1d4ed8" if e.seg == "N" else "#059669" for e in np.array(EPS)[m][ok]]); ax.axhline(0, color="k", lw=.5)
ax.set_xlabel("footprint slip sum (mm)"); ax.set_ylabel("interface contrast ON−PRE (%)"); ax.set_title(f"contrast vs episode size (Spearman {rs:+.2f})", fontsize=9)
ax = fig.add_subplot(gs[1, 2]); t12 = tab.head(12); ax.barh(range(12), t12.z.values[::-1], color=["#dc2626" if v < 0.05 else "#94a3b8" for v in t12.p_fwer.values[::-1]])
ax.set_yticks(range(12)); ax.set_yticklabels([f"{r['product'][:14]}|{r['cells'][:9]}|{r['window'][:10]}|{r['subset']}" for _, r in t12.iloc[::-1].iterrows()], fontsize=6); ax.set_xlabel("z"); ax.set_title(f"top-12 of {S_obs.size} hypotheses (red = FWER p<0.05); best FWER p = {tab.p_fwer.min():.2f}", fontsize=9)
fig.suptitle("Pattern search: coda dv/v vs modeled slow slip (max-statistic family-wise control)", fontsize=11)
fig.savefig("figures/sse_pattern_search.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_pattern_search.png, pattern_search.csv/.npz")
