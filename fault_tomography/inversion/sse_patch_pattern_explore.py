#!/usr/bin/env python
"""EXTENSIVE per-patch pattern exploration (user 2026-09-10): does the dv/v response inside a slip patch depend on (a) the
SEQUENCE the patch sits in (time since the previous event on the same site, prior slip history, concurrent/neighbouring
activity, chain position, previous response at the site, position in time) and (b) the UPPER PLATE (terrane class by
latitude, with Siletzia's core separated), plus patch properties (Mw, duration, peak rate, depth, migration), PRE level and
overlap contamination? Inputs: slip_patches_5mm.csv/.npz (unbridged), patch_response_table.csv (+ patch_prepost_null.npz
per-patch null draws), patch_dvv_maps_summary(_Mw6).csv, components_at_cells_trailing30.npz (prior cumulative slip).
Responses (interface, corrected, deseasoned): d_median_dur (DUR-PRE, PRE 120-31 d), d_median_post, d_median_dur_minus_ring,
std_ratio_dur (all with per-patch shift-null draws) and the 60-d-PRE DUR_minus_PRE from the map tables (permutation null).
Tests per set (Mw>=6.2 / Mw>=6.0 / all 217): Spearman for numeric features and Kruskal-Wallis for classes, p from the null
draws (response replaced by its null draw), FWER = max-|z| over features; multivariate CV R^2 (ridge + gradient boosting)
vs null draws; site consistency (sign agreement + ICC) vs null; lag-1 autocorrelation of responses at a site; terrane x
sequence interactions; overlap-controlled repeats. Terrane classes (centroid latitude proxy; APPROXIMATE):
  <43.2 Klamath/Franciscan | 43.2-46.8 Siletzia core (OR/SW WA Coast Range, thick mafic, low permeability) |
  46.8-48.6 Olympic accretionary complex + Crescent rim | >48.6 Vancouver Island (Wrangellia / Metchosin).
Outputs: res_catalog_sse/patch_pattern_features.csv, patch_pattern_tests.csv, patch_pattern_sites.csv,
figures/sse_patch_pattern_explore.png
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy.stats import spearmanr, kruskal, mannwhitneyu
SD = "fault_tomography/inversion/res_catalog_sse"; rng = np.random.RandomState(11); NDRAW_ML = 40
cat = pd.read_csv(f"{SD}/slip_patches_5mm.csv", parse_dates=["start", "end"]); zl = np.load(f"{SD}/slip_patches_5mm.npz", allow_pickle=True); L = zl["label"].T; days = pd.to_datetime(zl["days"]); G = len(days); gi = pd.Series(np.arange(G), index=days)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); ip = zc["patch"]; cum = zc["slip_thrust_cum"].astype(float) * 1e3
zn = np.load(f"{SD}/patch_prepost_null.npz", allow_pickle=True); pids = [int(x) for x in zn["pids"]]; INS = list(zn["inside"]); TS = zn["ts"]; TE = zn["te"]; PN = "interface δβ/β (corrected)"
NULL = {k: zn[f"{PN}|{k}"] for k in ("d_median_dur", "d_median_post", "d_median_dur_minus_ring", "std_ratio_dur")}
R = pd.read_csv(f"{SD}/patch_response_table.csv"); R = R[R["product"] == PN].set_index("id").reindex(pids)
MAPS = pd.concat([pd.read_csv(f"{SD}/patch_dvv_maps_summary_Mw6.csv")]).drop_duplicates("id").set_index("id")
C = cat.set_index("id"); sub = cat[~cat.minor & ~cat.broad_pulse].sort_values("start")          # substantive, all years (sequence context)
FP = {}
for _, r in sub.iterrows():
    ts, te = int(gi[r.start]), int(gi[r.end]); FP[r.id] = set(np.unique(np.where(L[:, ts:te + 1] == r.id)[0]))
# sites: agglomerative clustering of patch centroids (complete linkage, 150 km) -> repeating-footprint groups without the
# "largest footprint swallows the segment" failure of greedy overlap clustering
from scipy.cluster.hierarchy import linkage, fcluster
kmf = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XYs = kmf(sub.lat_centroid.values, sub.lon_centroid.values)
lab = fcluster(linkage(XYs, method="complete"), t=150, criterion="distance"); site_of = dict(zip(sub.id, lab - 1)); sites = [dict(ids=list(sub.id[lab - 1 == k])) for k in range(lab.max())]
site_events = {si: sub[sub.id.isin(s["ids"])].sort_values("start") for si, s in enumerate(sites)}
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XYc = km(sub.lat_centroid.values, sub.lon_centroid.values); sub_xy = dict(zip(sub.id, XYc))
def terrane(la): return "Klamath/Franciscan (<43.2)" if la < 43.2 else ("Siletzia core (43.2-46.8)" if la < 46.8 else ("Olympic/Crescent (46.8-48.6)" if la < 48.6 else "Vancouver Is. (>48.6)"))
TERR_ORDER = ["Klamath/Franciscan (<43.2)", "Siletzia core (43.2-46.8)", "Olympic/Crescent (46.8-48.6)", "Vancouver Is. (>48.6)"]
rows = []
for j, pid in enumerate(pids):
    c = C.loc[pid]; ts, te = int(TS[j]), int(TE[j]); si = site_of.get(pid); ev = site_events[si] if si is not None else None
    f = dict(id=pid, start=c.start, year=c.start.year, month=c.start.month, season=["DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF"][c.start.month - 1],
             Mw=c.Mw, duration_d=c.duration_d, area_km2=c.area_km2, peak_rate=c.peak_mm_per30d, mean_slip_mm=c.mean_slip_mm, migration=c.migration_km_per_d, strike_extent_km=c.strike_extent_km,
             depth_min=c.depth_min_km, depth_max=c.depth_max_km, depth_mid=(c.depth_min_km + c.depth_max_km) / 2, lat=c.lat_centroid, lon=c.lon_centroid, terrane=terrane(c.lat_centroid), seg=R.seg[pid], n_in=int(R.n_in[pid]),
             PRE_level=R.PRE_median[pid], PRE_sd=R.PRE_std[pid], site=si, n_events_site=len(ev) if ev is not None else 1)
    # sequence at the site
    fp_this = FP.get(pid, set()); ovl = sub[[pid2 != pid and len(FP[pid2] & fp_this) >= 0.5 * max(len(fp_this), 1) for pid2 in sub.id]]     # patches overlapping >= 50 % of THIS footprint
    prev = ovl[ovl.end < c.start]; nxt = ovl[ovl.start > c.end]; f["n_overlapping_patches"] = len(ovl)
    f["t_since_prev_site_d"] = (c.start - prev.end.iloc[-1]).days if len(prev) else np.nan; f["t_to_next_site_d"] = (nxt.start.iloc[0] - c.end).days if len(nxt) else np.nan
    f["n_prev_site_3yr"] = int(((c.start - prev.start).dt.days <= 1095).sum()) if len(prev) else 0; f["prev_Mw_site"] = prev.Mw.iloc[-1] if len(prev) else np.nan
    f["prev_response_site"] = R.d_median_dur.get(int(prev.id.iloc[-1]), np.nan) if len(prev) else np.nan; f["cycle_index_site"] = len(prev)
    # prior slip inside the footprint (trailing-frame cumulative thrust slip, mm)
    for yrs in (1, 2, 3):
        t0 = ts - 365 * yrs; f[f"prior_slip_{yrs}yr_mm"] = np.nanmean(cum[ts, INS[j]] - cum[t0, INS[j]]) if t0 >= 0 else np.nan
    # concurrent / neighbouring activity (substantive patches, not same site)
    oth = sub[(sub.id != pid) & (sub.id.map(site_of) != si)]; dstart = (oth.start - c.start).dt.days; dist = np.sqrt(((XYc[oth.index.map(lambda i: sub.index.get_loc(i))] - sub_xy[pid]) ** 2).sum(1))
    f["n_concurrent_60d"] = int((dstart.abs() <= 60).sum()); f["neighbor_prev90"] = bool(((dstart >= -90) & (dstart <= -1) & (dist <= 400)).any()); f["neighbor_next90"] = bool((((oth.start - c.end).dt.days >= 1) & ((oth.start - c.end).dt.days <= 90) & (dist <= 400)).any())
    f["north_prev90"] = bool(((dstart >= -90) & (dstart <= -1) & (oth.lat_centroid > c.lat_centroid + 1)).any()); f["south_prev90"] = bool(((dstart >= -90) & (dstart <= -1) & (oth.lat_centroid < c.lat_centroid - 1)).any())
    f["chain_first"] = not f["neighbor_prev90"]
    # overlap contamination (another substantive patch on the footprint during PRE 120-31 / POST 31-120)
    fp = np.array(sorted(FP.get(pid, set()))); nonminor = set(sub.id)
    def olap(a, b): 
        w = L[fp][:, max(a, 0):min(b, G - 1) + 1]; ids = set(np.unique(w[(w > 0) & (w != pid)])); return bool(ids & nonminor)
    f["pre_overlap"] = olap(ts - 120, ts - 31); f["post_overlap"] = olap(te + 31, te + 120)
    for k in ("d_median_dur", "d_median_post", "d_median_dur_minus_ring", "std_ratio_dur"): f[k] = R[k][pid]; f[k + "_nullsd"] = R[k + "_nullsd"][pid]
    f["DUR_minus_PRE60"] = MAPS.DUR_minus_PRE.get(pid, np.nan); f["POST_minus_PRE60"] = MAPS.POST_minus_PRE.get(pid, np.nan)
    rows.append(f)
F = pd.DataFrame(rows); F.to_csv(f"{SD}/patch_pattern_features.csv", index=False)
NUM = ["n_overlapping_patches", "Mw", "duration_d", "area_km2", "peak_rate", "mean_slip_mm", "migration", "strike_extent_km", "depth_min", "depth_max", "depth_mid", "lat", "lon", "year", "n_in", "PRE_level", "PRE_sd",
       "t_since_prev_site_d", "t_to_next_site_d", "n_prev_site_3yr", "prev_Mw_site", "prev_response_site", "cycle_index_site", "prior_slip_1yr_mm", "prior_slip_2yr_mm", "prior_slip_3yr_mm", "n_concurrent_60d", "n_events_site"]
CATF = ["terrane", "seg", "season", "neighbor_prev90", "neighbor_next90", "north_prev90", "south_prev90", "chain_first", "pre_overlap", "post_overlap"]
RESP = ["d_median_dur", "d_median_post", "d_median_dur_minus_ring", "std_ratio_dur", "DUR_minus_PRE60"]
SETS = {"Mw>=6.2": F.Mw >= 6.2, "Mw>=6.0": F.Mw >= 6.0, "all": np.ones(len(F), bool), "Mw>=6.0 no PRE overlap": (F.Mw >= 6.0) & ~F.pre_overlap}; SETS = {k: np.asarray(v, bool) for k, v in SETS.items()}
def null_draws(resp, sel):
    if resp in NULL: return NULL[resp][:, sel]
    v = F.loc[sel, resp].values; return np.array([rng.permutation(v) for _ in range(300)])                  # permutation null for the map-table metric
tests = []
for sname, sel in SETS.items():
    for resp in RESP:
        y = F.loc[sel, resp].values; ND = null_draws(resp, sel); ok_y = np.isfinite(y)
        if ok_y.sum() < 8: continue
        res = []
        for feat in NUM:
            x = F.loc[sel, feat].values; ok = ok_y & np.isfinite(x)
            if ok.sum() < 8 or np.nanstd(x[ok]) == 0: continue
            rho = spearmanr(x[ok], y[ok]).correlation; nr = np.array([spearmanr(x[ok], ND[d][ok], nan_policy="omit").correlation for d in range(ND.shape[0])]); nr = nr[np.isfinite(nr)]
            if len(nr) < 50: continue
            res.append(dict(set=sname, response=resp, feature=feat, kind="spearman", n=int(ok.sum()), stat=rho, null_sd=nr.std(), z=(rho - nr.mean()) / nr.std(), p=float(np.mean(np.abs(nr - nr.mean()) >= abs(rho - nr.mean()))), _nz=(nr - nr.mean()) / nr.std(), detail=""))
        for feat in CATF:
            g = F.loc[sel, feat].astype(str).values; ok = ok_y; cls = [c_ for c_ in pd.unique(g[ok]) if (g[ok] == c_).sum() >= 3]
            if len(cls) < 2: continue
            def H(v): return kruskal(*[v[ok & (g == c_)] for c_ in cls]).statistic
            h = H(np.nan_to_num(y, nan=np.nanmedian(y))); nh = np.array([H(np.nan_to_num(ND[d], nan=np.nanmedian(ND[d]))) for d in range(ND.shape[0])])
            means = "; ".join(f"{c_}: {np.nanmean(y[ok & (g == c_)]):+.2f} (n{(ok & (g == c_)).sum()})" for c_ in cls)
            res.append(dict(set=sname, response=resp, feature=feat, kind="kruskal", n=int(ok.sum()), stat=h, null_sd=nh.std(), z=(h - nh.mean()) / nh.std(), p=float(np.mean(nh >= h)), _nz=(nh - nh.mean()) / nh.std(), detail=means))
        nmin = min(len(r["_nz"]) for r in res); Z = np.array([r["_nz"][:nmin] for r in res]); zmax = np.max(np.abs(Z), 0)
        for r in res: r["p_fwer"] = float(np.mean(zmax >= abs(r["z"]))); r.pop("_nz"); tests.append(r)
T = pd.DataFrame(tests); T.to_csv(f"{SD}/patch_pattern_tests.csv", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 120)
print(f"features: {len(F)} patches, {len(NUM)} numeric + {len(CATF)} categorical; sites: {len(sites)} (sites with >=3 substantive events: {sum(len(s['ids']) >= 3 for s in sites)})")
print("\nNOTE: PRE_level (and PRE_sd for std_ratio) correlate with the CHANGE by construction (regression to the mean: change = DUR - PRE); the null pairs the OBSERVED PRE with null changes, so those rows are NOT evidence.")
print("\n=== hits p_fwer < 0.10 (any set / response) ==="); hits = T[T.p_fwer < 0.10].sort_values("p_fwer"); print(hits[["set", "response", "feature", "kind", "n", "stat", "z", "p", "p_fwer", "detail"]].round(3).to_string(index=False) if len(hits) else "none")
print("\n=== top 3 per set x response (uncorrected p) ==="); print(T.sort_values("p").groupby(["set", "response"]).head(3).sort_values(["set", "response", "p"])[["set", "response", "feature", "n", "stat", "z", "p", "p_fwer"]].round(3).to_string(index=False))
print("\n=== terrane class means (d_median_dur, Mw>=6.0 and Mw>=6.2) ===")
for sname in ("Mw>=6.2", "Mw>=6.0"):
    s = F[SETS[sname]]; print(sname); print(s.groupby("terrane").agg(n=("id", "size"), dur_minus_pre=("d_median_dur", "mean"), se=("d_median_dur", lambda x: x.std() / np.sqrt(len(x))), post_minus_pre=("d_median_post", "mean"), dur60=("DUR_minus_PRE60", "mean"), frac_pos=("d_median_dur", lambda x: (x > 0).mean()), t_since_prev_med=("t_since_prev_site_d", "median"), prior_slip_2yr=("prior_slip_2yr_mm", "mean")).reindex(TERR_ORDER).round(3).to_string())
# ---------- multivariate: CV R^2 vs null ----------
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import cross_val_score, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
ML = ["Mw", "duration_d", "area_km2", "peak_rate", "depth_mid", "lat", "year", "t_since_prev_site_d", "n_prev_site_3yr", "prior_slip_2yr_mm", "n_concurrent_60d", "cycle_index_site", "n_overlapping_patches"]   # PRE_level / PRE_sd EXCLUDED: they predict the change by regression to the mean
Xd = pd.get_dummies(F[ML + ["terrane"]], columns=["terrane"]).astype(float); Xd["neighbor_prev90"] = F.neighbor_prev90.astype(float); Xd["pre_overlap"] = F.pre_overlap.astype(float)
print("\n=== multivariate CV R^2 (5-fold, 3 repeats) vs null draws — features EXCLUDE PRE level/sd ===")
mlrows = []
for sname in ("Mw>=6.0", "all"):
    sel = SETS[sname]
    for resp in ("d_median_dur", "d_median_post"):
        y = F.loc[sel, resp].values; ok = np.isfinite(y); X = Xd.values[sel][ok]; yy = y[ok]; ND = NULL[resp][:, sel][:, ok]
        for mname, mk in (("ridge", lambda: make_pipeline(SimpleImputer(), StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 12)))), ("gbr", lambda: make_pipeline(SimpleImputer(), GradientBoostingRegressor(n_estimators=150, max_depth=2, learning_rate=0.05, subsample=0.8, random_state=0)))):
            def cvr2(t):
                return np.mean([cross_val_score(mk(), X, t, cv=KFold(5, shuffle=True, random_state=s_), scoring="r2").mean() for s_ in range(3)])
            obs = cvr2(yy); nul = np.array([cvr2(np.nan_to_num(ND[d], nan=np.nanmedian(ND[d]))) for d in range(NDRAW_ML)])
            mlrows.append(dict(set=sname, response=resp, model=mname, n=int(ok.sum()), cv_r2=obs, null_r2_mean=nul.mean(), null_r2_p95=np.percentile(nul, 95), p=float(np.mean(nul >= obs))))
            print(f"{sname:8s} {resp:14s} {mname:5s} n={ok.sum():3d}  CV R2 {obs:+.3f}  (null mean {nul.mean():+.3f}, p95 {np.percentile(nul, 95):+.3f}, p {np.mean(nul >= obs):.2f})", flush=True)
pd.DataFrame(mlrows).to_csv(f"{SD}/patch_pattern_ml.csv", index=False)
# ---------- site consistency + lag-1 autocorrelation ----------
print("\n=== site consistency (sites with >=3 patches in set): sign agreement and ICC vs null ===")
srows = []
for sname in ("Mw>=6.0", "all"):
    sel = SETS[sname]; Fs = F[sel]; ND = NULL["d_median_dur"][:, sel]
    big = [s for s, n in Fs.site.value_counts().items() if n >= 3]
    def agree(v, g): return np.mean([max(np.mean(v[g == s] > 0), np.mean(v[g == s] < 0)) for s in big])
    def icc(v, g):
        gm = np.array([np.nanmean(v[g == s]) for s in big]); wv = np.nanmean([np.nanvar(v[g == s]) for s in big]); bv = np.nanvar(gm); return bv / (bv + wv) if bv + wv > 0 else np.nan
    v = Fs.d_median_dur.values; g = Fs.site.values; a_obs, i_obs = agree(v, g), icc(v, g)
    a_nul = np.array([agree(np.nan_to_num(ND[d], nan=0), g) for d in range(ND.shape[0])]); i_nul = np.array([icc(ND[d], g) for d in range(ND.shape[0])])
    print(f"{sname}: {len(big)} sites | sign agreement {a_obs:.3f} (null {a_nul.mean():.3f} ± {a_nul.std():.3f}, p {np.mean(a_nul >= a_obs):.2f}) | ICC {i_obs:.3f} (null {np.nanmean(i_nul):.3f}, p {np.nanmean(i_nul >= i_obs):.2f})")
    for s in big:
        e = Fs[Fs.site == s].sort_values("start"); srows.append(dict(set=sname, site=s, n=len(e), lat=e.lat.mean(), terrane=e.terrane.iloc[0], mean_dur=e.d_median_dur.mean(), frac_pos=(e.d_median_dur > 0).mean(), mean_post=e.d_median_post.mean(), events=";".join(str(d.date()) for d in e.start)))
    # lag-1 autocorrelation of the response within sites
    pairs = [(i, k) for s in Fs.site.unique() for e in [Fs[Fs.site == s].sort_values("start")] for i, k in zip(e.index[:-1], e.index[1:])]
    if len(pairs) >= 8:
        a = np.array([F.d_median_dur[i] for i, k in pairs]); b = np.array([F.d_median_dur[k] for i, k in pairs]); rho = spearmanr(a, b).correlation
        pos = {pid: n for n, pid in enumerate(pids)}; NDf = NULL["d_median_dur"]
        nr = np.array([spearmanr(NDf[d][[pos[F.id[i]] for i, k in pairs]], NDf[d][[pos[F.id[k]] for i, k in pairs]], nan_policy="omit").correlation for d in range(NDf.shape[0])]); nr = nr[np.isfinite(nr)]
        print(f"   lag-1 autocorrelation of DUR-PRE between consecutive events at a site: rho {rho:+.3f} over {len(pairs)} pairs (null sd {nr.std():.3f}, p {np.mean(np.abs(nr) >= abs(rho)):.2f})")
pd.DataFrame(srows).to_csv(f"{SD}/patch_pattern_sites.csv", index=False)
# ---------- terrane x sequence interactions ----------
print("\n=== within-terrane Spearman of DUR-PRE with sequence features (Mw>=6.0; p from null draws) ===")
sel = SETS["Mw>=6.0"]; ND = NULL["d_median_dur"][:, sel]; Fs = F[sel].reset_index(drop=True)
for terr in TERR_ORDER:
    m = (Fs.terrane == terr).values
    if m.sum() < 8: print(f"{terr}: n={m.sum()} (too few)"); continue
    out = []
    for feat in ("t_since_prev_site_d", "prior_slip_2yr_mm", "n_concurrent_60d", "year", "PRE_level", "Mw", "prev_response_site"):
        x = Fs.loc[m, feat].values; y = Fs.loc[m, "d_median_dur"].values; ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 8: continue
        rho = spearmanr(x[ok], y[ok]).correlation; nr = np.array([spearmanr(x[ok], ND[d][m][ok], nan_policy="omit").correlation for d in range(ND.shape[0])]); nr = nr[np.isfinite(nr)]
        out.append(f"{feat} rho {rho:+.2f} (p {np.mean(np.abs(nr) >= abs(rho)):.2f}, n {ok.sum()})")
    print(f"{terr}: " + " | ".join(out))
# ---------- figure ----------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
COL = dict(zip(TERR_ORDER, ["#ea580c", "#059669", "#1d4ed8", "#7c3aed"])); Fs = F[F.Mw >= 6.0]
fig, ax = plt.subplots(3, 3, figsize=(19, 14), constrained_layout=True)
a = ax[0, 0]; data = [Fs[Fs.terrane == t].d_median_dur.dropna().values for t in TERR_ORDER]; a.boxplot(data, showfliers=False, widths=.5); a.set_xticks(range(1, 5)); a.set_xticklabels([t.split(" (")[0] for t in TERR_ORDER], fontsize=8)
for k, t in enumerate(TERR_ORDER):
    s = Fs[Fs.terrane == t]; a.plot(np.full(len(s), k + 1) + rng.uniform(-.15, .15, len(s)), s.d_median_dur, "o", color=COL[t], ms=4, alpha=.6); b = s[s.Mw >= 6.2]; a.plot(np.full(len(b), k + 1) + rng.uniform(-.15, .15, len(b)), b.d_median_dur, "o", mfc="none", color="k", ms=7)
a.axhline(0, color="k", lw=.5); a.set_ylabel("DUR − PRE (patch median dv/v, %)"); a.set_title("by upper-plate terrane class (Mw≥6.0; black rings = Mw≥6.2)", fontsize=9)
def sc(a, xk, yk, title, xlabel):
    for t in TERR_ORDER:
        s = Fs[Fs.terrane == t]; a.plot(s[xk], s[yk], "o", color=COL[t], ms=4, alpha=.7, label=t.split(" (")[0])
    ok = np.isfinite(Fs[xk]) & np.isfinite(Fs[yk]); r_ = spearmanr(Fs[xk][ok], Fs[yk][ok]).correlation if ok.sum() > 5 else np.nan; a.axhline(0, color="k", lw=.5); a.set_xlabel(xlabel); a.set_ylabel(yk); a.set_title(f"{title} (Spearman {r_:+.2f}, n {ok.sum()})", fontsize=9)
sc(ax[0, 1], "t_since_prev_site_d", "d_median_dur", "response vs time since previous event at the site", "days since previous substantive patch on this site"); ax[0, 1].legend(fontsize=7)
sc(ax[0, 2], "prior_slip_2yr_mm", "d_median_dur", "response vs prior 2-yr cumulative slip in footprint", "mm")
sc(ax[1, 0], "prev_response_site", "d_median_dur", "response vs previous event's response at the same site", "previous DUR − PRE (%)")
sc(ax[1, 1], "n_concurrent_60d", "d_median_dur", "response vs concurrent activity", "other substantive patches starting within ±60 d")
sc(ax[1, 2], "PRE_level", "d_median_dur", "response vs PRE level (regression-to-mean check)", "PRE median dv/v (%)")
a = ax[2, 0]
for t in TERR_ORDER:
    s = Fs[Fs.terrane == t]; a.plot(pd.to_datetime(s.start), s.d_median_dur, "o-", color=COL[t], ms=4, lw=.6, alpha=.8, label=t.split(" (")[0])
a.axhline(0, color="k", lw=.5); a.set_ylabel("DUR − PRE (%)"); a.set_title("response through time by terrane", fontsize=9); a.legend(fontsize=7)
a = ax[2, 1]
for flag, col, lab in (("chain_first", "#1d4ed8", "first in chain (no neighbour in prior 90 d)"), ("neighbor_prev90", "#ea580c", "neighbour patch started in prior 90 d")):
    s = Fs[Fs[flag]]; a.hist(s.d_median_dur.dropna(), bins=np.arange(-2.5, 2.6, .25), alpha=.5, color=col, label=f"{lab} (n {len(s)}, mean {s.d_median_dur.mean():+.2f})")
a.axvline(0, color="k", lw=.5); a.set_xlabel("DUR − PRE (%)"); a.set_title("chain position", fontsize=9); a.legend(fontsize=7)
a = ax[2, 2]; ss = pd.DataFrame(srows); ss = ss[ss.set == "Mw>=6.0"] if len(ss) else ss
if len(ss): a.scatter(ss.lat, ss.mean_dur, s=30 * ss.n, c=[COL[t] for t in ss.terrane], edgecolor="k"); a.axhline(0, color="k", lw=.5); a.set_xlabel("site latitude"); a.set_ylabel("site mean DUR − PRE (%)"); a.set_title(f"sites with ≥3 events (size = n events): mean response", fontsize=9)
fig.suptitle("Per-patch pattern exploration: sequence context and upper-plate terrane vs dv/v response inside the patch (interface δβ/β; unbridged catalog, post-2011)", fontsize=11)
fig.savefig("figures/sse_patch_pattern_explore.png", dpi=130); print("\nwrote figures/sse_patch_pattern_explore.png, patch_pattern_{features,tests,ml,sites}.csv")
