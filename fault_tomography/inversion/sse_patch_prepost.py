#!/usr/bin/env python
"""PRE / DURING / POST dv/v INSIDE the SSE slip patches (user 2026-09-09; unbridged reference catalog slip_patches_5mm).
Patch = union footprint of a space-time component (mesh triangles); inside cells = our 0.1-deg interface cells whose nearest
triangle is in the footprint. Products (on the sse daily grid): interface delta-beta/beta (production daily inversion, per-cell
patch-mapped transmission correction as in the explorer; T<=0.2 uncorrected) and the family 2-4 s product. Both deseasoned per
cell with a month-of-year climatology computed on days when the cell is NOT inside an active patch.
Windows in the common trailing-30-d frame (catalog start ts / end te are the first/last days the 30-d slip window exceeds 5 mm):
  PRE  = [ts-120, ts-31]  (dv/v windows entirely before the first above-threshold slip window)
  DUR  = [ts, te]
  POST = [te+31, te+120],  POST2 = [te+121, te+240]
Per patch x window (pooled over inside cell-days): mean, median, std, MAD, n; spatial std of per-cell means; temporal std of the
patch-mean series. Changes: DUR-PRE, POST-PRE (mean and median); std ratio DUR/PRE. Ring control = cells 100-300 km from the
footprint (not inside). Null = rigid circular time shift of the dv/v field (>= 1 yr, same shift for every patch in a draw,
NDRAW draws): p-values for the set-mean change and per-patch 95 % bands. Aligned stacks (patch-mean series vs days since ts).
Sets: post-2011, not minor, not broad_pulse: Mw>=6.2 all / north (centroid >= 46N) / south (< 46N); Mw>=6.0 all.
Outputs: res_catalog_sse/patch_prepost_patches.csv, patch_prepost_summary.csv, patch_prepost_null.npz (per-patch null draws + geometry, for
sse_patch_prepost_model.py), figures/sse_patch_prepost.png
Env: CATSUF ("" reference; e.g. _gap7_r50), NDRAW (300).
"""
import os, warnings, numpy as np, pandas as pd; warnings.filterwarnings("ignore")
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; SUF = os.environ.get("CATSUF", ""); NDRAW = int(os.environ.get("NDRAW", "300")); rng = np.random.RandomState(7)
cat = pd.read_csv(f"{SD}/slip_patches_5mm{SUF}.csv", parse_dates=["start", "end", "peak_date"])
zl = np.load(f"{SD}/slip_patches_5mm{SUF}.npz", allow_pickle=True); L = zl["label"].T; days = pd.to_datetime(zl["days"]); G = len(days); gi = pd.Series(np.arange(G), index=days)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); cells = zc["cell"]; clat, clon, dep, ip, valid = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"], zc["patch"], zc["valid"]; M = len(cells)
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); assert list(zp["cell"]) == list(cells)
Tc = pd.read_csv(f"{RD}/cell_T_patchmapped.csv").set_index("cell").T_cell.reindex([str(c) for c in cells]).values; Tc = np.where(np.isfinite(Tc) & (Tc > 0.2), Tc, 1.0)
Y = pd.DataFrame(np.asarray(zp["MODEL"], float) / Tc[None, :], index=pd.to_datetime(zp["days"])).reindex(days).values
zf = np.load(f"{SD}/product_family24.npz", allow_pickle=True); assert (pd.to_datetime(zf["days"]) == days).all(); F24 = np.asarray(zf["MODEL"], float)
active = (L[ip] > 0).T                                             # (G, M): cell inside an active patch that day
month = days.month.values
def deseason(X):
    Xm = np.where(active, np.nan, X); out = X.copy()
    for m in range(1, 13): s = month == m; out[s] -= np.nanmean(Xm[s], 0)[None]
    return out
PROD = {"interface δβ/β (corrected)": deseason(Y), "family 2-4 s": deseason(F24)}
sel0 = cat[(cat.start >= "2011-01-01") & ~cat.minor & ~cat.broad_pulse]
SETS = {"Mw>=6.2 all": sel0[sel0.Mw >= 6.2], "Mw>=6.2 north": sel0[(sel0.Mw >= 6.2) & (sel0.lat_centroid >= 46)], "Mw>=6.2 south": sel0[(sel0.Mw >= 6.2) & (sel0.lat_centroid < 46)], "Mw>=6.0 all": sel0[sel0.Mw >= 6.0]}
WIN = {"PRE": (-120, -31, "s"), "DUR": (0, 0, "se"), "POST": (31, 120, "e"), "POST2": (121, 240, "e")}
def wrange(ts, te, w):
    a, b, ref = WIN[w]
    if ref == "se": return ts, te
    r = ts if ref == "s" else te; return max(r + a, 0), min(r + b, G - 1)
def pstats(X, cs, idx):
    blk = X[idx][:, cs]; v = blk[np.isfinite(blk)]
    if len(v) < 20: return dict(mean=np.nan, median=np.nan, std=np.nan, mad=np.nan, n=len(v), sp_std=np.nan, t_std=np.nan)
    cm = np.nanmean(blk, 0); tm = np.nanmean(blk, 1)
    return dict(mean=v.mean(), median=np.median(v), std=v.std(), mad=np.median(np.abs(v - np.median(v))) * 1.4826, n=len(v), sp_std=np.nanstd(cm), t_std=np.nanstd(tm))
# patch geometry
P = {}
for _, r in sel0.iterrows():
    ts, te = int(gi[r.start]), int(gi[r.end]); fp = np.unique(np.where(L[:, ts:te + 1] == r.id)[0]); inside = np.where(valid & np.isin(ip, fp))[0]
    if len(inside) < 5: continue
    d = np.sqrt(((XY[:, None, :] - XY[inside][None, :, :]) ** 2).sum(-1)).min(1); ring = np.where(valid & (d >= 100) & (d <= 300))[0]
    P[r.id] = dict(ts=ts, te=te, inside=inside, ring=ring, seg="north" if r.lat_centroid >= 46 else "south", start=r.start, Mw=r.Mw, n_in=len(inside))
print(f"patches with >=5 inside cells: {len(P)} of {len(sel0)} (post-2011, not minor/broad)")
def changes(X, p, shift=0):
    """per-patch window stats; shift = circular time shift of the field (null); windows wrap via modulo indexing"""
    ts, te = p["ts"], p["te"]; out = {}
    for w in WIN:
        a, b = wrange(ts, te, w); idx = (np.arange(a, b + 1) + shift) % G
        out[w] = pstats(X, p["inside"], idx); out[w + "_ring"] = pstats(X, p["ring"], idx)
    return out
rows = []; NULL = {}
for pname, X in PROD.items():
    for pid, p in P.items():
        o = changes(X, p)
        rec = dict(product=pname, id=pid, start=p["start"].date(), seg=p["seg"], Mw=p["Mw"], n_in=p["n_in"], dur_d=p["te"] - p["ts"] + 1)
        for w in WIN:
            for k in ("mean", "median", "std", "mad", "n", "sp_std", "t_std"): rec[f"{w}_{k}"] = o[w][k]
            rec[f"{w}_ring_median"] = o[w + "_ring"]["median"]; rec[f"{w}_ring_mean"] = o[w + "_ring"]["mean"]
        rec["d_mean_dur"] = o["DUR"]["mean"] - o["PRE"]["mean"]; rec["d_median_dur"] = o["DUR"]["median"] - o["PRE"]["median"]
        rec["d_mean_post"] = o["POST"]["mean"] - o["PRE"]["mean"]; rec["d_median_post"] = o["POST"]["median"] - o["PRE"]["median"]; rec["d_median_post2"] = o["POST2"]["median"] - o["PRE"]["median"]
        rec["std_ratio_dur"] = o["DUR"]["std"] / o["PRE"]["std"]; rec["std_ratio_post"] = o["POST"]["std"] / o["PRE"]["std"]
        rec["d_median_dur_minus_ring"] = rec["d_median_dur"] - (o["DUR_ring"]["median"] - o["PRE_ring"]["median"]); rec["d_median_post_minus_ring"] = rec["d_median_post"] - (o["POST_ring"]["median"] - o["PRE_ring"]["median"])
        rows.append(rec)
    # null: same shift for all patches per draw
    nul = {k: np.full((NDRAW, len(P)), np.nan) for k in ("d_mean_dur", "d_median_dur", "d_mean_post", "d_median_post", "d_median_post2", "std_ratio_dur", "d_median_dur_minus_ring", "d_median_post_minus_ring")}
    for dI in range(NDRAW):
        sh = int(rng.randint(365, G - 365))
        for j, (pid, p) in enumerate(P.items()):
            o = changes(X, p, sh)
            nul["d_mean_dur"][dI, j] = o["DUR"]["mean"] - o["PRE"]["mean"]; nul["d_median_dur"][dI, j] = o["DUR"]["median"] - o["PRE"]["median"]
            nul["d_mean_post"][dI, j] = o["POST"]["mean"] - o["PRE"]["mean"]; nul["d_median_post"][dI, j] = o["POST"]["median"] - o["PRE"]["median"]; nul["d_median_post2"][dI, j] = o["POST2"]["median"] - o["PRE"]["median"]
            nul["std_ratio_dur"][dI, j] = o["DUR"]["std"] / o["PRE"]["std"]
            nul["d_median_dur_minus_ring"][dI, j] = nul["d_median_dur"][dI, j] - (o["DUR_ring"]["median"] - o["PRE_ring"]["median"]); nul["d_median_post_minus_ring"][dI, j] = nul["d_median_post"][dI, j] - (o["POST_ring"]["median"] - o["PRE_ring"]["median"])
    NULL[pname] = nul
T = pd.DataFrame(rows); T.to_csv(f"{SD}/patch_prepost_patches.csv", index=False)
np.savez_compressed(f"{SD}/patch_prepost_null.npz", pids=np.array(list(P)), products=np.array(list(PROD)), **{f"{pn}|{k}": v for pn, nul in NULL.items() for k, v in nul.items()},
                    inside=np.array([P[pid]["inside"] for pid in P], dtype=object), ring=np.array([P[pid]["ring"] for pid in P], dtype=object), ts=np.array([P[pid]["ts"] for pid in P]), te=np.array([P[pid]["te"] for pid in P]))
pids = list(P); summ = []
for pname in PROD:
    Tp = T[T["product"] == pname].set_index("id").reindex(pids)
    for sname, S in SETS.items():
        j = np.array([k for k, pid in enumerate(pids) if pid in set(S.id)]);
        if len(j) == 0: continue
        rec = dict(product=pname, set=sname, n_patches=len(j))
        for w in WIN:
            rec[f"{w}_median_of_medians"] = np.nanmedian(Tp[f"{w}_median"].values[j]); rec[f"{w}_mean_of_means"] = np.nanmean(Tp[f"{w}_mean"].values[j]); rec[f"{w}_std_pooled_median"] = np.nanmedian(Tp[f"{w}_std"].values[j])
        for k in NULL[pname]:
            obs = np.nanmean(Tp[k].values[j]); nd = np.nanmean(NULL[pname][k][:, j], 1); rec[k] = obs; rec[k + "_null_sd"] = np.nanstd(nd); rec[k + "_p"] = float(np.mean(np.abs(nd - np.nanmean(nd)) >= abs(obs - np.nanmean(nd))))
            if k.startswith("d_"): rec[k + "_frac_neg"] = float(np.nanmean(Tp[k].values[j] < 0))
        summ.append(rec)
SM = pd.DataFrame(summ); SM.to_csv(f"{SD}/patch_prepost_summary.csv", index=False)
pd.set_option("display.width", 250)
for pname in PROD:
    print(f"\n=== {pname} ===")
    s = SM[SM["product"] == pname]
    print(s[["set", "n_patches", "PRE_median_of_medians", "DUR_median_of_medians", "POST_median_of_medians", "POST2_median_of_medians", "PRE_std_pooled_median", "DUR_std_pooled_median"]].round(3).to_string(index=False))
    print(s[["set", "d_median_dur", "d_median_dur_null_sd", "d_median_dur_p", "d_median_dur_frac_neg", "d_mean_dur", "d_mean_dur_p", "d_median_post", "d_median_post_p", "d_median_post2", "d_median_post2_p", "std_ratio_dur", "std_ratio_dur_p", "d_median_dur_minus_ring", "d_median_dur_minus_ring_p", "d_median_post_minus_ring_p"]].round(3).to_string(index=False))
# ---------- figure ----------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
REL = np.arange(-150, 271); fig, ax = plt.subplots(2, 3, figsize=(19, 9.5), constrained_layout=True)
for i, (pname, X) in enumerate(PROD.items()):
    Tp = T[T["product"] == pname].set_index("id").reindex(pids)
    a = ax[i, 0]
    for sname, col in (("Mw>=6.2 north", "#1d4ed8"), ("Mw>=6.2 south", "#ea580c")):
        ids = [pid for pid in pids if pid in set(SETS[sname].id)]
        def stack(shift=0):
            S_ = []
            for pid in ids:
                p = P[pid]; idx = (p["ts"] + REL + shift) % G; s_ = np.nanmean(X[idx][:, p["inside"]], 1); S_.append(s_ - np.nanmean(s_[(REL >= -120) & (REL <= -31)]))
            return np.nanmean(S_, 0)
        st = stack(); nl = np.array([stack(int(rng.randint(365, G - 365))) for _ in range(100)])
        a.fill_between(REL, np.nanpercentile(nl, 2.5, 0), np.nanpercentile(nl, 97.5, 0), color=col, alpha=.12); a.plot(REL, st, color=col, lw=1.4, label=f"{sname} ({len(ids)})")
        dm = np.median([P[pid]["te"] - P[pid]["ts"] for pid in ids]); a.axvline(dm, color=col, ls=":", lw=1)
    a.axvspan(-120, -31, color="#e5e7eb", alpha=.5); a.axvline(0, color="k", lw=.8); a.axhline(0, color="k", lw=.5); a.set_xlabel("days since patch start (first day slip window > 5 mm)"); a.set_ylabel("patch-mean dv/v rel. PRE (%)"); a.set_title(f"{pname}: aligned stack inside the patch (band = shift null 95 %; dotted = median end)", fontsize=9); a.legend(fontsize=8)
    a = ax[i, 1]
    for seg, col in (("north", "#1d4ed8"), ("south", "#ea580c")):
        s_ = Tp[(Tp.seg == seg) & (Tp.Mw >= 6.2)]; j = np.array([pids.index(pid) for pid in s_.index]); lo = np.nanpercentile(NULL[pname]["d_median_dur"][:, j], 2.5, 0); hi = np.nanpercentile(NULL[pname]["d_median_dur"][:, j], 97.5, 0)
        x = pd.to_datetime(s_.start); a.vlines(x, lo, hi, color=col, alpha=.35, lw=3); a.plot(x, s_.d_median_dur, "o", color=col, ms=5, label=f"{seg}: DUR−PRE median")
        a.plot(x, s_.d_median_post, "s", mfc="none", color=col, ms=5, label=f"{seg}: POST−PRE median")
    a.axhline(0, color="k", lw=.5); a.set_ylabel("change of patch median dv/v (%)"); a.set_title("per patch (Mw ≥ 6.2): change inside the patch; bars = per-patch null 95 %", fontsize=9); a.legend(fontsize=7, ncol=2)
    a = ax[i, 2]
    s_ = Tp[(Tp.Mw >= 6.2)]; data = [s_[f"{w}_median"].dropna().values for w in WIN]
    a.boxplot(data, showfliers=False, widths=.5); a.set_xticks(range(1, len(WIN) + 1)); a.set_xticklabels(list(WIN))
    for k, w in enumerate(WIN): a.plot(np.full(len(data[k]), k + 1) + rng.uniform(-.12, .12, len(data[k])), data[k], ".", color="#475569", ms=4, alpha=.6)
    a.axhline(0, color="k", lw=.5); a.set_ylabel("patch median dv/v in window (%)"); a.set_title(f"window medians across patches (Mw ≥ 6.2, n={len(s_)}); std ratio DUR/PRE median = {np.nanmedian(s_.std_ratio_dur):.2f}", fontsize=9)
fig.suptitle("dv/v inside SSE slip patches (unbridged catalog, > 5 mm/30 d), PRE [−120,−31] · DUR [start,end] · POST [+31,+120] · POST2 [+121,+240] in the common trailing-30-d frame", fontsize=10)
fig.savefig("figures/sse_patch_prepost.png", dpi=130); print("wrote figures/sse_patch_prepost.png, patch_prepost_patches.csv, patch_prepost_summary.csv")
