#!/usr/bin/env python
"""PER-PATCH response table (user 2026-09-09: "the goal is to determine changes; the response may differ between patches, in time,
and between north / central / south"). Reads sse_patch_prepost.py + sse_patch_prepost_model.py outputs.
Per patch: dv/v change inside the patch (DUR-PRE and POST-PRE of the patch median; also minus ring), its OWN null sd / z / p from the
per-patch shift-null draws, and the model changes inside the same patch (slip in window, Coulomb, shear, normal, dilatation, von
Mises; DUR-PRE and POST-PRE), plus a "sensitivity" = dv/v change per kPa of shear-stress change.
Segments: north >= 46N, central 44-46N, south < 44N (centroid). Heterogeneity: per-segment counts of significant patches by sign,
Kruskal-Wallis of the change across segments, Spearman trend of the change with time within each segment, and the spread of
per-patch changes vs the spread expected from the null (over-dispersion = genuinely different responses).
Outputs: res_catalog_sse/patch_response_table.csv, patch_response_heterogeneity.csv, figures/sse_patch_response.png
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy.stats import kruskal, spearmanr
SD = "fault_tomography/inversion/res_catalog_sse"
T = pd.read_csv(f"{SD}/patch_prepost_patches.csv", parse_dates=["start"]); zn = np.load(f"{SD}/patch_prepost_null.npz", allow_pickle=True); pids = list(zn["pids"]); PRODS = list(zn["products"])
MD = pd.read_csv(f"{SD}/patch_prepost_model.csv"); cat = pd.read_csv(f"{SD}/slip_patches_5mm.csv").set_index("id")
seg = lambda la: "north" if la >= 46 else ("central" if la >= 44 else "south")
KEYM = {"slip_mm30_dur": ("slip_thrust_inc30", "DUR"), "slip_cum_mm_post": ("slip_thrust_cum", "d_post"), "coulomb_kPa_dur": ("coulomb", "d_dur"), "shear_kPa_dur": ("shear_stress", "d_dur"), "normal_kPa_dur": ("normal_stress", "d_dur"), "dilat_nstr_dur": ("dilatation", "d_dur"), "vonmises_kPa_dur": ("von_mises_stress", "d_dur"), "coulomb_kPa_post": ("coulomb", "d_post"), "shear_kPa_post": ("shear_stress", "d_post")}
rows = []
for pn in PRODS:
    Tp = T[T["product"] == pn].set_index("id").reindex(pids)
    for j, pid in enumerate(pids):
        r = dict(product=pn, id=pid, start=Tp.start[pid].date(), year=Tp.start[pid].year, seg=seg(cat.lat_centroid[pid]), lat=cat.lat_centroid[pid], Mw=cat.Mw[pid], dur_d=int(Tp.dur_d[pid]), n_in=int(Tp.n_in[pid]),
                 PRE_median=Tp.PRE_median[pid], DUR_median=Tp.DUR_median[pid], POST_median=Tp.POST_median[pid], PRE_std=Tp.PRE_std[pid], DUR_std=Tp.DUR_std[pid])
        for k in ("d_median_dur", "d_median_post", "d_median_post2", "d_median_dur_minus_ring", "d_median_post_minus_ring", "std_ratio_dur"):
            nd = zn[f"{pn}|{k}"][:, j]; nd = nd[np.isfinite(nd)]; o = Tp[k][pid]; r[k] = o; r[k + "_nnull"] = len(nd)
            r[k + "_nullsd"] = nd.std() if len(nd) > 20 else np.nan; r[k + "_z"] = (o - nd.mean()) / nd.std() if len(nd) > 20 and nd.std() > 0 else np.nan; r[k + "_p"] = float(np.mean(np.abs(nd - nd.mean()) >= abs(o - nd.mean()))) if len(nd) > 20 and np.isfinite(o) else np.nan
        for nm, (f, ch) in KEYM.items():
            m = MD[(MD.field == f) & (MD.id == pid)]; r[nm] = float(m[ch].iloc[0]) if len(m) else np.nan
        r["sens_dvv_per_kPa_shear"] = r["d_median_dur"] / r["shear_kPa_dur"] if np.isfinite(r["shear_kPa_dur"]) and abs(r["shear_kPa_dur"]) > 0.5 else np.nan
        rows.append(r)
R = pd.DataFrame(rows).sort_values(["product", "seg", "start"]); R.to_csv(f"{SD}/patch_response_table.csv", index=False)
# heterogeneity
H = []
for pn in PRODS:
    for mwmin in (6.2, 6.0):
        Rp = R[(R["product"] == pn) & (R.Mw >= mwmin)]
        for k in ("d_median_dur", "d_median_post", "d_median_dur_minus_ring"):
            groups = [Rp[Rp.seg == s][k].dropna().values for s in ("north", "central", "south")]
            kw = kruskal(*[g for g in groups if len(g) >= 3]).pvalue if sum(len(g) >= 3 for g in groups) >= 2 else np.nan
            for s, g in zip(("north", "central", "south", "all"), groups + [Rp[k].dropna().values]):
                sub = Rp if s == "all" else Rp[Rp.seg == s]
                if len(sub) < 3: continue
                sig = sub[sub[k + "_p"] < 0.05]; tr = spearmanr(pd.to_datetime(sub.start).map(lambda d_: d_.toordinal()), sub[k]) if len(sub) >= 6 else None
                H.append(dict(product=pn, Mw_min=mwmin, metric=k, seg=s, n=len(sub), mean=sub[k].mean(), median=sub[k].median(), sd_across_patches=sub[k].std(), null_sd_median=sub[k + "_nullsd"].median(),
                              overdispersion=sub[k].std() / sub[k + "_nullsd"].median(), n_sig_neg=int((sig[k] < 0).sum()), n_sig_pos=int((sig[k] > 0).sum()), frac_neg=float((sub[k] < 0).mean()),
                              trend_rho=tr.correlation if tr else np.nan, trend_p=tr.pvalue if tr else np.nan, kruskal_p_segments=kw))
HT = pd.DataFrame(H); HT.to_csv(f"{SD}/patch_response_heterogeneity.csv", index=False)
pd.set_option("display.width", 260)
for pn in PRODS:
    print(f"\n=== {pn}: per-patch responses, Mw>=6.2 (z/p from each patch's own shift null) ===")
    Rp = R[(R["product"] == pn) & (R.Mw >= 6.2)]
    print(Rp[["id", "start", "seg", "Mw", "dur_d", "n_in", "PRE_median", "DUR_median", "POST_median", "d_median_dur", "d_median_dur_z", "d_median_dur_p", "d_median_post", "d_median_post_p", "d_median_dur_minus_ring", "std_ratio_dur", "slip_mm30_dur", "coulomb_kPa_dur", "shear_kPa_dur", "normal_kPa_dur", "dilat_nstr_dur", "sens_dvv_per_kPa_shear"]].round(3).to_string(index=False))
    print(f"\n--- heterogeneity, {pn} ---"); print(HT[HT["product"] == pn][["Mw_min", "metric", "seg", "n", "mean", "median", "sd_across_patches", "null_sd_median", "overdispersion", "n_sig_neg", "n_sig_pos", "frac_neg", "trend_rho", "trend_p", "kruskal_p_segments"]].round(3).to_string(index=False))
# figure: per patch by segment and time
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
COL = {"north": "#1d4ed8", "central": "#059669", "south": "#ea580c"}
fig, ax = plt.subplots(3, 2, figsize=(18, 12), constrained_layout=True, sharex="col")
for c, pn in enumerate(PRODS[:2]):
    Rp = R[(R["product"] == pn) & (R.Mw >= 6.0)]
    for i, (k, lab) in enumerate((("d_median_dur", "DUR − PRE (patch median dv/v, %)"), ("d_median_post", "POST − PRE (%)"), ("d_median_dur_minus_ring", "DUR − PRE, patch minus ring (%)"))):
        a = ax[i, c]
        for s in ("north", "central", "south"):
            g = Rp[Rp.seg == s]; x = pd.to_datetime(g.start); big = g.Mw >= 6.2
            a.errorbar(x, g[k], yerr=1.96 * g[k + "_nullsd"], fmt="none", ecolor=COL[s], alpha=.3, lw=1)
            a.scatter(x[big], g[k][big], s=40 + 60 * (g.Mw[big] - 6.2), color=COL[s], edgecolor="k" if True else None, linewidth=[1.2 if p_ < 0.05 else 0.2 for p_ in g[k + "_p"][big]], label=f"{s} (Mw≥6.2: {big.sum()})", zorder=3)
            a.scatter(x[~big], g[k][~big], s=12, color=COL[s], alpha=.5, marker="x", zorder=2)
        a.axhline(0, color="k", lw=.6); a.set_ylabel(lab); a.set_title(f"{pn}: per-patch change inside the patch (bars = own null 95 %; bold edge = p<0.05; x = Mw 6.0-6.2)", fontsize=9)
        if i == 0: a.legend(fontsize=8, ncol=3)
fig.suptitle("Per-patch dv/v response inside SSE slip patches by segment and time (unbridged catalog, post-2011)", fontsize=11)
fig.savefig("figures/sse_patch_response.png", dpi=130); print("wrote figures/sse_patch_response.png, patch_response_table.csv, patch_response_heterogeneity.csv")
