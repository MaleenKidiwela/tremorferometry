"""Q8e analysis: rank WA-NW variants on backward preference; tides for WA-NW; family location table; figure."""
import json, sys, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd, q8_rtr_lib as L
OUT = f"{L.OUT}/tables"; FIG = f"{L.OUT}/figures"; RNG = np.random.default_rng(12)
S = pd.read_parquet(f"{OUT}/streaks_all_q8e_d2.parquet"); M = pd.read_parquet(f"{OUT}/members_data_q8e_d2.parquet")
if os.path.exists(f"{OUT}/streaks_all_q8e_posthoc_d2.parquet"):
    S = pd.concat([S, pd.read_parquet(f"{OUT}/streaks_all_q8e_posthoc_d2.parquet")]); M = pd.concat([M, pd.read_parquet(f"{OUT}/members_data_q8e_posthoc_d2.parquet")])
VARS = ["V0", "Pa", "Pb", "Pc", "F1", "F2", "F3", "A1", "Pd", "Pe"]
LAB = {"V0": "baseline", "Pa": "one family per source", "Pb": "quality-A sources", "Pc": "one family per quality-A source",
       "F1": "tremor front (direction + onset)", "F2": "single-front episodes", "F3": "tremor front, single-front episodes", "A1": "slab-contour axis 345.2°", "Pd": "POST-HOC counted once (dup_D)", "Pe": "POST-HOC counted once + quality-A sources"}
TIDE = f"{OUT}/tidal_stress_nw_sources_hourly.parquet"
geo = pd.read_csv(f"{OUT}/nw_fault_geometry.csv"); src2fid = geo.set_index("source").fid
RES = {}

# ---- 1. family location table (diagnostic 1)
F = pd.read_csv(f"{qd.ROOT}/data/catalogue/families_v1_3.csv"); F = F[F.region == "Royer_NW_WA"]
setup = pd.read_csv(f"{qd.ROOT}/lfe_extension_wa/results/wa_prod/setup/families.csv"); setup = setup[setup.reg == "NW"]
SG = pd.read_csv(f"{qd.ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3.csv")
AC = pd.read_csv(f"{qd.ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3_array_confirmed.csv").rename(columns={"group": "source"})
f2s = SG.merge(AC, on="fid").set_index("fid").source
Sx = pd.read_csv(f"{qd.R}/s3/sources_nw.csv").set_index("source")
tab = setup[["fid", "family_uid", "lat_use", "lon_use", "location_quality", "depth_km"]].merge(
    F[["family_uid", "location_method", "erh", "erz", "located", "representative", "dup_group", "nw_copy_of", "nw_link_type", "nw_source_group", "nw_array_group"]],
    on="family_uid", how="left")
tab["source"] = tab.fid.map(f2s); tab["source_quality"] = tab.source.map(Sx.quality); tab["is_best_fid"] = tab.fid.isin(Sx.best_fid.astype(int))
tab["source_as_km"] = tab.source.map(Sx.as_km); tab["source_dd_km"] = tab.source.map(Sx.dd_km)
tab.to_csv(f"{OUT}/q8e_nw_family_locations.csv", index=False)
RES["locations"] = dict(n_setup_families=len(tab), with_source=int(tab.source.notna().sum()), n_sources=int(tab.source.nunique()),
                        quality=tab.location_quality.value_counts().to_dict(), erh_median=float(tab.erh.median()), erh_gt10=int((tab.erh > 10).sum()),
                        copies=int(tab.nw_copy_of.notna().sum()), link=tab.nw_link_type.value_counts(dropna=False).to_dict(),
                        method=tab.location_method.value_counts(dropna=False).to_dict(),
                        source_quality=Sx.quality.value_counts().to_dict(), source_dd_range=[float(Sx.dd_km.min()), float(Sx.dd_km.max())],
                        sources_dd_abs_gt30=int((Sx.dd_km.abs() > 30).sum()))
print(json.dumps(RES["locations"], indent=0, default=str))

# ---- tides helper
def tides(D, ev, E):
    if len(D) == 0 or not os.path.exists(TIDE):
        return dict(n=len(D))
    Ei = E.set_index("ep"); res = []
    for r in D.itertuples():
        w = ev.loc[M[M.uid == r.uid].row.values].unit.map(src2fid).value_counts()
        e = Ei.loc[r.ep]
        T = pd.read_parquet(TIDE, columns=["fid", "time", "tau_updip_Pa"], filters=[("fid", "in", [int(f) for f in w.index]), ("time", ">=", e.w0), ("time", "<", e.w1 + pd.Timedelta(hours=2))])
        Pv = T.pivot(index="time", columns="fid", values="tau_updip_Pa"); Pv = Pv[[f for f in w.index if f in Pv.columns]]
        ww = w.reindex(Pv.columns).values.astype(float); ser = pd.Series(Pv.values @ ww / ww.sum(), index=Pv.index)
        grid = pd.date_range(e.w0, e.w1, freq="10min"); tsv = np.interp(grid.astype("int64"), ser.index.astype("int64"), ser.values)
        start = e.w0 + pd.Timedelta(hours=r.t_first); n = max(int(round((r.t_last - r.t_first) * 6)), 1); i0 = int(np.searchsorted(grid, start))
        g = ev[ev.ep == r.ep]; days = g.time.dt.floor("D").value_counts(); act = days[days >= 0.01 * len(g)].index
        ok = np.isin(grid.floor("D"), act) & (np.arange(len(grid)) + n < len(grid))
        mins = pd.Series(tsv).rolling(n + 1).min().shift(-n).values
        res.append(dict(all_pos=bool((tsv[i0:i0 + n + 1] > 0).all()), p_all=float(np.nanmean(mins[ok] > 0))))
    T = pd.DataFrame(res); k = len(T) / 5
    return dict(n_pooled=len(T), frac_all=float(T.all_pos.mean()), exp_all=float(T.p_all.mean()),
                p_eff=float(stats.binomtest(int(round(T.all_pos.mean() * k)), max(int(round(k)), 1), T.p_all.mean(), alternative="greater").pvalue))

# ---- 2. ranking table
rows = []
for v in VARS:
    cat = f"q8e_{v}"; ev, E = qd.q8e_nw(v, "data"); evc, Ec = qd.q8e_nw(v, "c1")
    d = S[(S.cat == cat) & (S["var"] == "raw|data")]; sc = S[(S.cat == cat) & (S["var"] == "raw|c1") & (S["mode"] == "pw")]
    pw = d[d["mode"] == "pw"]; pc2 = d[d["mode"] == "pw_c2"]; pc1 = d[d["mode"] == "pw_c1"]
    R, Fw = (pw.cls == "RTR").sum(), (pw.cls == "forward").sum()
    R2, F2 = (pc2.cls == "RTR").sum(), (pc2.cls == "forward").sum()
    share, share_c = R / max(R + Fw, 1), R2 / max(R2 + F2, 1)
    k, n = int(round(R / 5)), int(round((R + Fw) / 5))
    p_bin = float(stats.binomtest(k, n, share_c, alternative="greater").pvalue) if n else 1.0
    per = pw[pw.cls == "RTR"].groupby("ep").size().reindex(E.ep, fill_value=0) / 5
    perc = sc[sc.cls == "RTR"].groupby("ep").size().reindex(E.ep, fill_value=0) / 5
    diff = (per - perc).values; ratio = per.sum() / max(perc.sum(), 0.1)
    wp = float(stats.wilcoxon(diff, alternative="greater").pvalue) if np.any(diff != 0) else 1.0
    row = dict(variant=v, label=LAB[v], n_eps=len(E), n_families=int(ev.fam.nunique()), n_sources=int(ev.unit.nunique()),
               rtr_per_run=R / 5, forward_per_run=Fw / 5, backward_share=share, posperm_rtr=R2 / 5, posperm_fwd=F2 / 5, posperm_share=share_c,
               p_backward_vs_posperm=p_bin, scr_rtr_per_run=float(perc.sum()), ratio_data_scr=float(ratio), wilcoxon_p=wp,
               time_shuffled_rtr=(pc1.cls == "RTR").sum() / 5,
               span_data=float(pw[pw.cls.isin(["RTR", "forward"])].span_ev.median()), span_posperm=float(pc2[pc2.cls.isin(["RTR", "forward"])].span_ev.median()),
               speed_data=float(pw[pw.cls.isin(["RTR", "forward"])].v_ev.abs().median()), speed_posperm=float(pc2[pc2.cls.isin(["RTR", "forward"])].v_ev.abs().median()))
    row["cond1_share_gt_ctrl"] = bool(share > share_c); row["cond2_p_lt_001"] = bool(p_bin < 0.01)
    row["cond3_data_gt_scr"] = bool(ratio >= 2 and wp < 0.05)
    row["FIXED"] = bool(row["cond1_share_gt_ctrl"] and row["cond2_p_lt_001"] and row["cond3_data_gt_scr"])
    tt = {cls: tides(pw[pw.cls == cls], ev, E) for cls in ["RTR", "forward", "surge_only"]}
    tt["scrambled_RTR"] = tides(sc[sc.cls == "RTR"], evc, Ec)
    row["tides"] = tt
    rows.append(row); print(v, {k_: row[k_] for k_ in ["rtr_per_run", "forward_per_run", "backward_share", "posperm_share", "p_backward_vs_posperm", "ratio_data_scr", "wilcoxon_p", "FIXED"]}, flush=True)
A = pd.DataFrame(rows).sort_values("p_backward_vs_posperm")
A.drop(columns="tides").to_csv(f"{OUT}/q8e_variants_ranked.csv", index=False)
RES["variants"] = rows
json.dump(RES, open(f"{L.OUT}/results_q8e.json", "w"), indent=1, default=str)

# ---- figure
fig, axs = plt.subplots(1, 2, figsize=(13, 4.3))
a = axs[0]; xs = np.arange(len(A))
a.bar(xs - 0.2, A.backward_share, 0.38, color="#2a78d6", label="data (P-weighted)")
a.bar(xs + 0.2, A.posperm_share, 0.38, color="#b5b4ae", label="position-permuted control")
for i, r in enumerate(A.itertuples()):
    a.text(i - 0.2, r.backward_share + 0.02, f"{r.rtr_per_run:.0f}/{r.forward_per_run:.0f}", ha="center", fontsize=7)
a.axhline(0.5, color="k", lw=0.8, ls=":"); a.set_ylim(0, 1); a.set_xticks(xs); a.set_xticklabels(A.variant, fontsize=8)
a.set_ylabel("backward share = RTR / (RTR + forward)"); a.set_title("WA-NW variants (labels: RTR / forward per run)", loc="left"); a.legend(frameon=False, fontsize=8)
a = axs[1]
lab = [("RTR", "#2a78d6"), ("forward", "#eda100"), ("surge_only", "#1baf7a"), ("scrambled_RTR", "#eb6834")]
w = 0.2
for j, (cls, col) in enumerate(lab):
    obs = [r["tides"][cls].get("frac_all", np.nan) for r in A.to_dict("records")]
    exp = [r["tides"][cls].get("exp_all", np.nan) for r in A.to_dict("records")]
    a.bar(xs + (j - 1.5) * w, obs, w, color=col, label=f"{cls} observed")
    a.scatter(xs + (j - 1.5) * w, exp, color="k", s=10, zorder=3, label="expected" if j == 0 else None)
a.set_xticks(xs); a.set_xticklabels(A.variant, fontsize=8); a.set_ylim(0, 1)
a.set_ylabel("fraction wholly in positive up-dip shear"); a.set_title("WA-NW tides (black dots = expected)", loc="left"); a.legend(frameon=False, fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(f"{FIG}/q8e_wa_nw_variants.png", dpi=130); plt.close(fig)
print(A.drop(columns="tides").T.to_string())
