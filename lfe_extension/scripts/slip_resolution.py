#!/usr/bin/env python
"""slip_resolution (2026-10-07) -- how well does the Gualandi (2025) GNSS slip model resolve slip under the
southern-Vancouver-Island LFE families?  Reads slip_build.py outputs; read-only use of cascadia-sse-strain.

    python lfe_extension/scripts/slip_resolution.py

1. GNSS station geometry around the SVI family cloud (X.mat llh + data_mask availability).
2. Model resolution of the SAME linear inversion the repo reproduces (per IC k:
   R_k = Cm0_k G^T (G Cm0_k G^T + Cd_k)^-1 G, exponential prior lambda = 21.09 km, sigma0_k from the CV curves),
   IC-weighted by each IC's share of SVI slip variance.  Spike tests (= columns of R) at every patch holding an
   SVI family; averaging-kernel rows; distinguishability of neighbouring family patches.
3. The model's own spatial structure under the families: effective rank of the family x day slip matrix,
   correlation of family slip series vs separation.
4. Temporal: noise level of daily slip vs episode slip, spectral content, rise time of major episodes.
Outputs results/slip_compare/: resolution_metrics.json, resolution_patch_table.csv, resolution_pairs.csv,
svi_gnss_stations.csv, svi_rep_families.json, fig_resolution_spikes.png
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sse_strain import load_fault_ascii, load_solution          # noqa: E402
from sse_strain.geodesy import llh2localxy                       # noqa: E402
import slip_build as SB                                          # noqa: E402

OUT = SB.OUT


def ic_weights(sol, Ls, Ld, A, meta, Wr25):
    """Share of SVI along-rake slip variance carried by each retained IC."""
    th = np.deg2rad(meta.rake_deg.to_numpy())[:, None]
    pal = (Wr25 @ Ls) * np.cos(th) + (-(Wr25 @ Ld)) * np.sin(th)       # (n_fam, K)
    v = (pal ** 2).mean(0) * A.var(0)
    return v / v.sum()


def resolution_mats(sol, G):
    """Per-IC resolution matrices and posterior/prior dip-slip variance ratios (mirrors Solution.invert)."""
    opt = sol.options["inversion"]
    lam = float(np.ravel(opt["lambda_dip"])[0]); lam0 = float(np.ravel(opt["lambda0"])[0])
    sig = np.ravel(np.asarray(opt["sigma"], float))
    mis = np.asarray(sol.misfit_comps, float)
    ind = np.nanargmin(mis, axis=0)
    xyz = sol.fault.centroid_xyz_km
    D = np.sqrt(((xyz[:, None] - xyz[None]) ** 2).sum(-1))
    n_p = sol.fault.n_patch
    E = np.exp(-D / lam)
    out = []
    for k, ic in enumerate(sol.component_index):
        Em = sig[ind[k]] * lam0 / lam
        Cm = np.zeros((2 * n_p, 2 * n_p)); Cm[:n_p, :n_p] = Em ** 2 * E; Cm[n_p:, n_p:] = Em ** 2 * E
        Cd = np.diag(sol.var_U[:, ic])
        GC = G @ Cm
        Ainv_G = np.linalg.solve(GC @ G.T + Cd, G)
        R = Cm @ G.T @ Ainv_G
        Cp_diag = np.diag(Cm) - np.einsum("ij,ji->i", GC.T, np.linalg.solve(GC @ G.T + Cd, GC))
        out.append(dict(R=R, var_ratio=Cp_diag / np.diag(Cm), sigma0=sig[ind[k]]))
    return out, D


def kernel_stats(vec, p, xy, strike, area, D3):
    """vec: recovered dip slip (n_p,) for a unit spike at p (or an averaging row). Stats in km."""
    pos = np.clip(vec, 0, None) * area
    tot = pos.sum()
    phi = np.deg2rad(strike[p])
    us, ud = np.array([np.sin(phi), np.cos(phi)]), np.array([np.cos(phi), -np.sin(phi)])   # along strike, down-dip (horiz)
    dxy = xy - xy[p]
    a_s, a_d = dxy @ us, dxy @ ud
    c_s, c_d = (pos * a_s).sum() / tot, (pos * a_d).sum() / tot
    rms_s = np.sqrt((pos * (a_s - c_s) ** 2).sum() / tot)
    rms_d = np.sqrt((pos * (a_d - c_d) ** 2).sum() / tot)
    o = np.argsort(D3[p]); cs = np.cumsum(pos[o]) / tot
    r50 = float(D3[p][o][np.searchsorted(cs, 0.5)])
    return dict(self_frac=float(pos[p] / tot), centroid_shift_strike_km=float(c_s), centroid_shift_downdip_km=float(c_d),
                rms_strike_km=float(rms_s), rms_downdip_km=float(rms_d), r50_km=r50)


def main():
    C = np.load(f"{OUT}/slip_cache.npz")
    G, Ls, Ld, A = C["G"], C["Ls"], C["Ld"], C["A"]
    dates = pd.to_datetime(C["dates"])
    fault = load_fault_ascii(SB.MESH)
    sol = load_solution(SB.SOLDIR, fault)
    sol.L = C["L"]
    meta = pd.read_csv(f"{OUT}/slip_svi_family_meta.csv")
    Wz = np.load(f"{OUT}/slip_svi_weights.npz")
    assert (Wz["fid"] == meta.fid.to_numpy()).all()
    n_p = fault.n_patch
    M = {}

    # ---------------------------------------------------------------- 1. stations
    st = C["station_llh"]; names = C["station_names"]
    lat0, lon0 = meta.lat.mean(), meta.lon.mean()
    sx, sy = llh2localxy(st[:, 1], st[:, 0], (lon0, lat0))
    fx, fy = llh2localxy(meta.lat.to_numpy(), meta.lon.to_numpy(), (lon0, lat0))
    dist = np.hypot(sx, sy)
    from sse_strain.matio import load_mat
    X = load_mat(f"{SB.SOLDIR}/X.mat")["Xd"]
    mask = np.asarray(X["data_mask"], float)[::3]           # east rows -> per station availability
    yrs = dates.year.to_numpy()
    S = pd.DataFrame(dict(station=names, lon=st[:, 0], lat=st[:, 1], dist_to_svi_centroid_km=dist.round(1),
                          avail_all=mask.mean(1).round(3),
                          avail_2007_2013=mask[:, yrs <= 2013].mean(1).round(3),
                          avail_2014_2026=mask[:, yrs >= 2014].mean(1).round(3)))
    nearest_fam = np.sqrt((sx[:, None] - fx[None]) ** 2 + (sy[:, None] - fy[None]) ** 2)
    S["min_dist_to_any_svi_family_km"] = nearest_fam.min(1).round(1)
    S = S.sort_values("dist_to_svi_centroid_km")
    S.to_csv(f"{OUT}/svi_gnss_stations.csv", index=False)
    M["stations"] = dict(svi_family_centroid=[round(lat0, 3), round(lon0, 3)], n_total=int(len(S)),
                         **{f"n_within_{r}km": int((S.dist_to_svi_centroid_km <= r).sum()) for r in (25, 50, 75, 100, 150)},
                         within_75km=S[S.dist_to_svi_centroid_km <= 75][["station", "dist_to_svi_centroid_km",
                                                                         "avail_2007_2013", "avail_2014_2026"]].to_dict("records"),
                         nearest_station_per_family_km=dict(pd.Series(nearest_fam.min(0)).describe().round(1)))

    # cloud extent along strike / down dip (local frame of the mean strike)
    phi = np.deg2rad(np.median(fault.strike[meta.patch]))
    us, ud = np.array([np.sin(phi), np.cos(phi)]), np.array([np.cos(phi), -np.sin(phi)])
    a_s, a_d = np.column_stack([fx, fy]) @ us, np.column_stack([fx, fy]) @ ud
    M["cloud"] = dict(n_families=int(len(meta)), n_distinct_containing_patches=int(meta.patch.nunique()),
                      strike_deg=round(float(np.degrees(phi)), 1),
                      along_strike_extent_km=round(float(np.ptp(a_s)), 1), downdip_extent_km=round(float(np.ptp(a_d)), 1),
                      along_strike_5_95_km=round(float(np.subtract(*np.percentile(a_s, [95, 5]))), 1),
                      downdip_5_95_km=round(float(np.subtract(*np.percentile(a_d, [95, 5]))), 1),
                      bostock_depth_range_km=[float(meta.bostock_depth_km.min()), float(meta.bostock_depth_km.max())],
                      mesh_depth_range_km=[float(meta.mesh_depth_at_epi_km.min()), float(meta.mesh_depth_at_epi_km.max())],
                      mesh_minus_bostock_depth_km=dict(meta.mesh_minus_bostock_km.describe().round(2)),
                      patch_edge_km_median=round(float(np.median(np.linalg.norm(
                          fault.vertices_m()[:, [1, 2, 0]] - fault.vertices_m(), axis=2)) / 1e3), 1),
                      patch_area_km2_median=round(float(np.median(fault.area)), 1))

    # ---------------------------------------------------------------- 2. resolution
    w = ic_weights(sol, Ls, Ld, A, meta, Wz["r25"])
    RR, D3 = resolution_mats(sol, G)
    Rbar = sum(wk * r["R"] for wk, r in zip(w, RR))
    vr = sum(wk * r["var_ratio"] for wk, r in zip(w, RR))
    M["ic_weights_svi"] = {f"IC{int(i) + 1}": round(float(x), 3) for i, x in zip(sol.component_index, w)}
    M["prior"] = dict(type="exponential", lambda_km=21.09, m0=0.0, sigma0_per_IC={f"IC{int(i) + 1}": float(r["sigma0"])
                      for i, r in zip(sol.component_index, RR)})
    xy = fault.centroid_xyz_km[:, :2]
    rows = []
    for p in sorted(meta.patch.unique()):
        col = Rbar[n_p:, n_p + p]             # dip-slip response to a unit dip-slip spike at p
        row = Rbar[n_p + p, n_p:]             # averaging kernel of the dip-slip estimate at p
        ks = kernel_stats(col, p, xy, fault.strike, fault.area, D3)
        kr = kernel_stats(row, p, xy, fault.strike, fault.area, D3)
        rows.append(dict(patch=p, n_families=int((meta.patch == p).sum()),
                         centroid_depth_km=round(float(-fault.centroid_llh[p, 2]), 1),
                         diag=round(float(Rbar[n_p + p, n_p + p]), 4),
                         spike_potency_recovered=round(float((col * fault.area).sum() / fault.area[p]), 3),
                         posterior_over_prior_var=round(float(vr[n_p + p]), 3),
                         **{f"spike_{k}": round(v, 2) for k, v in ks.items()},
                         **{f"avgker_{k}": round(v, 2) for k, v in kr.items()}))
    P = pd.DataFrame(rows)
    P.to_csv(f"{OUT}/resolution_patch_table.csv", index=False)
    # distinguishability: correlation of averaging kernels between family patches
    pats = P.patch.to_numpy()
    K = Rbar[n_p:, n_p:][pats][:, :]              # rows: averaging kernels
    Kc = np.corrcoef(K)
    pr = []
    for i in range(len(pats)):
        for j in range(i + 1, len(pats)):
            pr.append(dict(p1=pats[i], p2=pats[j], dist_km=round(float(D3[pats[i], pats[j]]), 1),
                           avgker_corr=round(float(Kc[i, j]), 3)))
    PR = pd.DataFrame(pr)
    PR.to_csv(f"{OUT}/resolution_pairs.csv", index=False)
    bins = [0, 15, 25, 35, 50, 75, 120]
    PR["dbin"] = pd.cut(PR.dist_km, bins)
    M["resolution"] = dict(
        patch_summary=dict(P[["diag", "spike_potency_recovered", "posterior_over_prior_var", "spike_self_frac",
                              "spike_r50_km", "spike_rms_strike_km", "spike_rms_downdip_km", "avgker_r50_km",
                              "avgker_rms_strike_km", "avgker_rms_downdip_km"]].describe().loc[["min", "50%", "max"]]
                          .round(3).to_dict()),
        avgker_corr_by_distance=PR.groupby("dbin", observed=True).avgker_corr.agg(["count", "median", "min", "max"])
        .round(3).reset_index().astype({"dbin": str}).to_dict("records"))

    # ---------------------------------------------------------------- 3. model's own spatial structure
    D = pd.read_parquet(f"{OUT}/slip_svi_family_daily.parquet", columns=["fid", "date", "slip_mm", "slip_mm_r25", "slip_sd_mm"])
    Y = D.pivot(index="date", columns="fid", values="slip_mm")[meta.fid]
    Yc = Y - Y.mean()
    sv = np.linalg.svd(Yc.to_numpy(), compute_uv=False)
    ev = np.cumsum(sv ** 2) / (sv ** 2).sum()
    cc = np.corrcoef(Y.to_numpy().T)
    fd = np.sqrt((fx[:, None] - fx[None]) ** 2 + (fy[:, None] - fy[None]) ** 2)
    iu = np.triu_indices(len(meta), 1)
    cd = pd.DataFrame(dict(d=fd[iu], c=cc[iu]))
    cd["dbin"] = pd.cut(cd.d, [0, 5, 10, 20, 30, 50, 80])
    # updip vs downdip episode amplitude ratio: regress each family's slip on the cloud mean
    ym = Yc.to_numpy().mean(1)
    beta = (Yc.to_numpy() * ym[:, None]).sum(0) / (ym ** 2).sum()
    M["model_structure"] = dict(
        n_components_total=int(A.shape[1]),
        family_slip_matrix_modes_for_90_99_99p9pct=[int(np.searchsorted(ev, q) + 1) for q in (0.90, 0.99, 0.999)],
        variance_first_mode=round(float(ev[0]), 4),
        corr_by_separation=cd.groupby("dbin", observed=True).c.agg(["count", "median", "min"]).round(3).reset_index()
        .astype({"dbin": str}).to_dict("records"),
        loading_on_cloud_mean_vs_bostock_depth_spearman=round(float(pd.Series(beta).corr(meta.bostock_depth_km,
                                                                                         method="spearman")), 3),
        loading_on_cloud_mean_range=[round(float(beta.min()), 2), round(float(beta.max()), 2)])

    # ---------------------------------------------------------------- 4. temporal
    t = np.asarray(Y.index)
    mid = meta.iloc[(meta.bostock_depth_km - meta.bostock_depth_km.median()).abs().argsort().iloc[0]].fid
    y = Y[mid].to_numpy()
    smooth = pd.Series(y).rolling(31, center=True, min_periods=10).median().to_numpy()
    resid = y - smooth
    mad = lambda x: 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))
    d1 = np.diff(y)
    f = np.fft.rfftfreq(len(y), 1.0)
    Pw = np.abs(np.fft.rfft(y - y.mean())) ** 2
    hp = Pw[f > 1 / 7].sum() / Pw[1:].sum()
    ac = [float(pd.Series(resid).autocorr(l)) for l in (1, 2, 3, 5, 7)]
    sd_post = D[D.fid == mid].slip_sd_mm.median()
    M["temporal"] = dict(
        family_used=int(mid), temporal_prior="none: vbICA sources are i.i.d. mixtures of 4 Gaussians per IC "
        "(init_parameters.HMM = 0), so daily V values carry no smoothing prior; the repo's slip RATE low-passes "
        "the amplitudes with a 365-tap Hamming FIR at 1/7 per day then a 7-day linear fit",
        daily_noise_mm_resid_vs_31d_median=round(float(mad(resid)), 3),
        daily_noise_mm_first_diff_over_sqrt2=round(float(mad(d1) / np.sqrt(2)), 3),
        posterior_sd_mm_from_var_V_median=round(float(sd_post), 4),
        resid_autocorr_lags_1_2_3_5_7=[round(a, 3) for a in ac],
        power_fraction_periods_shorter_than_7d=round(float(hp), 3))
    # rise time for major SVI episodes at that family (on a 7-day running-mean series)
    E = pd.read_csv(SB.ROOT + "/lfe_extension_wa/results/episode_rates/episodes_svi_final.csv")
    E = E[E.major]
    ys = pd.Series(y, index=pd.DatetimeIndex(t)).rolling(7, center=True, min_periods=4).mean()
    rr = []
    for e in E.itertuples():
        a, b = pd.Timestamp(e.start), pd.Timestamp(e.end)
        seg = ys[a - pd.Timedelta(days=20): b + pd.Timedelta(days=20)]
        if len(seg) < 20:
            continue
        pre, post = ys[a - pd.Timedelta(days=12): a - pd.Timedelta(days=3)].mean(), ys[b + pd.Timedelta(days=3): b + pd.Timedelta(days=12)].mean()
        dlt = post - pre
        ns = mad(resid) * np.sqrt(2 / 9)
        r = dict(episode=e.episode, lfe_span_days=int(e.span_days), slip_increment_mm=round(float(dlt), 2),
                 snr=round(float(dlt / ns), 1))
        if dlt > 10 * ns:                       # only well-measured episodes; 10/90 % crossings are noise-prone
            fr = ((seg - pre) / dlt).to_numpy()
            for lo, hi in ((0.1, 0.9), (0.25, 0.75)):
                t_lo = seg.index[np.argmax(fr >= lo)]; t_hi = seg.index[np.argmax(fr >= hi)]
                r[f"rise_{int(lo * 100)}_{int(hi * 100)}_days"] = int((t_hi - t_lo).days)
        rr.append(r)
    RT = pd.DataFrame(rr)
    ok = RT.dropna(subset=["rise_25_75_days"]) if "rise_25_75_days" in RT else RT.iloc[:0]
    M["temporal"]["major_episodes"] = RT.to_dict("records")
    M["temporal"]["snr_ge_10_episodes"] = int(len(ok))
    M["temporal"]["rise_10_90_days_median_snr_ge_10"] = float(ok.rise_10_90_days.median())
    M["temporal"]["rise_25_75_days_median_snr_ge_10"] = float(ok.rise_25_75_days.median())
    M["temporal"]["lfe_span_days_median_snr_ge_10"] = float(ok.lfe_span_days.median())
    M["temporal"]["negative_increment_episodes"] = dict(n=int((RT.slip_increment_mm < 0).sum()),
                                                        most_negative_mm=float(RT.slip_increment_mm.min()))

    # ---------------------------------------------------------------- representative families (fig use)
    cnt = pd.read_parquet(SB.COUNTS_SVI)
    cnt = cnt[cnt.run == "prod"]
    net = cnt.pivot_table(index="fid", columns="series", values="n075", aggfunc="sum").fillna(0)
    meta["net075"] = (net["data"] - net.get("c1", 0)).reindex(meta.fid).fillna(0).to_numpy()
    q = meta.bostock_depth_km.quantile([1 / 3, 2 / 3]).to_numpy()
    picks = {}
    for name, m in (("updip", meta.bostock_depth_km <= q[0]),
                    ("middle", (meta.bostock_depth_km > q[0]) & (meta.bostock_depth_km <= q[1])),
                    ("downdip", meta.bostock_depth_km > q[1])):
        sub = meta[m].sort_values("net075", ascending=False)
        # choose the busiest family whose patch differs from the ones already picked
        for r in sub.itertuples():
            if r.patch not in [v["patch"] for v in picks.values()]:
                picks[name] = dict(fid=int(r.fid), patch=int(r.patch), bostock_depth_km=float(r.bostock_depth_km),
                                   mesh_depth_km=float(r.mesh_depth_at_epi_km), lat=float(r.lat), lon=float(r.lon),
                                   net075_total=float(r.net075))
                break
    json.dump(picks, open(f"{OUT}/svi_rep_families.json", "w"), indent=1)
    ip = {v["patch"]: k for k, v in picks.items()}
    pp = [picks[k]["patch"] for k in ("updip", "middle", "downdip")]
    M["rep_families"] = picks
    M["rep_pairs_avgker_corr"] = {f"{ip[a]}-{ip[b]}": dict(dist_km=round(float(D3[a, b]), 1),
                                  avgker_corr=round(float(np.corrcoef(Rbar[n_p + a, n_p:], Rbar[n_p + b, n_p:])[0, 1]), 3),
                                  slip_series_corr=round(float(np.corrcoef(
                                      Y[picks[ip[a]]["fid"]], Y[picks[ip[b]]["fid"]])[0, 1]), 3))
                                  for a, b in ((pp[0], pp[1]), (pp[1], pp[2]), (pp[0], pp[2]))}

    # ---------------------------------------------------------------- figure: spike responses
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    verts = np.column_stack([fault.lon.ravel(), fault.lat.ravel()])
    uniq, inv = np.unique(np.round(verts, 6), axis=0, return_inverse=True)
    trg = mtri.Triangulation(uniq[:, 0], uniq[:, 1], inv.reshape(n_p, 3))
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.6), sharey=True)
    for ax, (name, pk) in zip(axs, picks.items()):
        p = pk["patch"]
        col = Rbar[n_p:, n_p + p]
        tp = ax.tripcolor(trg, facecolors=col, cmap="RdBu_r", vmin=-0.07, vmax=0.07)
        ax.triplot(trg, lw=0.2, color="0.6")
        ax.plot(meta.lon, meta.lat, ".", ms=3, color="k")
        ax.plot(fault.centroid_llh[p, 0], fault.centroid_llh[p, 1], "*", ms=14, mfc="gold", mec="k")
        ax.plot(st[:, 0], st[:, 1], "^", ms=6, mfc="none", mec="g")
        ax.set_xlim(-125.6, -121.8); ax.set_ylim(47.4, 49.6)
        ax.set_aspect(1 / np.cos(np.deg2rad(48.5)))
        rP = P.set_index("patch").loc[p]
        ax.set_title(f"{name}: unit thrust spike,\npatch {p} ({rP.centroid_depth_km:.0f} km)\n"
                     f"recovered at spike {rP.diag:.2f};\nhalf-potency radius {rP.spike_r50_km:.0f} km", fontsize=10)
        ax.set_xlabel("longitude")
    axs[0].set_ylabel("latitude")
    fig.colorbar(tp, ax=axs, shrink=0.8, label="recovered thrust slip per unit input (IC-weighted R)")
    fig.suptitle("GNSS resolution under SVI LFE families (Gualandi 2025 inversion; black dots = families, "
                 "green triangles = GNSS)", fontsize=11)
    fig.savefig(f"{OUT}/fig_resolution_spikes.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    json.dump(M, open(f"{OUT}/resolution_metrics.json", "w"), indent=1, default=str)
    print(json.dumps(M, indent=1, default=str)[:12000])


if __name__ == "__main__":
    main()
