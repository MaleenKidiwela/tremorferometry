#!/usr/bin/env python
"""slip_figures (2026-10-07) -- figures + per-episode slip table for the SVI slip-vs-LFE comparison.

    python lfe_extension/scripts/slip_figures.py

Reads slip_build.py / slip_resolution.py outputs.  Dates are UTC calendar days (GNSS daily solutions and the
LFE day bins are both UTC days); no times of day are shown.
Writes results/slip_compare/:
  slip_svi_family_episode.csv      every SVI episode (episodes_svi_final.csv) x family: slip increment for
                                   patch / r10 / r25, its 1-sigma from the family's daily noise, net LFE count
  fig_slip_vs_lfe_rep_families.png slip (patch, r10, r25) + daily LFE counts P>=0.75 (data minus scrambled control)
  fig_svi_families_mesh_map.png    families on the mesh with GNSS stations
  fig_episode_slip_heatmap.png     slip increment per major episode at each family (+ vs LFE counts)
  fig_sanity_fig07_strip.png       repo fig07 (signed magnitude) vs this reconstruction in the 47.5-49.5 N strip
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sse_strain import load_fault_ascii, load_solution          # noqa: E402
import slip_build as SB                                          # noqa: E402

OUT = SB.OUT
COAST = f"{SB.REPO}/data/cascadia_coastline_ne10m.json"
EPIS = f"{SB.ROOT}/lfe_extension_wa/results/episode_rates/episodes_svi_final.csv"
VAR = (("patch", "slip_mm", "containing patch"), ("r10", "slip_mm_r10", "Gaussian 10 km"),
       ("r25", "slip_mm_r25", "Gaussian 25 km"))
COL = {"patch": "#1f4e9c", "r10": "#d9822b", "r25": "#3a9a5b"}


def mad(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def lfe_net():
    C = pd.read_parquet(SB.COUNTS_SVI)
    C = C[C.run == "prod"]
    N = C.pivot_table(index=["fid", "date"], columns="series", values="n075", aggfunc="sum").fillna(0)
    N["net075"] = N["data"] - N["c1"]
    N = N.reset_index()
    N["date"] = pd.to_datetime(N.date)
    return N[["fid", "date", "data", "c1", "net075"]]


def coast(ax):
    for seg in json.load(open(COAST)):
        a = np.asarray(seg)
        ax.plot(a[:, 0], a[:, 1], lw=0.6, color="0.35")


def main():
    meta = pd.read_csv(f"{OUT}/slip_svi_family_meta.csv")
    picks = json.load(open(f"{OUT}/svi_rep_families.json"))
    D = pd.read_parquet(f"{OUT}/slip_svi_family_daily.parquet")
    N = lfe_net()
    cache = np.load(f"{OUT}/slip_cache.npz")
    fault = load_fault_ascii(SB.MESH)

    # ------------------------------------------------- per-family daily noise and episode increments
    piv = {k: D.pivot(index="date", columns="fid", values=c) for k, c, _ in VAR}
    noise = {k: pd.Series({f: mad((piv[k][f] - piv[k][f].rolling(31, center=True, min_periods=10).median()).to_numpy())
                           for f in meta.fid}) for k, _, _ in VAR}
    E = pd.read_csv(EPIS)
    Nn = N.set_index(["fid", "date"]).net075
    rows = []
    for e in E.itertuples():
        a, b = pd.Timestamp(e.start), pd.Timestamp(e.end)
        pre = (piv["patch"].index >= a - pd.Timedelta(days=5)) & (piv["patch"].index < a)
        post = (piv["patch"].index > b) & (piv["patch"].index <= b + pd.Timedelta(days=5))
        lf = N[(N.date >= a) & (N.date <= b)].groupby("fid").net075.sum()
        for f in meta.fid:
            r = dict(episode=e.episode, start=e.start, end=e.end, major=bool(e.major), fid=int(f),
                     lfe_net075=float(lf.get(f, 0.0)))
            for k, _, _ in VAR:
                y = piv[k][f].to_numpy()
                r[f"dslip_mm_{k}"] = round(float(y[post].mean() - y[pre].mean()), 3)
                r[f"dslip_sd_mm_{k}"] = round(float(noise[k][f] * np.sqrt(2 / 5)), 3)
            rows.append(r)
    EP = pd.DataFrame(rows)
    EP.to_csv(f"{OUT}/slip_svi_family_episode.csv", index=False)

    # ------------------------------------------------- fig: slip vs LFE at 3 representative families
    zoom = (pd.Timestamp("2024-06-15"), pd.Timestamp("2024-11-15"))
    fig, axs = plt.subplots(3, 2, figsize=(16, 10.5), gridspec_kw=dict(width_ratios=[2.6, 1]))
    for i, name in enumerate(("updip", "middle", "downdip")):
        pk = picks[name]; f = pk["fid"]
        d = D[D.fid == f].set_index("date")
        n = N[N.fid == f].set_index("date").net075
        for j, (t0, t1) in enumerate(((pd.Timestamp("2014-01-01"), d.index.max()), zoom)):
            ax = axs[i, j]
            ax2 = ax.twinx()
            nn = n[(n.index >= t0) & (n.index <= t1)]
            if j == 0:                                   # 12 years: 7-day sums so the bursts are visible
                nn = nn.reindex(pd.date_range(t0, t1), fill_value=0).resample("7D").sum()
                ax2.fill_between(nn.index, 0, nn.clip(lower=0).to_numpy(), step="post", color="0.55", alpha=0.55, lw=0)
                lab2 = "LFEs/week, P≥0.75\n(data − scrambled control)"
            else:
                ax2.bar(nn.index, nn.clip(lower=0).to_numpy(), width=1.0, color="0.55", alpha=0.55, lw=0)
                lab2 = "LFEs/day, P≥0.75\n(data − scrambled control)"
            ax2.set_ylim(0, max(5, nn.max()) * 2.6)
            ax2.set_ylabel(lab2, fontsize=8, color="0.4")
            dd = d[(d.index >= t0) & (d.index <= t1)]
            for k, c, lab in VAR:
                y = dd[c] - dd[c].iloc[:30].mean()
                ax.plot(dd.index, y, lw=0.7 if j == 0 else 1.1, color=COL[k], label=lab)
            if j == 1:
                y = dd["slip_mm"] - dd["slip_mm"].iloc[:30].mean()
                s = dd["slip_sd_mm"]
                ax.fill_between(dd.index, y - s, y + s, color=COL["patch"], alpha=0.15, lw=0,
                                label="±1σ (vbICA var_V only)")
            ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
            ax.set_xlim(t0, t1)
            ax.set_ylabel("slip along rake (mm)\nrel. to secular trend")
            if j == 0:
                ax.set_title(f"{name}: family {f}  (Bostock depth {pk['bostock_depth_km']:.1f} km, mesh "
                             f"{pk['mesh_depth_km']:.1f} km, {pk['lat']:.2f}N {-pk['lon']:.2f}W)", fontsize=10, loc="left")
            else:
                ax.set_title("zoom: 2024-08 episode", fontsize=10, loc="left")
            if i == 0 and j == 1:
                ax.legend(fontsize=7, loc="upper left")
    fig.suptitle("Modelled slow slip (Gualandi 2025 vbICA, solution 2026-08-18) at three SVI LFE families vs their "
                 "daily LFE counts; slip offset to the mean of the first 30 days shown", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_slip_vs_lfe_rep_families.png", dpi=130)
    plt.close(fig)

    # ------------------------------------------------- fig: families on the mesh with GNSS
    verts = np.column_stack([fault.lon.ravel(), fault.lat.ravel()])
    uniq, inv = np.unique(np.round(verts, 6), axis=0, return_inverse=True)
    trg = mtri.Triangulation(uniq[:, 0], uniq[:, 1], inv.reshape(fault.n_patch, 3))
    st = cache["station_llh"]; names = [s[:-1] for s in cache["station_names"]]
    fig, ax = plt.subplots(figsize=(11, 9))
    tp = ax.tripcolor(trg, facecolors=-fault.centroid_llh[:, 2], cmap="Greys", vmin=10, vmax=60, alpha=0.55)
    ax.triplot(trg, lw=0.35, color="0.45")
    used = meta.patch.unique()
    for p in used:
        ax.fill(fault.lon[p], fault.lat[p], fc="none", ec="#7a3fb0", lw=1.2)
    coast(ax)
    sc = ax.scatter(meta.lon, meta.lat, c=meta.bostock_depth_km, cmap="viridis_r", s=22, ec="k", lw=0.3, zorder=5)
    for name, pk in picks.items():
        ax.plot(pk["lon"], pk["lat"], "*", ms=17, mfc="gold", mec="k", zorder=6)
        ax.annotate(f"{name} (fid {pk['fid']})", (pk["lon"], pk["lat"]), xytext=(6, 6), textcoords="offset points",
                    fontsize=9, weight="bold")
    ax.plot(st[:, 0], st[:, 1], "^", ms=8, mfc="#2ca02c", mec="k", zorder=6)
    lo, la = (-125.3, -122.3), (47.6, 49.4)
    for (x, y), nm in zip(st[:, :2], names):
        if lo[0] < x < lo[1] and la[0] < y < la[1]:
            ax.annotate(nm, (x, y), xytext=(3, -9), textcoords="offset points", fontsize=7, color="#1d6b1d")
    ax.set_xlim(*lo); ax.set_ylim(*la)
    ax.set_aspect(1 / np.cos(np.deg2rad(48.5)))
    fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.01, label="family depth, Bostock (km)")
    fig.colorbar(tp, ax=ax, shrink=0.6, pad=0.06, label="mesh patch centroid depth (km)")
    ax.set_title("SVI LFE families (115) on the Slab2 fallback mesh used by the repo; purple = patches containing a "
                 "family;\ngreen triangles = GNSS stations in the Gualandi solution", fontsize=10)
    ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    fig.savefig(f"{OUT}/fig_svi_families_mesh_map.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------------------------- fig: per-episode slip at each family
    M = EP[EP.major].copy()
    order = meta.sort_values("bostock_depth_km").fid.to_numpy()
    eps = M.episode.drop_duplicates().tolist()
    fig, axs = plt.subplots(1, 3, figsize=(19, 9), gridspec_kw=dict(width_ratios=[1, 1, 0.9]))
    for ax, k in zip(axs[:2], ("patch", "r25")):
        H = M.pivot(index="fid", columns="episode", values=f"dslip_mm_{k}").reindex(index=order, columns=eps)
        im = ax.imshow(H.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-15, vmax=15, interpolation="nearest")
        ax.set_xticks(range(len(eps))); ax.set_xticklabels(eps, rotation=90, fontsize=7)
        yt = np.arange(0, len(order), 10)
        ax.set_yticks(yt); ax.set_yticklabels([f"{meta.set_index('fid').bostock_depth_km[order[i]]:.0f}" for i in yt],
                                              fontsize=7)
        ax.set_ylabel("families sorted by Bostock depth (km) →")
        ax.set_title(f"thrust slip increment per major SVI episode ({dict((a, c) for a, _, c in VAR)[k]}), mm",
                     fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.7, pad=0.01)
    ax = axs[2]
    x = M.lfe_net075.to_numpy(); y = M.dslip_mm_r10.to_numpy()
    ax.scatter(x, y, s=5, alpha=0.35, color="#1f4e9c")
    ax.axhline(0, color="k", lw=0.5)
    sd = float(M.dslip_sd_mm_r10.median())
    ax.axhspan(-3 * sd, 3 * sd, color="0.8", alpha=0.5, label=f"±3σ daily-noise band ({3 * sd:.1f} mm)")
    rho = pd.Series(x).corr(pd.Series(y), method="spearman")
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel("net LFEs in episode (P≥0.75, data − scrambled control)")
    ax.set_ylabel("thrust slip increment (mm), Gaussian 10 km")
    ax.set_title(f"family x major-episode pairs (n={len(M)}), Spearman {rho:.2f}", fontsize=9)
    ax.legend(fontsize=8)
    fig.suptitle("Modelled slip per major SVI episode at each LFE family (increment = mean of 5 days after end − "
                 "5 days before start)", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_episode_slip_heatmap.png", dpi=130)
    plt.close(fig)

    # ------------------------------------------------- sanity: repo fig07 strip vs this reconstruction
    sol = load_solution(SB.SOLDIR, fault)
    sol.L = cache["L"]
    mag = sol.slip_magnitude(apply_rake_flip=True)          # repo / Gualandi signed magnitude (fig07 bottom)
    Ls, Ld, A = cache["Ls"], cache["Ld"], cache["A"]
    thr = -(Ld @ A.T)
    lat = fault.centroid_llh[:, 1]; dep = -fault.centroid_llh[:, 2]
    band = (dep >= 25) & (dep <= 50) & (lat >= 47.5) & (lat <= 49.5)
    w = fault.area[band] / fault.area[band].sum()
    t_repo = sol.dates; t_true = pd.to_datetime(cache["dates"])
    r_mag = w @ mag[band]; r_thr = w @ thr[band]
    rep = dict(strip="47.5-49.5N, 25-50 km depth, area-weighted mean",
               corr_repo_signed_magnitude_vs_thrust=round(float(np.corrcoef(r_mag, r_thr)[0, 1]), 4))
    fig, ax = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    ax[0].plot(t_repo, 1e3 * r_mag, lw=0.6, color="0.3", label="repo slip_magnitude(apply_rake_flip=True) (fig07 quantity)")
    ax[0].plot(t_true, -1e3 * r_thr, lw=0.6, color="#d62728", alpha=0.7, label="this work: −(thrust dip slip)")
    ax[0].set_ylabel("mm"); ax[0].legend(fontsize=8)
    ax[0].set_title(f"Strip mean; corr(repo signed magnitude, this-work thrust) = "
                    f"{rep['corr_repo_signed_magnitude_vs_thrust']:+.3f}  → repo sign makes slow slip NEGATIVE on this mesh",
                    fontsize=9)
    ax[1].plot(t_true, 1e3 * r_thr, lw=0.6, color="#1f4e9c")
    ax[1].set_ylabel("thrust slip (mm)"); ax[1].set_title("this work: thrust-positive dip slip, same strip", fontsize=9)
    for a in ax:
        a.set_xlim(pd.Timestamp("2009-01-01"), pd.Timestamp("2018-01-01"))
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_sanity_fig07_strip.png", dpi=130)
    plt.close(fig)
    S = json.load(open(f"{OUT}/sanity_repo_reproduction.json"))
    S["fig07_strip_check"] = rep
    json.dump(S, open(f"{OUT}/sanity_repo_reproduction.json", "w"), indent=1, default=str)
    print(json.dumps(rep))
    # summary numbers for RESOLUTION.md
    summ = dict(
        major_pairs=int(len(M)),
        spearman_dslip_r10_vs_lfe=round(float(rho), 3),
        median_dslip_sd_mm=dict(patch=float(M.dslip_sd_mm_patch.median()), r10=float(M.dslip_sd_mm_r10.median()),
                                r25=float(M.dslip_sd_mm_r25.median())),
        frac_major_pairs_above_3sd_r10=round(float((M.dslip_mm_r10 > 3 * M.dslip_sd_mm_r10).mean()), 3),
        frac_major_pairs_below_minus_3sd_r10=round(float((M.dslip_mm_r10 < -3 * M.dslip_sd_mm_r10).mean()), 3),
        per_episode_spread_across_families=M.groupby("episode").dslip_mm_r10.agg(["median", "min", "max"]).round(2)
        .reset_index().to_dict("records"))
    json.dump(summ, open(f"{OUT}/episode_slip_summary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != "per_episode_spread_across_families"}))


if __name__ == "__main__":
    main()
