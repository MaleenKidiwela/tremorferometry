"""S3 figures: per-episode along-strike vs time (LFEs coloured by Slab2 depth, PNSN tremor grey), per-region summary
(map of start points + migration arrows; start->end vs date), patch-strength heatmap."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s3"
FIG = f"{ROOT}/svi_science/figures"
TZ = "America/Los_Angeles"
C_NW, C_SE, C_NONE = "#2a78d6", "#eb6834", "#8a8984"     # categorical slots 1/2 + neutral
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
DEPTH_CMAP = LinearSegmentedColormap.from_list("depth", ["#a9cdf3", "#2a78d6", "#0d2c5a"])
DIV = LinearSegmentedColormap.from_list("div", ["#1f5fae", "#8fb6e6", "#ededeb", "#f2a07c", "#c4471a"])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 9.5, "figure.facecolor": "white"})


def headline(g):
    d = g[g.series == "data"]
    sel = d.selected_v2 if "selected_v2" in d else True
    return d[(d.P_same >= 0.75) & sel & d.counted_once]


def to_pt(t):
    return pd.to_datetime(t).dt.tz_localize("UTC").dt.tz_convert(TZ).dt.tz_localize(None)


def frames():
    info = json.load(open(f"{OUT}/build_info.json"))
    return {"SVI": info["svi_frame"], "NW": info["nw_frame"]}


def inv(fr, a, d):
    s = np.radians(fr["strike_az"])
    x = a * np.sin(s) + d * np.cos(s); y = a * np.cos(s) - d * np.sin(s)
    lat = fr["origin_lat"] + y / 111.195; lon = fr["origin_lon"] + x / (111.195 * np.cos(np.radians(fr["origin_lat"])))
    return lon, lat


def dirc(v, migrates, coherent=True):
    if not migrates or not coherent:
        return C_NONE
    return C_NW if v > 0 else C_SE


def episode_panel(ax, reg, r, h, tr, bl, vmin, vmax, ylim):
    w0 = pd.Timestamp(r.window_start_PT)
    th = (to_pt(tr.time) - w0) / pd.Timedelta(days=1) if len(tr) else []
    if len(tr):
        ax.scatter(th, tr.as_km, s=3, c="#c9c8c3", lw=0, zorder=1, rasterized=True, label="PNSN tremor")
    tl = (to_pt(h.time) - w0) / pd.Timedelta(days=1)
    jit = np.random.default_rng(0).normal(0, 0.6, len(h))
    sc = ax.scatter(tl, h.as_km + jit, s=5, c=h.slab2_depth_km, cmap=DEPTH_CMAP, vmin=vmin, vmax=vmax, lw=0, zorder=2, rasterized=True)
    # Theil-Sen line through the episode's median point
    tm = np.median(tl); am = np.median(h.as_km)
    xx = np.array([tl.min(), tl.max()])
    col = dirc(r.lfe_speed, r.migrates, r.coherent)
    ax.plot(xx, am + r.lfe_speed * (xx - tm), color=col, lw=2, zorder=3)
    t1 = r.lfe_t_first - (pd.Timestamp(r.window_start_PT) - pd.Timestamp(r.window_start_PT)).days
    ax.scatter([np.quantile(tl, 0.05)], [r.lfe_start10], marker="*", s=160, c=col, edgecolor="white", lw=1, zorder=4)
    for _, b in bl.iterrows():
        ax.scatter([b.return_day_idx + 0.5], [b.bin_center_km], marker="o", s=70, facecolor="none", edgecolor=INK, lw=1.2, zorder=5)
    ax.set_ylim(*ylim); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_xlabel(f"days after {w0:%Y-%m-%d} (PT)")
    sig = "" if (r.migrates and r.coherent) else (" (not significant)" if not r.migrates else " (no single front)")
    ax.set_title(f"{reg} {r.episode}: {r.n_headline} LFEs, {abs(r.lfe_speed):.1f} km/day "
                 f"{'toward NW' if r.lfe_speed > 0 else 'toward SE'}{sig}", loc="left", color=INK)
    return sc


def per_episode(reg, K, E, T, B, frs):
    k = K[K.region == reg].reset_index(drop=True)
    H = headline(E[E.ep.isin(k.episode)])
    vmin, vmax = np.nanquantile(H.slab2_depth_km, [0.02, 0.98])
    ylim = (np.floor(H.as_km.quantile(0.0) / 5) * 5 - 5, np.ceil(H.as_km.quantile(1.0) / 5) * 5 + 5)
    n = len(k); nc = 2; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(12, 2.6 * nr), constrained_layout=True)
    axs = axs.ravel()
    for i, r in k.iterrows():
        h = headline(E[E.ep == r.episode]); tr = T[T.ep == r.episode]; bl = B[(B.region == reg) & (B.episode == r.episode) & (B.series == "LFE")]
        sc = episode_panel(axs[i], reg, r, h, tr, bl, vmin, vmax, ylim)
        axs[i].set_ylabel("along-strike km\n(SE −  /  + NW)")
        # single-episode file
        f1, a1 = plt.subplots(figsize=(8, 4), constrained_layout=True)
        s1 = episode_panel(a1, reg, r, h, tr, bl, vmin, vmax, ylim); a1.set_ylabel("along-strike km (SE −  /  + NW)")
        f1.colorbar(s1, ax=a1, label="plate-interface depth at the source (km, Slab2)", shrink=0.8)
        a1.text(0.99, 0.02, "grey = PNSN tremor · star = start (first 10 %) · line = migration fit · open circle = activity returns after ≥ 3 quiet days",
                transform=a1.transAxes, ha="right", va="bottom", fontsize=7, color=INK2)
        f1.savefig(f"{FIG}/s3_ep_{reg}_{r.episode}.png", dpi=130); plt.close(f1)
    for j in range(n, len(axs)):
        axs[j].axis("off")
    fig.colorbar(sc, ax=axs[:n], label="plate-interface depth at the source (km, Slab2)", shrink=0.4, location="right")
    fig.suptitle(f"{reg}: every big ETS episode, LFE position along strike vs time\n(grey dots = PNSN tremor; star = start; "
                 f"line = migration fit: blue = toward NW, orange = toward SE, grey = no clear migration; open circle = activity returns after ≥ 3 quiet days)",
                 fontsize=9, color=INK)
    fig.savefig(f"{FIG}/s3_{reg.lower()}_episodes.png", dpi=130); plt.close(fig)


def summary(reg, K, E, frs):
    k = K[K.region == reg].reset_index(drop=True)
    fr = frs[reg]
    fig = plt.figure(figsize=(13, 5.6), constrained_layout=True)
    ax = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
    if reg == "SVI":
        U = pd.read_csv(f"{OUT}/families_svi_frame.csv"); ulon, ulat = U.lon, U.lat
        ext = [-124.6, -122.7, 47.9, 49.05]
    else:
        U = pd.read_csv(f"{OUT}/sources_nw.csv"); U = U[U.quality.isin(["A", "B"])]; ulon, ulat = U.lon, U.lat
        ext = [-124.2, -122.0, 47.2, 48.75]
    ax.set_extent(ext, ccrs.PlateCarree())
    ax.add_feature(cfeature.NaturalEarthFeature("physical", "coastline", "10m", edgecolor=INK2, facecolor="none", lw=0.6))
    ax.scatter(ulon, ulat, s=10, c="#c9c8c3", lw=0, transform=ccrs.PlateCarree(), zorder=2)
    gl = ax.gridlines(draw_labels=True, lw=0.3, color=GRID); gl.top_labels = gl.right_labels = False; gl.left_labels = True
    for _, r in k.iterrows():
        lo, la = inv(fr, r.lfe_start10, r.lfe_start10_dd)
        vx, vy = r.lfe_speed, (r.lfe_dd_speed if r.lfe_dd_p < 0.01 else 0.0)
        L = np.hypot(vx, vy); scale = 3.0     # km of arrow per km/day
        lo2, la2 = inv(fr, r.lfe_start10 + vx * scale, r.lfe_start10_dd + vy * scale)
        col = dirc(r.lfe_speed, r.migrates, r.coherent)
        ax.annotate("", xy=(lo2, la2), xytext=(lo, la), arrowprops=dict(arrowstyle="-|>", color=col, lw=2, mutation_scale=14),
                    transform=ccrs.PlateCarree(), zorder=4)
        ax.scatter([lo], [la], s=45, c=col, edgecolor="white", lw=1, zorder=5, transform=ccrs.PlateCarree())
        ax.text(lo + 0.02, la + 0.012, r.episode, fontsize=7, color=INK2, transform=ccrs.PlateCarree(), zorder=6)
    ax.set_title(f"{reg}: where each big episode started (dot) and which way it moved (arrow; length = 3 days of travel)", loc="left")
    # right: start -> end vs date
    a2 = fig.add_subplot(1, 2, 2)
    dates = pd.to_datetime(k.window_start_PT)
    for i, r in k.iterrows():
        col = dirc(r.lfe_speed, r.migrates, r.coherent)
        a2.annotate("", xy=(dates[i], r.lfe_end10), xytext=(dates[i], r.lfe_start10),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2, mutation_scale=14))
        a2.scatter([dates[i]], [r.lfe_start10], s=40, c=col, edgecolor="white", zorder=3)
        if np.isfinite(r.get("trem_start10", np.nan)):
            a2.scatter([dates[i] + pd.Timedelta(days=40)], [r.trem_start10], marker="x", s=30, c=INK2, zorder=3)
        a2.text(dates[i], max(r.lfe_start10, r.lfe_end10) + 3, f"{abs(r.lfe_speed):.0f}", fontsize=7, ha="center", color=INK2)
    a2.axhline(0, color=GRID, lw=1)
    a2.set_ylabel("along-strike km  (SE −  /  + NW)")
    a2.set_title("start (dot) → end (arrow head) of each big episode; × = PNSN tremor start; number = km/day", loc="left")
    a2.grid(color=GRID, lw=0.6); a2.set_axisbelow(True)
    from matplotlib.lines import Line2D
    a2.legend(handles=[Line2D([], [], color=C_NW, lw=2, label="moved toward NW"), Line2D([], [], color=C_SE, lw=2, label="moved toward SE"),
                       Line2D([], [], color=C_NONE, lw=2, label="no clear / no single front"),
                       Line2D([], [], marker="x", ls="", color=INK2, label="tremor start")], loc="lower left", fontsize=8, frameon=False)
    fig.savefig(f"{FIG}/s3_{reg.lower()}_summary.png", dpi=130); plt.close(fig)


def heatmap(K):
    P = pd.read_csv(f"{OUT}/patch_strength.csv")
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.2), constrained_layout=True, gridspec_kw=dict(width_ratios=[10, 14]))
    for ax, reg in zip(axs, ["SVI", "NW"]):
        p = P[P.region == reg].copy()
        k = K[K.region == reg].set_index("episode")
        eps = k.index.tolist()
        if reg == "SVI":
            G = pd.read_csv(f"{ROOT}/svi_science/results/s0/groups.csv").set_index("group")
            order = G.sort_values("along_strike_km").index.tolist()
            labels = [f"group {g} ({G.loc[g, 'along_strike_km']:+.0f} km, {G.loc[g, 'depth_km']:.0f} km deep)" for g in order]
            p["patch"] = p.patch.astype(int)
        else:
            order = sorted(p.patch.unique(), key=lambda s: float(s.split("..")[0])); labels = [f"{o} km" for o in order]
        M = p.pivot(index="patch", columns="episode", values="ratio_to_longterm").loc[order, eps]
        N = p.pivot(index="patch", columns="episode", values="n").loc[order, eps]
        L = np.log2(M.clip(lower=1 / 8).fillna(1 / 8))
        im = ax.imshow(L[::-1], cmap=DIV, norm=TwoSlopeNorm(0, -2.5, 2.5), aspect="auto")
        for i in range(L.shape[0]):
            for j in range(L.shape[1]):
                ax.text(j, L.shape[0] - 1 - i, int(N.iloc[i, j]), ha="center", va="center", fontsize=6.5,
                        color="white" if abs(L.iloc[i, j]) > 1.6 else INK)
        ax.set_yticks(range(len(order))); ax.set_yticklabels(labels[::-1], fontsize=7.5)
        ax.set_xticks(range(len(eps)))
        ax.set_xticklabels([f"{e}\n{'·' if not (k.loc[e, 'migrates'] and k.loc[e, 'coherent']) else ('↑NW' if k.loc[e, 'lfe_speed'] > 0 else '↓SE')}" for e in eps], fontsize=7, rotation=90)
        ax.set_title(f"{reg}: patch strength per big episode (top = NW end; number = LFEs)", loc="left")
        ax.spines[:].set_visible(False)
    fig.colorbar(im, ax=axs, shrink=0.7, label="log2(share in episode / long-term share)\n+1 = twice its usual share, −1 = half", ticks=[-2, -1, 0, 1, 2])
    fig.savefig(f"{FIG}/s3_patch_heatmap.png", dpi=130); plt.close(fig)


def main():
    K = pd.read_csv(f"{OUT}/episodes_kinematics.csv")
    B = pd.read_csv(f"{OUT}/bursts.csv")
    frs = frames()
    for reg in ["SVI", "NW"]:
        E = pd.read_parquet(f"{OUT}/events_{reg.lower()}.parquet")
        if reg == "SVI":
            F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")[["fid", "slab2_depth_km"]].rename(columns={"fid": "unit"})
        else:
            F = pd.read_csv(f"{OUT}/sources_nw.csv")[["source", "slab2_depth_km"]].rename(columns={"source": "unit"})
        E = E.merge(F, on="unit", how="left")
        if "selected_v2" not in E:
            pass
        T = pd.read_parquet(f"{OUT}/tremor_{reg.lower()}.parquet")
        per_episode(reg, K, E, T, B, frs)
        summary(reg, K, E, frs)
    heatmap(K)


if __name__ == "__main__":
    main()


def cross_region():
    """Overlapping SVI / NW big episodes on one along-margin axis (azimuth 330 deg from 48.3N 123.4W)."""
    K = pd.read_csv(f"{OUT}/episodes_kinematics.csv")
    ks, kn = K[K.region == "SVI"], K[K.region == "NW"]
    pairs = []
    for _, a in ks.iterrows():
        for _, b in kn.iterrows():
            a0, a1 = pd.Timestamp(a.window_start_PT), pd.Timestamp(a.window_end_PT)
            b0, b1 = pd.Timestamp(b.window_start_PT), pd.Timestamp(b.window_end_PT)
            if min(a1, b1) > max(a0, b0):
                pairs.append((a.episode, b.episode, min(a0, b0), max(a1, b1)))
    if not pairs:
        return
    lat0, lon0, az = 48.3, -123.4, np.radians(330)

    def along(lat, lon):
        x = (lon - lon0) * 111.195 * np.cos(np.radians(lat0)); y = (lat - lat0) * 111.195
        return x * np.sin(az) + y * np.cos(az)
    S = pd.read_parquet(f"{OUT}/events_svi.parquet"); N = pd.read_parquet(f"{OUT}/events_nw.parquet")
    TS = pd.read_parquet(f"{OUT}/tremor_svi.parquet"); TN = pd.read_parquet(f"{OUT}/tremor_nw.parquet")
    fig, axs = plt.subplots(len(pairs), 1, figsize=(11, 2.9 * len(pairs)), constrained_layout=True)
    axs = np.atleast_1d(axs)
    rows = []
    for ax, (es, en, t0, t1) in zip(axs, pairs):
        hs = headline(S[S.ep == es]); hn = headline(N[N.ep == en])
        tr = pd.concat([TS[TS.ep == es], TN[TN.ep == en]]).drop_duplicates(subset=["time", "lat", "lon"])
        ax.scatter((to_pt(tr.time) - t0) / pd.Timedelta(days=1), along(tr.lat, tr.lon), s=3, c="#d6d5d0", lw=0, rasterized=True)
        for h, c, lab in [(hs, "#2a78d6", f"SVI LFEs ({es})"), (hn, "#1baf7a", f"Washington-NW LFEs ({en})")]:
            ax.scatter((to_pt(h.time) - t0) / pd.Timedelta(days=1), along(h.lat, h.lon) + np.random.default_rng(1).normal(0, 0.6, len(h)),
                       s=5, c=c, lw=0, label=lab, rasterized=True)
        ax.set_title(f"SVI {es} + Washington-NW {en} on one along-margin axis (km toward N30W from 48.3°N 123.4°W)", loc="left")
        ax.set_xlabel(f"days after {t0:%Y-%m-%d} (PT)"); ax.set_ylabel("along-margin km\n(SE − / + NW)")
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True); ax.legend(loc="upper left", fontsize=7.5, frameon=False, markerscale=3)
        fs = to_pt(hs.time).quantile(0.05); fn = to_pt(hn.time).quantile(0.05)
        rows.append(dict(svi_episode=es, nw_episode=en, svi_first5pct_PT=f"{fs:%Y-%m-%d %H:%M}", nw_first5pct_PT=f"{fn:%Y-%m-%d %H:%M}",
                         nw_lead_days=(fs - fn) / pd.Timedelta(days=1)))
    fig.savefig(f"{FIG}/s3_cross_region.png", dpi=130); plt.close(fig)
    pd.DataFrame(rows).to_csv(f"{OUT}/cross_region_overlaps.csv", index=False)
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    cross_region()
