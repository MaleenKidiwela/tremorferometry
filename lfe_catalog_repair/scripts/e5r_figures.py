#!/usr/bin/env python
"""E5 step 5 STEP G figures (figures/e5r_*.png; colorblind-safe Okabe-Ito palette, labelled, units)."""
import json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"; RES = f"{ROOT}/results"; FIG = f"{ROOT}/figures"
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00", purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", black="#000000", grey="#7f7f7f")
T0 = pd.Timestamp("2012-01-01")
def dt(doy): return T0 + pd.to_timedelta(np.asarray(doy) - 1, unit="D")

# 1 download coverage
D = pd.read_csv(f"{RES}/e5r_download_table.csv")
fig, ax = plt.subplots(figsize=(11, 4))
col = [C["orange"] if s.startswith("PB") else C["blue"] for s in D.station]
ax.bar(D.station, D.days_data, color=col); ax.bar(D.station, D.days_nodata, bottom=D.days_data, color=C["red"], label="no-data days")
for i, r in enumerate(D.itertuples()):
    ax.text(i, r.days_data + 0.5, f"{r.size_MB/1e3:.1f} GB\n{r.fs_hz} Hz\nmin cov {r.coverage_min:.2f}", ha="center", va="bottom", fontsize=7)
ax.axhline(52, color=C["grey"], ls="--", lw=0.8); ax.set_ylim(0, 62); ax.set_ylabel("station-days with data (of 52)")
ax.set_title(f"E5 step 5 download: 15 stations, 2012-08-15..10-05, 3-C native rate; {D.size_MB.sum()/1e3:.1f} GB (blue = CN/UW surface BH?/HH?, orange = PB borehole EH?)")
ax.legend(loc="lower right", fontsize=8); ax.tick_params(axis="x", rotation=45)
fig.tight_layout(); fig.savefig(f"{FIG}/e5r_download_coverage.png", dpi=130); plt.close(fig)

# 2 template S-SNR vs N per station class
Q = pd.read_csv(f"{RES}/e5r_template_qc_2012.csv")
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for k, c in (("surface_bb", C["blue"]), ("borehole", C["orange"])):
    g = Q[Q.klass == k]
    ax[0].scatter(g.n_events, g.snr_db, s=9, alpha=0.5, color=c, label=f"{k} (stack S-SNR)")
    ax[0].scatter(g.n_events, g.snr_null_db, s=6, alpha=0.3, color=c, marker="x")
    u = g[g.usable & g.template]
    med = [u[f"snr_N{n}"].median() for n in (30, 60, 100, 200)]
    ax[1].plot([30, 60, 100, 200], med, "o-", color=c, label=f"{k}: median of nested sub-stacks (n={len(u)} usable)")
ax[0].set_xscale("log"); ax[0].set_xlabel("build-half events stacked N"); ax[0].set_ylabel("S-SNR (dB); x = same-N random-time null"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
ax[0].set_title("from-scratch 2012 templates: S-SNR vs N (all 15 stations x 116 families)")
for s, c in (("PGC", C["green"]), ("B011", C["purple"]), ("NLLB", C["sky"]), ("B003", C["red"])):
    u = Q[(Q.sta == s) & Q.usable & Q.template]
    if len(u):
        ax[1].plot([30, 60, 100, 200], [u[f"snr_N{n}"].median() for n in (30, 60, 100, 200)], "s--", ms=4, color=c, label=f"{s} (n={len(u)})")
ax[1].set_xlabel("N events in the sub-stack"); ax[1].set_ylabel("S-SNR (dB) at the full-stack S pick"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3); ax[1].set_title("S-SNR growth with N (usable templates)")
fig.tight_layout(); fig.savefig(f"{FIG}/e5r_template_snr_vs_N.png", dpi=130); plt.close(fig)

# 3 quiet-day FA vs threshold
fig, ax = plt.subplots(figsize=(7, 4.5))
for tag, ls in (("e5r_snr", "-"), ("e5r_eq", "--")):
    p = f"{RES}/e5r_score_{tag}_fa.csv"
    if not os.path.exists(p):
        continue
    F = pd.read_csv(p)
    for L, c in ((8.0, C["blue"]), (12.0, C["orange"])):
        g = F[F.L == L]
        ax.plot(g["T"], g.quiet_fa_per_tpl_day, ls, marker="o", color=c, label=f"{'SNR-weighted' if tag=='e5r_snr' else 'equal-weight'} {int(L)} s: quiet-day detections")
        ax.plot(g["T"], g.quiet_rolled_per_tpl_day, ls, marker="x", color=c, alpha=0.5, lw=0.8, label=f"   rolled control {int(L)} s")
ax.axhline(1.0, color=C["red"], lw=1, label="pre-registered bar: 1 / template / day"); ax.set_yscale("log")
ax.set_xlabel("network threshold (MAD)"); ax.set_ylabel("detections per template per quiet day (2012-08-15..08-29)"); ax.grid(alpha=0.3); ax.legend(fontsize=7)
ax.set_title("false-alarm rate on the 15 pre-SSE quiet days")
fig.tight_layout(); fig.savefig(f"{FIG}/e5r_quietday_fa.png", dpi=130); plt.close(fig)

# 4 recall by magnitude: measured vs predicted vs Bostock density
p = f"{RES}/e5r_score_e5r_snr_recall_mag.csv"
if os.path.exists(p):
    R = pd.read_csv(p); P = pd.read_csv(f"{RES}/e5_prediction_2012_recall.csv"); P = P[(P.set == "set15") & (P["T"] == 8.0) & (P.g == "g0.86")]
    gt = pd.read_parquet(f"{RES}/e5r_gt_bostock2012.parquet"); tst = gt[(gt.half == "test") & ~gt.guard]
    bins = ["M<1.6", "M1.6-1.8", "M1.8-2.0", "M>=2.0"]; edges = [(-9, 1.6), (1.6, 1.8), (1.8, 2.0), (2.0, 9)]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for i, L in enumerate((8.0, 12.0)):
        g = R[(R.L == L) & (R["T"] == 8.0)].set_index("mbin")
        x = np.arange(4)
        ax[i].bar(x - 0.3, [g.loc[b, "recall_own_2s"] for b in bins], 0.2, color=C["blue"], label="measured: own-family template, +-2 s")
        ax[i].bar(x - 0.1, [g.loc[b, "recall_cat_2s"] for b in bins], 0.2, color=C["sky"], label="measured: catalog (any family, declustered)")
        for alpha, c, m in ((0.38, C["orange"], "o"), (0.2, C["green"], "s")):
            q = P[(P.L == L) & (P["mode"].astype(str).str.contains("snr")) & np.isclose(P.alpha, alpha)]
            if len(q):
                ax[i].plot(x + 0.1, [float(q.iloc[0][f"recall_{b}"]) for b in bins], m, color=c, ms=7, label=f"E5 prediction, SNR-weighted, alpha={alpha}")
        ax[i].plot(x, [0.0, 0.20, 0.35, 0.60], "_", color=C["red"], ms=25, mew=2, label="pre-registered bar")
        ax[i].axhline(0.25, color=C["red"], ls=":", lw=1)
        ax2 = ax[i].twinx(); ax2.bar(x + 0.3, [int(((tst.mag >= lo) & (tst.mag < hi)).sum()) for lo, hi in edges], 0.2, color=C["grey"], alpha=0.5, label="Bostock test-half events")
        ax2.set_ylabel("Bostock test-half events per bin", color=C["grey"])
        ax[i].set_xticks(x); ax[i].set_xticklabels(bins); ax[i].set_ylim(0, 1); ax[i].set_ylabel("recall of Bostock test-half events, 8 MAD"); ax[i].set_title(f"{int(L)}-s templates, 8 MAD, 15 stations")
        ax[i].grid(alpha=0.3, axis="y"); ax[i].legend(fontsize=7, loc="upper left")
    fig.tight_layout(); fig.savefig(f"{FIG}/e5r_recall_by_mag.png", dpi=130); plt.close(fig)

# 5 daily counts ours vs his
p = f"{RES}/e5r_score_e5r_snr_daily.csv"
if os.path.exists(p):
    Dd = pd.read_csv(p); g = Dd[(Dd.L == 8.0) & (Dd["T"] == 8.0)]
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.plot(dt(g.doy), g.bostock, "-", color=C["black"], lw=1.5, label="Bostock 2012 catalog (all events)")
    ax.plot(dt(g.doy), g.ours, "-", color=C["blue"], lw=1.5, label="ours: 8 s / 8 MAD, cross-family declustered")
    ax.plot(dt(g.doy), g.ours_matched, "--", color=C["sky"], lw=1.2, label="ours matched to a Bostock event (+-2 s)")
    ax.plot(dt(g.doy), g.rolled, ":", color=C["red"], lw=1.2, label="rolled control (summed over families)")
    g12 = Dd[(Dd.L == 12.0) & (Dd["T"] == 8.0)]
    ax.plot(dt(g12.doy), g12.ours, "-", color=C["orange"], lw=1, alpha=0.8, label="ours: 12 s / 8 MAD")
    ax.axvspan(dt(228), dt(243), color=C["grey"], alpha=0.12, label="pre-SSE quiet days (FA reference)")
    ax.set_yscale("symlog", linthresh=10); ax.set_ylabel("events per UTC day"); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper left")
    ax.set_title("daily counts through the Sept-2012 SSE: our from-scratch catalog vs Bostock")
    fig.tight_layout(); fig.savefig(f"{FIG}/e5r_daily_counts.png", dpi=130); plt.close(fig)

# 6/7 front tracks and secondary fronts, ours vs his
def load_h(tag):
    f = {}
    for k in ("leadtrack", "firstburst", "secondary", "families"):
        # the harness writes families per RUN tag (no _mcall suffix), the others per Mc cut
        p = f"{RES}/e5_thomas_{k}_{tag if k != 'families' else tag.rsplit('_mcall', 1)[0]}.csv"
        f[k] = pd.read_csv(p) if os.path.exists(p) else None
    p = f"{RES}/e5_thomas_summary_{tag.rsplit('_mcall', 1)[0]}.csv"
    f["summary"] = pd.read_csv(p) if os.path.exists(p) else None
    return f
H = load_h("bostock2012_mcall"); O = load_h("ours2012_L8_mcall")
if O["leadtrack"] is not None and H["leadtrack"] is not None:
    # common along-strike origin: shift ours by the mean difference over shared families
    fh, fo = H["families"], O["families"]
    m = fh.merge(fo, on=["lat", "lon", "dep"], suffixes=("_h", "_o"))
    dx = float((m.x_strike_h - m.x_strike_o).mean()) if len(m) else 0.0
    sh, so = H["summary"], O["summary"]
    sh = sh[sh.mc == "none"].iloc[0]; so = so[so.mc == "none"].iloc[0]
    fig, ax = plt.subplots(figsize=(12, 5))
    for f_, c, name, s_, off in ((H, C["black"], "Bostock", sh, 0.0), (O, C["blue"], "ours (8 s / 8 MAD)", so, dx)):
        lt = f_["leadtrack"]; t = pd.to_datetime(lt.t_mid, unit="s")
        ax.plot(t, lt.x_lead + off, "-", color=c, lw=1.5, label=f"{name}: 6-h leading edge (10th pct), {s_.front_lead_kmday:+.2f} km/day in window")
        ax.plot(t, lt.x_mean + off, ":", color=c, lw=1, label=f"{name}: 6-h centroid")
        fb = f_["firstburst"]
        ax.scatter(pd.to_datetime(fb.t_fb, unit="s"), fb.x + off, s=14, color=c, alpha=0.6, marker="o" if name == "Bostock" else "^",
                   label=f"{name}: family first burst; Theil-Sen {s_.front_firstburst_kmday:+.2f} km/day (CI {s_.front_firstburst_ci_lo:+.2f}..{s_.front_firstburst_ci_hi:+.2f})")
        ax.axvspan(pd.Timestamp(s_.win_start), pd.Timestamp(s_.win_end), color=c, alpha=0.06)
    ax.set_ylabel("along-strike position (km, N320E positive = NW; Bostock's origin)"); ax.set_xlabel("2012 (UTC)"); ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="lower left")
    ax.set_title("main-front tracks: ours vs Bostock (shaded = each run's propagation window; reference 7.3 km/day, accept 5-10)")
    fig.tight_layout(); fig.savefig(f"{FIG}/e5r_front_tracks.png", dpi=130); plt.close(fig)
    # secondary fronts
    fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True, sharey=True)
    for i, (f_, c, name, s_, off, tag) in enumerate(((H, C["black"], "Bostock", sh, 0.0, "bostock2012_mcall"), (O, C["blue"], "ours (8 s / 8 MAD)", so, dx, "ours2012_L8_mcall"))):
        catp = f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv" if i == 0 else f"{RES}/e5_catalog_2012_sVI.parquet"
        fam = f_["families"].set_index("family_id")
        if i == 0:
            cc_ = pd.read_csv(catp); cc_["time"] = pd.to_datetime(dict(year=cc_.year, month=cc_.month, day=cc_.day, hour=cc_.hour, minute=cc_["min"], second=cc_.sec.astype(int)))
            cc_ = cc_[(cc_.time >= "2012-08-30") & (cc_.time < "2012-09-27")]
            key = fam.reset_index().set_index(["lat", "lon", "dep"]).x_strike
            xs = pd.Series(list(zip(cc_.lat, cc_.lon, cc_.dep))).map(key.to_dict()).to_numpy()
        else:
            cc_ = pd.read_parquet(catp); xs = cc_.family_id.map(fam.x_strike).to_numpy() + off
        ax[i].scatter(cc_.time, xs, s=1, color=C["grey"], alpha=0.35, label=f"{name}: events ({len(cc_)})")
        lt = f_["leadtrack"]; ax[i].plot(pd.to_datetime(lt.t_mid, unit="s"), lt.x_lead + off, "-", color=c, lw=1, label="leading edge")
        sf = f_["secondary"]
        if sf is not None and len(sf):
            for r in sf.itertuples():
                ax[i].plot([pd.to_datetime(r.t_start, unit="s"), pd.to_datetime(r.t_end, unit="s")], [r.core_x_min + off, r.core_x_max + off] if r.speed_kmh > 0 else [r.core_x_max + off, r.core_x_min + off],
                           "-", color=C["red"], lw=2.5, alpha=0.9)
            ax[i].plot([], [], "-", color=C["red"], lw=2.5, label=f"secondary fronts detected: {len(sf)} (null circ {s_.n_secondary_null_circ:.1f}, day-shuffle {s_.n_secondary_null_dayshuf:.1f}); core extent, direction = sign of speed")
        ax[i].set_ylabel("along-strike (km)"); ax[i].grid(alpha=0.3); ax[i].legend(fontsize=7, loc="lower left"); ax[i].set_title(name)
    ax[1].set_xlabel("2012 (UTC)"); ax[0].set_xlim(pd.Timestamp("2012-09-04"), pd.Timestamp("2012-09-24"))
    fig.tight_layout(); fig.savefig(f"{FIG}/e5r_secondary_fronts.png", dpi=130); plt.close(fig)
print("[figures] written:", sorted(f for f in os.listdir(FIG) if f.startswith("e5r_")))
