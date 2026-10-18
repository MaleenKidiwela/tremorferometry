#!/usr/bin/env python
"""E6 figures (figures/e6_*.png; Okabe-Ito colorblind-safe palette, labelled axes, units).  Reads whatever E6 results exist."""
import glob, json, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"; RES = f"{ROOT}/results"; FIG = f"{ROOT}/figures"
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00", purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", black="#000000", grey="#7f7f7f")
OP_TAG = sys.argv[1] if len(sys.argv) > 1 else "e6_mergeB"; OP_T = float(sys.argv[2]) if len(sys.argv) > 2 else 8.5
RUNS = [("e6_sites", "unmerged 15-station (E5r network)", C["black"]), ("e6_mergeA", "site-merged A (max w)", C["blue"]), ("e6_mergeB", "site-merged B (sqrt sum w^2)", C["orange"]),
        ("e6_mergeB_cap25", "B + PGC site <= 25 % weight", C["green"]), ("e6_mergeB_coinc", "B + 2-site coincidence (>= 3 MAD)", C["purple"]),
        ("e6_lev_a", "B + lever (a) far-station picks", C["sky"]), ("e6_lev_b", "B + lever (b) >= 5 dB stations", C["yellow"]), ("e6_lev_c", "B + lever (c) 2010+2012 templates", C["red"])]
FA = {t: pd.read_csv(f"{RES}/e6_fa_table_{t}.csv") for t, _, _ in RUNS if os.path.exists(f"{RES}/e6_fa_table_{t}.csv")}

# 1 data vs station-rolled vs site-rolled by threshold (quiet, post-SSE, ETS) for the unmerged network and the two variants
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
for i, (win, title) in enumerate((("quiet", "pre-SSE quiet days 08-15..08-29"), ("post", "post-SSE 09-26..10-05"), ("ets", "ETS days (>= 100 Bostock events)"))):
    for tag, lab, c in RUNS[:3]:
        if tag not in FA:
            continue
        F = FA[tag]
        ax[i].plot(F["T"], F[f"{win}_data"], "-o", color=c, label=f"{lab}: data")
        ax[i].plot(F["T"], F[f"{win}_roll_site"], "--s", color=c, alpha=0.8, label="   site-rolled control")
        ax[i].plot(F["T"], F[f"{win}_roll_sta"], ":x", color=c, alpha=0.6, label="   station-rolled control")
    ax[i].axhline(1.0, color=C["red"], lw=1, ls="-.", label="bar: 1 / template / day" if i < 2 else None)
    ax[i].set_yscale("log"); ax[i].set_xlabel("network threshold (MAD)"); ax[i].grid(alpha=0.3); ax[i].set_title(title, fontsize=10)
ax[0].set_ylabel("declustered detections per template per day"); ax[0].legend(fontsize=6.5, loc="lower left")
fig.suptitle("E6 step 1-2: the SITE-rolled control (one shift per co-located site) tracks the data on quiet and post-SSE days; the station-rolled control does not", fontsize=10)
fig.tight_layout(); fig.savefig(f"{FIG}/e6_fa_rates.png", dpi=130); plt.close(fig)

# 2 PGC-site share of the network sum at detections (unmerged network)
p = f"{RES}/e5r_detections_e6_sites.parquet"
if os.path.exists(p):
    d = pd.read_parquet(p); d = d[(d.L == 8) & (d.net_mad >= 8)]
    gt = pd.read_parquet(f"{RES}/e5r_gt_bostock2012.parquet"); gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    pdd = gt.groupby("doy").size(); ETS = set(pdd[pdd >= 100].index)
    sh = d.con_PGC / d.net_cc
    fig, ax = plt.subplots(figsize=(8, 4.5)); bins = np.linspace(-0.2, 1.4, 65)
    for k, (lab, m, c) in enumerate((("quiet-day detections (08-15..08-29)", d.doy.isin(range(228, 243)), C["red"]), ("post-SSE detections (09-26..10-05)", d.doy.isin(range(270, 280)), C["orange"]),
                                     ("ETS-day detections matched to an own-family Bostock event", d.doy.isin(ETS) & d.match_own, C["blue"]))):
        ax.hist(sh[m], bins=bins, density=True, histtype="step", lw=1.8, color=c, label=f"{lab}: n={int(m.sum())}, median {sh[m].median():.2f}")
    ax.hist(d.wsh_PGC, bins=bins, density=True, histtype="stepfilled", alpha=0.25, color=C["grey"], label=f"PGC-site WEIGHT share at the same detections (median {d.wsh_PGC.median():.2f})")
    ax.set_xlabel("PGC-site (PGC+B009+B010+B011) share of the weighted network sum at the detection sample"); ax.set_ylabel("density"); ax.grid(alpha=0.3); ax.legend(fontsize=7)
    ax.set_title("E6 step 1: 8 s / >= 8 MAD, unmerged 15-station SNR-weighted network", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/e6_pgc_share.png", dpi=130); plt.close(fig)

# 3 stack test at the excluded stations
p = f"{RES}/e6_stack_test_1b.csv"
if os.path.exists(p):
    S = pd.read_csv(p)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for i, sta in enumerate(("VGZ", "SNB")):
        s = S[S.sta == sta]
        for j, r in enumerate(s[s.set == "quiet"].itertuples()):
            Ns = [30, 100, 300, r.N]; v = [r.snr_N30, r.snr_N100, r.snr_N300, r.snr_db]
            ax[i].plot(Ns, v, "-o", color=list(C.values())[j], ms=4, label=f"family {int(r.fam)} quiet-day detections (N={r.N}): {r.verdict}")
        for r in s[s.set == "pooled_quiet"].itertuples():
            ax[i].plot([30, 100, 300, r.N], [r.snr_N30, r.snr_N100, r.snr_N300, r.snr_db], "-D", color=C["black"], lw=2.2, ms=6, label=f"POOLED quiet (N={r.N}): {r.verdict}")
            ax[i].axhspan(r.null_mean_db - 2 * r.null_sd_db, r.null_mean_db + 2 * r.null_sd_db, color=C["grey"], alpha=0.25, label=f"same-N random-time null (mean +- 2 sd, pooled)")
            ax[i].axhline(r.null_mean_db + 3, color=C["red"], ls="--", lw=1, label="REAL bar: null + 3 dB (peak at S)"); ax[i].axhline(r.null_mean_db + 1, color=C["orange"], ls=":", lw=1, label="NOISE bar: <= null + 1 dB at N >= 300")
        for r in s[s.set == "pooled_ets_control"].itertuples():
            ax[i].plot([30, 100, r.N], [r.snr_N30, r.snr_N100, r.snr_db], "-^", color=C["green"], lw=2, ms=7, label=f"POSITIVE CONTROL: ETS own-matched detections (N={r.N}), peak at {r.peak_offset_s:+.2f} s")
        ax[i].set_xscale("log"); ax[i].set_xlabel("N windows stacked (family moveout, no alignment refinement)"); ax[i].set_title(f"CN.{sta} (excluded from the network sum)", fontsize=10); ax[i].grid(alpha=0.3); ax[i].legend(fontsize=6.3, loc="upper left")
    ax[0].set_ylabel("stack S-SNR (dB): energy S-2..S+2 s over the pre-window S-20..S-16 s")
    fig.suptitle("E6 step 1b: are the quiet-day detections of the five worst families LFEs?  3-C 20-Hz stacks at the two excluded stations", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/e6_stack_test.png", dpi=130); plt.close(fig)

# 4 site-merge variants and extras: own recall vs quiet-day / post-SSE FA (points = thresholds)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for tag, lab, c in RUNS:
    if tag not in FA:
        continue
    F = FA[tag]; F = F[F.recall_own_all.notna()]
    for i, win in enumerate(("quiet", "post")):
        ax[i].plot(F[f"{win}_data"], F.recall_own_all, "-o", color=c, ms=5, label=lab)
        for r in F.itertuples():
            ax[i].annotate(f"{r.T:g}", (getattr(r, f"{win}_data"), r.recall_own_all), fontsize=6, xytext=(3, 2), textcoords="offset points", color=c)
for i, win in enumerate(("quiet days 08-15..08-29", "post-SSE 09-26..10-05")):
    ax[i].axvline(1.0, color=C["red"], ls="-.", lw=1, label="FA bar: <= 1 / template / day"); ax[i].axhline(0.25, color=C["red"], ls=":", lw=1, label="recall bar 0.25 (own, overall)"); ax[i].axhline(0.20, color=C["grey"], ls=":", lw=1, label="'merge kept sensitivity' 0.20")
    ax[i].set_xscale("log"); ax[i].set_xlabel(f"detections per template per day, {win}"); ax[i].set_ylabel("own-family recall of Bostock test-half events (+-2 s)"); ax[i].grid(alpha=0.3); ax[i].set_ylim(0, 0.42)
ax[0].legend(fontsize=6.5, loc="upper left")
fig.suptitle("E6 step 2-3: own recall vs false-alarm rate; numbers = network threshold (MAD)", fontsize=10)
fig.tight_layout(); fig.savefig(f"{FIG}/e6_recall_vs_fa.png", dpi=130); plt.close(fig)

# 5 lever table at the operating point
rows = []
for tag, lab, c in RUNS:
    if tag in FA:
        F = FA[tag]; r = F[F["T"] == OP_T]
        if len(r):
            r = r.iloc[0]; rows.append(dict(tag=tag, lab=lab, recall=r.recall_own_all, r16=r["recall_own_M1.6-1.8"], r20=r["recall_own_M>=2.0"], quiet=r.quiet_data, post=r.post_data, ets_site=r.ets_frac_roll_site, c=c))
if rows:
    Lv = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(12, 4.8)); x = np.arange(len(Lv))
    ax.bar(x - 0.2, Lv.recall, 0.2, color=[r for r in Lv.c], label="own recall, all test events")
    ax.bar(x, Lv.r16, 0.2, color=[r for r in Lv.c], alpha=0.6, label="own recall, M1.6-1.8")
    ax.bar(x + 0.2, Lv.r20, 0.2, color=[r for r in Lv.c], alpha=0.35, label="own recall, M>=2.0")
    for i, r in enumerate(Lv.itertuples()):
        ax.text(i, max(r.recall, r.r20) + 0.01, f"quiet {r.quiet:.2f}\npost {r.post:.2f}\nETS site-rolled {r.ets_site:.0%}", ha="center", fontsize=7)
    ax.axhline(0.25, color=C["red"], ls=":", lw=1); ax.axhline(0.20, color=C["grey"], ls=":", lw=1)
    ax.set_xticks(x); ax.set_xticklabels([l.replace(" + ", "\n+ ") for l in Lv.lab], fontsize=7); ax.set_ylabel(f"own-family recall at {OP_T:g} MAD (8 s)"); ax.set_ylim(0, 0.45); ax.grid(alpha=0.3, axis="y"); ax.legend(fontsize=7)
    ax.set_title(f"E6 step 3: recall levers and step-2 extras at {OP_T:g} MAD (labels: quiet / post-SSE detections per template per day; ETS-day site-rolled fraction)", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/e6_lever_table.png", dpi=130); plt.close(fig)

# 6 corrected Thomas comparison with the thinned-reference spread
p = glob.glob(f"{RES}/e6_thinned_reference_{OP_TAG}_T{OP_T:g}_runs.csv")
if p:
    R = pd.read_csv(p[0]); b = pd.read_csv(f"{RES}/e5_thomas_summary_bostock2012.csv"); b = b[b.mc.astype(str) == "none"].iloc[0]
    ours = R[R.kind == "ours"]; ref = R[R.kind == "ref_floor"]
    fig, ax = plt.subplots(1, 4, figsize=(15, 4.3))
    for i, (k, lab) in enumerate((("front_firstburst_kmday", "main front, first-burst Theil-Sen (km/day)"), ("reactivation_interval_d_p50", "re-rupture interval median (d)"),
                                  ("frac_fam_reactivated", "fraction of families re-activated"), ("n_secondary", "secondary-front COUNT (F6: not density-monotonic)"))):
        v = pd.to_numeric(ref[k], errors="coerce").dropna()
        ax[i].violinplot([v], positions=[1], showmedians=True)
        ax[i].scatter(np.full(len(v), 1) + np.random.default_rng(0).normal(0, 0.04, len(v)), v, s=10, color=C["grey"], label=f"thinned reference, {len(v)} seeds")
        ax[i].axhline(float(b[k]), color=C["black"], lw=1.5, label=f"Bostock full: {float(b[k]):.2f}")
        if len(ours):
            ax[i].axhline(float(ours.iloc[0][k]), color=C["blue"], lw=1.5, ls="--", label=f"ours ({OP_TAG} {OP_T:g} MAD): {float(ours.iloc[0][k]):.2f}")
        ax[i].set_title(lab, fontsize=9); ax[i].set_xticks([]); ax[i].grid(alpha=0.3, axis="y"); ax[i].legend(fontsize=6.5)
    fig.suptitle("E6 step 4: Thomas metrics -- ours vs Bostock full vs Bostock thinned per family at our own recall + uniform floor (20 seeds)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/e6_thomas_reference.png", dpi=130); plt.close(fig)
print("[figures] written:", sorted(f for f in os.listdir(FIG) if f.startswith("e6_")))
