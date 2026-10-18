#!/usr/bin/env python
"""E5 figures (Okabe-Ito colourblind-safe; every axis labelled with units).  Each panel is skipped if its input
is not on disk yet, so the script can be re-run as the steps finish."""
from __future__ import annotations
import os, sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
FIG = "/home/jovyan/tremorferometry/lfe_catalog_repair/figures"
OI = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": .25, "figure.dpi": 140, "savefig.bbox": "tight"})
ex = lambda f: os.path.exists(f"{RES}/{f}")


def fig_fa():
    tags = [("e5_budget_op", "operating subset (PGC/VGZ/SNB + B011/B926/LRIV/NLLB)"),
            ("e5_budget_ring", "full ring + transfers (10 stations)"), ("e5_budget_ring6", "E4 6-station ring")]
    tags = [t for t in tags if ex(f"e5_budget_table_{t[0]}.csv")]
    if not tags:
        return
    fig, ax = plt.subplots(2, len(tags), figsize=(4.2 * len(tags), 7.0), squeeze=False)
    fig.subplots_adjust(hspace=0.42, wspace=0.28)
    for j, (tag, lab) in enumerate(tags):
        T = pd.read_csv(f"{RES}/e5_budget_table_{tag}.csv")
        a = ax[0, j]
        for i, L in enumerate(sorted(T.L.unique())):
            g = T[T.L == L].sort_values("T")
            a.semilogy(g["T"], g.quiet_fa_mean.clip(1e-3), "o-", color=OI[i + 1], label=f"{L:.0f} s quiet mean")
            a.semilogy(g["T"], g.quiet_fa_p90.clip(1e-3), "s--", color=OI[i + 1], ms=3, lw=.8, label=f"{L:.0f} s quiet p90")
            a.semilogy(g["T"], g.ets_det_mean.clip(1e-3), "^:", color=OI[i + 1], ms=3, lw=.8, label=f"{L:.0f} s ETS-day total")
            a.semilogy(g["T"], g.quiet_roll_mean.clip(1e-3), "x", color=OI[i + 1], ms=4, label=f"{L:.0f} s quiet rolled ctrl")
        a.axhline(1.0, color="k", lw=.8); a.axhline(10.0, color="k", lw=.5, ls=":")
        a.set_xlabel("network threshold (MAD)"); a.set_ylabel("detections / template / day")
        a.set_title(f"(a{j+1}) {lab}\nquiet days (<= 50 Bostock ev/day) vs ETS days", fontsize=7.5)
        if j == 0:
            a.legend(fontsize=5, ncol=2)
        a = ax[1, j]
        for i, L in enumerate(sorted(T.L.unique())):
            g = T[T.L == L].sort_values("T")
            a.plot(g["T"], g.ets_precision_own, "o-", color=OI[i + 1], label=f"{L:.0f} s precision (own family, +-2 s)")
            a.plot(g["T"], g.ets_precision_any, "s--", color=OI[i + 1], ms=3, lw=.8, label=f"{L:.0f} s precision (any family)")
            a.plot(g["T"], g.recall_all_events, "^:", color=OI[i + 1], ms=3, label=f"{L:.0f} s recall (all Bostock events)")
        a.axhline(0.85, color="k", lw=.8)
        a.set_xlabel("network threshold (MAD)"); a.set_ylabel("fraction"); a.set_ylim(0, 1)
        a.set_title(f"(b{j+1}) ETS-day precision and recall", fontsize=7.5)
        if j == 0:
            a.legend(fontsize=5)
    fig.suptitle("E5 step 1 -- detection budget on the 2010 data (Bostock templates + build-half transfers)", fontsize=9)
    fig.savefig(f"{FIG}/e5_fa_budget.png"); plt.close(fig)


def fig_recall_curves():
    if not ex("e5_prediction_measured_curve.csv"):
        return
    C = pd.read_csv(f"{RES}/e5_prediction_measured_curve.csv")
    SP = pd.read_csv(f"{RES}/e5_prediction_spread.csv"); INS = pd.read_csv(f"{RES}/e5_prediction_insample.csv")
    NET = pd.read_csv(f"{RES}/e4_step1_network.csv"); bg = float(NET[(NET.L == 8.0) & (NET.magbin == "all")].net_rand_bg.iloc[0])
    sig = float(SP[SP.L == 8.0].sigma_log_pooled.iloc[0]); sigw = float(SP[SP.L == 8.0].sigma_log_within_family_robust.iloc[0])
    e_ins = float(INS[(INS.source == "station_measured") & (INS.sigma == "pooled")].pred_raw_median.iloc[0]) - bg
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.8))
    M = np.linspace(1.2, 2.5, 100)
    a = ax[0]
    for i, T in enumerate((6, 8, 10)):
        g = C[C.L == 8.0]
        a.plot(g.m_mid, g[f"rec{T}"], "o", color=OI[i + 1], ms=4, label=f"measured 2010 ring, >= {T} MAD")
        for al, ls in ((0.38, "-"), (0.46, "--"), (0.58, ":")):
            rm = bg + e_ins * 10 ** (al * (M - 1.95))
            a.plot(M, 1 - norm.cdf((np.log10(T) - np.log10(rm)) / sig), ls, color=OI[i + 1], lw=1,
                   label=f"model alpha {al}" if i == 1 else None)
    a.set_xlabel("Bostock magnitude"); a.set_ylabel("recall of catalogued events"); a.set_ylim(0, 1)
    a.set_title("(a) measured vs modelled recall, 8 s, 2010 ring\n(raw = 2.73 + e_net 10^(alpha dM); pooled log-spread)", fontsize=7.5)
    a.legend(fontsize=5.5)
    a = ax[1]
    for i, L in enumerate((8.0, 12.0)):
        g = C[C.L == L]
        for T, mk in ((8, "o"), (10, "s")):
            a.plot(g.m_mid, g[f"rec{T}"], mk + "-", color=OI[i + 1], ms=3, label=f"{L:.0f} s, >= {T} MAD")
    a.set_xlabel("Bostock magnitude"); a.set_ylabel("recall"); a.set_ylim(0, 1); a.legend(fontsize=6)
    a.set_title("(b) measured recall by template length (2010 ring, +-1 s)", fontsize=7.5)
    a = ax[3]
    if ex("e5_prediction_measured_curve_within.csv"):
        CW = pd.read_csv(f"{RES}/e5_prediction_measured_curve_within.csv")
        for i, L in enumerate((8.0, 12.0)):
            g = CW[CW.L == L]
            for T, mk in ((6, "^"), (8, "o"), (10, "s")):
                a.plot(g.dM_mid, g[f"rec{T}"], mk + "-", color=OI[i + 1], ms=3, label=f"{L:.0f} s, >= {T} MAD")
        for al, ls in ((0.38, "-"), (0.58, ":")):
            dm = np.linspace(-0.6, 0.6, 50); rm = bg + e_ins * 10 ** (al * dm)
            a.plot(dm, 1 - norm.cdf((np.log10(8) - np.log10(rm)) / sig), ls, color="k", lw=.8, label=f"model >= 8, alpha {al}")
        a.set_xlabel("magnitude minus the family's median magnitude"); a.set_ylabel("recall"); a.set_ylim(0, 1); a.legend(fontsize=5.5)
        a.set_title("(d) WITHIN-family measured recall (distance confound removed)", fontsize=7.5)
    a = ax[2]
    if ex("e5_prediction_2012_recall.csv"):
        R = pd.read_csv(f"{RES}/e5_prediction_2012_recall.csv")
        bins = ["M<1.6", "M1.6-1.8", "M1.8-2.0", "M>=2.0"]
        k = 0
        for sname, ls in (("set15", "-"), ("set18_with_LZB_YOUB_PFB", "--")):
            for L in (8.0, 12.0):
                for mode in ("equal", "snrw"):
                    r = R[(R["set"] == sname) & (R.L == L) & (R["mode"] == mode) & (R["g"] == "g0.86") & (R.alpha == 0.38) & (R["T"] == 8.0)]
                    if not len(r):
                        continue
                    r = r.iloc[0]
                    a.plot(range(4), [r[f"recall_{b}"] for b in bins], ls, marker="o", ms=3, color=OI[(k % 7) + 1],
                           label=f"{'15 stations' if sname=='set15' else '+LZB/YOUB/PFB'}, {L:.0f} s, {mode} (all {r.recall_all:.2f})")
                    k += 1
        a.set_xticks(range(4)); a.set_xticklabels(bins); a.set_ylim(0, 1); a.set_ylabel("predicted recall of Bostock 2012 events at 8 MAD")
        a.set_title("(c) 2012 prediction, alpha 0.38, g 0.86, pooled spread", fontsize=7.5); a.legend(fontsize=5)
    fig.savefig(f"{FIG}/e5_recall_curves.png"); plt.close(fig)


def fig_weighting():
    runs = [("e5_w_op_eq", "op subset, equal"), ("e5_w_op_snr", "op subset, SNR-weighted"),
            ("e5_w_ring_eq", "ring+transfers, equal"), ("e5_w_ring_snr", "ring+transfers, SNR-weighted")]
    if not all(ex(f"e4_tt_events_{t}.parquet") for t, _ in runs):
        return
    if not ex("e5_weighting_table.csv"):
        return
    W = pd.read_csv(f"{RES}/e5_weighting_table.csv")
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
    order = ["M<1.5", "M1.5-1.8", "M1.8-2.1", "M>=2.1"]
    for i, (t, lab) in enumerate(runs):
        g = W[(W.tag == t) & (W.L == 8.0) & (W.fams == "no_1_15")].set_index("magbin")
        ax[0].plot(range(4), [g.excess.get(m, np.nan) for m in order], "o-", color=OI[i + 1], label=lab)
        ax[1].plot(range(4), [g.frac_ge8_2s.get(m, np.nan) for m in order], "o-", color=OI[i + 1], label=lab)
        ax[1].plot(range(4), [g.shift_frac_ge8.get(m, np.nan) for m in order], "x:", color=OI[i + 1], lw=.7)
    for a, yl, tt in ((ax[0], "held-out network excess over random (MAD)", "(a) test-half excess, 8 s, families 1/15 excluded"),
                      (ax[1], "fraction >= 8 MAD (x: rolled control)", "(b) held-out recall at 8 MAD")):
        a.set_xticks(range(4)); a.set_xticklabels(order); a.set_ylabel(yl); a.set_title(tt, fontsize=7.5)
    ax[0].axhline(3.3, color="k", lw=.7); ax[1].axhline(0.22, color="k", lw=.7); ax[0].legend(fontsize=6)
    if ex("e5_weighting_corr.csv"):
        Cc = pd.read_csv(f"{RES}/e5_weighting_corr.csv")
        ax[2].scatter(Cc.snr_db, Cc.excess_test, s=8, c=[OI[1 + i % 7] for i in pd.factorize(Cc.sta)[0]])
        for s_, g in Cc.groupby("sta"):
            ax[2].annotate(s_, (g.snr_db.median(), g.excess_test.median()), fontsize=6)
        ax[2].set_xlabel("build-half stack S-SNR (dB)"); ax[2].set_ylabel("test-half single-station excess (MAD)")
        ax[2].set_title(f"(c) do build weights predict test registration?\nSpearman rho = {Cc.attrs.get('rho', float(Cc.rho.iloc[0])) if 'rho' in Cc else np.nan:.2f}", fontsize=7.5)
    fig.savefig(f"{FIG}/e5_weighting.png"); plt.close(fig)


def fig_front():
    if not ex("e5_thomas_summary_bostock2012.csv"):
        return
    S = pd.read_csv(f"{RES}/e5_thomas_summary_bostock2012.csv")
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    a = ax[0]
    t0 = pd.Timestamp("2012-08-30").value / 1e9
    for i, tag in enumerate(["mcall", "mc1.8", "mc2.0"]):
        trk = pd.read_csv(f"{RES}/e5_thomas_leadtrack_bostock2012_{tag}.csv")
        a.plot((trk.t_mid - t0) / 86400, trk.x_lead, "-", color=OI[i + 1], lw=1, label=f"leading edge (p10), {tag}")
        a.plot((trk.t_mid - t0) / 86400, trk.x_mean, ":", color=OI[i + 1], lw=.8, label=f"centroid, {tag}")
    FB = pd.read_csv(f"{RES}/e5_thomas_firstburst_bostock2012_mcall.csv")
    a.plot(FB.day, FB.x, "k.", ms=3, label="first burst per family (>= 4 ev, 10-min gaps)")
    r = S[S.mc == "none"].iloc[0]
    ws = (pd.Timestamp(r.win_start).value / 1e9 - t0) / 86400; we = (pd.Timestamp(r.win_end).value / 1e9 - t0) / 86400
    a.axvspan(ws, we, color=OI[3], alpha=.12, label=f"propagation window {r.win_start[5:10]}..{r.win_end[5:10]}")
    xx = np.array([ws, we]); a.plot(xx, r.front_firstburst_kmday * (xx - ws) + FB[(FB.day >= ws) & (FB.day <= we)].x.max(), "k--", lw=1,
                                   label=f"first-burst Theil-Sen {r.front_firstburst_kmday:.1f} km/day (CI {r.front_firstburst_ci_lo:.1f}..{r.front_firstburst_ci_hi:.1f})")
    a.set_xlabel("days after 2012-08-30 00:00 UTC"); a.set_ylabel("along-strike position (km; + = NW)")
    a.set_title("(a) main-front tracks, Bostock 2012\n(reference 7.3 km/day = 80 km / 11 d; accept 5-10)", fontsize=7.5); a.legend(fontsize=5)
    a = ax[1]
    sf = pd.read_csv(f"{RES}/e5_thomas_secondary_bostock2012_mcall.csv")
    ev = pd.read_parquet(f"{RES}/e4_thomas_events_bostock2012_mcall.parquet") if ex("e4_thomas_events_bostock2012_mcall.parquet") else None
    if ev is not None:
        a.scatter((ev.t - t0) / 86400, ev.x_strike, s=.3, c="#bbbbbb", rasterized=True)
    for i, rr in enumerate(sf.itertuples()):
        a.plot([(rr.t_start - t0) / 86400, (rr.t_end - t0) / 86400], [rr.core_x_min, rr.core_x_max] if hasattr(rr, "core_x_min") else [rr.x_min, rr.x_max],
               "-", color=OI[6], lw=1.5)
        a.plot((rr.t_mid - t0) / 86400, rr.x_mean, "o", color=OI[6], ms=3)
    fam = pd.read_csv(f"{RES}/e5_thomas_families_bostock2012.csv")
    x141 = fam.loc[fam.bostock_id == 141, "x_strike"]
    if len(x141):
        a.axhline(float(x141.iloc[0]), color=OI[5], lw=.8, ls="--", label="family 141 (Thomas: six fronts 12-16 Sept)")
    a.axvspan(15, 16, color=OI[1], alpha=.15, label="14 Sept (Thomas Fig. 2C, ~12 km/h)")
    a.set_xlim(6, 26); a.set_xlabel("days after 2012-08-30"); a.set_ylabel("along-strike position (km)")
    a.set_title(f"(b) secondary fronts: {len(sf)} detected (nulls: circular {r.n_secondary_null_circ:.1f}, day-shuffle {r.n_secondary_null_dayshuf:.1f})\n"
                f"bars = core extent (events within 5 km of the fit)", fontsize=7.5); a.legend(fontsize=5.5)
    a = ax[2]
    mcs = S.mc.tolist(); x = np.arange(len(mcs))
    a.plot(x, -S.front_firstburst_kmday, "o-", color=OI[1], label="main front, first-burst (km/day)")
    a.fill_between(x, -S.front_firstburst_ci_hi, -S.front_firstburst_ci_lo, color=OI[1], alpha=.15)
    a.plot(x, -S.front_lead_kmday, "s--", color=OI[2], label="main front, leading edge (km/day)")
    a.axhspan(5, 10, color="k", alpha=.06); a.axhline(7.3, color="k", lw=.6)
    a2 = a.twinx()
    a2.plot(x, S.n_secondary, "^-", color=OI[6], label="secondary fronts detected")
    a2.plot(x, S.n_secondary_null_circ_p95, "^:", color=OI[6], lw=.8, label="circular-shift null p95")
    a2.plot(x, S.n_secondary_null_dayshuf_p95, "v:", color=OI[7], lw=.8, label="day-shuffle null p95")
    a2.plot(x, 100 * S.frac_fam_reactivated, "d-.", color=OI[3], label="families re-activated (>= 12 h quiet), %")
    a2.set_ylabel("count / percent")
    a.set_xticks(x); a.set_xticklabels([f"{m}\n{k*100:.0f}%" for m, k in zip(mcs, S.kept_frac)]); a.set_xlabel("completeness cut Mc (events kept)")
    a.set_ylabel("front speed (km/day, southward +)"); a.set_ylim(0, 12)
    a.set_title("(c) decimation: what each metric needs", fontsize=7.5)
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, fontsize=5.5, loc="upper left")
    fig.savefig(f"{FIG}/e5_thomas_harness.png"); plt.close(fig)


if __name__ == "__main__":
    for f in (fig_fa, fig_recall_curves, fig_weighting, fig_front):
        try:
            f(); print(f"[fig] {f.__name__} ok")
        except Exception as e:
            print(f"[fig] {f.__name__} FAILED: {e!r}")
