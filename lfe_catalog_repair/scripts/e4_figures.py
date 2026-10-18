#!/usr/bin/env python
"""E4 figures. Okabe-Ito colourblind-safe palette; every panel labelled with units."""
from __future__ import annotations
import numpy as np, pandas as pd, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
FIG = "/home/jovyan/tremorferometry/lfe_catalog_repair/figures"
OI = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": .25,
                     "figure.dpi": 140, "savefig.bbox": "tight"})


def f1_mad():
    reg = pd.read_csv(f"{RES}/e4_step1_registration.csv")
    net = pd.read_csv(f"{RES}/e4_step1_network.csv")
    ms = pd.read_csv(f"{RES}/e4_step1_magslope.csv")
    ch = pd.read_csv(f"{RES}/e4_step1_channels.csv")
    order = ["M<1.5", "M1.5-1.8", "M1.8-2.1", "M>=2.1"]
    fig, ax = plt.subplots(2, 2, figsize=(9, 6.6))
    a = ax[0, 0]
    for i, s in enumerate(["PGC", "LZB", "VGZ", "SNB", "YOUB", "PFB"]):
        g = reg[(reg.sta == s) & (reg.L == 8.0) & reg.magbin.isin(order)].set_index("magbin")
        a.plot(range(len(order)), [g.excess_meas.get(m, np.nan) for m in order], "o-",
               color=OI[i + 1], label=s, ms=4)
    a.set_xticks(range(len(order))); a.set_xticklabels(order, rotation=20)
    a.set_ylabel("single-station excess (MAD)\nevents minus random times")
    a.set_title("(a) per-station registration, 8-s template, 2010 outer ring")
    a.legend(fontsize=6, ncol=2)
    a = ax[0, 1]
    for i, L in enumerate(sorted(net.L.unique())):
        g = net[(net.L == L) & net.magbin.isin(order)].set_index("magbin")
        a.plot(range(len(order)), [g.excess.get(m, np.nan) for m in order], "o-",
               color=OI[i + 1], label=f"{L:.0f} s", ms=4)
    a.set_xticks(range(len(order))); a.set_xticklabels(order, rotation=20)
    a.set_ylabel("network excess (MAD)")
    a.set_title("(b) network (6 stations, median 3 with a template)")
    a.legend(title="template", fontsize=6)
    a = ax[1, 0]
    for i, L in enumerate(sorted(net.L.unique())):
        g = net[(net.L == L) & net.magbin.isin(order)].set_index("magbin")
        a.plot(range(len(order)), [g.frac_ge8.get(m, np.nan) for m in order], "o-",
               color=OI[i + 1], label=f"{L:.0f} s", ms=4)
        a.plot(range(len(order)), [g.shift_frac_ge8.get(m, np.nan) for m in order], "--",
               color=OI[i + 1], lw=.8)
    a.set_xticks(range(len(order))); a.set_xticklabels(order, rotation=20)
    a.set_ylabel("fraction of Bostock events >= 8 MAD")
    a.set_title("(c) recall at an 8-MAD network threshold\n(dashed = per-station time-shifted control)")
    a.legend(title="template", fontsize=6)
    a = ax[1, 1]
    a.plot(ch.L, ch.e_1chan, "o-", color=OI[1], label="1 channel (Z)")
    a.plot(ch.L, ch.e_3chan, "s-", color=OI[2], label="3 channels (Z,N,E)")
    a.plot(ch.L, ch.e_network, "^-", color=OI[3], label="network (~9 channels)")
    a.set_xlabel("template length (s)"); a.set_ylabel("median excess (MAD)")
    a.set_title("(d) channel scaling: the LFE is on the horizontals")
    a.legend(fontsize=6)
    fig.suptitle("E4 step 1 — corrected Bostock-template test, 2010 outer ring "
                 "(PGC LZB VGZ SNB YOUB PFB), 15,596 catalogued LFEs", y=1.01, fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/e4_step1_registration.png"); plt.close(fig)


def f2_transfer():
    qc = pd.read_csv(f"{RES}/e4_transfer_stackqc_t2.csv")
    Q = pd.read_csv(f"{RES}/e4_step2_q.csv")
    N = pd.read_csv(f"{RES}/e4_step2_network.csv")
    S = pd.read_csv(f"{RES}/e4_step2_registration.csv")
    fig, ax = plt.subplots(2, 2, figsize=(9, 6.6))
    a = ax[0, 0]
    for i, s in enumerate(sorted(qc.sta.unique())):
        g = qc[qc.sta == s]
        a.scatter(g.n_events, g.snr_db, s=6, color=OI[i % 8], label=s, alpha=.7)
    a.axhline(0, color="k", lw=.7)
    a.set_xscale("log"); a.set_xlabel("events stacked (build half)")
    a.set_ylabel("stack S-SNR (dB)")
    a.set_title("(a) stack quality vs number of events\n(0 dB = the same-N random-time null)")
    a.legend(fontsize=5, ncol=2)
    a = ax[0, 1]
    g = qc[qc.has_seisbench]
    a.scatter(g.n_events, g.cc_vs_seisbench, s=6, color=OI[5], alpha=.6)
    a.axhline(0.6, color=OI[6], ls="--", lw=1, label="Merlin bar 0.6")
    a.set_xscale("log"); a.set_xlabel("events stacked (build half)")
    a.set_ylabel("cc(our 2010 re-stack, SeisBench stack)")
    a.set_title("(b) q-control: our stacking reproduces Bostock's template")
    a.legend(fontsize=6)
    a = ax[1, 0]
    b = S[(S.leg == "C_test") & (S.L == 8.0) & (S.magbin == "M1.8-2.1")]
    d = qc.groupby("sta").epi_km.median()
    for i, r in enumerate(b.itertuples()):
        c = OI[3] if r.sta in ("B011", "B926") else (OI[6] if r.sta in ("LRIV", "STW") else OI[2])
        a.scatter(d.get(r.sta, np.nan), r.excess3, s=40, color=c)
        a.annotate(r.sta, (d.get(r.sta, np.nan), r.excess3), fontsize=6,
                   xytext=(3, 3), textcoords="offset points")
    a.set_xlabel("median family epicentral distance (km)")
    a.set_ylabel("held-out single-station excess (MAD)")
    a.set_title("(c) transferred (green/orange) vs ring (blue), 8-s template")
    a = ax[1, 1]
    for i, leg in enumerate(["A_test", "B_test", "C_test"]):
        g = N[(N.leg == leg) & (N.magbin == "all")]
        lab = {"A_test": "ring, published stacks", "B_test": "ring, our 2010 re-stacks",
               "C_test": "ring + B011/B926/LRIV/STW"}[leg]
        a.plot(g.L, g.excess, "o-", color=OI[i + 1], label=lab, ms=4)
    a.set_xlabel("template length (s)"); a.set_ylabel("network excess (MAD), held out")
    a.set_title("(d) network gain: adding weak stations to an\nequal-weight sum does NOT help")
    a.legend(fontsize=6)
    fig.suptitle("E4 step 2 — q-control and template transfer, held out on alternate 2-h blocks",
                 y=1.01, fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/e4_step2_transfer.png"); plt.close(fig)


def f3_inject():
    G = pd.read_csv(f"{RES}/e4_inject_recall.csv")
    SN = pd.read_csv(f"{RES}/e4_inject_snr_incoherent.csv")
    real = pd.read_csv(f"{RES}/e3_inject_vs_real.csv")
    real = real[(real.kind == "real") & (real.threshold == 0.1)]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    a = ax[0]
    for i, (k, n) in enumerate([("seisbench", -1), ("own", 0), ("own", 20)]):
        g = G[(G.kind == k) & (G.nstack == n) & (G.threshold == 0.1)].sort_values("snr_rms")
        lab = "SeisBench family stack" if k == "seisbench" else \
            ("own stack, all events" if n == 0 else f"own stack, N={n}")
        a.plot(g.snr_rms, g.recall, "o-", color=OI[i + 1], label=lab, ms=4)
    # re-anchor each real magnitude bin by interpolating the MEASURED incoherent s/n curve
    sn = np.interp(real.mag_median, SN.mag_median, SN.snr_inc)
    a.plot(sn, real.recall, "s-", color=OI[6], ms=5,
           label="REAL Bostock events,\nre-anchored (measured incoherent s/n)")
    a.plot(real.snr_rms, real.recall, "s--", color=OI[0], ms=4, alpha=.5,
           label="REAL, E3 anchor (coherent 0.078)")
    a.set_xscale("log"); a.set_xlabel("s/n (rms over 6 s from S-1 s)")
    a.set_ylabel("S-pick recall within +-1 s (thr 0.1)")
    a.set_title("(a) injection vs the re-anchored real curve, PGC")
    a.legend(fontsize=5.5, loc="upper left")
    a = ax[1]
    a.errorbar(SN.mag_median, SN.snr_inc, yerr=SN.snr_inc_sem, fmt="o-", color=OI[6],
               ms=4, label="measured incoherent s/n")
    a.axhline(0.078, color=OI[0], ls="--", lw=1, label="D3 coherent anchor 0.078")
    iso = pd.read_csv(f"{RES}/e4_inject_snr_isolated.csv") if \
        os.path.exists(f"{RES}/e4_inject_snr_isolated.csv") else None
    if iso is not None:
        a.plot(iso.mag_med, iso.snr_inc, "^-", color=OI[3], ms=4,
               label="isolated events only (<=1 neighbour in +-30 s)")
    a.set_yscale("log"); a.set_xlabel("Bostock magnitude"); a.set_ylabel("s/n at PGC")
    a.set_title("(b) the x-axis anchor, measured (no stacking)")
    a.legend(fontsize=6)
    fig.suptitle("E4 step 3 — corrected injection and the re-anchored real recall curve",
                 y=1.04, fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/e4_step3_injection.png"); plt.close(fig)


def f4_null():
    T = pd.read_csv(f"{RES}/e4_assoc_null_table.csv")
    R = pd.read_csv(f"{RES}/e4_assoc_bostock_recall.csv")
    old = pd.read_csv(f"{RES}/e3_nullevents_rates.csv")
    old = old[(old.set == "vi2010") & (old.model == "cascadia")]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    a = ax[0]
    w = T[T["shift"] == "station"].reset_index(drop=True)   # `shift` shadows DataFrame.shift
    x = np.arange(len(w))
    a.bar(x - .2, w.null_over_obs, .4, color=OI[2], label="E4: tol 1.2 s, std 0.5 s, z 20-50 km")
    o = []
    for r in w.itertuples():
        m = old[(old.threshold == r.thr) & (old.n_min == r.n_stations)]
        o.append(float(m.ratio_all.iloc[0]) if len(m) else np.nan)
    a.bar(x + .2, o, .4, color=OI[1], label="E3: tol 2.5 s, std 1.5 s, z 0-70 km")
    a.axhline(0.10, color=OI[6], ls="--", lw=1, label="acceptance 0.10")
    a.set_xticks(x); a.set_xticklabels([f"thr {r.thr}\nN>={r.n_stations}" for r in w.itertuples()],
                                       fontsize=6)
    a.set_ylabel("shifted-pick null / observed")
    a.set_title("(a) chance fraction of the associated catalogue")
    a.legend(fontsize=6)
    for i, r in enumerate(w.itertuples()):
        a.annotate(f"{r.obs_per_day:.0f}/d", (x[i] - .2, r.null_over_obs),
                   ha="center", va="bottom", fontsize=5.5)
    a = ax[1]
    g = R[(R.tol_s == 2.0) & (R.thr == 0.3) & (R.n_stations == 4)]
    ordr = ["M<1.5", "M1.5-1.8", "M1.8-2.1", "M>=2.1"]
    g = g.set_index("magbin")
    a.bar(range(4), [g.recall.get(m, np.nan) for m in ordr], color=OI[2], label="recall")
    a.bar(range(4), [g.chance.get(m, np.nan) for m in ordr], color=OI[6], label="chance")
    a.set_xticks(range(4)); a.set_xticklabels(ordr, rotation=20)
    a.set_ylabel("fraction of Bostock events matched (+-2 s)")
    a.set_title("(b) Bostock recall, thr 0.3 / N>=4 (71 ev/day)")
    a.legend(fontsize=6)
    fig.suptitle("E4 step 4 — association re-null at tolerances matched to the measured residual",
                 y=1.04, fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/e4_step4_assoc_null.png"); plt.close(fig)


def f5_thomas():
    S = pd.read_csv(f"{RES}/e4_thomas_summary_bostock2012.csv")
    ev = pd.read_parquet(f"{RES}/e4_thomas_events_bostock2012_mcall.parquet")
    trk = pd.read_csv(f"{RES}/e4_thomas_fronttrack_bostock2012_mcall.csv")
    dec = pd.read_csv(f"{RES}/e4_thomas_decay_bostock2012_mcall.csv")
    fig, ax = plt.subplots(2, 2, figsize=(9.5, 6.6))
    a = ax[0, 0]
    t0 = ev.t.min()
    a.scatter((ev.t - t0) / 86400, ev.x_strike, s=.6, c=OI[2], alpha=.35, lw=0)
    a.plot(trk.day, trk.x_mean, "-", color=OI[6], lw=1.5, label="6-h activity centroid")
    a.plot(trk.day, trk.x_lead, "-", color=OI[1], lw=1, label="6-h leading edge (p10)")
    v = float(S[S.mc == "none"].front_speed_centroid_kmday.iloc[0])
    a.set_xlabel("days since 2012-08-30"); a.set_ylabel("along-strike (km, + = NW)")
    a.set_title(f"(a) Bostock 2012 SSE space-time; centroid speed {v:.1f} km/day\n"
                f"(Thomas: ~5 km/day southward = negative here)")
    a.legend(fontsize=6)
    a = ax[0, 1]
    a.plot(S.kept_frac * 100, -S.front_speed_centroid_kmday, "o-", color=OI[6],
           label="6-h centroid")
    a.plot(S.kept_frac * 100, -S.front_speed_lead_kmday, "s-", color=OI[1],
           label="6-h leading edge")
    a.plot(S.kept_frac * 100, -S.front_speed_kmday, "^-", color=OI[2],
           label="first-activation fit")
    a.axhspan(5 / 1.5, 5 * 1.5, color=OI[4], alpha=.35, label="Thomas ~5 km/day +-1.5x")
    a.set_xscale("log"); a.set_xlabel("% of Bostock's events kept (magnitude cut)")
    a.set_ylabel("southward front speed (km/day)")
    a.set_title("(b) DECIMATION: the main front survives to 3 % completeness")
    a.legend(fontsize=6)
    a = ax[1, 0]
    a.plot(S.kept_frac * 100, S.rerupture_frac_le5d, "o-", color=OI[3],
           label="fraction of bursts re-rupturing within 5 d")
    a.plot(S.kept_frac * 100, S.rerupture_interval_d_median * 5, "s-", color=OI[1],
           label="median inter-burst interval (d) x5")
    a.plot(S.kept_frac * 100, S.n_secondary / max(S.n_secondary.max(), 1), "^-", color=OI[6],
           label="secondary fronts (normalised)")
    a.plot(S.kept_frac * 100, S.n_secondary_null / max(S.n_secondary.max(), 1), "^--",
           color=OI[0], label="secondary-front NULL")
    a.set_xscale("log"); a.set_xlabel("% of Bostock's events kept")
    a.set_title("(c) re-rupture survives; secondary fronts never\nexceed the circular-shift null")
    a.legend(fontsize=5.5)
    a = ax[1, 1]
    for i, c in enumerate(["tide_R_third1_PLACEHOLDER", "tide_R_third2_PLACEHOLDER",
                           "tide_R_third3_PLACEHOLDER"]):
        a.plot(S.kept_frac * 100, S[c], "o-", color=OI[i + 1], label=f"third {i+1} of the SSE")
    a.set_xscale("log"); a.set_xlabel("% of Bostock's events kept")
    a.set_ylabel("Schuster R (mean resultant length)")
    a.set_title("(d) tidal sensitivity grows through the SSE\n*** PLACEHOLDER equilibrium tide ***")
    a.legend(fontsize=6)
    fig.suptitle("E4 step 5 — Thomas-2025 harness on Bostock's own Sept-2012 catalogue "
                 "(28,149 events, 130 families) + decimation", y=1.01, fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/e4_step5_thomas.png"); plt.close(fig)


def f6_pred():
    T = pd.read_csv(f"{RES}/e4_prediction_table.csv")
    T = T[T.roster == "all"]
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for i, (col, lab) in enumerate([("Mc_8MAD_median", "8-MAD threshold (chance-clean)"),
                                    ("Mc_4MAD_median", "4-MAD threshold (Bostock-equivalent)")]):
        a = ax[i]
        for j, (reg, g) in enumerate(T.groupby("region")):
            g = g.sort_values("year")
            a.plot(g.year, g[col], "o-", color=OI[j % 8], ms=3, label=reg)
        a.axhline(1.70, color="k", ls="--", lw=1, label="Bostock 2012 median M1.70")
        a.set_xlabel("year"); a.set_ylabel("predicted completeness magnitude")
        a.set_title(f"({'ab'[i]}) {lab}")
        if i == 1:
            a.legend(fontsize=5.5, ncol=2)
    fig.suptitle("E4 step 7 — predicted network-template-matching completeness per region-year",
                 y=1.04, fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/e4_step7_prediction.png"); plt.close(fig)


if __name__ == "__main__":
    for f in (f1_mad, f2_transfer, f3_inject, f4_null, f5_thomas, f6_pred):
        try:
            f(); print("[fig]", f.__name__, "ok")
        except Exception as e:
            print("[fig]", f.__name__, "FAILED:", type(e).__name__, e)
