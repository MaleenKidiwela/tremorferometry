#!/usr/bin/env python
"""E3 figures.  Colour-blind-safe (Okabe-Ito), labelled panels, explicit units.

  figures/e3_null_rates.png        step 1 -- observed vs shifted-pick null event rates
  figures/e3_relocation.png        step 2 -- configs A/B (S-P, real picks) and C (synthetic S-only)
  figures/e3_injection.png         step 3 -- injected recall(SNR) + the real curve + F vs threshold
  figures/e3_classifier.png        step 4 -- best-family cc, the null excess and the MAD scale
  figures/e3_family_offsets.png    step 5 -- the per-family Bostock time base
  figures/e3_denoise.png           step 3b -- denoised vs raw recall, F, and LFE fidelity
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
FIG = "lfe_catalog_repair/figures"
# Okabe-Ito
OB = "#0072B2"; OR = "#D55E00"; OG = "#009E73"; OY = "#E69F00"
OP = "#CC79A7"; OS = "#56B4E9"; OK = "#000000"; OGY = "#999999"
plt.rcParams.update({"figure.dpi": 130, "font.size": 8.5, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.axisbelow": True})


def panel(ax, letter, title):
    ax.set_title(f"({letter}) {title}", loc="left", fontsize=9, fontweight="bold")


def fig_null():
    R = pd.read_csv(f"{RES}/e3_nullevents_rates.csv")
    try:
        H = pd.read_csv(f"{RES}/e3_nullhr_rates.csv")
    except FileNotFoundError:
        H = None
    f, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
    g = R[(R.set == "vi2010") & (R.model == "cascadia")]
    w = 0.26
    for k, (n, gg) in enumerate(g.groupby("n_min")):
        x = np.arange(len(gg))
        ax[0].bar(x + (k - 1) * w, gg.obs_per_day, w, color=[OB, OG, OY][k],
                  label=f"observed N$\\geq${n}")
        ax[0].errorbar(x + (k - 1) * w, gg.null_all_per_day, yerr=gg.null_all_sd, fmt="_",
                       ms=11, mew=2, color=OR, ecolor=OK, capsize=2,
                       label="shifted-pick null" if k == 0 else None)
        ax[0].set_xticks(np.arange(len(gg)))
        ax[0].set_xticklabels([f"{t:g}" for t in gg.threshold])
    ax[0].set_yscale("log"); ax[0].set_xlabel("LFEDetect pick threshold")
    ax[0].set_ylabel("events per day")
    ax[0].legend(fontsize=7, ncol=2)
    panel(ax[0], "a", "vi2010 cascadia: observed vs chance")

    for i, (n, gg) in enumerate(g.groupby("n_min")):
        ax[1].plot(gg.threshold, gg.ratio_all, "o-", color=[OB, OG, OY][i], label=f"N$\\geq${n}")
        if H is not None:
            hh = H[(H.set == "vi2010") & (H.model == "cascadia") & (H.n_min == n)]
            ax[1].plot(hh.threshold, hh.ratio_all, "s--", ms=4, color=[OB, OG, OY][i], alpha=0.6)
    ax[1].axhspan(0.5, 1.05, color=OR, alpha=0.12)
    ax[1].axhspan(0.0, 0.15, color=OG, alpha=0.12)
    ax[1].text(0.21, 0.72, "chance-dominated (>50 %)", fontsize=7, color=OR)
    ax[1].text(0.21, 0.055, "activity PASS (<15 %)", fontsize=7, color=OG)
    ax[1].set_ylim(0, 1.05); ax[1].set_xlabel("LFEDetect pick threshold")
    ax[1].set_ylabel("null / observed event rate")
    ax[1].legend(fontsize=7, title="solid: per-station shift\ndashed: per-station-per-hour",
                 title_fontsize=6)
    panel(ax[1], "b", "chance fraction of the catalogue")

    v = R[(R.set == "vi2017") & (R.model == "cascadia")]
    x = np.arange(len(v)); w = 0.2
    ax[2].bar(x - 1.5 * w, v.obs_ets_per_day, w, color=OB, label="observed, ETS days")
    ax[2].bar(x - 0.5 * w, v.null_ets_per_day, w, color=OS, label="null, ETS days")
    ax[2].bar(x + 0.5 * w, v.obs_quiet_per_day, w, color=OR, label="observed, quiet days")
    ax[2].bar(x + 1.5 * w, v.null_quiet_per_day, w, color=OY, label="null, quiet days")
    for i, r in enumerate(v.itertuples()):
        ax[2].text(i, max(r.obs_ets_per_day, 1) * 1.25,
                   f"obs {r.obs_ets_quiet:.1f}$\\times$\nnull {r.null_ets_quiet:.1f}$\\times$",
                   ha="center", fontsize=6.5)
    ax[2].set_yscale("log"); ax[2].set_xticks(x)
    ax[2].set_xticklabels([f"N$\\geq${int(n)}" for n in v.n_min])
    ax[2].set_ylabel("events per day"); ax[2].legend(fontsize=6.5)
    ax[2].set_ylim(0.05, 900)
    panel(ax[2], "c", "vi2017: chance reproduces the ETS/quiet contrast")
    f.suptitle("E3 step 1 — chance-association null (5 shifted-pick seeds; error bars = 1 s.d.)",
               fontsize=10)
    f.tight_layout(rect=(0, 0, 1, 0.94))
    f.savefig(f"{FIG}/e3_null_rates.png"); plt.close(f)
    print("[fig] e3_null_rates.png")


def fig_reloc():
    A = pd.read_csv(f"{RES}/e3_family_reloc.csv")
    frames = []
    for t, lab in (("s005", "0.05 s"), ("s030", "0.30 s"), ("s050", "0.50 s"),
                   ("sheavy", "heavy tail")):
        p = f"{RES}/e3_geom_synth_{t}.csv"
        if os.path.exists(p):
            d = pd.read_csv(p); d["lab"] = lab; frames.append(d)
    C = pd.concat(frames, ignore_index=True) if frames else None
    f, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
    parts = [A[A.config == "A"].dh_km.dropna(), A[A.config == "B"].dh_km.dropna()]
    labs = ["A: S$-$P, all stns", "B: S$-$P, 6 stns"]
    if C is not None:
        for n in (3, 4, 5, 6):
            s = C[(C.lab == "0.30 s") & (C.n_sta == n)].dh_km.dropna()
            if len(s):
                parts.append(s); labs.append(f"C: S-only, N={n}\n($\\sigma$=0.3 s)")
    bp = ax[0].boxplot(parts, showfliers=False, widths=0.6, patch_artist=True)
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor(OG if i < 2 else OB); b.set_alpha(0.65)
    for m in bp["medians"]:
        m.set_color(OK)
    ax[0].set_xticklabels(labs, fontsize=6.5, rotation=20, ha="right")
    ax[0].set_ylabel("horizontal offset from the published\nfamily hypocentre (km)")
    ax[0].axhline(5, color=OR, ls="--", lw=1)
    ax[0].text(0.4, 5.6, "audit pass bar 5 km", color=OR, fontsize=6.5)
    panel(ax[0], "a", "model+method (A,B) vs geometry (C)")
    ax[0].set_ylim(0, 70)

    parts, labs = [A[A.config == "A"].dz_km.dropna(), A[A.config == "B"].dz_km.dropna()], labs[:2]
    if C is not None:
        for n in (3, 4, 5, 6):
            s = C[(C.lab == "0.30 s") & (C.n_sta == n)].dz_km.dropna()
            if len(s):
                parts.append(s); labs.append(f"C: N={n}")
    bp = ax[1].boxplot(parts, showfliers=False, widths=0.6, patch_artist=True)
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor(OG if i < 2 else OB); b.set_alpha(0.65)
    for m in bp["medians"]:
        m.set_color(OK)
    ax[1].axhline(0, color=OK, lw=0.8)
    ax[1].set_xticklabels(labs, fontsize=6.5, rotation=20, ha="right")
    ax[1].set_ylabel("depth offset (km; negative = too shallow)")
    panel(ax[1], "b", "depth: S-only collapses, S$-$P does not")

    if C is not None:
        for i, (lab, g) in enumerate(C.groupby("lab")):
            s = g.groupby("n_sta").dh_km.median()
            ax[2].plot(s.index, s.values, "o-", color=[OB, OG, OY, OP][i % 4], label=lab)
        e2 = {3: 33.0, 4: 25.0, 5: 15.5, 6: 11.0}
        ax[2].plot(list(e2), list(e2.values()), "k^--", lw=1.6,
                   label="E2 associated events\n(Bostock time-matched)")
        for thr, mk in ((0.1, "*"), (0.3, "P")):
            fp = f"{RES}/e3_realpick_summary_t{thr}.csv"
            if os.path.exists(fp):
                D = pd.read_csv(fp)
                D = D[D.n >= 5]
                if len(D):
                    ax[2].plot(D.n_sta, D.dh_med, mk, ms=12, color=OR, ls="none",
                               label=f"D: REAL picks, known\nevents, no PyOcto (thr {thr})")
        ax[2].set_xlabel("stations used"); ax[2].set_ylabel("median horizontal offset (km)")
        ax[2].legend(fontsize=6.5, title="synthetic S-pick noise", title_fontsize=6)
        ax[2].set_xticks([3, 4, 5, 6])
    panel(ax[2], "c", "geometry explains N=3; N$\\geq$4 needs association error")
    f.suptitle("E3 step 2 — family relocation control (129 Bostock families)", fontsize=10)
    f.tight_layout(rect=(0, 0, 1, 0.94))
    f.savefig(f"{FIG}/e3_relocation.png"); plt.close(f)
    print("[fig] e3_relocation.png")


def fig_inject():
    V = pd.read_csv(f"{RES}/e3_inject_vs_real.csv")
    F = pd.read_csv(f"{RES}/e3_inject_falsealarm.csv")
    B = pd.read_csv(f"{RES}/e3_real_recall_by_burst.csv")
    f, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
    for i, (m, c) in enumerate((("cascadia", OB), ("jcms", OG))):
        for stn, ls in (("PGC", "-"), ("LZB", "--"), ("YOUB", ":")):
            g = V[(V.kind == "injected") & (V.model == m) & (V.station == stn)
                  & (V.threshold == 0.1)].sort_values("snr_rms")
            if len(g):
                ax[0].plot(g.snr_rms, g.recall, ls, marker="o", ms=3, color=c,
                           label=f"{m} @ {stn}")
    r = V[(V.kind == "real") & (V.threshold == 0.1)].sort_values("snr_rms")
    ax[0].plot(r.snr_rms, r.recall, "s-", color=OR, lw=2, label="REAL, PGC vs Bostock\n(by magnitude)")
    ax[0].axvline(0.078, color=OK, ls=":", lw=1)
    ax[0].axvspan(0.032, 0.078, color=OK, alpha=0.07)
    ax[0].text(0.081, 0.58, "D3 s/n anchor\n0.078 (0.032 lo)", fontsize=6.5)
    ax[0].set_xscale("log"); ax[0].set_xlabel("s/n  (rms of signal over 6 s / rms of 1-8 Hz noise)")
    ax[0].set_ylabel("recall within $\\pm$1 s (chance-corrected)")
    ax[0].legend(fontsize=6, ncol=2)
    panel(ax[0], "a", "injected vs real recall, pick threshold 0.1")

    for i, (lab, g) in enumerate(B.groupby("group")):
        ax[1].plot(g.snr_rms, g.recall, "o-", color=[OB, OG, OY][i % 3], label=lab)
    ax[1].set_xscale("log"); ax[1].set_xlabel("s/n (same conversion)")
    ax[1].set_ylabel("real recall within $\\pm$1 s")
    ax[1].legend(fontsize=6.5, title="Bostock events within $\\pm$30 s", title_fontsize=6)
    ax[1].set_ylim(0, 0.35)
    panel(ax[1], "b", "burst context does NOT explain the gap")

    fb = F[F.run == "base"].groupby(["station", "model", "threshold"]).n_S.mean().reset_index()
    for i, (m, c) in enumerate((("cascadia", OB), ("jcms", OG))):
        for stn, ls in (("PGC", "-"), ("LZB", "--"), ("YOUB", ":")):
            g = fb[(fb.model == m) & (fb.station == stn)].sort_values("threshold")
            ax[2].plot(g.threshold, g.n_S, ls, marker="o", ms=3, color=c, label=f"{m} @ {stn}")
    S = pd.read_csv(f"{RES}/e3_F_scaling.csv")
    s4 = S[S.n_min == 4]
    ax2 = ax[2].twinx()
    ax2.plot(s4.threshold, s4.F_reduction_needed, "d-", color=OR, lw=2)
    ax2.set_ylabel("F reduction needed for chance $\\leq$10 % of real\n(N$\\geq$4)", color=OR)
    ax2.tick_params(axis="y", colors=OR); ax2.grid(False)
    ax2.axhline(1, color=OR, ls=":", lw=1)
    ax[2].set_yscale("log"); ax[2].set_xlabel("LFEDetect pick threshold")
    ax[2].set_ylabel("false alarms F (S picks per quiet station-day)")
    ax[2].legend(fontsize=6, ncol=2, loc="upper right")
    panel(ax[2], "c", "false-alarm rate and how far it must fall")
    f.suptitle("E3 step 3 — injection recall vs SNR and the false-alarm rate on real noise",
               fontsize=10)
    f.tight_layout(rect=(0, 0, 1, 0.94))
    f.savefig(f"{FIG}/e3_injection.png"); plt.close(f)
    print("[fig] e3_injection.png")


def fig_classifier():
    p = f"{RES}/e3_classifier_events.parquet"
    if not os.path.exists(p):
        print("[fig] classifier data missing -- skipped"); return
    E = pd.read_parquet(p)
    # keep_default_na=False: the set NAMED "null" is otherwise read back as NaN by pandas
    S = pd.read_csv(f"{RES}/e3_classifier_summary.csv", keep_default_na=False)
    for c in S.columns:
        if c != "set":
            S[c] = pd.to_numeric(S[c], errors="coerce")
    med = float(S.random_median.iloc[0]); mad = float(S.random_mad.iloc[0])
    p99 = float(S.null_p99.iloc[0])
    f, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
    cols = {"random": OGY, "null": OR, "assoc": OB, "bostock": OG}
    bins = np.linspace(min(E.best_cc.quantile(0.001), 0), E.best_cc.quantile(0.9995), 70)
    for k, c in cols.items():
        v = E[E.set == k].best_cc.dropna()
        if len(v):
            ax[0].hist(v, bins=bins, density=True, histtype="step", lw=1.6, color=c,
                       label=f"{k} (n={len(v)})")
    ax[0].axvline(p99, color=OK, ls="--", lw=1)
    ax[0].text(p99, ax[0].get_ylim()[1] * 0.85, " null 99th pct", fontsize=6.5)
    ax[0].set_xlabel("best-family network cc (6 stations, 3-C, $\\pm$3 s)")
    ax[0].set_ylabel("density"); ax[0].legend(fontsize=6.5)
    panel(ax[0], "a", "best-family cc by event set")

    x = np.arange(len(S));
    ax[1].bar(x, S.frac_above_null_p99, color=[cols.get(k, OP) for k in S.set])
    ax[1].axhline(0.01, color=OK, ls="--", lw=1)
    ax[1].text(0.1, 0.013, "1 % = the null by construction", fontsize=6.5)
    ax[1].set_xticks(x); ax[1].set_xticklabels(S.set, rotation=25, ha="right", fontsize=6.5)
    ax[1].set_ylabel("fraction with best cc > null 99th pct")
    panel(ax[1], "b", "excess over the chance-event null")

    # panel (c) uses the PER-TEMPLATE MAD (each family's own random-time background), not the
    # best-of-129 MAD -- the latter is a maximum over 129 draws and is not the scale D3 quoted.
    M = pd.read_csv(f"{RES}/e3_classifier_mad.csv", keep_default_na=False)
    for c in M.columns:
        if c != "set":
            M[c] = pd.to_numeric(M[c], errors="coerce")
    Z = np.load(f"{RES}/e3_classifier_perfam.npz")
    Rm = Z["cc_random"]
    fmed = np.nanmedian(Rm, axis=0)
    fmad = np.nanmedian(np.abs(Rm - fmed[None, :]), axis=0)
    for k, c, lab in (("random", OGY, "random times, best family"),
                      ("null", OR, "chance events, best family"),
                      ("assoc", OB, "associated events, best family")):
        A = Z.get(f"cc_{k}")
        if A is None:
            continue
        bz = np.nanmax((A - fmed[None, :]) / np.clip(fmad[None, :], 1e-9, None), axis=1)
        ax[2].hist(bz, bins=np.linspace(-2, 16, 70), density=True, histtype="step", lw=1.5,
                   color=c, label=lab)
    mb = M[M.set == "bostock_own_template"]
    ob = E[(E.set == "bostock") & E.own_cc.notna()]
    if len(mb):
        # rebuild the own-template z from the saved matrix ordering is not needed: use the CSV
        ax[2].axvline(float(mb.best_mad_median.iloc[0]), color=OG, lw=2.5)
        ax[2].text(float(mb.best_mad_median.iloc[0]) + 0.3, 0.18,
                   f"Bostock's OWN events vs\ntheir OWN template:\n"
                   f"median {float(mb.best_mad_median.iloc[0]):.1f} MAD, "
                   f"p90 {float(mb.best_mad_p90.iloc[0]):.1f}",
                   fontsize=7, color=OG)
    for m, lab, ls in ((3.5, "D3's 3.5 MAD", ":"), (8.0, "8 MAD", "--")):
        ax[2].axvline(m, color=OK, ls=ls, lw=1)
        ax[2].text(m + 0.1, 0.02, lab, fontsize=6.5, rotation=90)
    ax[2].set_xlabel("network cc, in MAD of that template's own random-time background")
    ax[2].set_ylabel("density"); ax[2].legend(fontsize=6.5, loc="upper right")
    panel(ax[2], "c", "Bostock's own events do not stand out here")
    f.suptitle("E3 step 4 — Bostock family templates as a classifier", fontsize=10)
    f.tight_layout(rect=(0, 0, 1, 0.94))
    f.savefig(f"{FIG}/e3_classifier.png"); plt.close(f)
    print("[fig] e3_classifier.png")


def fig_offsets():
    F = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    FS = pd.read_csv(f"{RES}/e3_bostock_famsta_offsets.csv")
    acc = F[F.accept.astype(bool)]
    f, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
    ax[0].hist(acc.offset_s, bins=np.arange(0, 22, 0.5), color=OB, alpha=0.85)
    ax[0].axvline(12.46, color=OR, lw=2, ls="--", label="E2's single constant +12.46 s")
    ax[0].axvline(acc.offset_s.median(), color=OK, lw=1.5,
                  label=f"median {acc.offset_s.median():+.2f} s")
    ax[0].set_xlabel("per-family offset: equivalent origin $-$ DB time (s)")
    ax[0].set_ylabel("families"); ax[0].legend(fontsize=6.5)
    panel(ax[0], "a", f"{len(acc)} of {len(F)} families, s.d. {acc.offset_s.std():.2f} s")

    ax[1].scatter(acc.epi_ref, acc.offset_s, s=14, c=acc.dep, cmap="viridis")
    cb = plt.colorbar(ax[1].collections[0], ax=ax[1]); cb.set_label("family depth (km)", fontsize=7)
    ax[1].set_xlabel("epicentral distance to PGC (km)")
    ax[1].set_ylabel("offset (s)")
    panel(ax[1], "b", "weak distance trend (Spearman $-$0.39)")

    fsa = FS[(FS.bkey >= 0) & (FS.accept.astype("boolean").fillna(False))]
    g = fsa[fsa.station != "PGC"].groupby("bkey").offset_s.median()
    j = acc.set_index("bkey").join(g.rename("other"), how="inner").dropna(subset=["other"])
    ax[2].plot([0, 22], [0, 22], color=OGY, lw=1)
    ax[2].scatter(j.offset_s, j.other, s=14, color=OG)
    r = np.corrcoef(j.offset_s, j.other)[0, 1]
    ax[2].set_xlabel("offset at PGC, 2010-05..11 (s)")
    ax[2].set_ylabel("median offset at $\\geq$3 other stations,\n2010-08..09 (s)")
    ax[2].text(0.05, 0.92, f"r = {r:.3f}, n = {len(j)}\nresidual s.d. "
               f"{(j.offset_s-j.other).std():.2f} s", transform=ax[2].transAxes, fontsize=7)
    panel(ax[2], "c", "cross-validated on independent stations and window")
    f.suptitle("E3 step 5 — the Bostock time base is a per-FAMILY property", fontsize=10)
    f.tight_layout(rect=(0, 0, 1, 0.94))
    f.savefig(f"{FIG}/e3_family_offsets.png"); plt.close(f)
    print("[fig] e3_family_offsets.png")


def fig_denoise():
    p = f"{RES}/e3_denoise_recall.csv"
    if not os.path.exists(p):
        print("[fig] denoise data missing -- skipped"); return
    R = pd.read_csv(p); F = pd.read_csv(f"{RES}/e3_denoise_falsealarm.csv")
    FI = pd.read_csv(f"{RES}/e3_denoise_fidelity.csv")
    RD = pd.read_csv(f"{RES}/e3_denoise_realdays.csv")
    cols = {"raw": OB, "dn": OR, "wiener": OG}
    f, ax = plt.subplots(1, 4, figsize=(14, 3.5))
    q = R[(R.tol == 1.0) & (R.threshold == 0.3) & (R.model == "cascadia")]
    for arm, c in cols.items():
        for stn, ls in (("PGC", "-"), ("LZB", "--"), ("YOUB", ":")):
            g = q[(q.arm == arm) & (q.station == stn)].groupby("level").recall.mean()
            if len(g):
                ax[0].plot(g.index, g.values, ls, marker="o", ms=3.5, color=c,
                           label=f"{arm} @ {stn}")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("injection level (peak amplitude / host 1-8 Hz rms)")
    ax[0].set_ylabel("injected recall $\\pm$1 s (chance-corrected)")
    ax[0].legend(fontsize=5.5, ncol=3); panel(ax[0], "a", "injected recall (cascadia, thr 0.3)")

    for arm, c in cols.items():
        g = F[(F.arm == arm) & (F.model == "cascadia")].groupby("threshold").n_S_base.mean()
        ax[1].plot(g.index, g.values, "o-", color=c, label=arm)
    ax[1].set_yscale("log"); ax[1].set_xlabel("pick threshold")
    ax[1].set_ylabel("F (S picks per quiet station-day)")
    ax[1].legend(fontsize=7); panel(ax[1], "b", "false-alarm rate")

    for arm, c in cols.items():
        g = FI[FI.arm == arm].groupby("level").amp_ratio.mean()
        ax[2].plot(g.index, g.values, "o-", color=c, label=f"{arm}: amplitude kept")
        g2 = FI[FI.arm == arm].groupby("level")["corr"].mean()
        ax[2].plot(g2.index, g2.values, "s--", ms=3, color=c, alpha=0.6)
    ax[2].axhline(1, color=OK, lw=0.8)
    ax[2].set_xscale("log"); ax[2].set_xlabel("injection level (peak / host rms)")
    ax[2].set_ylabel("solid: rms kept;  dashed: waveform corr.")
    ax[2].legend(fontsize=6.5); panel(ax[2], "c", "what each arm does to the injected LFE")

    if len(RD):
        w = RD[(RD.threshold == 0.1)].groupby(["gt", "arm"]).agg(
            excess=("excess", "sum"), n_gt=("n_gt", "sum")).reset_index()
        w["recall"] = w.excess / w.n_gt
        # `gt` shadows DataFrame.gt, so bracket access is required here
        gts = sorted(w["gt"].unique())
        x = np.arange(len(gts))
        bw = 0.26
        for i, (arm, c) in enumerate(cols.items()):
            v = [float(w[(w["gt"] == g) & (w.arm == arm)].recall.iloc[0])
                 if len(w[(w["gt"] == g) & (w.arm == arm)]) else 0 for g in gts]
            ax[3].bar(x + (i - 1) * bw, v, bw, color=c, label=arm)
        ax[3].set_xticks(x); ax[3].set_xticklabels(gts)
        ax[3].set_ylabel("$\\pm$1 s peak-excess recall, thr 0.1")
        ax[3].legend(fontsize=7)
        panel(ax[3], "d", "REAL PGC 2010 ground-truth days")
    f.suptitle("E3 step 3b — does denoising before picking help?  DeepDenoiser and a "
               "spectral-Wiener baseline vs raw, all arms at 100 Hz", fontsize=10)
    f.tight_layout(rect=(0, 0, 1, 0.93))
    f.savefig(f"{FIG}/e3_denoise.png"); plt.close(f)
    print("[fig] e3_denoise.png")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    for fn in (fig_null, fig_offsets, fig_inject, fig_reloc, fig_classifier, fig_denoise):
        try:
            fn()
        except Exception as e:
            import traceback
            print(f"[fig] {fn.__name__} FAILED: {e}")
            traceback.print_exc()
