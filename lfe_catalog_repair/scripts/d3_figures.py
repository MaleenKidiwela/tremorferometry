#!/usr/bin/env python
"""D3 stage 3d — figures/d3_pilot_validation.png (the before/after panel set)."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

R = Path("lfe_catalog_repair/results")
F = Path("lfe_catalog_repair/figures")
PNSN = "catalogs/pnsn_tremor_cascadia_full.csv"
KS = (8, 9, 10)
COL = {"old": "0.45", "new8": "#1b7f4b", "new9": "#c8531b", "new10": "#7b3ea8"}
LBL = {"old": "OLD cc>=0.80", "new8": "NEW MAD 8", "new9": "NEW MAD 9", "new10": "NEW MAD 10"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", required=True, help="comma-separated station tags")
    ap.add_argument("--examples", default="", help="station:template,station:template,... for panel D")
    ap.add_argument("--coinc", default="", help="coincidence csv (optional)")
    ap.add_argument("--lin", default=str(R / "d3_lin_groundtruth.csv"))
    ap.add_argument("--out", default=str(F / "d3_pilot_validation.png"))
    a = ap.parse_args()
    stas = [s.strip() for s in a.stations.split(",") if s.strip()]
    val = {s: pd.read_csv(R / f"d3_validation_{s}.csv") for s in stas}
    F.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(17, 13.5))
    gs = fig.add_gridspec(4, 3, hspace=0.42, wspace=0.28)

    # ---- A: ETS/quiet ratio, old vs new, per template (ECDF) ----
    ax = fig.add_subplot(gs[0, 0])
    for rule in ["old"] + [f"new{k}" for k in KS]:
        v = np.concatenate([val[s][f"ratio_{rule}"].values for s in stas])
        v = v[np.isfinite(v) & (v > 0)]
        if not len(v):
            continue
        x = np.sort(v)
        ax.plot(x, 1 - np.arange(len(x)) / len(x), color=COL[rule], lw=1.8, label=f"{LBL[rule]} (n={len(v)})")
    ax.axvline(5, color="k", ls=":", lw=1)
    ax.text(5.2, 0.9, "pass bar\nETS/quiet=5", fontsize=7)
    ax.set_xscale("log")
    ax.set_xlabel("ETS / quiet detection-rate ratio (per template)")
    ax.set_ylabel("fraction of templates above")
    ax.set_title("A  ETS responsiveness: old vs new rule", fontsize=10, loc="left")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # ---- B: per-station bar of frac(ratio>=5) and frac(ratio>=2) ----
    ax = fig.add_subplot(gs[0, 1])
    w = 0.2
    rules = ["old"] + [f"new{k}" for k in KS]
    xs = np.arange(len(stas))
    for i, rule in enumerate(rules):
        f5 = [float((val[s][f"ratio_{rule}"] >= 5).mean()) for s in stas]
        ax.bar(xs + (i - 1.5) * w, f5, w, color=COL[rule], label=LBL[rule])
    ax.set_xticks(xs)
    ax.set_xticklabels(stas, fontsize=8)
    ax.set_ylabel("fraction of templates with ETS/quiet >= 5")
    ax.set_title("B  fraction passing the ETS/quiet >= 5 bar", fontsize=10, loc="left")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, axis="y")

    # ---- C: Spearman(daily count, local tremor) distribution ----
    ax = fig.add_subplot(gs[0, 2])
    for rule in rules:
        v = np.concatenate([val[s][f"rho_{rule}"].values for s in stas])
        v = v[np.isfinite(v)]
        if not len(v):
            continue
        ax.hist(v, bins=np.linspace(-1, 1, 41), histtype="step", color=COL[rule], lw=1.6,
                density=True, label=f"{LBL[rule]} med {np.median(v):+.3f}")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel(r"Spearman $\rho$(daily detections, local 30-d tremor)")
    ax.set_ylabel("density")
    ax.set_title("C  tremor correlation per template", fontsize=10, loc="left")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # ---- D-F: daily counts vs tremor for 3 example families ----
    ex = []
    if a.examples:
        for tok in a.examples.split(","):
            s, t = tok.split(":", 1)
            ex.append((s.strip(), t.strip()))
    else:
        for s in stas:
            v = val[s].copy()
            v = v[np.isfinite(v["ratio_new9"])]
            if len(v):
                ex.append((s, v.sort_values("ratio_new9", ascending=False).template.iloc[0]))
        ex = ex[:3]
    tre = pd.read_csv(PNSN, usecols=["time", "lat", "lon"], parse_dates=["time"])
    tre["day"] = tre.time.dt.floor("D")
    for j, (s, tmpl) in enumerate(ex[:3]):
        ax = fig.add_subplot(gs[1, j])
        det = pd.read_parquet(R / f"d3_detections_{s}.parquet",
                              columns=["template", "t", "old", "new9"])
        det = det[det["template"].astype(str) == tmpl]
        if not len(det):
            ax.set_title(f"D{j+1}  {s} {tmpl}: no detections", fontsize=9, loc="left")
            continue
        det["day"] = pd.to_datetime(det["t"], unit="s").dt.floor("D")
        lat, lon = [float(x) for x in tmpl.split("__")[0].split("_")]
        d = np.hypot((tre.lat - lat) * 111.2, (tre.lon - lon) * 111.2 * np.cos(np.radians(lat)))
        ts = tre[d <= 40].groupby("day").size()
        days = pd.date_range(det.day.min(), det.day.max(), freq="D")
        roll = ts.reindex(days, fill_value=0).rolling(30, min_periods=1).sum()
        ax2 = ax.twinx()
        ax2.fill_between(days, roll.values, color="0.85", zorder=0)
        ax2.set_ylabel("PNSN tremor, 30-d", color="0.5", fontsize=8)
        ax2.tick_params(labelsize=7, colors="0.5")
        for rule in ["old", "new9"]:
            c = det[det[rule]].groupby("day").size().reindex(days, fill_value=0)
            rho = spearmanr(c.values, roll.values)[0]
            ax.plot(days, c.values, color=COL[rule], lw=0.9,
                    label=f"{LBL[rule]} rho={rho:+.2f}", zorder=2)
        ax.set_ylabel("detections / day", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(f"D{j+1}  {s}  {tmpl}", fontsize=9, loc="left")
        ax.legend(fontsize=6.5, loc="upper left")

    # ---- G: Lin recall/precision vs threshold (POOLED -- per template the rates are too low) ----
    ax = fig.add_subplot(gs[2, 0])
    lf = Path(str(a.lin).replace(".csv", "_pooled.csv"))
    if lf.exists():
        L = pd.read_csv(lf)
        order = [r for r in ["old", "new8", "new9", "new10"] if r in set(L.rule)]
        x = np.arange(len(order))
        for key, c, mk, lab in [("recall", "#1b6ca8", "o", "recall"),
                                ("precision", "#c0392b", "s", "precision"),
                                ("precision_shift", "0.5", "^", "precision, +600 s control")]:
            v = [float(L.loc[L.rule == r, key].iloc[0]) for r in order]
            ax.plot(x, v, mk + "-", color=c, label=lab)
        ax.set_xticks(x)
        ax.set_xticklabels([LBL[r] for r in order], fontsize=7, rotation=20)
        ax.set_yscale("log")
        ax.set_ylabel("rate vs Lin/Bostock", fontsize=8)
        ax.legend(fontsize=7)
        ax.text(0.03, 0.06, "precision sits ON the shifted\ncontrol at every rule = chance",
                transform=ax.transAxes, fontsize=7.5, color="#c0392b")
    else:
        ax.text(.5, .5, "no pooled Lin file", ha="center")
    ax.set_title("G  Lin/Bostock recall & precision (pooled)", fontsize=10, loc="left")
    ax.grid(alpha=0.3)

    # ---- H: coincidence bars ----
    ax = fig.add_subplot(gs[2, 1])
    if a.coinc and Path(a.coinc).exists():
        C = pd.read_csv(a.coinc)
        order = ["old", "new8", "new9", "new10"]
        xs = np.arange(len(order))
        for i, key in enumerate(["frac_ets", "frac_quiet", "frac_shift"]):
            v = [C.loc[C.rule == r, key].mean() for r in order]
            ax.bar(xs + (i - 1) * 0.26, v, 0.26,
                   label={"frac_ets": "ETS", "frac_quiet": "quiet",
                          "frac_shift": "+600 s control"}[key],
                   color=["#1b7f4b", "#4a90c2", "0.6"][i])
        ax.set_xticks(xs)
        ax.set_xticklabels([LBL[r] for r in order], fontsize=7, rotation=20)
        ax.set_yscale("log")
        ax.set_ylabel("coincident fraction (+-3 s)", fontsize=8)
        ax.legend(fontsize=7)
    else:
        ax.text(.5, .5, "no coincidence file", ha="center")
    ax.set_title("H  network coincidence (surface x borehole)", fontsize=10, loc="left")
    ax.grid(alpha=0.3, axis="y")

    # ---- I: detection rate per template-day, old vs new ----
    ax = fig.add_subplot(gs[2, 2])
    summ = pd.read_csv(R / "d3_station_summary.csv")
    xs = np.arange(len(stas))
    for i, rule in enumerate(rules):
        v = [float(summ[(summ.station == s) & (summ.rule == rule)].det_per_tmpl_day.iloc[0])
             if len(summ[(summ.station == s) & (summ.rule == rule)]) else np.nan for s in stas]
        ax.bar(xs + (i - 1.5) * w, v, w, color=COL[rule], label=LBL[rule])
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(stas, fontsize=8)
    ax.set_ylabel("detections / template / day")
    ax.set_title("I  detection rate (the false-alarm load)", fontsize=10, loc="left")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, axis="y")

    # ---- J: threshold distributions / attainability ----
    ax = fig.add_subplot(gs[3, 0])
    for s in stas:
        f = R / f"d3_thresholds_{s}.csv"
        if not f.exists():
            continue
        T = pd.read_csv(f)
        for k, ls in zip(KS, ["-", "--", ":"]):
            v = T[f"thr{k}"].dropna().values
            if not len(v):
                continue
            x = np.sort(v)
            ax.plot(x, np.arange(len(x)) / len(x), ls, lw=1.3, label=f"{s} k={k}")
    ax.axvline(1.0, color="r", lw=1)
    ax.axvline(0.80, color="k", ls=":", lw=1)
    ax.text(0.805, 0.05, "old 0.80", fontsize=7, rotation=90)
    ax.text(1.005, 0.05, "unattainable", fontsize=7, rotation=90, color="r")
    ax.set_xlabel("threshold = median + k*MAD (per template x 30-d window)")
    ax.set_ylabel("cumulative fraction")
    ax.set_title("J  where the MAD thresholds land", fontsize=10, loc="left")
    ax.legend(fontsize=6, ncol=2)
    ax.grid(alpha=0.3)

    # ---- K: SNR distributions ----
    ax = fig.add_subplot(gs[3, 1])
    for s in stas:
        det = pd.read_parquet(R / f"d3_detections_{s}.parquet", columns=["snr", "old", "new9"])
        for rule, ls in [("old", "--"), ("new9", "-")]:
            v = det.loc[det[rule], "snr"].values
            v = v[np.isfinite(v) & (v > 0)]
            if len(v) < 10:
                continue
            ax.hist(v, bins=np.logspace(0, 2, 40), histtype="step", ls=ls, lw=1.4,
                    density=True, label=f"{s} {LBL[rule]} (med {np.median(v):.1f})")
    ax.set_xscale("log")
    ax.set_xlabel("detection SNR (peak |Z| in window / pre-window RMS)")
    ax.set_ylabel("density")
    ax.set_title("K  detection SNR", fontsize=10, loc="left")
    ax.legend(fontsize=6)
    ax.grid(alpha=0.3)

    # ---- L: cc_stack of new vs old detections ----
    ax = fig.add_subplot(gs[3, 2])
    for s in stas:
        det = pd.read_parquet(R / f"d3_detections_{s}.parquet", columns=["cc_stack", "old", "new9"])
        for rule, ls in [("old", "--"), ("new9", "-")]:
            v = det.loc[det[rule], "cc_stack"].values
            if len(v) < 10:
                continue
            ax.hist(v, bins=np.linspace(0.5, 1.0, 60), histtype="step", ls=ls, lw=1.4,
                    density=True, label=f"{s} {LBL[rule]}")
    ax.set_xlabel("3-C stacked cc")
    ax.set_ylabel("density")
    ax.set_title("L  stacked-cc distribution", fontsize=10, loc="left")
    ax.legend(fontsize=6)
    ax.grid(alpha=0.3)

    fig.suptitle("D3 pilot — waveform re-detection with 3-C stacking and per-template MAD thresholds "
                 "vs the stored fixed-cc rule", fontsize=12)
    fig.savefig(a.out, dpi=125, bbox_inches="tight")
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
