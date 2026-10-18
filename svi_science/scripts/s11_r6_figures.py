#!/usr/bin/env python
"""s11_r6_figures -- S11 round 6 figures (PLAN.md §10). Points coloured by weight, no error bars.
fig_r6_t1.png   : per-major deep-pulse recovery (waveform injection) vs constraint leverage h_ee, E3 / E4 / adopted, L-SEL.
fig_r6_size.png : ΔD(e) vs episode size per region x measure (headline ΔD source), coloured by 1/σ², thin line in size
                  order, WLS line."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
R6 = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/round6"
FIG = f"{R6}/figures"; os.makedirs(FIG, exist_ok=True)
REG = [("svi", "SVI"), ("wa_nw", "WA-NW")]
NAMES = {"G": "GNSS offset D (mm)", "M0": "Tier-B moment (1e18 N·m)", "G_WA": "WA GNSS offset (mm)", "T": "tremor (1,000 epicentres)",
         "L": "LFE count (1,000)", "Du": "duration (10 d)", "Pc": "patch coverage", "Fs": "family-set sumP (1,000)",
         "Fc": "family-set coverage"}


def main():
    ad = json.load(open(f"{R6}/adopted.json")); a = ad["task2_source"]
    # --- T1 cause figure
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (sub, lab) in zip(axs, REG):
        P = pd.read_csv(f"{R6}/{sub}/tables/r6_t1_permajor.csv")
        pa = f"{R6}/{sub}/tables/r6a_t1_permajor.csv"
        if os.path.exists(pa):
            P = pd.concat([P, pd.read_csv(pa)], ignore_index=True)
        P = P[P.setup == "two3_20"]
        for f, c, m in (("E3", "#888888", "s"), ("E4", "#d95f02", "o"), ("E4y", "#7570b3", "v"), ("E4lu", "#1b9e77", "^")):
            if f not in set(P.form):
                continue
            q = P[P.form == f]
            ax.scatter(q.h_ee, q.R_wave, s=36, c=c, marker=m, label=f"{f} (median {q.R_wave.median():.2f})", edgecolor="white", lw=0.6)
        ax.axhspan(0.8, 1.1, color="#eeeeee", lw=0, zorder=0)
        ax.set_xlabel("leverage h_ee of the episode in the D constraint basis"); ax.set_ylabel("deep-pulse recovery (waveform injection)")
        ax.set_title(f"{lab}, L-SEL: per-major recovery"); ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r6_t1.png", dpi=150); plt.close(fig)
    # --- size figure
    for src, pre, fname, flab in (("headline", "r6", "fig_r6_size.png", f"{a}, L-SEL, pre-registered headline"),
                                  ("amend", "r6p", "fig_r6_size_E4lu.png", "E4lu, L-SEL, post-hoc descriptive column")):
        if not os.path.exists(f"{R6}/svi/tables/{pre}_size_regress.csv"):
            continue
        panels = []
        for sub, lab in REG:
            R = pd.read_csv(f"{R6}/{sub}/tables/{pre}_size_regress.csv"); R = R[(R.source == src) & ~R.majors_only]
            S = pd.read_csv(f"{R6}/{sub}/tables/r6_sizes.csv").set_index("episode")
            D = pd.read_csv(f"{R6}/{sub}/tables/{pre}_size_dD.csv"); D = D[D.source == src].set_index("episode")
            for r in R.itertuples():
                panels.append((lab, r, S, D))
        nc = 4; nr = int(np.ceil(len(panels) / nc))
        fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.4 * nr), squeeze=False)
        for ax, (lab, r, S, D) in zip(axs.ravel(), panels):
            ee = r.episodes.split(";")
            x = S.loc[ee, r.measure].to_numpy(float); y = D.loc[ee, "dD"].to_numpy(float); w = 1 / D.loc[ee, "sig"].to_numpy(float) ** 2
            o = np.argsort(x)
            ax.plot(x[o], y[o], color="#bbbbbb", lw=0.8, zorder=1)
            sc = ax.scatter(x, y, c=w, cmap="viridis", s=34, zorder=2, edgecolor="white", lw=0.5)
            for i, e in enumerate(ee):
                if bool(S.loc[e, "major"]) and (x[i] >= np.percentile(x, 85) or abs(y[i]) >= np.percentile(np.abs(y), 90)):
                    ax.annotate(e, (x[i], y[i]), fontsize=6.5, color="#444444", xytext=(3, 3), textcoords="offset points")
            xx = np.linspace(x.min(), x.max(), 10); ax.plot(xx, r.a + r.b * xx, color="#d95f02", lw=1.4, ls="--", zorder=3)
            ax.axhline(0, color="#999999", lw=0.6)
            ax.set_title(f"{lab}: {NAMES[r.measure]}\nb {r.b:+.3f} ± {r.sig_b:.3f} %/unit, scrambled p {r.p_scr:.2f}", fontsize=8.5)
            ax.set_xlabel(NAMES[r.measure], fontsize=8); ax.set_ylabel("ΔD (%)", fontsize=8); ax.tick_params(labelsize=7)
            cb = fig.colorbar(sc, ax=ax, pad=0.01); cb.ax.tick_params(labelsize=6); cb.set_label("weight 1/σ²", fontsize=6.5)
            ax.grid(alpha=0.2)
        for ax in axs.ravel()[len(panels):]:
            ax.axis("off")
        fig.suptitle(f"S11 round 6: per-episode deep-weighted change ΔD vs episode size ({flab}; no error bars)", fontsize=10)
        fig.tight_layout(); fig.savefig(f"{FIG}/{fname}", dpi=150); plt.close(fig)
    open(f"{R6}/FIGURES_DONE", "w").write("done\n")


if __name__ == "__main__":
    main()
