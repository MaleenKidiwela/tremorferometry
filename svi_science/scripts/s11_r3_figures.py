#!/usr/bin/env python
"""s11_r3_figures -- round-3 figures (user rules: no error bars on dv/v, points coloured by weight, line graphs).
    python s11_r3_figures.py     -> s11_depth/round3/figures/fig_r3_{kerr,q3,zones}.png
"""
from __future__ import annotations
import os, json
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

B = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/round3"
FIG = f"{B}/figures"; os.makedirs(FIG, exist_ok=True)
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # S11 figure palette (fixed order)
REG = [("svi", "SVI"), ("wa_nw", "WA-NW")]
GEN = ["ell10", "ell20", "ell40", "ell80", "two5_40", "two3_20"]
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6",
                     "grid.linewidth": 0.5, "font.size": 9})


def fig_kerr():
    fig, ax = plt.subplots(2, 4, figsize=(17, 7.5))
    rules = [("R1_Q2rec", "R1: deep-step recovery (Q2)", (0.7, 1.3)), ("R3_Bleak", "R3: shallow annual → deep annual leak", (0, 0.25)),
             ("R5_Psleak", "R5: one-episode shallow pulse → ΔD leak", (0, 0.35))]
    for i, (sub, name) in enumerate(REG):
        T = pd.read_csv(f"{B}/{sub}/tables/kerr_rules.csv")
        x = np.arange(len(GEN))
        for j, (col, title, band) in enumerate(rules):
            a = ax[i, j]
            a.axhspan(*band, color="#dfeee4", zorder=0, label="pass band")
            for k, (s, lab) in enumerate((("NOM", "nominal ell40"), ("AVG6", "average of 6"), ("L-SEL", "lapse-chosen (two3_20)"))):
                q = T[T.setup == s].set_index("gen").reindex(GEN)
                a.plot(x, q[col], "-o", color=CAT[k], lw=1.5, ms=5, label=lab)
            lo = [float(T[(T.setup == f"LOKO_{g}") & (T.gen == g)][col].iloc[0]) for g in GEN]
            a.plot(x, lo, "D", color=CAT[3], ms=8, mec="k", mew=0.6, label="leave-one-kernel-out (left-out kernel)")
            a.set_xticks(x); a.set_xticklabels(GEN, rotation=30, fontsize=8)
            a.set(title=f"{name}: {title}", xlabel="kernel used to make the injected signal")
            if i == 0 and j == 0:
                h, l = a.get_legend_handles_labels()
        L = pd.read_csv(f"{B}/{sub}/tables/kerr_lapse_choice.csv").set_index("variant").reindex(GEN)
        a = ax[i, 3]
        a.plot(x, L.dchi2, "-o", color="#52514e", lw=1.5, ms=6)
        a.set_xticks(x); a.set_xticklabels(GEN, rotation=30, fontsize=8)
        a.set(title=f"{name}: which kernel fits the lapse dependence", ylabel="Δχ² vs best (lower = better)")
    fig.legend(h, l, loc="upper center", ncol=5, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(f"{FIG}/fig_r3_kerr.png", dpi=130); plt.close(fig)


def fig_q3():
    fig, ax = plt.subplots(3, 2, figsize=(16, 10), gridspec_kw=dict(height_ratios=[1.6, 1.6, 1]))
    for i, (sub, name) in enumerate(REG):
        T = pd.read_csv(f"{B}/{sub}/tables/r3q3_dD.csv"); T = T[T.major]; T["mid"] = pd.to_datetime(T.mid)
        for r, (cfg, lab) in enumerate((("B0", "B0 raw path, P weights (round-2 headline)"), ("B2", "B2 cc-weighted all-pairs graph (≥ 50-event nodes)"))):
            a = ax[r, i]; q = T[T.config == cfg].sort_values("mid")
            a.plot(q.mid, q.dD, "-", color="#9a9a96", lw=0.9, zorder=1)
            w = 1 / q.sig ** 2
            sc = a.scatter(q.mid, q.dD, c=w / w.max(), cmap="viridis", vmin=0, vmax=1, s=40, edgecolor="k", lw=0.4, zorder=2)
            for rr in q.itertuples():
                if abs(rr.dD / rr.sig) >= 2.5:
                    a.annotate(rr.episode, (rr.mid, rr.dD), xytext=(0, 6), textcoords="offset points", fontsize=7, ha="center")
            a.axhline(0, color="grey", lw=0.6, ls=":")
            a.set(ylabel="ΔD (%) of each major vs its ±3-yr neighbours", title=f"{name} — {lab}; σ in tables, not drawn")
            fig.colorbar(sc, ax=a, pad=0.01).set_label("relative weight 1/σ²")
        S = pd.read_csv(f"{B}/{sub}/tables/r3q3_summary.csv").set_index("config").reindex(["B0", "B0n50", "B1", "B2", "B2p"])
        a = ax[2, i]; x = np.arange(len(S))
        a.axhspan(0.8, 1.1, color="#dfeee4", zorder=0)
        a.axhspan(0, 0.35, color="#f3e6da", zorder=0)
        a.plot(x, S.deep_median, "-o", color=CAT[0], lw=1.5, label="deep pulse recovery (pass 0.8–1.1, green)")
        a.plot(x, S.shal_median, "-o", color=CAT[1], lw=1.5, label="shallow pulse leak (pass ≤ 0.35, tan)")
        a.plot(x, S.sig_major_median * 5, "--s", color=CAT[2], lw=1.2, label="median σ per major × 5 (%)")
        a.set_xticks(x); a.set_xticklabels(S.index); a.set(title=f"{name}: pulse calibration per estimator", ylim=(0, 1.6))
        a.legend(fontsize=7, loc="upper center", ncol=3)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r3_q3.png", dpi=130); plt.close(fig)


def fig_zones():
    fig, ax = plt.subplots(1, 2, figsize=(16, 4.8))
    for i, (sub, name) in enumerate(REG):
        P = pd.read_csv(f"{B}/{sub}/tables/zones_peryear.csv"); R_ = pd.read_csv(f"{B}/{sub}/tables/zones_result.csv")
        d = json.load(open(f"{B}/{sub}/tables/zones_def.json")); nz = d["run_zones"]; bnd = d["bounds"][str(nz)]
        a = ax[i]
        for k, z in enumerate(sorted(P.zone.unique())):
            q = P[P.zone == z].sort_values("year")
            a.plot(q.year, q.dD, "-", color=CAT[k], lw=1.4, label=f"zone {z}")
            w = 1 / q.sig ** 2
            a.scatter(q.year, q.dD, c=CAT[k], s=20 + 80 * w / w.max(), edgecolor="k", lw=0.4, zorder=3)
            c = R_[(R_["item"] == "data") & (R_.zone == z)]
            if len(c):
                a.axhline(float(c.composite.iloc[0]), color=CAT[k], lw=1, ls="--")
        a.axhline(0, color="grey", lw=0.6, ls=":")
        p = R_[R_["item"] == "data_equal_zones_test"].p.iloc[0]
        a.set(xlabel="year", ylabel="zone ΔD (%) of that year's majors",
              title=f"{name}: {nz} along-strike deep zones (bounds {np.round(bnd, 0).tolist()} km, SE→NW)\ndashed = all-major composite; "
                    f"equal-zones p = {p:.2f}; marker size = weight", )
        a.title.set_fontsize(9); a.legend(fontsize=7); a.xaxis.get_major_locator().set_params(integer=True)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r3_zones.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    fig_kerr(); fig_q3(); fig_zones(); print("figures written", FIG)
