#!/usr/bin/env python
"""s10_svd_figures -- figures for the SVD-Wiener step of the S10 pilot (rules: results/s10_dvv/pilot/svd/SVD.md §6).
(1) packed day-stack image raw | SVD side by side (B011, PGC); (2) daily dv/v (= -eps) vs calendar time, raw and SVD,
+- bootstrap sigma (SVD re-run inside each draw), ETS shading, epoch boundaries."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S

TAG = sys.argv[1] if len(sys.argv) > 1 else ""
OUTD = f"{S.OUT}/svd"; T, F = f"{OUTD}/tables", f"{OUTD}/figures"
os.makedirs(F, exist_ok=True)
BLUE, ORANGE, GRAY, INK, MUTED = "#2a78d6", "#eb6834", "#8a8984", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
EPI = pd.read_csv(f"{S.ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")
D = pd.read_parquet(f"{T}/daily_svd{TAG}.parquet")
Z = np.load(f"{T}/stacks_svd{TAG}.npz")
KEYS = list(Z["keys"]); TAU = Z["tau"]
STAS = ["PB.B011", "PB.B009", "PB.B010", "PB.B926", "CN.PGC"]


def series_mats(fam, ns):
    out = []
    for i, k in enumerate(KEYS):
        f, n_, ep = k.split("|")
        if int(f) == fam and n_ == ns:
            out.append((ep, Z[f"St_{i}"], Z[f"Rf_{i}"], pd.to_datetime(Z[f"days_{i}"]), float(Z["tS"][i])))
    return sorted(out, key=lambda x: x[0])


def stacks_fig(fam, stas=("PB.B011", "CN.PGC")):
    fig, axes = plt.subplots(2, 2 * len(stas), figsize=(6.2 * len(stas) * 1.6, 13),
                             gridspec_kw=dict(height_ratios=[1, 4]), squeeze=False)
    for j, ns in enumerate(stas):
        ser = series_mats(fam, ns)
        if not ser:
            continue
        for kk, (lab, which) in enumerate((("raw (realigned day stacks)", 1), ("SVD-Wiener", 2))):
            ax0, ax = axes[0, 2 * j + kk], axes[1, 2 * j + kk]
            rows, days, bounds = [], [], []
            for s in ser:
                A = s[which].astype(float)
                # shift each epoch to its own S pick so all rows share tau = time after picked S
                sh = int(round(s[4] * S.FS))
                A = np.roll(A, -sh, axis=1)
                bounds.append((len(rows), s[0])); rows.extend(list(A)); days.extend(list(s[3]))
            M_ = np.array(rows)
            g = M_.mean(0)
            ax0.plot(TAU, g / np.abs(g).max(), color=INK, lw=0.8)
            ax0.axvspan(2, 4, color=ORANGE, alpha=0.18, lw=0, label="2–4 s coda (dv/v window)")
            ax0.axvline(0, color="r", lw=0.8, ls="--"); ax0.set_xlim(-14, 10)
            ax0.set_title(f"fam {fam} at {ns}, {lab}: mean of {len(days)} day rows", fontsize=9, loc="left")
            norm = M_ / (np.abs(M_).max(1, keepdims=True) + 1e-12)
            ax.imshow(norm, aspect="auto", cmap="seismic", vmin=-1, vmax=1, interpolation="nearest",
                      extent=[TAU[0], TAU[-1], len(days) - 0.5, -0.5])
            for x in (0, 2, 4):
                ax.axvline(x, color="k", lw=0.5, ls="--" if x == 0 else "-")
            yt = [i for i in range(len(days)) if i == 0 or days[i].year != days[i - 1].year]
            ax.set_yticks(yt); ax.set_yticklabels([f"{days[i]:%Y-%m-%d}" for i in yt], fontsize=7)
            for k0, ep in bounds:
                ax.axhline(k0 - 0.5, color="lime", lw=1.2); ax.text(-13.8, k0 + 1.5, f"epoch {ep}", color="lime", fontsize=8)
            ax.set_xlim(-14, 10); ax.set_xlabel("time after picked S (s)")
        axes[0, 2 * j].legend(loc="upper left", fontsize=7, frameon=False)
    fig.suptitle(f"Family {fam}: one row per active day (ΣP ≥ 3), each row normalised; left of each pair raw, right after "
                 f"the original pipeline's SVD-Wiener filter (per station-epoch)", y=0.995, fontsize=10)
    fig.tight_layout(); fig.savefig(f"{F}/s10_svd_stacks_fam{fam}.png", dpi=120); plt.close(fig)


def shade(ax):
    for r in EPI.itertuples():
        ax.axvspan(pd.Timestamp(r.start), pd.Timestamp(r.end) + pd.Timedelta(days=1),
                   color="#d9d8d3" if r.big else "#eeede9", lw=0, zorder=0)


def daily_fig(fam):
    fig, axs = plt.subplots(len(STAS), 1, figsize=(13, 12), sharex=True)
    for ax, ns in zip(axs, STAS):
        shade(ax)
        for eid, a, b in S.EPOCHS.get(ns, []):
            for d in (a, b):
                if "2002" < d < "2027":
                    ax.axvline(pd.Timestamp(d), color=INK, lw=0.9, ls="--", zorder=1)
        g = D[(D.fam == fam) & (D.netsta == ns)]
        t = pd.to_datetime(g.date)
        txt = []
        for tag, col, mk, lab, dx in (("raw", BLUE, "o", "raw day stack", -12), ("svd", ORANGE, "s", "SVD-Wiener", 12)):
            ax.errorbar(t + pd.Timedelta(days=dx), -100 * g[f"{tag}_eps24"], yerr=100 * g[f"{tag}_sig24"], fmt=mk, ms=3,
                        color=col, ecolor=col, elinewidth=0.6, alpha=0.85, label=f"{lab} ± bootstrap σ", zorder=3)
            txt.append(f"{tag}: median σ {100 * g[f'{tag}_sig24'].median():.2f} %, measured {g[f'{tag}_eps24'].notna().sum()}/{len(g)}")
        ax.axhline(0, color=MUTED, lw=0.5); ax.set_ylim(-1.15, 1.15); ax.set_ylabel(f"{ns}\ndv/v (%)")
        ax.text(0.005, 0.88, "   ".join(txt), transform=ax.transAxes, fontsize=8, color=MUTED)
    axs[0].legend(loc="upper right", fontsize=7, ncol=2, frameon=False)
    axs[0].set_title(f"Family {fam}: daily dv/v (= −ε) in 2–4 s after S, raw vs SVD-Wiener; each epoch its own zero "
                     f"(dashed = epoch boundary; grey = ETS episodes, darker = big; points offset ±12 d for visibility)",
                     fontsize=9, loc="left")
    axs[-1].set_xlim(pd.Timestamp("2003-01-01"), pd.Timestamp("2026-10-01"))
    fig.tight_layout(); fig.savefig(f"{F}/s10_svd_daily_fam{fam}.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    fams = sorted(D.fam.unique())
    for f in fams:
        stacks_fig(f); daily_fig(f)
    print("figures written", fams)
