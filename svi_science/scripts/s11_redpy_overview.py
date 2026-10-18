#!/usr/bin/env python
"""s11_redpy_overview -- cross-family overview figures: (1) episode-to-episode SD by estimator/node type vs the existing
family-10 PGC baselines; (2) co-located PGC vs B011 node dv/v for identical event sets; (3) graph error checks."""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_redpy_common as R
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
FAMS = [10, 87, 17]; STA = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
# (1) scatter
fig, axs = plt.subplots(1, 3, figsize=(17, 5), sharey=True)
mk = {"day": "o", "N20": "s", "N50": "D", "episode": "^"}
for ax, f in zip(axs, FAMS):
    T4 = pd.read_csv(f"{R.TB}/T4_scatter_fam{f}.csv"); T4 = T4[(T4.subset == "all") & (T4.weight == "P")]
    for i, st in enumerate(STA):
        for j, (est, col) in enumerate((("graph_boot", "k"), ("reference", "tab:red"))):
            for k_, kind in enumerate(mk):
                q = T4[(T4.station == st) & (T4.est == est) & (T4.kind == kind)]
                if len(q):
                    x = i + (j - 0.5) * 0.3 + (k_ - 1.5) * 0.06
                    ax.plot(x, 100 * q.ep_sd.iloc[0], mk[kind], color=col, ms=6, mfc="none" if j else col,
                            label=f"{'graph' if j == 0 else 'reference'} {kind}" if i == 0 else None)
    if f == 10:
        ax.axhline(0.20, color="grey", ls="--", lw=1); ax.text(3.4, 0.205, "raw Z 2–4 s episode stacks 0.20 %", fontsize=7, ha="right")
        ax.axhline(0.122, color="grey", ls=":", lw=1); ax.text(3.4, 0.126, "raw 3C 2–5 s day→episode (CC>0.5) 0.122 %", fontsize=7, ha="right")
    ax.set_xticks(range(4)); ax.set_xticklabels(STA); ax.set_title(f"Family {f}: episode-to-episode SD of dv/v", fontsize=10)
    ax.set_ylim(0, 0.35)
axs[0].set_ylabel("SD of episode means (%)"); axs[0].legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(f"{R.FIG}/overview_episode_scatter.png", dpi=110); plt.close(fig)
# (2) co-located
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, f in zip(axs, FAMS):
    V = pd.read_parquet(f"{R.TB}/coloc_nodes_fam{f}.parquet"); q = V[(V.kind == "N20") & (V.subset == "all")]
    P = q.pivot_table(index="node", columns="netsta", values="v") * -100
    rho = pd.read_csv(f"{R.TB}/coloc_rho_fam{f}.csv"); rho = rho[(rho.kind == "N20") & (rho.subset == "all")].set_index(["sta1", "sta2"])
    for st, c in (("PB.B009", "tab:blue"), ("PB.B010", "tab:green"), ("PB.B011", "tab:orange")):
        r = rho.loc[("CN.PGC", st)]
        ax.plot(P["CN.PGC"], P[st], "o", ms=4, color=c, alpha=0.7, label=f"{st}: ρ {r.rho:.2f} (p {r.p:.3f})")
    lim = 0.7; ax.plot([-lim, lim], [-lim, lim], "k:", lw=0.7); ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("PGC dv/v (%), N=20 chunk"); ax.set_ylabel("borehole dv/v (%), same events"); ax.legend(fontsize=7)
    ax.set_title(f"Family {f}: co-located instruments, identical events (2005-11 → 2015-01)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{R.FIG}/overview_colocated.png", dpi=110); plt.close(fig)
# (3) error checks
fig, axs = plt.subplots(1, 2, figsize=(13, 4.4))
for f, c in zip(FAMS, ["k", "tab:blue", "tab:orange"]):
    T2 = pd.read_csv(f"{R.TB}/T2_misclosure_fam{f}.csv"); T2 = T2[T2.subset == "all"]
    for k_, kind in enumerate(["day", "N20", "N50", "episode"]):
        q = T2[T2.kind == kind]; axs[0].plot(np.full(len(q), k_) + FAMS.index(f) * 0.12 - 0.12, q.ratio_med, "o", color=c, label=f"family {f}" if k_ == 0 else None)
    T3b = pd.read_csv(f"{R.TB}/T3b_shuffled_fam{f}.csv")
    for k_, kind in enumerate(["N20", "N50"]):
        q = T3b[T3b.kind == kind]
        axs[1].plot(np.full(len(q), k_) + FAMS.index(f) * 0.12 - 0.12, q.chi2_red, "o", color=c, mfc="none", label=f"family {f}, misclosure σ" if k_ == 0 else None)
        axs[1].plot(np.full(len(q), k_ + 2) + FAMS.index(f) * 0.12 - 0.12, q.chi2_red_boot, "o", color=c, label=f"family {f}, bootstrap σ" if k_ == 0 else None)
axs[0].axhspan(0.8, 1.25, color="g", alpha=0.1); axs[0].set_xticks(range(4)); axs[0].set_xticklabels(["day", "N20", "N50", "episode"])
axs[0].set_ylabel("misclosure σ / bootstrap σ (median per series)"); axs[0].set_title("T2: graph misclosure error vs event bootstrap", fontsize=9); axs[0].legend(fontsize=7)
axs[1].set_yscale("log"); axs[1].axhspan(0.5, 1.5, color="g", alpha=0.1); axs[1].set_xticks(range(4)); axs[1].set_xticklabels(["N20 misclosure", "N50 misclosure", "N20 bootstrap", "N50 bootstrap"])
axs[1].set_ylabel("χ²_red of episode means vs 0"); axs[1].set_title("T3b: events shuffled in time (no time structure)", fontsize=9); axs[1].legend(fontsize=6)
fig.tight_layout(); fig.savefig(f"{R.FIG}/overview_error_checks.png", dpi=110); plt.close(fig)
print("ok")
