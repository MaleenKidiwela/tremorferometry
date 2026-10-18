"""Q1 figures: residual share vs time since last active; example episode maps; what explains patch strength."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.colors import LinearSegmentedColormap

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q1_patch_budget"
FIG = f"{ROOT}/svi_science/figures"
BLUE, ORANGE, AQUA, GRAY, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#f0efec", "#0b0b0b", "#52514e"
DIV = LinearSegmentedColormap.from_list("div", [BLUE, GRAY, ORANGE])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white"})
T = json.load(open(f"{OUT}/tests.json"))


def demean(df, col):
    return df[col] - df.groupby("ti")[col].transform("mean")


def fig_resid():
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, name, title in [(axs[0], "SVI_primary", "S. Vancouver Island 2015–2026 (10 big ETS, 12 groups)"),
                            (axs[1], "NW_primary", "Washington-NW 2015–2026 (14 big ETS, 10-km bins)"),
                            (axs[2], "JOINT_2003_2026", "S. Vancouver Island 2003–2026 (Bostock 10 + ours 10)")]:
        df = pd.read_csv(f"{OUT}/rows_{name}.csv")
        df = df[~df.cens & np.isfinite(df.r)].copy()
        df["x"] = np.log10(df.T1)
        if name != "JOINT_2003_2026":
            K = np.load(f"{OUT}/kin_{name}.npz", allow_pickle=True)
            df["inr"] = K["inr"][df.p.to_numpy(), df.v.to_numpy()]
        else:
            df["inr"] = 1.0
        df["ry"] = demean(df, "r"); df["xx"] = demean(df, "x")
        m = df.inr > 0.5
        ax.scatter(10 ** df.x[m], df.r[m], s=18, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.6,
                   label="patch inside the tremor extent")
        if (~m).any():
            ax.scatter(10 ** df.x[~m], df.r[~m], s=18, facecolor="white", edgecolor=ORANGE, linewidth=1.0,
                       label="patch outside the tremor extent")
        # binned means (within-episode adjusted, added back to the overall mean)
        q = pd.qcut(df.x, 5, duplicates="drop")
        bm = df.groupby(q, observed=True).agg(x=("x", "median"), y=("ry", "mean"))
        ax.plot(10 ** bm.x, bm.y + df.r.mean(), color=INK, lw=2, marker="o", ms=5, label="binned mean (episode-adjusted)")
        ax.axhline(0, color=INK2, lw=0.6, ls=":")
        ax.set_xscale("log"); ax.set_xlabel("days since the patch was last active")
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_ylabel("residual share  log(obs / expected)")
        t = T[name]["T1"]
        ax.set_title(title, fontsize=9, color=INK)
        ax.text(0.02, 0.03, f"slope {t['beta']:+.2f} per SD, order-shuffle p = {t['p_order']:.3f}", transform=ax.transAxes,
                fontsize=8, color=INK2)
        ax.legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(f"{FIG}/q1_resid_vs_tsince.png", dpi=150); plt.close(fig)


def fig_examples(region, name, ncol=5):
    K = np.load(f"{OUT}/kin_{name}.npz", allow_pickle=True)
    df = pd.read_csv(f"{OUT}/rows_{name}.csv", dtype={"patch": str})
    PA = pd.read_csv(f"{OUT}/patches.csv", dtype={"patch": str})
    pdf = "group" if region == "SVI" else "bin10"
    pa = PA[(PA.catalogue == region) & (PA.patchdef == pdf)].set_index("patch")
    CL = pd.read_parquet(f"{OUT}/counts_long.parquet")
    TR = pd.read_csv(f"{OUT}/tremor_ep.csv"); TR = TR[TR.region == region].set_index("episode")
    eps = list(K["eps"]); n = len(eps); nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.7 * ncol, 2.6 * nrow + 0.6), squeeze=False, layout="constrained")
    for i, ep in enumerate(eps):
        ax = axs[i // ncol, i % ncol]
        d = df[df.episode == ep].set_index("patch")
        cnt = CL[(CL.catalogue == region) & (CL.variant == "once075") & (CL.patchdef == pdf) & (CL.episode == ep)].set_index("patch").n
        lat, lon = pa.loc[d.index, "lat"], pa.loc[d.index, "lon"]
        sc = ax.scatter(lon, lat, c=d.r.clip(-2.5, 2.5), cmap=DIV, vmin=-2.5, vmax=2.5, s=40 + 3 * np.sqrt(cnt.reindex(d.index).fillna(0)) * 8,
                        edgecolor=INK2, linewidth=0.5)
        inr = K["inr"][:, i]
        for (pt, row), ir in zip(d.iterrows(), inr):
            lab = "—" if row.cens else f"{row.T1:.0f}"
            ax.annotate(lab, (pa.loc[pt, "lon"], pa.loc[pt, "lat"]), xytext=(4, 4), textcoords="offset points", fontsize=6, color=INK2)
            if ir < 0.5:
                ax.scatter(pa.loc[pt, "lon"], pa.loc[pt, "lat"], marker="x", s=25, color=INK, linewidth=0.8)
        tr = TR.loc[ep]
        dirtxt = {1: "front → NW", -1: "front → SE", 0: "no front"}[int(tr.trem_dir)]
        ax.set_title(f"{ep}\n{int(cnt.sum())} LFEs, {dirtxt}", fontsize=7.5, color=INK)
        ax.set_aspect(1 / np.cos(np.radians(lat.mean())), adjustable="datalim")
        ax.tick_params(labelsize=6); ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3)); ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    for j in range(n, nrow * ncol):
        axs[j // ncol, j % ncol].axis("off")
    cb = fig.colorbar(sc, ax=axs, shrink=0.6, pad=0.02)
    cb.set_label("residual share log(obs/expected)\n orange = more than usual, blue = less", fontsize=7)
    fig.suptitle(f"{region}: each big ETS. Colour = patch residual share; number = days since the patch was last active "
                 f"(— = not active since record start)\nx = outside the PNSN tremor extent; circle size = LFE count", fontsize=8, color=INK)
    fig.savefig(f"{FIG}/q1_example_episodes_{region.lower()}.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def fig_explain():
    pr = pd.read_csv(f"{OUT}/partial_r2.csv")
    ev = json.load(open(f"{OUT}/episode_level_values.json"))
    labels = {"in_reach": "inside tremor extent", "km_outside": "km outside tremor extent", "u_along_path": "position along front path",
              "dir_mode": "same-direction pattern", "T1": "days since last active", "T2": "previous big-ETS residual",
              "T3": "cumulative past residual", "detectability": "station availability", "kinematics_set": "ALL front terms",
              "budget_set": "ALL budget terms", "all": "everything"}
    order = ["in_reach", "km_outside", "u_along_path", "dir_mode", "T1", "T2", "T3", "detectability", "kinematics_set", "budget_set", "all"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    ax = axs[0]; y = np.arange(len(order)); h = 0.38
    for k, (name, col, lab) in enumerate([("SVI_primary", BLUE, "S. Vancouver Island"), ("NW_primary", ORANGE, "Washington-NW")]):
        q = pr[pr.design == name].set_index("predictor").reindex(order)
        ax.barh(y + (k - 0.5) * h, q.r2_within, height=h * 0.9, color=col, label=lab)
    ax.set_yticks(y); ax.set_yticklabels([labels[o] for o in order], fontsize=8); ax.invert_yaxis()
    ax.set_xlabel("share of within-episode variance explained (R²)"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("What explains which patches are strong in a big ETS", fontsize=9)
    ax = axs[1]
    for name, col, lab in [("SVI_primary", BLUE, "S. Vancouver Island"), ("NW_primary", ORANGE, "Washington-NW")]:
        v = ev[name]
        ax.scatter(v["trem_width"], v["U"], color=col, s=36, edgecolor="white", linewidth=0.8, label=lab)
        for e, x, u in zip(v["eps"], v["trem_width"], v["U"]):
            ax.annotate(e, (x, u), xytext=(3, 2), textcoords="offset points", fontsize=6, color=INK2)
    ax.set_yscale("log"); ax.set_xlabel("along-strike extent of PNSN tremor in the episode (km, 5–95 %)")
    ax.set_ylabel("how uneven the episode is (χ²/df of patch counts)")
    s = T["SVI_primary|episode_level"]["trem_width"]; n_ = T["NW_primary|episode_level"]["trem_width"]
    ax.set_title(f"Episodes that cover less ground are more uneven (Spearman SVI {s['rho']:+.2f}, NW {n_['rho']:+.2f})", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIG}/q1_partial_r2.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    fig_resid()
    fig_examples("SVI", "SVI_primary")
    fig_examples("NW", "NW_primary")
    fig_explain()
    print("figures done")
