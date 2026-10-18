#!/usr/bin/env python
"""s11_figures -- S11 pilot figures (user rules: dv/v points coloured by CC, NO error bars, CC-weighted averages, lines).
dv/v = -eps, in %. Uncertainties are in tables/text only.
  fig_kernels.png        a_d and a_s vs family depth per station (ell 40, 2-5 s, 3C) + coefficient per layer
  fig_D_S_series.png     D(e) (deep 25-45 km) and S_site(e) per episode, raw and SVD paths
  fig_cells_by_site.png  cell dv/v (series-demeaned) per site, CC-coloured, CC-weighted episode means as a line
  fig_seasonal.png       D and S vs day of year
  fig_injection.png      injected-minus-natural D for runs A (deep step) and B (shallow annual)
  fig_gates_scatter.png  split-network and lapse comparisons of D
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]
SITECOL = dict(zip(C.SITES, CAT))
CMAP, VMIN, VMAX = "viridis", 0.8, 1.0
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6",
                     "grid.linewidth": 0.5, "font.size": 9})


def head():
    try:
        return json.load(open(f"{C.TB}/headline.json"))["headline"]
    except Exception:
        return "raw"


def fig_kernels():
    K = pd.read_parquet(f"{C.TB}/kernel_coef.parquet")
    k = K[(K.variant == "ell40") & (K.win == C.HEAD) & (K.comp == "3c")].groupby(["fam", "netsta"]).first().reset_index()
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for i, ns in enumerate(C.STATIONS):
        q = k[k.netsta == ns].sort_values("fam_depth")
        ax[0].plot(q.fam_depth, q.a_d, "o", ms=4, color=CAT[i], label=ns)
        ax[1].plot(q.fam_depth, q.a_s, "o", ms=4, color=CAT[i], label=ns)
    ax[0].set(xlabel="family depth (km)", ylabel="a_d (25–45 km)", title="deep coefficient, 2–5 s past S, 3C, ell 40 km")
    ax[1].set(xlabel="family depth (km)", ylabel="a_s (0–5 km)", title="shallow coefficient")
    ax[0].legend(fontsize=7, ncol=2)
    cols = [c for c in K.columns if c.startswith("c_")]
    lay = np.array([[k[k.netsta == ns][[f"c_{l}|{z}" for z in C.ZONES]].sum(1).median() for l in C.LAYER_NAMES] for ns in C.STATIONS])
    y = np.arange(len(C.LAYER_NAMES))
    for i, ns in enumerate(C.STATIONS):
        ax[2].plot(lay[i], y, "-o", ms=3, lw=1.5, color=CAT[i], label=ns)
    ax[2].set_yticks(y); ax[2].set_yticklabels(C.LAYER_NAMES); ax[2].invert_yaxis(); ax[2].axvline(0, color="grey", lw=0.6)
    ax[2].set(xlabel="coefficient per layer (median over families)", ylabel="layer (km)", title="where 2–5 s dv/v is sensitive")
    fig.tight_layout(); fig.savefig(f"{C.FIG}/fig_kernels.png", dpi=130); plt.close(fig)


def fig_series():
    fig, ax = plt.subplots(2, 2, figsize=(16, 7.5), sharex=True)
    for j, path in enumerate(("raw", "svd")):
        p = f"{C.TB}/inv_{path}_head.csv"
        if not os.path.exists(p):
            continue
        par = pd.read_csv(p, parse_dates=["mid"])
        D = par[par.kind == "D"].sort_values("mid")
        a = ax[0, j]
        a.plot(D.mid, D.value, "-", color="#52514e", lw=1, zorder=2)
        sc = a.scatter(D.mid, D.value, c=D.cc_w, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=12 + 2 * np.sqrt(D.n_cells), zorder=3,
                       edgecolor="k", lw=0.3)
        a.axhline(0, color="grey", lw=0.6, ls=":")
        a.set_title(f"{path} path{' (HEADLINE)' if path == head() else ''}: deep D(e), 25–45 km layer\n"
                    f"RMS {np.sqrt(np.mean(D.value ** 2)):.3f} %, median σ {D.sig.median():.3f} % (σ not drawn)", fontsize=9)
        a.set_ylabel("dv/v (%)")
        b = ax[1, j]
        for s in C.SITES:
            q = par[(par.kind == "S") & (par.site == s)].sort_values("mid")
            b.plot(q.mid, q.value, "-", color=SITECOL[s], lw=1.2, label=f"S {s}")
        b.axhline(0, color="grey", lw=0.6, ls=":"); b.set_ylabel("dv/v (%)"); b.legend(fontsize=7, ncol=4)
        b.set_title(f"{path} path: shallow S_site(e), 0–5 km, per site", fontsize=9)
    fig.colorbar(sc, ax=ax[0, :], shrink=0.8, label="CC-weighted mean cc of cells")
    fig.savefig(f"{C.FIG}/fig_D_S_series.png", dpi=125, bbox_inches="tight"); plt.close(fig)


def fig_cells():
    path = head()
    c = pd.read_csv(f"{C.TB}/cells_{path}_head.csv")
    c["mid"] = c.episode.map(C.episodes().mid)
    w = 1 / c.sig ** 2
    key = c.fam.astype(str) + "|" + c.netsta + "|" + c.epoch
    c["dm"] = c.d - (c.d * w).groupby(key).transform("sum") / w.groupby(key).transform("sum")
    fig, ax = plt.subplots(len(C.SITES), 1, figsize=(15, 11), sharex=True)
    for i, s in enumerate(C.SITES):
        q = c[c.site == s]
        sc = ax[i].scatter(q.mid, q.dm, c=q.cc, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=6, alpha=0.6, lw=0)
        e = q.groupby("mid").apply(lambda g: pd.Series(dict(v=np.average(g.dm, weights=g.cc), n=len(g))), include_groups=False)
        e = e[e.n >= 3]
        ax[i].plot(e.index, e.v, "-", color="k", lw=1.2)
        ax[i].axhline(0, color="grey", lw=0.6, ls=":"); ax[i].set_ylim(-0.5, 0.5)
        ax[i].set_title(f"{s} site — {len(q)} cells; line = CC-weighted episode mean (≥ 3 cells)", fontsize=9)
        ax[i].set_ylabel("dv/v (%)")
    fig.colorbar(sc, ax=ax, shrink=0.5, label="cell cc (3C, 2–5 s)")
    fig.savefig(f"{C.FIG}/fig_cells_by_site.png", dpi=120, bbox_inches="tight"); plt.close(fig)


def fig_seasonal():
    path = head()
    par = pd.read_csv(f"{C.TB}/inv_{path}_head.csv", parse_dates=["mid"])
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
    D = par[par.kind == "D"]
    sc = ax[0].scatter(D.mid.dt.dayofyear, D.value, c=D.cc_w, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=18, edgecolor="k", lw=0.3)
    f = C.annual_fit(D.mid, D.value.to_numpy(), D.sig.to_numpy()); x = np.arange(1, 366)
    ax[0].plot(x, f["c"] + f["a_cos"] * np.cos(2 * np.pi * (x - 1) / 365.25) + f["a_sin"] * np.sin(2 * np.pi * (x - 1) / 365.25),
               "-", color="k", lw=1.2)
    ax[0].set(xlabel="day of year (episode mid)", ylabel="dv/v (%)",
              title=f"deep D: annual amplitude {f['amp']:.3f} ± {f['amp_sig']:.3f} %")
    for s in C.SITES:
        q = par[(par.kind == "S") & (par.site == s)]
        fs = C.annual_fit(q.mid, q.value.to_numpy(), q.sig.to_numpy())
        ax[1].plot(q.mid.dt.dayofyear, q.value, "o", ms=3, color=SITECOL[s], alpha=0.6)
        if np.isfinite(fs["amp"]):
            ax[1].plot(x, fs["c"] + fs["a_cos"] * np.cos(2 * np.pi * (x - 1) / 365.25) + fs["a_sin"] * np.sin(2 * np.pi * (x - 1) / 365.25),
                       "-", color=SITECOL[s], lw=1.5, label=f"{s}: {fs['amp']:.3f} ± {fs['amp_sig']:.3f} %")
    ax[1].legend(fontsize=7); ax[1].set(xlabel="day of year", title="shallow S per site: annual fits")
    for a in ax:
        a.axhline(0, color="grey", lw=0.6, ls=":")
    fig.colorbar(sc, ax=ax[0], label="CC-weighted mean cc")
    fig.tight_layout(); fig.savefig(f"{C.FIG}/fig_seasonal.png", dpi=125); plt.close(fig)


def fig_injection():
    p = f"{C.TB}/inject_deltas.csv"
    if not os.path.exists(p):
        return
    d = pd.read_csv(p, parse_dates=["mid"])
    fig, ax = plt.subplots(2, 2, figsize=(15, 6.5), sharex=True)
    for j, path in enumerate(("raw", "svd")):
        for i, (run, truth) in enumerate((("A", lambda t: 0.3 * (t >= "2015-01-01")),
                                          ("B", lambda t: 0 * t.dt.year))):
            q = d[(d.path == path) & (d.run == run) & (d.kind == "D")].sort_values("mid")
            ax[i, j].plot(q.mid, q.value, "-o", ms=3, lw=1, color=CAT[0], label="recovered D (injected − natural)")
            ax[i, j].plot(q.mid, truth(q.mid), "-", lw=1.5, color=CAT[1], label="injected deep change")
            if run == "B":
                for s in C.SITES:
                    qs = d[(d.path == path) & (d.run == run) & (d.kind == "S") & (d.site == s)].sort_values("mid")
                    ax[i, j].plot(qs.mid, qs.value, "-", lw=0.8, color=SITECOL[s], alpha=0.6, label=f"S {s}")
            ax[i, j].set_title(f"{path}: injection {run} ({'deep step +0.3 %' if run == 'A' else 'shallow annual 0.1 % at all sites'})",
                               fontsize=9)
            ax[i, j].set_ylabel("dv/v (%)"); ax[i, j].legend(fontsize=7, ncol=3)
    fig.tight_layout(); fig.savefig(f"{C.FIG}/fig_injection.png", dpi=125); plt.close(fig)


def fig_gates():
    path = head()
    import s11_invert as I
    from s11_gates import Dser
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    try:
        a = Dser(pd.read_csv(f"{C.TB}/inv_{path}_w2035.csv")); b = Dser(pd.read_csv(f"{C.TB}/inv_{path}_w3560.csv"))
        j = pd.concat([a.value.rename("a"), b.value.rename("b"), a.cc_w.rename("cc")], axis=1).dropna()
        sc = ax[1].scatter(j.a, j.b, c=j.cc, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=18, edgecolor="k", lw=0.3)
        lim = np.nanmax(np.abs(j[["a", "b"]].to_numpy())) * 1.1
        ax[1].plot([-lim, lim], [-lim, lim], "-", color="grey", lw=0.8)
        ax[1].set(xlabel="D from 2–3.5 s (%)", ylabel="D from 3.5–6 s (%)", title="gate (f): lapse")
    except FileNotFoundError:
        pass
    c = pd.read_csv(f"{C.TB}/cells_{path}_head.csv")

    def sub(m):
        q = c[m].reset_index(drop=True); key = q.fam.astype(str) + "|" + q.netsta + "|" + q.epoch
        q = q[key.map(key.value_counts()) >= 3].reset_index(drop=True)
        return Dser(I.invert(q, extras=False)[0])
    A, B = sub(c.site == "PGC"), sub(c.site != "PGC")
    j = pd.concat([A.value.rename("a"), B.value.rename("b"), A.cc_w.rename("cc")], axis=1).dropna()
    ax[0].scatter(j.a, j.b, c=j.cc, cmap=CMAP, vmin=VMIN, vmax=VMAX, s=18, edgecolor="k", lw=0.3)
    lim = np.nanmax(np.abs(j[["a", "b"]].to_numpy())) * 1.1
    ax[0].plot([-lim, lim], [-lim, lim], "-", color="grey", lw=0.8)
    ax[0].set(xlabel="D from PGC site (%)", ylabel="D from B926 + VGZ + SNB (%)", title="gate (d): split network")
    fig.tight_layout(); fig.savefig(f"{C.FIG}/fig_gates_scatter.png", dpi=125); plt.close(fig)


if __name__ == "__main__":
    for f in (fig_kernels, fig_series, fig_cells, fig_seasonal, fig_injection, fig_gates):
        try:
            f(); print("ok", f.__name__)
        except Exception as e:
            import traceback; print("FAILED", f.__name__, e, traceback.format_exc()[-600:])
    open(f"{C.OUT}/FIGURES_DONE", "w").write(pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S\n"))
