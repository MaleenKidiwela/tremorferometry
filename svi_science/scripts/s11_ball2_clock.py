#!/usr/bin/env python
"""s11_ball2_clock -- S11 ballistic round 2, §1 clock / timing-glitch correction (round2/PLAN.md).
Quartet (PGC, B009, B010, B011) instrument terms per (family, bin) from the 6 common-row pair differences (LS, median
gauge), aggregated over families with leave-family-out; corrected co-located floor and far pairs (PGC ref, B011 ref,
site-mean ref); far-station terms; detector-frame drift; glitch calls; metadata comparison.

    python s11_ball2_clock.py [--tag _test]
"""
from __future__ import annotations
import os, sys, json, argparse, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_ball_analysis as A1
from s11_ball_analysis import demean, scatter, decomp

OUT = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic/round2"
TB, FIG = f"{OUT}/tables", f"{OUT}/figures"
QUART = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
FAR = ["PB.B926", "CN.SNB", "CN.VGZ"]
FLOOR_SETS = {"CN.PGC": ["CN.PGC-PB.B009", "CN.PGC-PB.B010", "CN.PGC-PB.B011"],
              "PB.B011": ["PB.B009-PB.B011", "PB.B010-PB.B011"]}
SEED = 20261012
SIGMIN = 0.5
PKEY = ["fam", "pair", "epA", "epB", "bin"]


def prep(pr, val="dS", sig="sig_dS", rail="railS"):
    d = pr[(pr.kind == "data")].copy()
    d.loc[d[rail], val] = np.nan
    d["episode"] = d["lab"]
    out = []
    for bk, g in d.groupby("bin"):
        out.append(demean(g, val, PKEY))
    D = pd.concat(out, ignore_index=True)
    D["sig"] = np.maximum(D[sig], SIGMIN)
    D["nodeA"] = D.A + "|" + D.epA; D["nodeB"] = D.B + "|" + D.epB
    return D


def quartet_terms(D):
    """per (fam, bin, lab): LS node lags from quartet pairs, median gauge -> frame fam, bin, lab, node, k, var"""
    Q = D[D.A.isin(QUART) & D.B.isin(QUART)]
    rows = []
    for (f, bk, lab), g in Q.groupby(["fam", "bin", "lab"]):
        nodes = sorted(set(g.nodeA) | set(g.nodeB))
        if len(nodes) < 2:
            continue
        ix = {n: i for i, n in enumerate(nodes)}
        Am = np.zeros((len(g), len(nodes))); w = 1 / g.sig.to_numpy()
        Am[np.arange(len(g)), [ix[n] for n in g.nodeB]] = 1; Am[np.arange(len(g)), [ix[n] for n in g.nodeA]] = -1
        k = np.linalg.lstsq(Am * w[:, None], g.r.to_numpy() * w, rcond=None)[0]
        k = k - np.median(k)
        prec = np.abs(Am).T @ (w ** 2)                       # sum of edge weights per node
        for n in nodes:
            rows.append(dict(fam=f, bin=bk, lab=lab, node=n, k=k[ix[n]], var=1 / max(prec[ix[n]], 1e-9),
                             date=g.date.min(), cc=float(np.nanmedian(np.r_[g.ccS_A, g.ccS_B]))))
    return pd.DataFrame(rows)


def aggregate(T):
    """K_inst per (bin, lab, node) over families + leave-family-out values per (fam, bin, lab, node)"""
    T = T.copy(); T["w"] = 1 / T["var"]
    g = T.groupby(["bin", "lab", "node"])
    T["W"] = g.w.transform("sum")
    T["WK"] = (T.w * T.k).groupby([T.bin, T.lab, T.node]).transform("sum")
    T["nf"] = g.k.transform("size")
    T["K_lfo"] = (T.WK - T.w * T.k) / (T.W - T.w)
    T.loc[T.nf - 1 < 3, "K_lfo"] = np.nan
    Kt = T.groupby(["bin", "lab", "node"]).apply(lambda x: pd.Series(dict(
        K=(x.w * x.k).sum() / x.w.sum(), n_fam=len(x), sd_fam=x.k.std(ddof=1) if len(x) > 1 else np.nan,
        date=x.date.min(), cc=x.cc.median()))).reset_index()
    Kt["sigK"] = Kt.sd_fam / np.sqrt(Kt.n_fam)
    sa = T.merge(Kt[["bin", "lab", "node", "K"]], on=["bin", "lab", "node"])
    agree = sa.groupby(["bin", "lab", "node"]).apply(lambda x: np.mean(np.sign(x.k) == np.sign(x.K.iloc[0])))
    Kt = Kt.merge(agree.rename("sign_agree").reset_index(), on=["bin", "lab", "node"])
    Kt["glitch"] = (Kt.n_fam >= 3) & (np.abs(Kt.K) >= 3 * Kt.sigK) & (Kt.sign_agree >= 0.7)
    Kt["station"] = Kt.node.str.split("|").str[0]; Kt["epoch"] = Kt.node.str.split("|").str[1]
    return T, Kt


def correct(D, T):
    lfo = T.set_index(["fam", "bin", "lab", "node"]).K_lfo
    def look(f, b, l, n):
        return lfo.get((f, b, l, n), np.nan)
    kA = np.array([look(f, b, l, n) for f, b, l, n in zip(D.fam, D.bin, D.lab, D.nodeA)])
    kB = np.array([look(f, b, l, n) for f, b, l, n in zip(D.fam, D.bin, D.lab, D.nodeB)])
    D = D.copy()
    isQ = D.A.isin(QUART) & D.B.isin(QUART)
    corr = np.where(isQ, -(np.nan_to_num(kB) - np.nan_to_num(kA)), np.nan_to_num(kA))
    has = np.where(isQ, np.isfinite(kA) & np.isfinite(kB), np.isfinite(kA))
    D["kcorr"] = corr; D["has_corr"] = has
    D["rc0"] = D.r + D["kcorr"]
    D["rc"] = D.rc0 - D.groupby(PKEY).rc0.transform("mean")     # re-demean per series
    return D


def site_mean(D):
    F = D[D.B.isin(FAR) & D.A.isin(QUART)]
    g = F.groupby(["fam", "B", "epB", "bin", "lab"])
    S = g.agg(rs=("rc", "mean"), r_unc=("r", "mean"), nA=("A", "size"), sig=("sig", "median"), date=("date", "min"),
              ccX=("ccS_B", "median"), ccA=("ccS_A", "median"), n=("n", "median")).reset_index()
    S = S.rename(columns={"B": "far", "epB": "epoch"})
    for c in ("rs", "r_unc"):
        S[c] = S[c] - S.groupby(["fam", "far", "epoch", "bin"])[c].transform("mean")
    S["K"] = S.groupby(["fam", "far", "epoch", "bin"]).rs.transform("size")
    S = S[S.K >= 3]
    S["sid"] = S.groupby(["fam", "far", "epoch", "bin"]).ngroup()
    return S


def sc_of(d, col):
    if len(d) == 0:
        return np.nan, np.nan
    e = d.assign(r=d[col])
    return scatter(e)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); a = ap.parse_args(); tag = a.tag
    os.makedirs(FIG, exist_ok=True)
    rng = np.random.default_rng(SEED)
    pr = pd.read_parquet(f"{TB}/lags_pair{tag}.parquet")
    st = pd.read_parquet(f"{TB}/lags_station{tag}.parquet")
    D = prep(pr)
    T = quartet_terms(D)
    T, Kt = aggregate(T)
    D = correct(D, T)
    D.to_parquet(f"{TB}/pair_corrected{tag}.parquet", index=False)
    Kt.to_csv(f"{TB}/clock_terms{tag}.csv", index=False)
    Sm = site_mean(D); Sm.to_parquet(f"{TB}/far_sitemean{tag}.parquet", index=False)
    res = []
    # ---- floors before / after (EP headline; YR, Y2 columns)
    for bk in ("EP", "YR", "Y2"):
        Db = D[D.bin == bk]
        for ref, pl in FLOOR_SETS.items():
            d = Db[Db.pair.isin(pl)]
            dc = d[d.has_corr]
            b0, b0r = sc_of(d, "r"); b1, b1r = sc_of(dc, "r"); a1, a1r = sc_of(dc, "rc")
            r = dict(bin=bk, item=f"floor {ref}", n_values=len(d), n_corr=len(dc), sd_before_all=b0, rob_before_all=b0r,
                     sd_before=b1, rob_before=b1r, sd_after=a1, rob_after=a1r)
            # decomposition (episode-common part) before / after on the corrected subset
            for col, lab_ in (("r", "before"), ("rc", "after")):
                o, _ = decomp(dc.assign(r=dc[col]), "sig", "x", rng)
                r[f"sd_ep_{lab_}"] = np.sqrt(max(o.get("var_ep", np.nan), 0)) if np.isfinite(o.get("var_ep", np.nan)) else np.nan
                r[f"p_{lab_}"] = o.get("perm_p", np.nan)
            r["pass_rob_le_3"] = bool(a1r <= 3.0)
            res.append(r)
        # far pairs
        for ref in ("CN.PGC", "PB.B011"):
            for far in FAR:
                d = Db[Db.pair == f"{ref}-{far}"]; dc = d[d.has_corr]
                b0, b0r = sc_of(d, "r"); b1, b1r = sc_of(dc, "r"); a1, a1r = sc_of(dc, "rc")
                r = dict(bin=bk, item=f"far {ref}-{far}", n_values=len(d), n_corr=len(dc), sd_before_all=b0,
                         rob_before_all=b0r, sd_before=b1, rob_before=b1r, sd_after=a1, rob_after=a1r)
                for col, lab_ in (("r", "before"), ("rc", "after")):
                    o, _ = decomp(dc.assign(r=dc[col]), "sig", "x", rng)
                    ve = o.get("var_ep", np.nan)
                    r[f"sd_ep_{lab_}"] = np.sqrt(max(ve, 0)) if np.isfinite(ve) else np.nan
                    r[f"p_{lab_}"] = o.get("perm_p", np.nan)
                r["reduced"] = bool((a1 < b1) and (a1r < b1r))
                res.append(r)
        for far in FAR:
            d = Sm[(Sm.bin == bk) & (Sm.far == far)]
            a1, a1r = sc_of(d, "rs"); b1, b1r = sc_of(d, "r_unc")
            r = dict(bin=bk, item=f"far SITE-{far}", n_values=len(d), n_corr=len(d), sd_before=b1, rob_before=b1r,
                     sd_after=a1, rob_after=a1r)
            for col, lab_ in (("r_unc", "before"), ("rs", "after")):
                o, _ = decomp(d.assign(r=d[col], episode=d.lab), "sig", "x", rng)
                ve = o.get("var_ep", np.nan)
                r[f"sd_ep_{lab_}"] = np.sqrt(max(ve, 0)) if np.isfinite(ve) else np.nan
                r[f"p_{lab_}"] = o.get("perm_p", np.nan)
            res.append(r)
    R = pd.DataFrame(res); R.to_csv(f"{TB}/clock_floor_far{tag}.csv", index=False)
    # ---- far-station terms (site-mean corrected) per bin label: clock + near-site path, inseparable
    if not len(Sm):                                   # Amendment 4 variants: quartet only
        print(pd.DataFrame(res).round(2).to_string()); return
    Ft = Sm.groupby(["bin", "far", "lab"]).apply(lambda x: pd.Series(dict(
        F=np.average(x.rs, weights=1 / x.sig ** 2), n_fam=len(x), sd_fam=x.rs.std(ddof=1) if len(x) > 1 else np.nan,
        sign_agree=np.mean(np.sign(x.rs) == np.sign(np.average(x.rs, weights=1 / x.sig ** 2))),
        date=x.date.min(), cc=x.ccX.median()))).reset_index()
    Ft["sigF"] = Ft.sd_fam / np.sqrt(Ft.n_fam)
    Ft["flag"] = (Ft.n_fam >= 3) & (np.abs(Ft.F) >= 3 * Ft.sigF) & (Ft.sign_agree >= 0.7)
    Ft.to_csv(f"{TB}/far_terms{tag}.csv", index=False)
    # ---- detector-frame drift: station-level all-row S lags, mean over stations and families per label
    s = st[(st.kind == "data") & ~st.railS].copy(); s["episode"] = s.lab
    dd = pd.concat([demean(g, "lagS", ["fam", "netsta", "epoch", "bin"]) for _, g in s.groupby("bin")])
    Od = dd.groupby(["bin", "lab"]).agg(O=("r", "mean"), n=("r", "size"), date=("date", "min"), sd=("r", "std")).reset_index()
    Od["sigO"] = Od.sd / np.sqrt(Od.n)
    Od.to_csv(f"{TB}/frame_drift{tag}.csv", index=False)
    # ---- summary print
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(R.round(2).to_string())
        print(Kt[Kt.glitch & (Kt.bin == "EP")].sort_values("date")[["station", "epoch", "lab", "date", "K", "sigK", "n_fam",
                                                                       "sign_agree"]].round(2).to_string())
        print(Ft[Ft.flag & (Ft.bin == "EP")].sort_values("date").round(2).to_string())
    figures(Kt, Ft, R, Od, tag)


def figures(Kt, Ft, R, Od, tag):
    from matplotlib.colors import Normalize
    norm = Normalize(0.5, 1.0); cmap = plt.get_cmap("viridis")
    K = Kt[Kt.bin == "EP"].copy(); K["date"] = pd.to_datetime(K.date)
    F = Ft[Ft.bin == "EP"].copy(); F["date"] = pd.to_datetime(F.date)
    O = Od[Od.bin == "EP"].copy(); O["date"] = pd.to_datetime(O.date)
    rows = QUART + FAR + ["frame"]
    fig, axs = plt.subplots(len(rows), 1, figsize=(13, 2.0 * len(rows)), sharex=True)
    for ax, nm in zip(axs, rows):
        if nm in QUART:
            d = K[K.station == nm].sort_values("date"); y = d.K; c = d.cc; flag = d.glitch
            lab = f"{nm.split('.')[1]} instrument term\n(ms, vs quartet median)"
        elif nm in FAR:
            d = F[F.far == nm].sort_values("date"); y = d.F; c = d.cc; flag = d.flag
            lab = f"{nm.split('.')[1]} − PGC site\n(clock + near-site path, ms)"
        else:
            d = O.sort_values("date"); y = d.O; c = pd.Series(np.full(len(d), 0.9)); flag = pd.Series(np.zeros(len(d), bool))
            lab = "detector frame /\norigin drift (ms)"
        ax.plot(d.date, y, color="0.6", lw=0.8, zorder=1)
        sc = ax.scatter(d.date, y, c=c, cmap=cmap, norm=norm, s=14, lw=0, zorder=2)
        ax.scatter(d.date[flag.to_numpy()], y[flag.to_numpy()], s=60, facecolors="none", edgecolors="r", lw=1.2, zorder=3)
        ax.axhline(0, color="0.5", lw=0.5); ax.set_ylabel(lab, fontsize=8)
        lim = np.nanpercentile(np.abs(y), 99) * 1.2 if len(d) else 5
        ax.set_ylim(-max(lim, 3), max(lim, 3))
    fig.colorbar(sc, ax=axs, shrink=0.3, label="median S cc")
    fig.suptitle("Episode timing terms (round 2). Red rings = 3σ, ≥70 % of families agree in sign", fontsize=10)
    fig.savefig(f"{FIG}/fig_clock_terms{tag}.png", dpi=110, bbox_inches="tight"); plt.close(fig)
    # floor / far before-after (line per item)
    E = R[R.bin == "EP"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (cb, ca, t) in zip(axs, (("sd_before", "sd_after", "dof-corrected SD"), ("rob_before", "rob_after", "robust (MAD)"))):
        for i, r in enumerate(E.itertuples()):
            ax.plot([0, 1], [getattr(r, cb), getattr(r, ca)], marker="o", lw=1.2,
                    label=r.item.replace("CN.", "").replace("PB.", ""))
        ax.set_xticks([0, 1]); ax.set_xticklabels(["before", "after clock correction"]); ax.set_ylabel(f"scatter, {t} (ms)")
        if cb == "rob_before":
            ax.axhline(3, color="r", lw=0.8, ls="--")
        ax.grid(color="0.9")
    axs[1].legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.suptitle("Episode S-difference scatter before and after the quartet clock correction (red: 3 ms rule)", fontsize=10)
    fig.savefig(f"{FIG}/fig_floor_before_after{tag}.png", dpi=110, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
