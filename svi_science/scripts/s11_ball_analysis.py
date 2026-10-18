#!/usr/bin/env python
"""s11_ball_analysis -- S11 ballistic test analysis (ballistic/PLAN.md §3-8 + Amendments 1-2).
Demeaned per series; floor (co-located dS), far-pair and S-P variance decomposition (episode-common vs family-specific
vs bootstrap noise), within-series episode-label permutation (10,000 draws), verdicts, dbeta/beta conversion, figures.

    python s11_ball_analysis.py [--tag _test]
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic"
TB, FIG = f"{OUT}/tables", f"{OUT}/figures"
NPERM, SEED = 10000, 20261012
FAR = ["PB.B926", "CN.SNB", "CN.VGZ"]
COLOC = {"CN.PGC": ["CN.PGC-PB.B009", "CN.PGC-PB.B010", "CN.PGC-PB.B011"], "PB.B011": ["PB.B011-PB.B009", "PB.B011-PB.B010"]}
STATIONS = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011", "PB.B926", "CN.VGZ", "CN.SNB"]


def demean(df, val, keys, kmin=3):
    d = df[np.isfinite(df[val])].copy()
    d["K"] = d.groupby(keys)[val].transform("size")
    d = d[d.K >= kmin].copy()
    d["r"] = d[val] - d.groupby(keys)[val].transform("mean")
    d["sid"] = d.groupby(keys).ngroup()
    return d


def scatter(d):
    if len(d) == 0:
        return np.nan, np.nan
    dof = (d.groupby("sid").size() - 1).sum()
    sd = np.sqrt((d.r ** 2).sum() / dof)
    fac = np.sqrt(np.mean(d.K / (d.K - 1)))
    rob = 1.4826 * np.median(np.abs(d.r)) * fac
    return float(sd), float(rob)


def anova(y, g):
    """one-way: y values, g integer group labels -> MSB, MSW, n0, F"""
    G = g.max() + 1
    n = np.bincount(g, minlength=G).astype(float); s = np.bincount(g, y, minlength=G)
    mu = s / np.maximum(n, 1); N = n.sum(); gm = y.mean()
    ok = n > 0; Gk = ok.sum()
    msb = (n[ok] * (mu[ok] - gm) ** 2).sum() / (Gk - 1)
    msw = ((y - mu[g]) ** 2).sum() / (N - Gk)
    n0 = (N - (n[ok] ** 2).sum() / N) / (Gk - 1)
    return msb, msw, n0, msb / msw


def decomp(d, sigcol, label, rng, min_fams=3):
    """d: demeaned frame (r, sid, episode, sigcol) -> dict + per-episode table"""
    out = dict(item=label, n_values=len(d), n_series=d.sid.nunique())
    out["scatter_ms"], out["scatter_robust_ms"] = scatter(d)
    e_n = d.groupby("episode").size()
    keep = e_n[e_n >= min_fams].index
    q = d[d.episode.isin(keep)].reset_index(drop=True)
    out.update(n_used=len(q), n_episodes=len(keep))
    if len(keep) < 3 or len(q) < 10:
        return out, None
    g = q.episode.astype("category").cat.codes.to_numpy()
    y = q.r.to_numpy(); msb, msw, n0, F = anova(y, g)
    shrink = np.mean((q.K - 1) / q.K)
    noise = float(np.nanmean(q[sigcol] ** 2) * shrink)
    out.update(MSB=msb, MSW=msw, n0=n0, F=F, var_ep=(msb - msw) / n0, var_within=msw, var_noise=noise,
               var_fs=max(0.0, msw - noise), chi2_noise=float(np.nanmean(q.r ** 2 / (q[sigcol] ** 2 * shrink))))
    # permutation within series
    sid = q.sid.to_numpy(); order = np.argsort(sid, kind="stable"); sid_s = sid[order]
    Fs = np.empty(NPERM)
    for i in range(NPERM):
        key = sid_s + rng.random(len(sid_s))
        perm = order[np.argsort(key, kind="stable")]           # within-series shuffle of positions
        yp = np.empty_like(y); yp[order] = y[perm]
        Fs[i] = anova(yp, g)[3]
    out["perm_p"] = float((1 + (Fs >= F).sum()) / (NPERM + 1)); out["F_perm95"] = float(np.quantile(Fs, 0.95))
    out["F_perm99"] = float(np.quantile(Fs, 0.99))
    ept = q.groupby("episode").agg(date=("date", "min"), n_fam=("r", "size"), common=("r", "mean"),
                                   sd_within=("r", "std")).reset_index()
    ept["item"] = label
    return out, ept


def verdict(row, floor, floor_ep_var=None):
    sc = row.get("scatter_ms", np.nan)
    if not np.isfinite(sc):
        return "no data"
    if sc <= 2 * floor:
        return "FLAT"
    p = row.get("perm_p", 1.0); ve = row.get("var_ep", np.nan); vf = row.get("var_fs", np.nan)
    if p < 0.01 and ve > vf:
        if floor_ep_var is not None and not (ve > floor_ep_var):
            return "PATH? (not above co-located common part)"
        return "PATH"
    if p < 0.01:
        return "MIXED"
    return "WANDER"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); a = ap.parse_args()
    tag = a.tag
    rng = np.random.default_rng(SEED)
    os.makedirs(FIG, exist_ok=True)
    st = pd.read_parquet(f"{TB}/lags_station{tag}.parquet")
    pr = pd.read_parquet(f"{TB}/lags_pair{tag}.parquet")
    geo = pd.read_csv(f"{TB}/series_geom{tag}.csv")
    spmap = geo.groupby(["fam", "netsta"]).SmP_pick.median()
    res, eps_tabs = [], []
    summ = {}
    # ---- basic quality
    for kind in ("data", "ctrl"):
        s = st[st.kind == kind]; p = pr[pr.kind == kind]
        summ[f"{kind}_station_rows"] = int(len(s)); summ[f"{kind}_pair_rows"] = int(len(p))
        summ[f"{kind}_railS_frac"] = float(s.railS.mean()); summ[f"{kind}_railP_frac"] = float(s.railP.mean())
        summ[f"{kind}_ccS_med"] = float(s.ccS.median()); summ[f"{kind}_ccP_med"] = float(s.ccP.median())
    d = st[st.kind == "data"]
    summ["inj5S_median"] = float(d.inj5S.median()); summ["inj5S_p05_p95"] = [float(d.inj5S.quantile(.05)), float(d.inj5S.quantile(.95))]
    q = d.groupby("netsta").agg(n=("lagS", "size"), ccS=("ccS", "median"), sigS=("sigS", "median"), ccP=("ccP", "median"),
                                ccPz=("ccPz", "median"), sigP=("sigP", "median"), sigPz=("sigPz", "median"),
                                sigSmP=("sigSmP", "median"), railS=("railS", "mean"), railP=("railP", "mean"),
                                inj5S=("inj5S", "median"))
    q.to_csv(f"{TB}/station_quality{tag}.csv")
    # ---- values
    pr = pr.copy(); pr["dSmP"] = pr.dS - pr.dP; pr["sig_dSmP"] = np.hypot(pr.sig_dS, pr.sig_dP)
    pr.loc[pr.railS, "dS"] = np.nan; pr.loc[pr.railS | pr.railP, "dSmP"] = np.nan
    st = st.copy(); st["SmP"] = st.lagS - st.lagP
    st.loc[st.railS | st.railP, "SmP"] = np.nan; st.loc[st.railS, "lagS"] = np.nan
    PK = ["fam", "pair", "epA", "epB"]; SK = ["fam", "netsta", "epoch"]
    floors = {}
    allfloor_SmP = []
    for subset, cls_ok in (("all", None), ("majors", "major")):
        for kind in ("data", "ctrl"):
            P_ = pr[(pr.kind == kind) & ((pr.cls == cls_ok) if cls_ok else True)]
            S_ = st[(st.kind == kind) & ((st.cls == cls_ok) if cls_ok else True)]
            # floors per reference
            for ref, pl in COLOC.items():
                dd = demean(P_[P_.pair.isin(pl)], "dS", PK)
                o, ept = decomp(dd, "sig_dS", f"floor dS {ref}-coloc", rng)
                o.update(subset=subset, kind=kind, ref=ref, kind_item="floor", obs="dS")
                res.append(o); floors[(subset, kind, ref)] = o
                if ept is not None:
                    ept["subset"], ept["kind"] = subset, kind; eps_tabs.append(ept)
                for pp in pl:                                         # per pair and per PGC epoch
                    for epA in sorted(P_[P_.pair == pp].epA.unique()):
                        dd1 = demean(P_[(P_.pair == pp) & (P_.epA == epA)], "dS", PK)
                        sc = scatter(dd1)
                        res.append(dict(item=f"floor dS {pp} {epA}", subset=subset, kind=kind, ref=ref, kind_item="floor_part",
                                        obs="dS", n_values=len(dd1), scatter_ms=sc[0], scatter_robust_ms=sc[1],
                                        sig_med=float(dd1.sig_dS.median()) if len(dd1) else np.nan))
            ddS = demean(P_[P_.pair.isin(sum(COLOC.values(), []))], "dSmP", PK)
            o, _ = decomp(ddS, "sig_dSmP", "floor d(S-P) coloc", rng)
            o.update(subset=subset, kind=kind, ref="both", kind_item="floor", obs="d(S-P)"); res.append(o)
            floors[(subset, kind, "SmP")] = o
            # far pairs
            for ref in ("CN.PGC", "PB.B011"):
                fl = floors[(subset, kind, ref)]
                fl_sc = fl.get("scatter_ms", np.nan); fl_ve = fl.get("var_ep", np.nan)
                for far in FAR:
                    pp = f"{ref}-{far}"
                    dd = demean(P_[P_.pair == pp], "dS", PK)
                    o, ept = decomp(dd, "sig_dS", f"dS {pp}", rng)
                    # conversion: common SD -> dbeta/beta
                    fams = dd.fam.unique()
                    dTS = np.nanmean([(spmap.get((f, far), np.nan) - spmap.get((f, ref), np.nan)) * 1.75 / 0.75 for f in fams]) if len(fams) else np.nan
                    TSX = np.nanmean([spmap.get((f, far), np.nan) * 1.75 / 0.75 for f in fams]) if len(fams) else np.nan
                    ve = o.get("var_ep", np.nan); sd_ep = np.sqrt(max(ve, 0)) if np.isfinite(ve) else np.nan
                    o.update(subset=subset, kind=kind, ref=ref, kind_item="far", obs="dS", floor_ms=fl_sc,
                             ratio_to_floor=o.get("scatter_ms", np.nan) / fl_sc, dTS_s=dTS, TS_X_s=TSX,
                             sd_ep_ms=sd_ep, dbb_ep_uniform_pct=sd_ep / 1e3 / abs(dTS) * 100 if np.isfinite(dTS) else np.nan,
                             dbb_ep_Xleg_pct=sd_ep / 1e3 / TSX * 100, dbb_floor_Xleg_pct=fl_sc / 1e3 / TSX * 100,
                             sig_med=float(dd.sig_dS.median()) if len(dd) else np.nan)
                    o["verdict"] = verdict(o, fl_sc, fl_ve)
                    res.append(o)
                    if ept is not None:
                        ept["subset"], ept["kind"] = subset, kind; eps_tabs.append(ept)
            # S-P per station
            fl_smp = floors[(subset, kind, "SmP")].get("scatter_ms", np.nan)
            for ns in STATIONS:
                dd = demean(S_[S_.netsta == ns], "SmP", SK)
                o, ept = decomp(dd, "sigSmP", f"S-P {ns}", rng)
                SPm = np.nanmean([spmap.get((f, ns), np.nan) for f in dd.fam.unique()]) if len(dd) else np.nan
                ve = o.get("var_ep", np.nan); sd_ep = np.sqrt(max(ve, 0)) if np.isfinite(ve) else np.nan
                o.update(subset=subset, kind=kind, ref=ns, kind_item="SmP", obs="S-P", floor_ms=fl_smp,
                         ratio_to_floor=o.get("scatter_ms", np.nan) / fl_smp, SP_s=SPm, sd_ep_ms=sd_ep,
                         dbb_ep_constVpVs_pct=sd_ep / 1e3 / SPm * 100, dbb_ep_VsOnly_pct=sd_ep / 1e3 / (SPm * 1.75 / 0.75) * 100,
                         dbb_scatter_constVpVs_pct=o.get("scatter_ms", np.nan) / 1e3 / SPm * 100,
                         sig_med=float(dd.sigSmP.median()) if len(dd) else np.nan)
                o["verdict"] = verdict(o, fl_smp)
                res.append(o)
                if ept is not None:
                    ept["subset"], ept["kind"] = subset, kind; eps_tabs.append(ept)
                # single-station S lag (contains the catalogue-time drift) for context
                dd2 = demean(S_[S_.netsta == ns], "lagS", SK)
                o2, _ = decomp(dd2, "sigS", f"S lag {ns} (with catalogue drift)", rng)
                o2.update(subset=subset, kind=kind, ref=ns, kind_item="Slag_context", obs="S lag"); res.append(o2)
    R = pd.DataFrame(res)
    for c in ("var_ep", "var_fs", "var_noise", "var_within"):
        if c in R:
            R["sd_" + c[4:]] = np.sqrt(R[c].clip(lower=0))
    R.to_csv(f"{TB}/decomposition{tag}.csv", index=False)
    E = pd.concat(eps_tabs, ignore_index=True) if eps_tabs else pd.DataFrame()
    E.to_csv(f"{TB}/episode_common{tag}.csv", index=False)
    json.dump(summ, open(f"{TB}/summary{tag}.json", "w"), indent=1, default=float)
    cols = ["item", "subset", "kind", "n_values", "n_series", "n_episodes", "scatter_ms", "scatter_robust_ms", "sd_ep",
            "sd_fs", "sd_noise", "F", "perm_p", "ratio_to_floor", "verdict"]
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(R[(R.kind_item != "floor_part")][[c for c in cols if c in R]].round(3).to_string())
        print(R[R.kind_item == "floor_part"][["item", "subset", "kind", "n_values", "scatter_ms", "scatter_robust_ms", "sig_med"]].round(2).to_string())
        print(q.round(2)); print(summ)
    figures(pr, st, R, tag)


def figures(pr, st, R, tag):
    from matplotlib.colors import Normalize
    norm = Normalize(0.3, 1.0); cmap = plt.get_cmap("viridis")
    # ---- dS series
    items = [f"{r}-{f}" for r in ("CN.PGC", "PB.B011") for f in FAR] + sum(COLOC.values(), [])
    fig, axs = plt.subplots(len(items), 1, figsize=(13, 2.3 * len(items)), sharex=True)
    fig.subplots_adjust(top=0.97, hspace=0.35)
    P_ = pr[pr.kind == "data"]
    for ax, pp in zip(axs, items):
        dd = demean(P_[P_.pair == pp], "dS", ["fam", "pair", "epA", "epB"])
        dd["date"] = pd.to_datetime(dd.date); dd["cc"] = np.minimum(dd.ccS_A, dd.ccS_B)
        for sid, g in dd.groupby("sid"):
            g = g.sort_values("date"); ax.plot(g.date, g.r, color="0.8", lw=0.5, zorder=1)
        sc = ax.scatter(dd.date, dd.r, c=dd.cc, cmap=cmap, norm=norm, s=9, lw=0, zorder=2)
        ec = dd.groupby("episode").agg(date=("date", "min"), m=("r", "mean"), n=("r", "size"))
        ec = ec[ec.n >= 3].sort_values("date")
        ax.plot(ec.date, ec.m, color="k", lw=1.6, zorder=3)
        lim = np.nanpercentile(np.abs(dd.r), 98) * 1.2 if len(dd) else 10
        ax.set_ylim(-lim, lim); ax.axhline(0, color="0.5", lw=0.5)
        row = R[(R.item == f"dS {pp}") & (R.subset == "all") & (R.kind == "data")]
        txt = ""
        if len(row):
            r0 = row.iloc[0]; txt = f"  scatter {r0.scatter_ms:.1f} ms, common SD {r0.get('sd_ep', np.nan):.1f} ms, p {r0.get('perm_p', np.nan):.4f}"
        ax.set_ylabel(f"{pp.replace('CN.', '').replace('PB.', '')}\nΔS (ms)", fontsize=8)
        ax.set_title(f"{pp}: S-time difference per episode (demeaned per family series){txt}", fontsize=8, loc="left")
    fig.colorbar(sc, ax=axs, shrink=0.3, label="S cc with reference (lower of the two stations)")
    fig.suptitle("Inter-station S differences per episode: grey = one family, black = episode-common mean (≥3 families)", fontsize=10)
    fig.savefig(f"{FIG}/fig_dS_series{tag}.png", dpi=110, bbox_inches="tight"); plt.close(fig)
    # ---- S-P series
    fig, axs = plt.subplots(len(STATIONS), 1, figsize=(13, 2.3 * len(STATIONS)), sharex=True)
    fig.subplots_adjust(top=0.94, hspace=0.35)
    S_ = st[st.kind == "data"]
    for ax, ns in zip(axs, STATIONS):
        dd = demean(S_[S_.netsta == ns], "SmP", ["fam", "netsta", "epoch"]); dd["date"] = pd.to_datetime(dd.date)
        for sid, g in dd.groupby("sid"):
            g = g.sort_values("date"); ax.plot(g.date, g.r, color="0.8", lw=0.5, zorder=1)
        sc = ax.scatter(dd.date, dd.r, c=dd.ccP, cmap=cmap, norm=norm, s=9, lw=0, zorder=2)
        ec = dd.groupby("episode").agg(date=("date", "min"), m=("r", "mean"), n=("r", "size"))
        ec = ec[ec.n >= 3].sort_values("date"); ax.plot(ec.date, ec.m, color="k", lw=1.6, zorder=3)
        lim = np.nanpercentile(np.abs(dd.r), 98) * 1.2 if len(dd) else 10
        ax.set_ylim(-lim, lim); ax.axhline(0, color="0.5", lw=0.5)
        row = R[(R.item == f"S-P {ns}") & (R.subset == "all") & (R.kind == "data")]
        txt = ""
        if len(row):
            r0 = row.iloc[0]; txt = f"  scatter {r0.scatter_ms:.1f} ms, common SD {r0.get('sd_ep', np.nan):.1f} ms, p {r0.get('perm_p', np.nan):.4f}"
        ax.set_ylabel(f"{ns.split('.')[1]}\nS−P (ms)", fontsize=8)
        ax.set_title(f"{ns}: S−P change per episode{txt}", fontsize=8, loc="left")
    fig.colorbar(sc, ax=axs, shrink=0.3, label="P cc with reference")
    fig.suptitle("S−P per episode: grey = one family, black = episode-common mean (≥3 families)", fontsize=10)
    fig.savefig(f"{FIG}/fig_SmP_series{tag}.png", dpi=110, bbox_inches="tight"); plt.close(fig)
    # ---- variance decomposition bars (SD equivalents, linear)
    cols = {"sd_ep": ("episode-common", "#2a6fdb"), "sd_fs": ("family-specific", "#e07b39"), "sd_noise": ("bootstrap noise", "#8c8c8c")}
    T0 = R[(R.subset == "all") & (R.kind == "data") & R.kind_item.isin(["floor", "far", "SmP"])]
    order = list(T0["item"])
    fig, axs = plt.subplots(1, 2, figsize=(14, 7.5))
    for ax, kind in zip(axs, ("data", "ctrl")):
        T = R[(R.subset == "all") & (R.kind == kind) & R.kind_item.isin(["floor", "far", "SmP"])].set_index("item").reindex(order)
        y = np.arange(len(T)); h = 0.26
        for k, (c, (lab, col)) in enumerate(cols.items()):
            ax.barh(y + (k - 1) * h, T[c].fillna(0), height=h - 0.03, color=col, label=lab)
        ax.set_yticks(y); ax.set_yticklabels([i.replace("CN.", "").replace("PB.", "") for i in order], fontsize=8)
        ax.set_ylim(len(T) - 0.5, -0.5)
        xm = np.nanmax(T[list(cols)].to_numpy()) if np.isfinite(T[list(cols)].to_numpy()).any() else 1
        ax.set_xlim(0, xm * 1.35); ax.set_xlabel("SD of component (ms)"); ax.grid(axis="x", color="0.9")
        ax.set_title("data" if kind == "data" else "scrambled control (few values; most railed)", fontsize=10)
        for yi, (_, r) in zip(y, T.iterrows()):
            txt = f"p={r.perm_p:.4f}" if np.isfinite(r.get("perm_p", np.nan)) else "too few episodes"
            ax.text(xm * 1.33, yi, txt, fontsize=7, color="0.3", ha="right", va="center")
    h_, l_ = axs[0].get_legend_handles_labels(); fig.legend(h_, l_, loc="lower center", ncol=3, fontsize=9)
    fig.subplots_adjust(wspace=0.35, bottom=0.12)
    fig.suptitle("Variance decomposition (all episode units, n ≥ 20): episode-common vs family-specific vs bootstrap noise; "
                 "p = episode-label permutation", fontsize=10)
    fig.savefig(f"{FIG}/fig_variance_decomp{tag}.png", dpi=110, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
