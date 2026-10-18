#!/usr/bin/env python
"""s10_tide_figures -- figures for the S10 tide / sub-daily follow-up. dv/v plots: points coloured by CC with the
reference, marker = station, no error bars (user rule). Times in PT.
    python s10_tide_figures.py [--tag _test]
"""
from __future__ import annotations
import os, sys, argparse, time
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_tide_common as T

MK = ["o", "s", "^", "v", "D", "P", "X", "*", "<", ">", "h", "p", "d", "8", "H", "1", "2", "3", "4", "+", "x"]
CMAP = plt.get_cmap("viridis")
NORM = matplotlib.colors.Normalize(0.2, 1.0)
DES = {"D1": "D1 within-day differential", "D2": "D2 day-matched pooled", "D3": "D3 pooled vs epoch stack"}


def markers(stations):
    return {s: MK[i % len(MK)] for i, s in enumerate(sorted(stations))}


def lim(v, pooled):
    v = np.abs(v[np.isfinite(v)])
    if not len(v):
        return 0.5
    return max(np.nanpercentile(v, 95) * 1.15, np.nanmax(np.abs(pooled)) * 1.3 if np.isfinite(pooled).any() else 0, 0.05)


def phase_fig(B, RU, fs, win, tag, mk):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True)
    for j, des in enumerate(("D1", "D2", "D3")):
        for i, var in enumerate(("th_udss", "udss")):
            a = ax[i, j]
            g = B[(B.famset == fs) & (B["var"] == var) & (B.design == des) & (B.win == win)]
            st = g[g.level != "pooled"]; po = g[g.level == "pooled"].sort_values("bin")
            xcol = "bin_center" if var == "th_udss" else "x_mean"
            if var == "udss" and fs == "all":
                a.text(0.5, 0.5, "value bins differ by family\n(see per-family figures)", ha="center", va="center", transform=a.transAxes)
                continue
            for ns, gg in st.groupby("level"):
                a.scatter(gg[xcol], gg.dvv_pct, c=gg.cc, cmap=CMAP, norm=NORM, marker=mk.get(ns, "o"), s=28, alpha=0.75,
                          linewidths=0.6 if mk.get(ns, "o") in "1234+x" else 0.2, edgecolors="none" if mk.get(ns, "o") in "1234+x" else "k")
            a.scatter(po[xcol], po.dvv_pct, c=po.cc, cmap=CMAP, norm=NORM, marker="o", s=150, edgecolors="k", linewidths=1.6, zorder=5)
            r = RU[(RU.famset == fs) & (RU["var"] == var) & (RU.design == des) & (RU["mode"] == "rep")]
            if var == "th_udss" and len(r):
                ww = po.w.to_numpy(); e = -po.dvv_pct.to_numpy() / 100
                xx = np.linspace(-22.5, 337.5, 200)
                # pooled cos fit (eps); dv/v = -eps
                X = np.c_[np.ones(8), np.cos(np.radians(po.bin_center)), np.sin(np.radians(po.bin_center))]
                okb = np.isfinite(e) & (ww > 0)
                if okb.sum() >= 4:
                    co = np.linalg.lstsq(X[okb] * np.sqrt(ww[okb, None]), e[okb] * np.sqrt(ww[okb]), rcond=None)[0]
                    a.plot(xx, -100 * (co[0] + co[1] * np.cos(np.radians(xx)) + co[2] * np.sin(np.radians(xx))), "k-", lw=1.5)
                rr = r.iloc[0]
                a.set_title(f"{DES[des]}\nA = {rr.A_eps_pct:.3f} %, control 99th pct {rr.A_ctrl99_pct:.3f} %, p = {rr.p_A:.3f}", fontsize=9)
                a.set_xlabel("M2 UDSS phase (deg; 0 = peak up-dip shear, −90 = peak rate)")
                a.set_xticks(np.arange(0, 360, 45))
            elif len(r):
                rr = r.iloc[0]
                a.set_title(f"{DES[des]}\ndv/v slope {rr.slope_dvv_pct_per_unit:+.4f} %/kPa, control SD {rr.slope_ctrl_sd_pct:.4f}, p = {rr.p_slope:.3f}", fontsize=9)
                a.set_xlabel("UDSS at repeat time (kPa), 8 quantile bins")
            vals = st.dvv_pct.to_numpy(); L = lim(vals, po.dvv_pct.to_numpy())
            a.set_ylim(-L, L); a.axhline(0, color="0.6", lw=0.8)
            nclip = int(np.sum(np.abs(vals[np.isfinite(vals)]) > L))
            if nclip:
                a.text(0.99, 0.02, f"{nclip} station points outside", ha="right", fontsize=7, transform=a.transAxes, color="0.4")
            a.set_ylabel(f"dv/v (%), {win[0]}–{win[1]} s")
    h = [Line2D([], [], marker=m, ls="", color="0.3", label=s) for s, m in mk.items()]
    h.append(Line2D([], [], marker="o", ls="", mfc="w", mec="k", mew=1.6, ms=11, label="pooled (all stations)"))
    fig.legend(handles=h, loc="outside right center", fontsize=7, title="station")
    fig.colorbar(matplotlib.cm.ScalarMappable(NORM, CMAP), ax=ax, label="CC with reference", shrink=0.6, location="bottom")
    fl = "families 10+87+17" if fs == "all" else f"family {fs}"
    fig.suptitle(f"S10 tidal-phase dv/v, {fl}, window {win[0]}–{win[1]} s (per-station bins and pooled; no error bars — "
                 "uncertainties in RESULTS tables)", fontsize=10)
    fig.savefig(f"{T.FIG}/s10_tide_phase_{'all' if fs == 'all' else 'fam' + fs}_{win}{tag}.png", dpi=130)
    plt.close(fig)


def control_fig(RU, tag):
    """observed amplitude vs the random-shift control, tidal and fake cycles (not a dv/v plot)."""
    vars_ = ["th_udss", "th_sn", "th_fake12", "th_fake13.5", "th_fake10.5", "c2:th_udss"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for j, des in enumerate(("D1", "D2", "D3")):
        a = ax[j]
        for i, fs in enumerate(("10", "87", "17", "all")):
            for k, v in enumerate(vars_):
                r = RU[(RU.famset == fs) & (RU["var"] == v) & (RU.design == des) & (RU["mode"] == "rep")]
                if not len(r):
                    continue
                r = r.iloc[0]; x = k + (i - 1.5) * 0.18
                a.plot([x, x], [r.A_ctrl_med_pct, r.A_ctrl99_pct], color="0.75", lw=5, solid_capstyle="butt")
                a.plot(x, r.A_eps_pct, marker=["o", "s", "^", "D"][i], color=["C0", "C1", "C2", "k"][i], ls="",
                       ms=7, label=("family " + fs if fs != "all" else "pooled families") if k == 0 else None)
        a.set_xticks(range(len(vars_))); a.set_xticklabels(["M2 UDSS", "M2 σn", "fake 12.0 h", "fake 13.5 h", "fake 10.5 h",
                                                             "scrambled\ncontrol"], fontsize=8)
        a.set_ylabel("fitted cosine amplitude |dv/v| (%)"); a.set_title(f"{DES[des]}, 2–4 s\ngrey bar = control median → 99th pct", fontsize=9)
        a.set_yscale("log")
    ax[0].legend(fontsize=8)
    fig.savefig(f"{T.FIG}/s10_tide_controls{tag}.png", dpi=130); plt.close(fig)


def subday_fig(B, FR, fam, win, tag, mk):
    eps_ = ["2004-07", "2005-09", "2003-02", "2013-09", "2008-05", "2011-08", "2024-08", "2023-04", "2019-08"]
    g = B[(B.fam == fam) & (B.res == "6h")].copy()
    fig, ax = plt.subplots(3, 3, figsize=(16, 11), constrained_layout=True)
    for a, ep in zip(ax.ravel(), eps_):
        ge = g[g.episode == ep].copy()
        ge = ge[ge.groupby(["netsta", "epoch"]).eps24.transform("count") >= 2]     # one-bin series carry no within-episode change
        fr = FR[(FR.fam == fam) & (FR.episode == ep)].iloc[0]
        pt = lambda x: pd.Timestamp(x).tz_localize("UTC").tz_convert("America/Los_Angeles").tz_localize(None)
        if len(ge):
            t_all = pd.to_datetime(ge.bin_mid).tolist() + ([pd.Timestamp(fr.front)] if pd.notna(fr.front) else [])
            a.set_xlim(pt(min(t_all) - pd.Timedelta("8h")), pt(max(t_all) + pd.Timedelta("8h")))
        elif pd.notna(fr.front):
            a.set_xlim(pt(pd.Timestamp(fr.front) - pd.Timedelta("2D")), pt(pd.Timestamp(fr.front) + pd.Timedelta("2D")))
        if len(ge):
            # per series demeaned (CC^2-weighted episode mean)
            def dm(x):
                e = x[f"eps{win}"]; w = x[f"cc{win}"] ** 2; ok = e.notna()
                return e - (np.sum(w[ok] * e[ok]) / np.sum(w[ok]) if ok.any() else np.nan)
            ge["d"] = -100 * ge.groupby(["netsta", "epoch"], group_keys=False).apply(dm)
            ge["tp"] = pd.to_datetime(ge.bin_mid).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.tz_localize(None)
            for ns, gg in ge.groupby("netsta"):
                a.scatter(gg.tp, gg.d, c=gg[f"cc{win}"], cmap=CMAP, norm=NORM, marker=mk.get(ns, "o"), s=30,
                          edgecolors="none" if mk.get(ns, "o") in "1234+x" else "k", linewidths=0.6 if mk.get(ns, "o") in "1234+x" else 0.2)
            L = lim(ge.d.to_numpy(), np.array([np.nan])); a.set_ylim(-L, L)
        tf = pd.Timestamp(fr.front).tz_localize("UTC").tz_convert("America/Los_Angeles").tz_localize(None) if pd.notna(fr.front) else None
        if tf is not None:
            a.axvline(tf, color="crimson", lw=1.5, ls="--")
            a.text(tf, a.get_ylim()[1] * 0.92, " front" + (" (S3 tremor)" if "S3" in fr.front_src else " (family onset)"),
                   color="crimson", fontsize=7)
        a.axhline(0, color="0.6", lw=0.8)
        n = len(ge); ns_ = ge.netsta.nunique() if n else 0
        a.set_title(f"{ep}: {n} 6-h bins (≥ 20 repeats) in series with ≥ 2 bins, {ns_} stations", fontsize=9)
        a.set_ylabel(f"dv/v − series episode mean (%), {win[0]}–{win[1]} s", fontsize=8)
        a.tick_params(axis="x", labelsize=7, rotation=30)
        a.set_xlabel("time (PT)", fontsize=8)
    h = [Line2D([], [], marker=m, ls="", color="0.3", label=s) for s, m in mk.items() if s in set(g.netsta)]
    fig.legend(handles=h, loc="outside right center", fontsize=7, title="station")
    fig.colorbar(matplotlib.cm.ScalarMappable(NORM, CMAP), ax=ax, label="CC with reference", shrink=0.5, location="bottom")
    fig.suptitle(f"S10 sub-daily (6-h) dv/v inside the biggest episodes, family {fam}, {win[0]}–{win[1]} s; reference = "
                 "leave-bin-out epoch stack; no error bars", fontsize=10)
    fig.savefig(f"{T.FIG}/s10_subday_fam{fam}_{win}{tag}.png", dpi=120); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); a = ap.parse_args()
    T.ensure()
    B = pd.read_csv(f"{T.TB}/phase_bins{a.tag}.csv", dtype={"famset": str, "win": str})
    RU = pd.read_csv(f"{T.TB}/phase_rules{a.tag}.csv", dtype={"famset": str})
    SB = pd.read_csv(f"{T.TB}/subday_bins.csv")
    FR = pd.read_csv(f"{T.TB}/subday_fronts.csv")
    mk = markers(set(B.level[B.level != "pooled"]) | set(SB.netsta))
    for fs in ("10", "87", "17", "all"):
        if (B.famset == fs).any():
            for win in ("24", "35"):
                phase_fig(B, RU, fs, win, a.tag, mk)
    control_fig(RU, a.tag)
    for fam in (10, 87, 17):
        for win in ("24", "35"):
            subday_fig(SB, FR, fam, win, a.tag, mk)
    if not a.tag:
        open(f"{T.OUT}/FIGURES_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
