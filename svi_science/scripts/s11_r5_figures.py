#!/usr/bin/env python
"""s11_r5_figures -- S11 round 5 figures (PLAN.md §9). dv/v points coloured by weight (1/σ²), no error bars, line graphs.
  fig_r5_forms.png : σ per major and the prior shift per form (both regions, NOM / L-SEL), the 0.1σ line
  fig_r5_dD.png    : ΔD per episode, selected form (+ adopted smooth Sc), L-SEL, both regions
  fig_r5_q1.png    : WA / SVI yearly cycle: per-station kernel-free amplitude 2-3.5 vs 3.5-6 s by group; Q1 band over kernels
  fig_r5_pgc.png   : 2019-11 per instrument: apparent episode values 2016-2023 and the lapse profile vs shapes
  fig_r5_band.png  : Q3 composite over fit kernels (band), both regions
    python s11_r5_figures.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_r5_common as R5

R = R5.R5
FIG = f"{R}/figures"; os.makedirs(FIG, exist_ok=True)
SUBS = {"SVI": "svi", "WA_NW": "wa_nw"}
CAT = ["#0072B2", "#D55E00", "#009E73"]          # validated (light): fixed order
INK, MUTED = "#222222", "#888888"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                     "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6})


def tb(reg, f):
    return f"{R}/{SUBS[reg]}/tables/{f}"


def fig_forms():
    fig, ax = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True)
    for j, reg in enumerate(SUBS):
        p = tb(reg, "r5_forms.csv")
        if not os.path.exists(p):
            continue
        D = pd.read_csv(p)
        forms = list(dict.fromkeys(D.form))
        for k, setup in enumerate(("ell40", "two3_20")):
            d = D[D.setup == setup].set_index("form").reindex(forms)
            ax[0, j].plot(forms, d.sig_major_median, "-o", color=CAT[k], lw=2, ms=6, label=f"{setup} ({'NOM' if setup == 'ell40' else 'L-SEL'})")
            ax[1, j].plot(forms, d.prior_max, "-o", color=CAT[k], lw=2, ms=6, label=setup)
        ax[0, j].axhline(0.10, color=MUTED, lw=1, ls="--"); ax[1, j].axhline(0.10, color=MUTED, lw=1, ls="--")
        ax[0, j].set_title(f"{reg}: σ per major ΔD (%)"); ax[1, j].set_title(f"{reg}: largest shift when a prior is added (σ)")
        ax[1, j].set_yscale("log"); ax[1, j].tick_params(axis="x", rotation=30)
        ax[0, j].legend(frameon=False)
    fig.suptitle("Round 5 task 1: one shallow term per unit, exact basis (rule: shift < 0.1σ, dashed)")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r5_forms.png", dpi=140); plt.close(fig)


def fig_dD():
    sel = json.load(open(f"{R}/selected.json"))["selected"] if os.path.exists(f"{R}/selected.json") else None
    ad = json.load(open(f"{R}/sc_adoption.json")).get("adopted_knot_yr") if os.path.exists(f"{R}/sc_adoption.json") else None
    fig, ax = plt.subplots(2, 1, figsize=(11, 6.5))
    for j, reg in enumerate(SUBS):
        p = tb(reg, "r5_forms_dD.csv")
        if sel is None or not os.path.exists(p):
            continue
        T = pd.read_csv(p); T = T[(T.form == sel) & (T.setup == "two3_20")].copy(); T["mid"] = pd.to_datetime(T.mid)
        T = T.sort_values("mid"); w = 1 / T.sig ** 2
        ax[j].plot(T.mid, T.dD, "-", color=MUTED, lw=0.8, zorder=1)
        sc = ax[j].scatter(T.mid, T.dD, c=w, cmap="Blues", vmin=0, s=np.where(T.major, 46, 18), edgecolor=INK, linewidth=0.3, zorder=2)
        if ad is not None and os.path.exists(tb(reg, "r5_sc_dD.csv")):
            S = pd.read_csv(tb(reg, "r5_sc_dD.csv")); S = S[(S.form == f"{sel}+Scs{int(ad)}") & (S.setup == "two3_20")].copy()
            S["mid"] = pd.to_datetime(S.mid); S = S.sort_values("mid")
            ax[j].plot(S.mid, S.dD, "-", color=CAT[1], lw=1.2, label=f"+ smooth regional shallow ({int(ad)}-yr knots)")
            ax[j].legend(frameon=False, loc="upper left")
        ax[j].axhline(0, color=MUTED, lw=0.8)
        cb = fig.colorbar(sc, ax=ax[j], pad=0.01); cb.set_label("weight 1/σ²")
        ax[j].set_ylabel("ΔD (%)"); ax[j].set_title(f"{reg}: deep-weighted episode change, form {sel}, two3_20 (large dots = majors)")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r5_dD.png", dpi=140); plt.close(fig)


def fig_q1():
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    names = {"P": "PGC site", "B": "other boreholes", "O": "other stations"}
    for j, reg in enumerate(SUBS):
        hs = "b" if reg == "WA_NW" else "a"
        p = tb(reg, f"q1_stations_{hs}_two3_20.csv")
        if os.path.exists(p):
            S = pd.read_csv(p)
            for k, g in enumerate(("P", "B", "O")):
                s = S[S.grp == g]
                w = 1 / (s.A2035_sig ** 2 + s.A3560_sig ** 2)
                ax[0, j].scatter(s.A2035, s.A3560, s=20 + 200 * w / max(w.max(), 1e-9), color=CAT[k], alpha=0.8, label=names[g], edgecolor="white", linewidth=0.5)
                for r in s.itertuples():
                    if g == "P":
                        ax[0, j].annotate(r.netsta.split(".")[1], (r.A2035, r.A3560), fontsize=7, color=INK, xytext=(3, 3), textcoords="offset points")
            lim = np.nanmax(np.abs(np.r_[S.A2035, S.A3560])) * 1.1 if len(S) else 0.3
            xx = np.array([-lim, lim])
            ax[0, j].plot(xx, xx, color=MUTED, lw=1, ls="--", label="ratio 1 (receiver)")
            ax[0, j].plot(xx, xx * float(np.nanmedian(S.a_d_ratio)), color=INK, lw=1, label="deep (median a_d ratio)")
            ax[0, j].plot(xx, xx * float(np.nanmedian(S.a_s_ratio)), color=CAT[1], lw=1, ls=":", label="shallow (median a_s ratio)")
            ax[0, j].set_xlabel("in-phase annual amplitude 2-3.5 s (%)"); ax[0, j].set_ylabel("3.5-6 s (%)")
            ax[0, j].set_title(f"{reg} (set {hs}): per-station yearly cycle, kernel-free (size = weight)"); ax[0, j].legend(frameon=False, fontsize=7)
        p = tb(reg, "q1_band.csv")
        if os.path.exists(p):
            B = pd.read_csv(p); B = B[B.set == hs]
            ax[1, j].plot(B.setup, B.Q1_r4, "-o", color=CAT[0], lw=2, label="Q1, round-4 form")
            ax[1, j].plot(B.setup, B.Q1_R, "-o", color=CAT[1], lw=2, label="Q1-R (per-station receiver term)")
            for r in B[B.plausible].itertuples():
                ax[1, j].axvspan(list(B.setup).index(r.setup) - 0.4, list(B.setup).index(r.setup) + 0.4, color="#eef3f8", zorder=0)
            ax[1, j].set_ylabel("deep-weighted yearly amplitude (%)"); ax[1, j].set_title(f"{reg}: Q1 over fit kernels (shaded = plausible)")
            ax[1, j].legend(frameon=False); ax[1, j].tick_params(axis="x", rotation=30)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r5_q1.png", dpi=140); plt.close(fig)


def fig_pgc():
    p = tb("SVI", "pgc_profile.csv")
    if not os.path.exists(p):
        return
    import s11_r4_combined as CB
    import s11_r4_common as R4
    H, _ = R5.load_H(lambda m: None)
    PR = pd.read_csv(p)
    inst = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
    fig, ax = plt.subplots(2, 4, figsize=(15, 7), sharey="row")
    for i, ns in enumerate(inst):
        c = H[(H.netsta == ns) & ~H.far].copy()
        c["dm"] = c.d - c.groupby("series").d.transform("mean")
        c["w"] = 1 / c.sig ** 2
        e = c.groupby("episode").apply(lambda g: pd.Series(dict(v=np.average(g.dm, weights=g.w), w=g.w.sum(), cc=g.cc.mean())), include_groups=False).reset_index()
        e["mid"] = e.episode.map(R4.U.mid); e = e.sort_values("mid"); e = e[(e.mid >= "2015-06-01") & (e.mid <= "2024-01-01")]
        ax[0, i].plot(e.mid, e.v, color=MUTED, lw=0.8, zorder=1)
        sc = ax[0, i].scatter(e.mid, e.v, c=e.cc, cmap="Blues", vmin=0.3, vmax=1.0, s=26, edgecolor=INK, linewidth=0.3, zorder=2)
        x = e[e.episode == "2019-11"]
        ax[0, i].scatter(x.mid, x.v, s=120, facecolor="none", edgecolor=CAT[1], linewidth=1.5, zorder=3)
        ax[0, i].axhline(0, color=MUTED, lw=0.8); ax[0, i].set_title(f"{ns}: episode dv/v − series mean\n(ring = 2019-11)")
        q = PR[(PR.netsta == ns) & PR.win.astype(str).isin(["23", "34", "45", "56"])].sort_values("t0")
        tc = (q.t0 + q.t1) / 2
        ax[1, i].plot(tc, q.dE, "-o", color=INK, lw=2, ms=7, label="data ΔE(2019-11)")
        for k, (col, lab) in enumerate((("a_d_two3_20", "deep shape"), ("a_s_two3_20", "shallow shape"), ("timing_f", "timing shape"))):
            f = q[col].to_numpy(); w = 1 / q.sig.to_numpy() ** 2
            if np.isfinite(f).all() and len(f):
                kk = np.sum(w * f * q.dE) / np.sum(w * f * f); ax[1, i].plot(tc, kk * f, "-", color=CAT[k], lw=1.5, label=lab)
        ax[1, i].axhline(0, color=MUTED, lw=0.8); ax[1, i].set_xlabel("lapse after S (s)")
        ax[1, i].set_title(f"{ns}: 2019-11 vs lapse window"); ax[1, i].legend(frameon=False, fontsize=7)
    ax[0, 0].set_ylabel("dv/v (%)"); ax[1, 0].set_ylabel("apparent ΔE (%)")
    fig.colorbar(sc, ax=ax[0, :].tolist(), pad=0.01).set_label("mean cc")
    fig.savefig(f"{FIG}/fig_r5_pgc.png", dpi=140, bbox_inches="tight"); plt.close(fig)


def fig_band():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for j, reg in enumerate(SUBS):
        p = tb(reg, "r5_band_q3.csv")
        if not os.path.exists(p):
            continue
        B = pd.read_csv(p)
        for k, f in enumerate(list(dict.fromkeys(B.form))[:3]):
            b = B[B.form == f]
            ax[j].plot(b.fit_setup, b.composite, "-o", color=CAT[k], lw=2, label=f"{f} composite")
            ax[j].plot(b.fit_setup, 2 * b.composite_sig, ":", color=CAT[k], lw=1.2, label=f"{f} 2σ")
            ax[j].plot(b.fit_setup, -2 * b.composite_sig, ":", color=CAT[k], lw=1.2)
        for i_, s in enumerate(B.fit_setup.unique()):
            if s in R5.PLAUS[reg]:
                ax[j].axvspan(i_ - 0.4, i_ + 0.4, color="#eef3f8", zorder=0)
        ax[j].axhline(0, color=MUTED, lw=0.8); ax[j].set_title(f"{reg}: mean slow-slip ΔD over fit kernels (shaded = plausible)")
        ax[j].set_ylabel("composite ΔD (%)"); ax[j].tick_params(axis="x", rotation=30); ax[j].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_r5_band.png", dpi=140); plt.close(fig)


def main():
    for f in (fig_forms, fig_dD, fig_q1, fig_pgc, fig_band):
        try:
            f(); print("ok", f.__name__)
        except Exception as e:
            import traceback
            print("FAILED", f.__name__, e, traceback.format_exc()[-600:])
    open(f"{R}/FIGURES_DONE", "w").write(R5.pt_now() + "\n")


if __name__ == "__main__":
    main()
