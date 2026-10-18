#!/usr/bin/env python
"""s12_figures -- figures for S12-LVL LFE side (forward calibration, reflection record sections + scan, superposed epoch).
Line graphs; dv/v points coloured by weight; no error bars on dv/v.
    python s12_figures.py [fwd] [refl] [epoch]
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, f"{ROOT}/svi_science/scripts"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s11_common as C
OUT = f"{ROOT}/svi_science/results/s12_lvl/lfe"; TB, FIG = f"{OUT}/tables", f"{OUT}/figures"
os.makedirs(FIG, exist_ok=True)
WINS = ["1-2", "2-3", "3-4", "4-6", "2-5"]
XPOS = {"pre": -5, "0": 0, "1": 1, "2": 2, "3": 3, "4-6": 5, "7-10": 8.5, "11-20": 15.5, "tail": 40}


def fig_fwd():
    S = pd.read_csv(f"{TB}/fwd_coef_summary.csv")
    R = pd.read_csv(f"{TB}/fwd_rays.csv")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    xs = np.arange(4); wl = ["1-2", "2-3", "3-4", "4-6"]
    for v, c in (("ell20", "C0"), ("ell40", "C1"), ("two3_20", "C2"), ("ell40_cvm", "0.5")):
        g = S[S.variant == v].set_index("win").loc[wl]
        ax[0].plot(xs, 2 * g.a_lvl_med, "-o", color=c, label=f"{v} median")
        ax[0].plot(xs, 2 * g.a_lvl_q10, ":", color=c); ax[0].plot(xs, 2 * g.a_lvl_q90, ":", color=c)
    ax[0].set_xticks(xs); ax[0].set_xticklabels([f"{w} s" for w in wl]); ax[0].set_ylabel("coda dv/v for +2 % in the LVL (%)")
    ax[0].set_title("Forward target: S-anchored coda stretch\n(dotted: pairs 10-90 %)"); ax[0].legend(fontsize=7)
    t = "H3.5_v3"
    ax[1].scatter(R.epi_km, R[f"lagS_SrS_b_{t}"], s=6, c="C3", label="SrS off LVL base (H 3.5 km, Vs 3.0)")
    ax[1].scatter(R.epi_km, R[f"lagS_SrS_m_{t}"], s=6, c="C4", label="SrS off oceanic Moho (7 km)")
    ax[1].scatter(R.epi_km, R[f"lagS_PrS_b_{t}"], s=6, c="C2", label="PrS (LVL base)")
    ax[1].set_xlabel("epicentral distance (km)"); ax[1].set_ylabel("time after direct S (s)")
    ax[1].set_title("Predicted later phases (flat layers)"); ax[1].legend(fontsize=7); ax[1].set_xlim(0, 110)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_fwd.png", dpi=130); plt.close(fig)


def fig_refl():
    z = np.load(f"{TB}/refl_traces.npz")
    T = z["TREL"]; st = z["netsta"]; epi = z["epi"]
    HG = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]; VG = [2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8]
    ia, ib = HG.index(3.5), VG.index(3.0)
    sites = {"PGC site": ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"], "B926": ["PB.B926"], "SNB": ["CN.SNB"], "VGZ": ["CN.VGZ"]}
    fig, axs = plt.subplots(2, 4, figsize=(17, 9), gridspec_kw=dict(height_ratios=[3, 1.3]))
    for k, (nm, sts) in enumerate(sites.items()):
        m = np.isin(st, sts); idx = np.flatnonzero(m)[np.argsort(epi[m])]
        R = z["RH"][idx]
        tm = (T >= -1) & (T <= 6)
        ax = axs[0, k]
        ax.imshow(np.clip(R[:, tm], -1, 1), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1,
                  extent=[-1, 6, len(idx), 0], interpolation="nearest")
        ax.plot(z["lag_SrS_b"][idx, ia, ib], np.arange(len(idx)) + 0.5, "k.", ms=2, label="SrS LVL base")
        ax.plot(z["lag_SrS_m"][idx, ia, ib], np.arange(len(idx)) + 0.5, "m.", ms=2, label="SrS Moho")
        ax.set_title(f"{nm}: {len(idx)} pairs, horizontal residual log-envelope\n(sorted by distance {epi[idx].min():.0f}-{epi[idx].max():.0f} km)", fontsize=9)
        ax.set_xlabel("time after S (s)"); ax.set_ylabel("pair (by distance)")
        if k == 0:
            ax.legend(fontsize=7, loc="lower right")
        ax = axs[1, k]
        bins = np.quantile(epi[idx], [0, 1 / 3, 2 / 3, 1]) if len(idx) > 6 else [epi[idx].min(), epi[idx].max()]
        for j in range(len(bins) - 1):
            mm = idx[(epi[idx] >= bins[j]) & (epi[idx] <= bins[j + 1])]
            if len(mm) < 2:
                continue
            ax.plot(T[tm], np.nanmean(z["RH"][mm][:, tm], 0), label=f"{bins[j]:.0f}-{bins[j + 1]:.0f} km (n {len(mm)})")
            ax.axvline(np.nanmedian(z["lag_SrS_b"][mm, ia, ib]), color=f"C{j}", ls=":", lw=1)
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("time after S (s)"); ax.set_ylabel("mean residual ln E")
        ax.legend(fontsize=6); ax.set_xlim(-1, 6)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_refl_records.png", dpi=120); plt.close(fig)
    # scan grids + moveout-corrected stacks
    sc = pd.read_csv(f"{TB}/refl_scan.csv").set_index("phase")
    phs = [p for p in ["SrS_b", "SrS_m", "SrP_b", "PrP_b", "SrS_b_mf"] if os.path.exists(f"{TB}/refl_grid_{p}.npy")]
    fig, axs = plt.subplots(2, len(phs), figsize=(4 * len(phs), 7.5))
    for k, ph in enumerate(phs):
        g = np.load(f"{TB}/refl_grid_{ph}.npy")
        ax = axs[0, k]
        im = ax.imshow(g, origin="lower", aspect="auto", cmap="viridis", extent=[VG[0] - 0.1, VG[-1] + 0.1, -0.5, len(HG) - 0.5])
        ax.set_yticks(range(len(HG))); ax.set_yticklabels(HG); ax.set_xlabel("LVL Vs (km/s)"); ax.set_ylabel("H (km)")
        r = sc.loc[ph]
        ax.set_title(f"{ph}: best {r.best_peak:.3f} (H {r.best_H}, Vs {r.best_V})\nscan p {r.p_scan:.3f}; scrambled max q95 {r.scr_max_q95:.3f}", fontsize=8)
        plt.colorbar(im, ax=ax, fraction=0.05)
        ax = axs[1, k]
        key = {"SrS_b": "RH", "SrS_m": "RH", "SrP_b": "RP", "PrP_b": "RP", "SrS_b_mf": "MF"}[ph]
        R = z[key]; L = z[f"lag_{ph.replace('_mf', '')}"]
        add = -z["sp"] if ph == "PrP_b" else (z["sp_pred"] - z["sp"] if ph == "SrP_b" else 0)
        w = np.sqrt(z["n"])
        tt = np.arange(-1.2, 1.2 + 1e-9, 0.025)
        for lab, (a_, b_), col in (("nominal H3.5 Vs3.0", (ia, ib), "k"), ("best", (HG.index(r.best_H), VG.index(r.best_V)), "C3")):
            S_ = np.zeros(len(tt)); W_ = np.zeros(len(tt))
            for i in range(len(R)):
                v = np.interp(tt + L[i, a_, b_] + (add[i] if np.ndim(add) else 0), T, R[i], left=np.nan, right=np.nan)
                o = np.isfinite(v); S_[o] += w[i] * v[o]; W_[o] += w[i]
            ax.plot(tt, S_ / np.maximum(W_, 1e-9), color=col, label=lab)
        rng = np.random.default_rng(1)
        for d in range(20):
            perm = np.arange(len(R))
            for s_ in np.unique(st):
                ii = np.flatnonzero(st == s_); perm[ii] = rng.permutation(ii)
            S_ = np.zeros(len(tt)); W_ = np.zeros(len(tt))
            for i in range(len(R)):
                v = np.interp(tt + L[perm[i], ia, ib] + (add[i] if np.ndim(add) else 0), T, R[i], left=np.nan, right=np.nan)
                o = np.isfinite(v); S_[o] += w[i] * v[o]; W_[o] += w[i]
            ax.plot(tt, S_ / np.maximum(W_, 1e-9), color="0.7", lw=0.6, zorder=0)
        ax.axvline(0, color="k", ls=":"); ax.set_xlabel("time relative to the predicted phase (s)")
        ax.set_title("moveout-corrected stack (grey: scrambled lags)", fontsize=8); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_refl_scan.png", dpi=120); plt.close(fig)


def fig_epoch(sub):
    p = f"{TB}/{sub}/epoch_pooled.csv"
    if not os.path.exists(p):
        return
    T = pd.read_csv(p)
    Cc = pd.read_csv(f"{TB}/{sub}/epoch_ctrl.csv").set_index("obs") if os.path.exists(f"{TB}/{sub}/epoch_ctrl.csv") else None
    fig, axs = plt.subplots(2, 4, figsize=(18, 8))
    obs = [f"dvv_{w}" for w in WINS] + ["day_slope", "lnPS_Z", "lnPS_H"]
    for k, q in enumerate(obs):
        ax = axs.flat[k]
        g = T[(T.obs == q) & T.cls.isin(XPOS) & (T.n_series >= 10)].copy(); g["x"] = g.cls.map(XPOS); g = g.sort_values("x")
        ax.plot(g.x, g.value, "-", color="0.4", lw=1)
        sc = ax.scatter(g.x, g.value, c=np.log10(g.weight), cmap="viridis", s=45, zorder=3)
        for cl, mk in (("betw_early", "s"), ("betw_late", "D"), ("G0..15", "^"), ("G-15..-1", "v"), ("G16..45", ">")):
            r = T[(T.obs == q) & (T.cls == cl) & (T.n_series >= 10)]
            if len(r):
                ax.scatter([60 if cl.startswith("betw") else {"G0..15": 7.5, "G-15..-1": -8, "G16..45": 30}[cl]], r.value,
                           marker=mk, facecolor="none", edgecolor="C3", s=50, label=cl, zorder=4)
        r020 = T[(T.obs == q) & (T.cls == "0..20")]
        if len(r020):
            ax.axhline(float(r020.value.iloc[0]), color="C1", lw=1, ls="--", label="0..20 aggregate")
        if q.startswith("dvv_") and "target_pct" in T and len(r020):
            tg = float(r020.target_pct.iloc[0])
            ax.axhline(tg, color="C3", lw=1.5, ls=":", label=f"forward target +{tg:.2f} % (2 % in LVL)")
        if Cc is not None and q in Cc.index:
            ax.scatter([10], [Cc.loc[q, "ctrl_mean"]], marker="x", color="C0", s=50, label="time-of-day control (0..20)")
        ax.axhline(0, color="k", lw=0.5)
        unit = {"day_slope": "Δ ln-energy slope (1/s)", "lnPS_Z": "Δ ln(P/S) Z (%)", "lnPS_H": "Δ ln(P/S) H (%)"}.get(q, "dv/v − between-ETS (%)")
        ax.set_title(q, fontsize=10); ax.set_ylabel(unit); ax.set_xlabel("days after the front reached the family")
        if k == 0:
            ax.legend(fontsize=6, loc="best")
    fig.suptitle(f"{sub.upper()}: superposed epoch, class − between-ETS (points coloured by log10 weight)", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_epoch_{sub}.png", dpi=120); plt.close(fig)
    # lapse profile of the headline contrast vs forward target
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    wl = ["1-2", "2-3", "3-4", "4-6"]; xs = np.arange(4)
    g = T[(T.cls == "0..20") & T.obs.isin([f"dvv_{w}" for w in wl])].set_index("obs").loc[[f"dvv_{w}" for w in wl]]
    ax.plot(xs, g.value, "-", color="0.4")
    ax.scatter(xs, g.value, c=np.log10(g.weight), cmap="viridis", s=60, zorder=3, label="data: ETS days 0..20 − between")
    if "target_pct" in g:
        ax.plot(xs, g.target_pct, "r:", label="forward target, +2 % in LVL")
        ax.plot(xs, g.target_pct / 2, "r--", lw=0.8, label="+1 % in LVL")
    if "pre" in set(T.cls):
        gp = T[(T.cls == "pre") & T.obs.isin([f"dvv_{w}" for w in wl])].set_index("obs").loc[[f"dvv_{w}" for w in wl]]
        ax.plot(xs, gp.value, "-s", color="C2", ms=4, label="pre-front − between")
    ax.axhline(0, color="k", lw=0.5); ax.set_xticks(xs); ax.set_xticklabels([f"{w} s" for w in wl])
    ax.set_ylabel("dv/v (%)"); ax.set_title(f"{sub.upper()}: lapse profile"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_epoch_lapse_{sub}.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    what = sys.argv[1:] or ["fwd", "refl", "epoch"]
    if "fwd" in what:
        fig_fwd()
    if "refl" in what and os.path.exists(f"{TB}/refl_scan.csv"):
        fig_refl()
    if "epoch" in what:
        for sub in ("svi", "wa_nw"):
            fig_epoch(sub)
    print("figures done")


def fig_novel():
    fig, axs = plt.subplots(2, 4, figsize=(19, 8.5))
    for row, sub in enumerate(("svi", "wa_nw")):
        p = f"{TB}/{sub}"
        if not os.path.exists(f"{p}/novel_A_curve.csv"):
            continue
        # A: daily curve vs front lag
        Ac = pd.read_csv(f"{p}/novel_A_curve.csv"); A = pd.read_csv(f"{p}/novel_A.csv")
        ax = axs[row, 0]
        for wn, col in (("2-5", "k"), ("4-6", "C1")):
            g = Ac[(Ac.win == wn) & (Ac.weight > 0)]
            ax.plot(g.day, g.value, "-", color=col, lw=0.8, label=f"data {wn} s")
            ax.scatter(g.day, g.value, c=np.log10(g.weight + 1e-9), cmap="viridis", s=18, zorder=3)
        g = Ac[Ac.win == "2-5"]
        ax.plot(g.day, g.target_curve, "r:", label="target: +2 % LVL switching on with the front")
        a = A[(A.win == "2-5") & (A.v_kmd == 10)].iloc[0]
        ax.set_title(f"{sub.upper()} A: dv/v − between-ETS vs day after front\nA = {a.A:+.3f} ± {a.A_sig:.3f} % (target +{a.A_target:.2f}), scrambled p {a.scr_p:.2f}", fontsize=8)
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("day after the front reached the family"); ax.set_ylabel("dv/v (%)")
        ax.legend(fontsize=6); ax.set_ylim(-0.6, 0.6)
        # B: recharge curve
        ax = axs[row, 1]
        if os.path.exists(f"{p}/novel_B_curve.csv"):
            Bc = pd.read_csv(f"{p}/novel_B_curve.csv"); Bk = pd.read_csv(f"{p}/novel_B_check.csv").set_index("obs")
            for q, col in (("dvv_2-5", "k"), ("dvv_4-6", "C1")):
                g = Bc[Bc.obs == q]
                ax.plot(g.mid, g.value, "-", color=col, lw=1, label=q)
                ax.scatter(g.mid, g.value, c=np.log10(g.weight), cmap="viridis", s=40, zorder=3)
            tm = np.linspace(2, 15, 50)
            ax.plot(tm, 0.45 * (np.exp(-tm * 30.44 / 45) - np.exp(-7 * 30.44 / 45)), "r:", lw=1, label="target shape: +2 % LVL decaying, τ 45 d")
            ax.plot(tm, -0.45 / 14 * (tm - 7.5), "r--", lw=0.8, label="target: 0.45 % linear fall over 14 mo")
            b = Bk.loc["dvv_2-5"]
            ax.set_title(f"{sub.upper()} B: between-ETS dv/v vs months since ETS (annual removed)\nslope {b.slope:+.3f} %/yr, scrambled p {b.scr_p:.3f}, LOSO {b.loso_min:+.2f}..{b.loso_max:+.2f}", fontsize=8)
            ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("months since the previous ETS front"); ax.set_ylabel("dv/v (%)"); ax.legend(fontsize=6)
        # C: source scalars by front class
        ax = axs[row, 2]
        Cc = pd.read_csv(f"{p}/novel_C.csv"); Cc = Cc[Cc.adjusted & (Cc.minus == "betw")]
        xs = {"pre": -1, "0": 0, "1": 1, "2-3": 2.5, "4-20": 8}
        for q, col in (("sr", "C0"), ("width", "C2"), ("lnPS_Z_rms", "C3"), ("slope", "C4")):
            g = Cc[Cc.obs == q].copy(); g["x"] = g.cls.map(xs); g = g.sort_values("x")
            ax.plot(g.x, g.value, "-o", color=col, ms=4, label=q)
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("days after front (class)"); ax.set_ylabel("class − between-ETS")
        ax.set_title(f"{sub.upper()} C: source/decay scalars (size-adjusted)", fontsize=8); ax.legend(fontsize=6)
        # D: per-major contrast vs duration
        ax = axs[row, 3]
        Mj = pd.read_csv(f"{p}/novel_D_majors.csv"); Dd = pd.read_csv(f"{p}/novel_D.csv")
        g = Mj[(Mj.win == "2-5") & np.isfinite(Mj.Du)]
        sc = ax.scatter(10 * g.Du, g.value, c=np.log10(1 / g.sigma ** 2), cmap="viridis", s=40)
        d = Dd[(Dd.win == "2-5") & (Dd["size"] == "Du") & (Dd.form == "size_only")].iloc[0]
        xx = np.linspace(10 * g.Du.min(), 10 * g.Du.max(), 10)
        ax.plot(xx, d.mean_contrast + d.slope * (xx / 10 - g.Du.mean()), "k-", lw=1, label=f"fit {d.slope:+.3f} ± {d.slope_sig:.3f} %/10 d")
        ax.plot(xx, d.target_slope * xx / 10, "r:", label=f"target ∝ duration ({d.target_slope:.2f} %/10 d)")
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("episode duration (d)"); ax.set_ylabel("ETS days 0..20 − between (%)")
        ax.set_title(f"{sub.upper()} D: per-episode contrast vs size (2-5 s)", fontsize=8); ax.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(f"{FIG}/fig_novel.png", dpi=115); plt.close(fig)


if __name__ == "__main__" and "novel" in sys.argv[1:]:
    fig_novel(); print("novel figure done")


def fig_recharge():
    for sub in ("wa_nw", "svi"):
        p = f"{TB}/{sub}"
        if not os.path.exists(f"{p}/r2_verdict.json"):
            continue
        R = json.load(open(f"{p}/r2_verdict.json"))
        fig, axs = plt.subplots(1, 4, figsize=(19, 4.4))
        cur = pd.read_csv(f"{p}/r2_curve.csv")
        ax = axs[0]
        for k, (lab, g) in enumerate(cur.groupby("set")):
            ax.plot(g.mid, g.value, "-", color=f"C{k}", lw=1, label=lab)
            ax.scatter(g.mid, g.value, c=np.log10(g.weight), cmap="viridis", s=30, zorder=3)
        tm = np.linspace(2, 30, 60)
        ax.plot(tm, 0.45 * (np.exp(-tm * 30.44 / 45) - np.exp(-7 * 30.44 / 45)), "r:", label="target shape (+2 % LVL, τ 45 d)")
        h = R["headline_all"]
        ax.set_title(f"{sub.upper()}: between-ETS dv/v 2-5 s vs months since ETS (season removed)\nall families {h['b']:+.3f} ± {h['sig']:.3f} %/yr, scrambled p {h['scr_p']:.3f}", fontsize=8)
        ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("months since the previous front"); ax.set_ylabel("dv/v (%)"); ax.legend(fontsize=6)
        ax = axs[1]
        PF = pd.read_csv(f"{p}/r2_per_family.csv")
        if len(PF):
            ax.scatter(PF.n, PF.slope, c=np.where(PF.has_front, "C3", "C0"), s=20)
            ax.axhline(0, color="k", lw=0.5); ax.set_xscale("log")
            c1 = R["C1_per_family"]
            ax.set_title(f"per-family slopes (red: 11 families of §6)\n{c1['n_neg']}/{c1['n']} negative, binomial p {c1['binom_p']:.3f}", fontsize=8)
            ax.set_xlabel("between-ETS day measurements"); ax.set_ylabel("slope (%/yr)"); ax.set_ylim(-5, 5)
        ax = axs[2]
        wl = ["1-2", "2-3", "3-4", "4-6"]; xs = np.arange(4)
        y = [R["lapse"][w]["b"] for w in wl]
        ax.plot(xs, y, "k-o", label="data")
        for nm, pp in R["lapse_patterns"].items():
            pat = {"LVL": [0.191, 0.223, 0.216, 0.201], "shallow 0-5 km": [0.166, 0.263, 0.304, 0.321], "25-45 km": [0.649, 0.543, 0.505, 0.445]}[nm]
            ax.plot(xs, pp["k"] * np.array(pat), ":", label=f"{nm} pattern (χ² {pp['chi2']:.1f})")
        ax.set_xticks(xs); ax.set_xticklabels([f"{w} s" for w in wl]); ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel("slope (%/yr)"); ax.set_title("lapse dependence of the slope vs forward patterns", fontsize=8); ax.legend(fontsize=6)
        ax = axs[3]
        SH = pd.read_csv(f"{p}/r2_shape.csv")
        for lab, g in SH[SH.form == "exp"].groupby("set"):
            li = SH[(SH.set == lab) & (SH.form == "linear")].aic.iloc[0]
            ax.plot(g.tau, g.aic - li, "-o", ms=3, label=lab)
        ax.axhline(0, color="k", lw=0.5); ax.set_xscale("log"); ax.set_xlabel("τ (d)"); ax.set_ylabel("AIC(exp) − AIC(linear)")
        ax.set_title("shape: exponential recovery vs linear", fontsize=8); ax.legend(fontsize=6)
        fig.tight_layout(); fig.savefig(f"{FIG}/fig_recharge_{sub}.png", dpi=120); plt.close(fig)


if __name__ == "__main__" and "recharge" in sys.argv[1:]:
    fig_recharge(); print("recharge figure done")
if __name__ == "__main__" and "lowband" in sys.argv[1:] and os.path.exists(f"{TB}/refl_scan_lb.csv"):
    print(pd.read_csv(f"{TB}/refl_scan_lb.csv").to_string())
