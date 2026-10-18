#!/usr/bin/env python
"""Figures for PLAN.md §6 Amendment 5 (s11_eq_netsplit.py outputs)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

EQ = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/tomo/eq"
TT = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/tomo/tables"
T = f"{EQ}/tables"
BLUE, ORANGE, AQUA, GREY, INK = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984", "#0b0b0b"
DIV = LinearSegmentedColormap.from_list("div", [BLUE, "#f0efec", ORANGE])
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#52514e",
                     "axes.labelcolor": INK, "xtick.color": "#52514e", "ytick.color": "#52514e"})
MK = {"UW/US": ("o", BLUE), "PB": ("s", AQUA), "temporary": ("^", ORANGE)}
EQK = {"WA_NW": 1.775, "WA_SW": 1.735, "SVI": 1.755}


def grp(net, cls):
    return "PB" if net == "PB" else ("UW/US" if cls == "permanent" else "temporary")


def fig_map():
    S = pd.read_csv(f"{T}/netsplit_station_residuals.csv")
    R = pd.read_csv(f"{T}/netsplit_residual_rows.csv")
    fig, ax = plt.subplots(2, 2, figsize=(11, 9.5))
    norm = TwoSlopeNorm(0, -0.5, 0.5)
    for j, (reg, ps) in enumerate((("WA_NW", "BV3"), ("WA_SW", "CV3"))):
        a = ax[0, j]
        P = pd.read_csv(f"{TT}/picks_{ps}.csv", low_memory=False).drop_duplicates("fam")
        a.scatter(P.lon, P.lat, s=6, c="#c3c2b7", lw=0, label="LFE families", zorder=1)
        s = S[(S.reg == reg)].copy()
        s["g"] = [grp(n, c) for n, c in zip(s.net, s.cls)]
        for g, (m, _) in MK.items():
            q = s[s.g == g]
            sc = a.scatter(q.lon, q.lat, c=q.med_r_rel_kept.clip(-0.5, 0.5), cmap=DIV, norm=norm, marker=m, s=110,
                           edgecolor=INK, lw=0.8, zorder=3, label=g)
        s["key"] = s.lat.round(2).astype(str) + s.lon.round(2).astype(str)
        for _, gg in s.groupby("key"):
            txt = "\n".join(f"{r.sta} {r.med_r_rel_kept:+.2f}" if np.isfinite(r.med_r_rel_kept) else f"{r.sta} n/a"
                             for r in gg.itertuples())
            a.annotate(txt, (gg.lon.iloc[0], gg.lat.iloc[0]), xytext=(6, 4), textcoords="offset points", fontsize=6,
                       color="#52514e")
        a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
        a.set_title(f"{reg.replace('_', ' ')}: per-station median tS−tP residual (kept rows, family-demeaned)\n"
                    "vs model (κ 1.773 × CVM v1.7 ray tP)", fontsize=8.5)
        a.set_xlabel("longitude"); a.set_ylabel("latitude")
        a.set_aspect(1 / np.cos(np.radians(47.3)))
        leg = a.legend(loc="lower left", fontsize=7, frameon=False)
        for h in leg.legend_handles[1:]:
            h.set_facecolor("#f0efec")
        cb = fig.colorbar(sc, ax=a, shrink=0.7, extend="both")
        cb.set_label("median residual (s); orange = more S−P than model")
        b = ax[1, j]
        d = R[(R.reg == reg) & R.kept_allfit].copy()
        d["g"] = [grp(n, c) for n, c in zip(d.net, d.cls)]
        for g, (m, col) in MK.items():
            q = d[d.g == g]
            b.scatter(q.epi_km, q.r_rel, s=12, marker=m, color=col, alpha=0.55, lw=0, label=f"{g} ({len(q)} rows)")
            if len(q) >= 15:
                bins = np.arange(0, 110, 10)
                mid = 0.5 * (bins[1:] + bins[:-1])
                med = [q[(q.epi_km >= lo) & (q.epi_km < hi)].r_rel.median() if ((q.epi_km >= lo) & (q.epi_km < hi)).sum() >= 8
                       else np.nan for lo, hi in zip(bins[:-1], bins[1:])]
                b.plot(mid, med, color=col, lw=2)
        b.axhline(0, color=GREY, lw=0.8)
        b.set_ylim(-1.0, 1.0)
        b.set_xlabel("epicentral distance (km)"); b.set_ylabel("tS−tP residual vs model, family-demeaned (s)")
        b.set_title(f"{reg.replace('_', ' ')}: rows kept by the clipped fit; lines = 10-km medians", fontsize=8.5)
        b.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(f"{EQ}/figures/fig_netsplit_residual_map.png", dpi=150)
    plt.close(fig)


def fig_kappa():
    A = pd.read_csv(f"{T}/netsplit_class_kappa.csv")
    C = pd.read_csv(f"{T}/netsplit_class_scrambled.csv")
    K = pd.read_csv(f"{T}/netsplit_scatter_kappa.csv")
    M = pd.read_csv(f"{T}/netsplit_synth_kappa.csv")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    # A: κ per class
    a = ax[0, 0]
    subs = ["all", "permanent", "temporary", "PB", "UW_US", "XU", "permanent_common", "temporary_common",
            "permanent_le50", "temporary_le50"]
    for j, (reg, col) in enumerate((("WA_NW", BLUE), ("WA_SW", ORANGE))):
        g = A[(A.reg == reg) & (A.nmin == A.groupby("subset").nmin.transform("max"))].set_index("subset")
        x = np.arange(len(subs)) + (j - 0.5) * 0.3
        y = [g.kappa.get(s, np.nan) for s in subs]; e = [g.sig.get(s, np.nan) for s in subs]
        ok = [bool(g["pass"].get(s, False)) if s in g.index else False for s in subs]
        a.errorbar(x, y, yerr=e, fmt="none", ecolor=col, lw=1)
        a.scatter(x, y, color=[col if o else "white" for o in ok], edgecolor=col, s=40, zorder=3,
                  label=f"{reg.replace('_', ' ')} (filled = passes gates)")
        yp = [g.kappa_profile.get(s, np.nan) for s in subs]
        a.scatter(x + 0.12, yp, marker="D", color="white", edgecolor=col, s=22, zorder=3, lw=1,
                  label=f"{reg.replace('_', ' ')} start-free profile κ")
        a.axhline(EQK[reg], color=col, ls="--", lw=1)
    a.axhspan(0.9, 1.4, color="#f0efec", zorder=0)
    a.text(0, 1.3, "κ < 1.4: collapsed clipped fit (mixture, rms ≥ 1.5 s)", fontsize=7, color="#52514e")
    a.set_xticks(range(len(subs))); a.set_xticklabels([s.replace("_", " ") for s in subs], rotation=40, ha="right", fontsize=7)
    a.set_ylabel("LFE κ (clipped tP-on-tS; diamonds = profile)"); a.set_ylim(1.0, 2.6)
    a.set_title("A  κ by station class; dashed = near-interface earthquakes", fontsize=9)
    a.legend(fontsize=7, frameon=False)
    # B: scrambled control
    b = ax[0, 1]
    for j, (reg, col) in enumerate((("WA_NW", BLUE), ("WA_SW", ORANGE))):
        for k, r in enumerate(C[C.reg == reg].itertuples()):
            yy = j * 2 + k
            b.plot([r.scr_p025, r.scr_p975], [yy, yy], color=col, lw=6, alpha=0.35, solid_capstyle="butt")
            b.plot(r.scr_med, yy, "|", color=col, ms=12)
            b.plot(r.d_real, yy, "o", color=INK, ms=6)
            b.annotate(f"real {r.d_real:+.3f} ({r.real_pct:.0f} pct)", (r.d_real, yy), xytext=(6, 5),
                       textcoords="offset points", fontsize=7)
    b.set_yticks(range(4))
    b.set_yticklabels(["NW, all stations", "NW, within distance terciles", "SW, all stations", "SW, within terciles"], fontsize=7)
    b.axvline(0, color=GREY, lw=0.8)
    b.set_xlabel("Δκ = κ(permanent) − κ(temporary)")
    b.set_title("B  real Δκ (dot) vs scrambled control 2.5–97.5 % (band)\n(class labels shuffled across stations)", fontsize=9)
    # C: near-station κ with estimators
    c = ax[1, 0]
    sets = [("all", 4), ("le50", 3), ("le40", 3)]
    est = [("kappa_E1", "E1 (as before)", "o"), ("kappa_IV", "IV (headline)", "s"), ("kappa_E3", "E3 (σS co-located)", "^"),
           ("kappa_SIMEX", "SIMEX", "v")]
    for j, (reg, col) in enumerate((("WA_NW", BLUE), ("WA_SW", ORANGE), ("SVI", AQUA))):
        g = K[K.reg == reg]
        for m, (ec, lab, mk) in enumerate(est):
            y = [g[(g.set == s) & (g.nmin == n)][ec].mean() for s, n in sets]
            x = np.arange(3) + (j - 1) * 0.27 + (m - 1.5) * 0.05
            c.scatter(x, y, marker=mk, color=col if ec == "kappa_IV" else "white", edgecolor=col, s=30, zorder=3,
                      label=lab if j == 0 else None)
            if ec == "kappa_IV":
                e = [g[(g.set == s) & (g.nmin == n)].sig_IV.mean() for s, n in sets]
                c.errorbar(x, y, yerr=e, fmt="none", ecolor=col, lw=1)
        c.axhline(EQK[reg], color=col, ls="--", lw=1)
        c.text(2.45, EQK[reg], reg.replace("_", " "), color=col, fontsize=7, va="center")
    c.set_xticks(range(3)); c.set_xticklabels(["all stations", "≤ 50 km", "≤ 40 km"])
    c.set_ylabel("LFE κ"); c.set_ylim(1.5, 2.5)
    c.set_title("C  near-station κ: uncorrected (E1) vs scatter-corrected; dashed = earthquakes", fontsize=9)
    c.legend(fontsize=7, frameon=False, loc="upper left")
    # D: synthetics
    d = ax[1, 1]
    for j, (reg, col) in enumerate((("WA_NW", BLUE), ("WA_SW", ORANGE))):
        for kt, ls in ((1.75, "-"), (1.90, ":")):
            g = M[(M.reg == reg) & (M.kappa_true == kt)]
            y1 = [g[(g.set == s) & (g.nmin == n)].bias_E1.mean() for s, n in sets]
            y2 = [g[(g.set == s) & (g.nmin == n)].bias_IV.mean() for s, n in sets]
            d.plot(range(3), y1, ls=ls, marker="o", color=col, mfc="white", label=f"{reg.replace('_', ' ')} E1, κ true {kt}")
            d.plot(range(3), y2, ls=ls, marker="s", color=col, label=f"{reg.replace('_', ' ')} IV, κ true {kt}")
    d.axhline(0, color=GREY, lw=0.8)
    d.set_xticks(range(3)); d.set_xticklabels(["all stations", "≤ 50 km", "≤ 40 km"])
    d.set_ylabel("recovered − true κ")
    d.set_title("D  synthetics at the real WA geometry and LFE noise\n(IV instrument from locations perturbed by 5 km)", fontsize=9)
    d.legend(fontsize=6.5, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(f"{EQ}/figures/fig_netsplit_kappa.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    fig_map()
    fig_kappa()
