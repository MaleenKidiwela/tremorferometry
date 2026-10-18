#!/usr/bin/env python
"""P3D figures (figures/s_*.png)."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

FIG = S.FIG
C_OLD, C_NEW, C_PQ = "#4c6ef5", "#f76707", "#12b886"


def fig_premise():
    PR = pd.read_csv(f"{S.P3D}/premise_family_sites_v2.csv")
    SI = pd.read_csv(f"{S.P3D}/sites_p3d_v2.csv")
    F = S.families()
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    # (a) map
    a = ax[0]
    a.scatter(F.lon, F.lat, s=12, c="0.6", label="130 Bostock families", zorder=2)
    for r in SI.itertuples():
        c = C_PQ if r.role == "new" and r.site in ("LHLYB", "JRDRB", "ALBHB", "DAOB",
                                                   "OHOH", "ODEER") else (
            C_NEW if r.role == "new" else C_OLD)
        a.plot(r.lon, r.lat, "^", ms=11, color=c, mec="k", mew=0.6, zorder=3)
        a.annotate(r.site, (r.lon, r.lat), fontsize=6.5, xytext=(2, 4),
                   textcoords="offset points")
    a.set_xlabel("longitude"); a.set_ylabel("latitude")
    a.set_title("(a) sites after the co-located merge\nblue = frozen, orange = P3D new, green = PQ (V2)")
    a.grid(alpha=.25)
    # (b) premise correction
    b = ax[1]
    w = 0.26
    x = np.arange(3)
    for i, (lab, cols, col) in enumerate([
            ("before", ["n45_before", "n55_before"], C_OLD),
            ("V1 after", ["n45_after_V1", "n55_after_V1"], C_NEW),
            ("V2 after (+PQ)", ["n45_after_V2", "n55_after_V2"], C_PQ)]):
        b.bar([0 + (i - 1) * w, 1 + (i - 1) * w], [PR[cols[0]].median(), PR[cols[1]].median()],
              width=w, color=col, label=lab, ec="k", lw=.4)
    b.set_xticks([0, 1]); b.set_xticklabels(["sites $\\leq$45 km", "sites $\\leq$55 km"])
    b.set_ylabel("median over the 130 families")
    b.set_title("(b) the premise correction\n(CLRS/OSQM add NO geometry)")
    b.legend(fontsize=8); b.grid(alpha=.25, axis="y")
    # (c) rho
    R = pd.read_csv(f"{S.P3D}/step2_rho.csv")
    pq = f"{S.P3D}/step2_rho_pqdip.csv"                 # PQ legs (dip-corrected), addendum 2026-09-19b
    if os.path.exists(pq):
        R = pd.concat([R, pd.read_csv(pq)], ignore_index=True)
    R = R[R.ref == "PGC"].sort_values("ratio")
    c = ax[2]
    ST = pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta")
    d = [float(np.median(S.hyp_km(F.lat.values, F.lon.values, F.dep.values,
                                  ST.lat[s], ST.lon[s]))) for s in R.sta]
    PQSET = set(S.V2_3C)
    c.axhline(1.0, color="k", lw=1, ls="--")
    c.axhline(3.0, color="crimson", lw=1, ls=":")
    c.text(78, 3.05, "Merlin's discovery bar ($\\rho\\geq3\\times$)", color="crimson", fontsize=8)
    c.scatter(d, R.ratio, s=60, ec="k", zorder=3,
              c=[C_PQ if ss in PQSET else C_NEW for ss in R.sta])
    for dd, rr, ss in zip(d, R.ratio, R.sta):
        c.annotate(ss, (dd, rr), fontsize=7.5, xytext=(3, 4), textcoords="offset points")
    c.set_yscale("log"); c.set_ylim(0.1, 5)
    c.set_xlabel("median hypocentral distance to the family field (km)")
    c.set_ylabel(r"$\rho_{\rm new}\,/\,\rho_{\rm PGC}$")
    c.set_title("(c) single-event registration, same events\norange = P3D new, green = PQ (addendum)")
    c.grid(alpha=.25, which="both")
    fig.tight_layout()
    fig.savefig(f"{FIG}/s_premise_and_rho.png", dpi=150)
    print(f"-> {FIG}/s_premise_and_rho.png")


def fig_cf():
    p = f"{S.P3D}/step4a_residuals.csv"
    if not os.path.exists(p):
        return
    D = pd.read_csv(p)
    SU = pd.read_csv(f"{S.P3D}/step4a_summary.csv")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for site in SU.sort_values("med_abs_resid_s").site.head(8):
        g = D[D.site == site]
        h, e = np.histogram(g.resid_s, bins=np.arange(-5, 5.01, 0.25), density=True)
        ax[0].plot(0.5 * (e[1:] + e[:-1]), h, lw=1.4,
                   label=f"{site} ({np.median(np.abs(g.resid_s)):.2f} s)")
    ax[0].axhline(1 / 10.0, color="k", ls="--", lw=1, label="uniform (no timing info)")
    ax[0].set_xlabel("CF peak time $-$ Crosson S (s)")
    ax[0].set_ylabel("density")
    ax[0].set_title("(a) CF timing residual at hp events, test episodes")
    ax[0].legend(fontsize=7.5); ax[0].grid(alpha=.25)
    k = SU.sort_values("med_abs_resid_s")
    ax[1].barh(range(len(k)), k.med_abs_resid_s,
               color=[C_NEW if x else "0.7" for x in k.le60], ec="k", lw=.4)
    ax[1].axvline(1.5, color="crimson", ls="--", lw=1.2, label="kill bar 1.5 s")
    ax[1].axvline(2.5, color="k", ls=":", lw=1.2, label="uniform chance 2.5 s")
    ax[1].set_yticks(range(len(k))); ax[1].set_yticklabels(k.site, fontsize=8)
    ax[1].set_xlabel("median |residual| (s)")
    ax[1].set_title("(b) per site (orange = inside the 60 km rule)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.25, axis="x")
    fig.tight_layout()
    fig.savefig(f"{FIG}/s_cf_timing.png", dpi=150)
    print(f"-> {FIG}/s_cf_timing.png")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_premise()
    fig_cf()
