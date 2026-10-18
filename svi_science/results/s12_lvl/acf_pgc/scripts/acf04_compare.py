#!/usr/bin/env python
"""S12c: summary table (pass rules) + ACF curve vs receiver-function curve figure.
Usage: acf04_compare.py --tag full|pilot2
Out: tables/summary_<tag>.csv, figures/compare_rf_<tag>_<STA>.png
Recovery fit Δt(k) = A exp(-k/τ) (k >= 0, ±3-d curve) is reported for every curve; it only means something where the
0-15 d bin passes P3.
"""
import sys, os, glob, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from acf_common import WS, ROOT

ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="full"); a = ap.parse_args()
TGT = json.load(open(f"{WS}/tables/predicted_targets.json"))
rows = []
for f in sorted(glob.glob(f"{WS}/tables/epoch_{a.tag}_*.csv")):
    parts = os.path.basename(f)[:-4].split("_")
    ns, band, mask = parts[-4], parts[-3] + "_" + parts[-2], parts[-1]
    T = pd.read_csv(f)
    for ph, pr in [("SS", "SS"), ("PS", "PSr"), ("PP", "ZZ")]:
        for bn in ["post0_15", "post15_30", "post30_60", "post60_120"]:
            g = T[(T.pair == pr) & (T.bin == bn)].set_index("metric")
            if "ddt" not in g.index:
                continue
            k = f"kappa2.3"; tg = {"SS": "dSS_ms", "PS": "dPS_ms", "PP": "dPP_ms_if_dVp_half"}[ph]
            rows.append(dict(station=ns, band=band, mask=mask, phase=ph, bin=bn,
                             ddt_ms=1e3 * g.loc["ddt", "obs"], scr_rsd_ms=1e3 * g.loc["ddt", "scr_rsd"],
                             pct_neg=g.loc["ddt", "pct_neg"], ddt_ccpk_ms=1e3 * g.loc["ddt_ccpk", "obs"] if "ddt_ccpk" in g.index else np.nan,
                             ddt_str_ms=1e3 * g.loc["ddt_str", "obs"], ctl_ms=1e3 * g.loc["ddt_ctl", "obs"],
                             ctl_pct=g.loc["ddt_ctl", "pct_neg"], dcc=g.loc["dcc", "obs"], dcc_pct=g.loc["dcc", "pct_neg"],
                             amp_lvl=g.loc["amp_lvl", "obs"], amp_pct=g.loc["amp_lvl", "pct_neg"],
                             target_ms=TGT.get(ns, {}).get("target_delays", {}).get(k, {}).get(tg, np.nan)))
S = pd.DataFrame(rows)
if len(S):
    S["z_vs_target"] = S.target_ms / S.scr_rsd_ms
    S["P3_pass"] = (S.ddt_ms < 0) & (S.pct_neg <= 0.05)
    S.to_csv(f"{WS}/tables/summary_{a.tag}.csv", index=False)
    pd.set_option("display.width", 250)
    print(S[S.bin == "post0_15"].round(3).to_string(index=False))

# RF curve vs ACF curve, per station; primary channel = b05_2 tremorfree SS
for ns in sorted(set(S.station)) if len(S) else []:
    rf = f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin/tables/full/recovery_{ns}.csv"
    fig, axs = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    if os.path.exists(rf):
        R = pd.read_csv(rf)
        axs[0].plot(R.centre, R.cc_sv_clean_w15, "k-", lw=1.2, label="RF ΔCC_SV (±7.5 d)")
        axs[0].plot(R.centre, R.cc_sv_clean_w15_lo95, color="0.5", lw=0.6, ls="--", label="RF scrambled 95 %")
        axs[0].plot(R.centre, R.cc_sv_clean_w15_hi95, color="0.5", lw=0.6, ls="--")
        axs[0].axhline(0, color="k", lw=0.4); axs[0].legend(fontsize=7); axs[0].set_ylabel("RF ΔCC_SV")
    fits = []
    for ax, (band, ph, mem) in zip(axs[1:], [("b05_2", "SS", ["RR", "TT"]), ("b1_4", "SS", ["RR", "TT"]), ("b1_4", "PS", ["ZR+", "ZR-"])]):
        for mask, ls in [("tremorfree", "-"), ("none", ":")]:
            cf = f"{WS}/tables/curve_{a.tag}_{ns}_{band}_{mask}.csv"
            if not os.path.exists(cf):
                continue
            C = pd.read_csv(cf)
            c = C[(C.hw == 3) & (C.pair.isin(mem))].groupby("k")[["ddt", "cc_lvl"]].mean().reset_index()
            ax.plot(c.k, 1e3 * c.ddt, color="0.5", lw=0.8, ls=ls, label=f"{mask}")
            if mask == "tremorfree":
                sc = ax.scatter(c.k, 1e3 * c.ddt, c=c.cc_lvl, cmap="viridis", s=12, vmin=0.5, vmax=1, zorder=3)
                plt.colorbar(sc, ax=ax, label="CC (LVL window)")
                # recovery fit
                d = c[c.k >= 0]; best = None
                for tau in np.arange(3, 121, 1.0):
                    e = np.exp(-d.k.values / tau); A = np.dot(e, d.ddt.values) / np.dot(e, e)
                    r = np.sum((d.ddt.values - A * e) ** 2)
                    if best is None or r < best[2]: best = (tau, A, r)
                fits.append(dict(station=ns, band=band, phase=ph, tau_d=best[0], A_ms=1e3 * best[1]))
                ax.plot(d.k, 1e3 * best[1] * np.exp(-d.k / best[0]), color="tab:red", lw=0.8,
                        label=f"fit A={1e3 * best[1]:.1f} ms, τ={best[0]:.0f} d")
            scf = f"{WS}/tables/curve_scr_{a.tag}_{ns}_{band}_{mask}.npy"
            if ph == "SS" and mask == "tremorfree" and os.path.exists(scf):
                sc_ = np.load(scf); kk = np.arange(-60, 121, 3)
                ax.plot(kk, 1e3 * np.nanpercentile(sc_, 5, 0), color="0.3", lw=0.5, ls="--", label="scrambled 5/95 %")
                ax.plot(kk, 1e3 * np.nanpercentile(sc_, 95, 0), color="0.3", lw=0.5, ls="--")
        tg = TGT.get(ns, {}).get("target_delays", {}).get("kappa2.3", {}).get({"SS": "dSS_ms", "PS": "dPS_ms"}[ph], np.nan)
        ax.axhline(tg, color="tab:red", ls=":", lw=0.8, label=f"target {tg:.0f} ms")
        ax.axhline(0, color="k", lw=0.4); ax.set_ylabel(f"{band} {ph}\nΔt LVL−pre (ms)"); ax.legend(fontsize=7, loc="lower right")
    for ax in axs:
        ax.axvline(0, color="k", lw=0.5); ax.axvspan(0, 15, color="tab:red", alpha=0.05)
    axs[-1].set_xlabel("days after ETS midpoint (ALBH GPS)")
    fig.suptitle(f"{ns} ({a.tag}): receiver-function change (top) vs ACF differential LVL delay (±3-d boxcar)", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{WS}/figures/compare_rf_{a.tag}_{ns}.png", dpi=110); plt.close(fig)
    if fits:
        pd.DataFrame(fits).to_csv(f"{WS}/tables/recovery_fit_{a.tag}_{ns}.csv", index=False)
print("done")
