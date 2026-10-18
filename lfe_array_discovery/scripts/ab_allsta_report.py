#!/usr/bin/env python
"""Compare hypo runs (headline GSC 1.73 R1.5) -- current 48-family result vs the all-station
templates (refined, and the common-time control) -- and draw the catalogue figure.

  python scripts/ab_allsta_report.py --runs old=results/ab_royer/NW/hypo,allsta=...,common=... \
      --fig results/ab_royer/NW/hypo_allsta/NW_catalogue_allsta.png --main allsta
-> <dir of fig>/allsta_comparison.csv (+ truth-score tables truth_<name>.csv)
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

TRUTH = f"{HERE}/../results/ab_royer/truth/royer2014_templates_digitized.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--fig", required=True)
    ap.add_argument("--main", default="allsta")
    ap.add_argument("--sub", default="NW")
    a = ap.parse_args()
    runs = dict(x.split("=", 1) for x in a.runs.split(","))
    od = os.path.dirname(a.fig)
    lo, hi = A.SUBARRAY[a.sub]["depth_pub"]
    rows = []
    for name, d in runs.items():
        P = pd.read_csv(f"{d}/picks.csv")
        H = pd.read_csv(f"{d}/loc_headline.csv")
        okP = P.tP.notna() & (P.rP >= 1.5)
        okS = P.tS.notna() & (P.rS >= 1.5) if "tS" in P else pd.Series(False, index=P.index)
        per = (okP.astype(int) + okS.astype(int)).groupby(P.grp).sum()
        nsta = (okP | okS).groupby(P.grp).sum()
        # families Hyp2000 was given: >= 2 stations with a pick (write_hyp_inputs)
        n_in = int((nsta >= 2).sum())
        tmp = f"{od}/_truth_in_{name}.csv"
        H.rename(columns={"dep": "depth"})[["lon", "lat", "depth"]].to_csv(tmp, index=False)
        sc = {}
        for ts in ("Hyp2000", "hypoDD"):
            o = f"{od}/truth_{name}_{ts}.csv"
            subprocess.run([sys.executable, f"{HERE}/ab_truth_score.py", tmp, "--truth-set", ts,
                            "--subarrays", a.sub, "--out", o], check=True, capture_output=True)
            sc[ts] = pd.read_csv(o).iloc[0]
        S = sc["Hyp2000"]
        rows.append(dict(
            run=name, n_families=int(P.grp.nunique()), n_with_picks=n_in,
            picks_per_fam_med=float(per.median()), sta_with_pick_med=float(nsta.median()),
            n_located=len(H), dep_med=float(H.dep.median()),
            frac_31_46=float(H.dep.between(lo, hi).mean()),
            erh_mean=float(H.erh.mean()), erh_med=float(H.erh.median()),
            erz_mean=float(H.erz.mean()), erz_med=float(H.erz.median()),
            n_erh_lt5=int((H.erh < 5).sum()), n_dd_kept=int(H.dd_kept.sum()),
            truth_n_in_region=int(S.n_ours_in_region), truth_match_med_km=float(S.median_match_km),
            truth_frac3=float(S.frac_3km), truth_n3=int(round(S.frac_3km * S.n_truth)),
            truth_frac5=float(S.frac_5km), truth_frac10=float(S.frac_10km),
            truth_p3=float(S.p_3km), truth_p5=float(S.p_5km), truth_p10=float(S.p_10km),
            truth_null95_3=float(S.null_p95_3km), dz_med=float(S.dz_median_km),
            truthDD_n3=int(round(sc["hypoDD"].frac_3km * sc["hypoDD"].n_truth)),
            truthDD_p3=float(sc["hypoDD"].p_3km)))
    C = pd.DataFrame(rows)
    C.to_csv(f"{od}/allsta_comparison.csv", index=False)
    pd.set_option("display.width", 250)
    print(C.T.to_string())
    # ---- figure (as NW_first_catalogue.png) ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    H = pd.read_csv(f"{runs[a.main]}/loc_headline.csv")
    T = pd.read_csv(TRUTH)
    T = T[(T.subarray == a.sub) & (T.locator == "Hyp2000")]
    fig, ax = plt.subplots(1, 2, figsize=(14.5, 6))
    ax[0].scatter(T.lon, T.lat, s=22, facecolors="none", edgecolors="0.5",
                  label="Royer & Bostock 2014 (Fig.4)")
    bad = H.erh >= 5
    ax[0].scatter(H.lon[bad], H.lat[bad], s=14, c="tab:orange", label=f"ours, erh>=5 km (n={bad.sum()})")
    sc_ = ax[0].scatter(H.lon[~bad], H.lat[~bad], c=H.dep[~bad], s=45, cmap="viridis", vmin=25, vmax=50,
                        edgecolors="k", label=f"ours, erh<5 km (n={(~bad).sum()})")
    plt.colorbar(sc_, ax=ax[0], label="depth km")
    la = [s for s in A.POOL[a.sub]]
    ax[0].scatter([A.COORD[s][1] for s in la], [A.COORD[s][0] for s in la], marker="^", s=30,
                  c="tab:red", label="NW 25-station pool")
    ax[0].set_xlim(-124.1, -122.2)
    ax[0].set_ylim(46.6, 48.4)
    ax[0].legend(loc="lower left", fontsize=8)
    ax[0].set_title(f"{a.sub} all-station templates, Hyp2000 (GSC 1-D, Vp/Vs 1.73)")
    bins = np.arange(0, 62, 2)
    ax[1].axvspan(lo, hi, color="tab:green", alpha=0.12)
    ax[1].hist(H.dep, bins=bins, alpha=0.7, label=f"ours ({len(H)} located)")
    ax[1].hist(T.depth_km.dropna(), bins=bins, histtype="step", color="k", lw=1.2,
               label=f"Royer {a.sub}")
    ax[1].set_xlabel("depth km")
    ax[1].legend(loc="upper left")
    r = C[C.run == a.main].iloc[0]
    ax[1].set_title(f"depth median {r.dep_med:.1f} km, {100 * r.frac_31_46:.0f}% in {lo:.0f}-{hi:.0f} km; "
                    f"truth <=3 km {r.truth_n3}/{len(T)}", fontsize=10)
    fig.tight_layout()
    fig.savefig(a.fig, dpi=110)
    print("figure", a.fig)


if __name__ == "__main__":
    main()
