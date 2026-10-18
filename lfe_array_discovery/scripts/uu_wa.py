#!/usr/bin/env python
"""Apply the rebuilt locator to Washington and redo the Chestler comparison.  Prefix `uu_`.

Three populations, the SAME chain and the SAME table for each (memory `rank-dont-discard`):
  real   results/o_oly2010/2010/        -- the arrays-as-arrays 2010 Olympic discovery run
  null   results/o_oly2010/2010_null/   -- its whole-chain site-shift null
  armB   results/o_oly2010/chestler_arm/-- Chestler & Creager's OWN detection chain, run literally

EVERY family keeps its row; every criterion is a COLUMN (n_sta, n_det, n_days, cost, edge flag,
bootstrap radius, |dz| to Slab2, distance to the nearest published Olympic coordinate).  No cut
is applied anywhere in this file.  The comparison with Chestler's coordinates is reported as
distance distributions and as counts within 5/10/20 km, real beside null beside Arm B.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uu_common as U                                                       # noqa: E402

NEAR = (5.0, 10.0, 20.0)


def hyp(la1, lo1, dp1, la2, lo2, dp2):
    return np.sqrt(U.epi_km(la1, lo1, la2, lo2) ** 2
                   + (np.asarray(dp1, float) - np.asarray(dp2, float)) ** 2)


def nearest(L, C):
    """For every located family, distance to the nearest published coordinate (3-D and epi)."""
    out3, oute, oid = [], [], []
    for r in L.itertuples():
        d3 = hyp(r.lat, r.lon, r.dep, C.latitude.values, C.longitude.values, C.depth.values)
        de = U.epi_km(r.lat, r.lon, C.latitude.values, C.longitude.values)
        i = int(np.argmin(d3))
        out3.append(float(d3[i]))
        oute.append(float(de.min()))
        oid.append(int(C.unidentified_source_row.iloc[i]))
    return np.array(out3), np.array(oute), np.array(oid)


def summarise(tag, L, C):
    ok = L[L.located == True].copy()                                        # noqa: E712
    d = dict(pop=tag, n_families=len(L), n_located=int(len(ok)))
    if not len(ok):
        return d, ok
    d3, de, cid = nearest(ok, C)
    ok["chest_km_3d"] = d3
    ok["chest_km_epi"] = de
    ok["chest_row"] = cid
    d.update(rail=float(ok.edge.mean()),
             med_nearest_3d=float(np.median(d3)), med_nearest_epi=float(np.median(de)),
             p10_nearest_3d=float(np.quantile(d3, 0.1)))
    if "boot_km" in ok.columns and ok.boot_km.notna().any():
        d["boot_km_med"] = float(np.nanmedian(ok.boot_km))
        d["boot_km_p90"] = float(np.nanquantile(ok.boot_km, 0.9))
    if "ctrl_pass" in ok.columns:
        d["ctrl_pass_frac"] = float(ok.ctrl_pass.fillna(False).astype(bool).mean())
    for k in NEAR:
        d[f"fam_within_{k:.0f}km"] = int((d3 <= k).sum())
        d[f"fam_within_{k:.0f}km_frac"] = float((d3 <= k).mean())
    # reverse: how many published coordinates have a family near them
    for k in NEAR:
        hit = 0
        for r in C.itertuples():
            if (hyp(ok.lat.values, ok.lon.values, ok.dep.values,
                    r.latitude, r.longitude, r.depth) <= k).any():
                hit += 1
        d[f"coords_hit_{k:.0f}km"] = hit
    d["n_coords"] = len(C)
    return d, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--null", required=True)
    ap.add_argument("--armb", default="")
    ap.add_argument("--det-real", default="")
    ap.add_argument("--det-null", default="")
    ap.add_argument("--det-armb", default="")
    ap.add_argument("--sf-real", default=f"{U.RES}/sf_final/table_oly2010_real.csv")
    ap.add_argument("--sf-null", default=f"{U.RES}/sf_final/table_oly2010_null.csv")
    ap.add_argument("--sf-col", default="dedupe_n_partners")
    ap.add_argument("--sf-min", type=float, default=20.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    C = pd.read_csv(U.CHESTLER)
    SF = {}
    for tag, sp in (("real", a.sf_real), ("null", a.sf_null)):
        if sp and os.path.exists(sp):
            T = pd.read_csv(sp)
            keep = [c for c in ("grp", a.sf_col, "n_det", "n_days", "med_cc", "max_nmad",
                                "rate_tremor", "tremor_quiet_ratio", "co_excess", "dedupe_kept",
                                "tierA", "tierB", "tierC", "lat", "lon", "dep",
                                "near_chestler_km") if c in T.columns]
            T = T[keep].rename(columns={"lat": "old_lat", "lon": "old_lon", "dep": "old_dep",
                                        "near_chestler_km": "old_near_chestler_km",
                                        "n_det": "sf_n_det", "n_days": "sf_n_days"})
            SF[tag] = T
    rows, tabs, full_tabs = [], {}, {}
    for tag, p, dp in (("real", a.real, a.det_real), ("null", a.null, a.det_null),
                       ("armB", a.armb, a.det_armb)):
        if not p or not os.path.exists(p):
            continue
        L = pd.read_csv(p)
        if tag in SF:
            # every sf_final ranking criterion becomes a COLUMN here; nothing is filtered
            L = L.merge(SF[tag], on="grp", how="left", suffixes=("", "_sf"))
            L["sf_top"] = L[a.sf_col].fillna(-1) >= a.sf_min
        if dp and os.path.exists(dp):
            D = pd.read_parquet(dp)
            n = D.groupby("grp").size().rename("n_det")
            nd = D.groupby("grp").date.nunique().rename("n_days")
            L = L.merge(n, on="grp", how="left").merge(nd, on="grp", how="left")
        d, ok = summarise(tag, L, C)
        rows.append(d)
        tabs[tag] = ok
        full = L.merge(ok[["grp", "chest_km_3d", "chest_km_epi", "chest_row"]],
                       on="grp", how="left") if len(ok) else L
        full.insert(0, "population", tag)
        full.to_csv(f"{a.out}/wa_{tag}_ranked.csv", index=False)
        full_tabs[tag] = full
    S = pd.DataFrame(rows)
    S.to_csv(f"{a.out}/wa_chestler_summary.csv", index=False)
    # separation curves: real-vs-null counts as a function of each criterion (no thresholds used)
    if "real" in tabs and "null" in tabs:
        cur = []
        for col, lo, hi, n in (("chest_km_3d", 0, 40, 41), ("cost", 0, 4, 41),
                               ("n_sta", 5, 22, 18), ("boot_km", 0, 40, 41)):
            if col not in tabs["real"].columns:
                continue
            for t in np.linspace(lo, hi, n):
                if col in ("n_sta",):
                    r = int((tabs["real"][col] >= t).sum())
                    u = int((tabs["null"][col] >= t).sum())
                else:
                    r = int((tabs["real"][col] <= t).sum())
                    u = int((tabs["null"][col] <= t).sum())
                cur.append(dict(criterion=col, threshold=float(t), real=r, null=u,
                                excess=r - u, ratio=(r / u) if u else np.nan))
        pd.DataFrame(cur).to_csv(f"{a.out}/wa_separation_curves.csv", index=False)
    # ---------------------------------------------------------------------------------------
    # THE PRIORITY QUESTION.  results/sf_final/table_oly2010_{real,null}.csv report 280 REAL
    # 2010-Olympic templates with >= 20 dedupe partners and ZERO null templates that rich
    # (the null maximum is 16).  Those 280 are the strongest candidate family population
    # Washington has produced.  Here they are RE-LOCATED with the rebuilt locator and placed
    # against Chestler's published coordinates, beside the null's best-ranked rows -- the
    # equal-size head of the null ordered by the SAME criterion, which is the matched
    # comparator, not a threshold.  This selects nothing: the 280 keep their rows in
    # wa_real_ranked.csv and carry `sf_top` as one more column.
    # ---------------------------------------------------------------------------------------
    if "real" in full_tabs and a.sf_col in full_tabs["real"].columns:
        srows, subs = [], {}
        R = full_tabs["real"]
        sel = R[R.sf_top == True]                                           # noqa: E712
        subs[f"real[{a.sf_col}>={a.sf_min:.0f}]"] = sel
        if "null" in full_tabs and a.sf_col in full_tabs["null"].columns:
            N = full_tabs["null"].sort_values(a.sf_col, ascending=False)
            subs[f"null[top {len(sel)} by {a.sf_col}]"] = N.head(len(sel))
            subs[f"null[{a.sf_col}>={a.sf_min:.0f}]"] = N[N[a.sf_col] >= a.sf_min]
        subs["real[all]"] = R
        if "null" in full_tabs:
            subs["null[all]"] = full_tabs["null"]
        for tag, T in subs.items():
            d, ok = summarise(tag, T, C)
            if a.sf_col in T.columns:
                d[f"{a.sf_col}_min"] = (float(T[a.sf_col].min()) if len(T) else np.nan)
                d[f"{a.sf_col}_med"] = (float(T[a.sf_col].median()) if len(T) else np.nan)
            srows.append(d)
            if tag.startswith(f"real[{a.sf_col}"):
                ok = ok.copy()
                ok["population"] = tag
                ok.to_csv(f"{a.out}/wa_sf_top_real_located.csv", index=False)
        SS = pd.DataFrame(srows)
        SS.to_csv(f"{a.out}/wa_sf_top_summary.csv", index=False)
        json.dump(srows, open(f"{a.out}/wa_sf_top_summary.json", "w"), indent=1)
        with pd.option_context("display.width", 240, "display.max_columns", 60):
            U.log(f"\nTHE {len(sel)} RICHEST REAL 2010-OLYMPIC TEMPLATES vs the null's best "
                  f"ranked rows, against Chestler's {len(C)} published coordinates\n"
                  + SS.round(2).to_string(index=False))

    json.dump(rows, open(f"{a.out}/wa_chestler_summary.json", "w"), indent=1)
    with pd.option_context("display.width", 240, "display.max_columns", 60):
        U.log("\n" + S.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
