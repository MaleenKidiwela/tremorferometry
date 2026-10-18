#!/usr/bin/env python
"""STAGE 5 -- SCORE the Puget / S-Olympic beds exactly as the 2010 Olympic bed was scored.
Prefix `y_`.

RANK, DO NOT DISCARD (memory `rank-dont-discard`).  Every grown family keeps its row, every
criterion is a COLUMN, the site-shift null gets the IDENTICAL table, and the evidence is the
real-vs-null SEPARATION CURVE on each criterion.  No gate is applied anywhere in this file.

THE HEADLINE IS AN EXISTENCE TEST, and it does not need anybody's event list.  Royer & Bostock
report that LFE families EXIST in Washington on the plate interface at **29-46 km** depth.  That
claim is falsifiable on our side without their catalogue: run their data and either we produce
locatable, repeating, interface-depth families above the whole-chain null in that band, or we
have failed to reproduce them.  Their catalogue would only say WHICH families we recovered, never
WHETHER any are there.  So the first thing this file computes is: how many located families, real
and null, fall inside the bed's own sub-region AND inside 29-46 km, as a SEPARATION CURVE over
solution quality (bootstrap radius, cost, station count) -- never a count past a threshold.  The
2010 N-Olympic bed, which worked, is carried in the SAME tables as the positive control.

Reported, all with the null's identical row:

 (0) THE EXISTENCE TEST -- in-region / in-band counts real vs null as a curve over quality, with
     bootstrap radii, beside the 2010 positive control (`y_existence_<bed>.csv`,
     `y_existence_curve_<bed>.csv`);
 (a) real vs null separation on EVERY column (`y_separation_<bed>.csv`), both directions, over
     the full range of each criterion;
 (b) distance of located families to Chestler & Creager's 69 published Olympic coordinates, which
     is only meaningful where the Puget beds OVERLAP the Olympics -- reported for all families
     and again restricted to lat >= 47.3 (`y_chestler_<bed>.csv`);
 (c) the ALONG-STRIKE and DEPTH distribution against the Slab2 interface and the Royer & Bostock
     band.  That band and Slab2 are REFERENCE LINES drawn on the distribution, never gates: the
     fraction of real and of null families inside them is one more column.
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
import y_common as Y                                                        # noqa: E402

RB_DEPTH_BAND = (29.0, 46.0)        # Royer & Bostock's stated WA family depths: REFERENCE LINE
OLY_LAT_MIN = 47.30                 # where these beds overlap Chestler's Olympic coverage
NEAR = (5.0, 10.0, 20.0)
# the bed's own sub-region (n_common.LAT0..LON1) and the 2010 N-Olympic control's
REGION = {"2008": (46.2, 47.8, -125.0, -121.5), "2007": (46.2, 47.8, -125.0, -121.5),
          "2010": (47.3, 48.7, -125.0, -122.2)}
# the 2010 N-Olympic bed that WORKED, carried through the identical code as the positive control
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
CONTROL = {"2010 real (positive control)": f"{RES}/u_locate/wa/wa_real_ranked.csv",
           "2010 null (positive control)": f"{RES}/u_locate/wa/wa_null_ranked.csv"}

# criteria ranked HIGH-is-better (>=) and LOW-is-better (<=); both directions are curved
HI = ["n_sta", "n_pair", "n_det", "n_days", "med_cc", "max_nmad", "boot_n"]
LO = ["cost", "boot_km", "boot_dz", "h_med", "chest_km_3d", "chest_km_epi", "abs_dz",
      "cost_shuf_p05"]


def epi_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = (np.asarray(x, float) for x in (la1, lo1, la2, lo2))
    return np.sqrt((111.320 * np.cos(np.radians((la1 + la2) / 2)) * (lo1 - lo2)) ** 2
                   + (110.574 * (la1 - la2)) ** 2)


def hyp(la1, lo1, d1, la2, lo2, d2):
    return np.sqrt(epi_km(la1, lo1, la2, lo2) ** 2
                   + (np.asarray(d1, float) - np.asarray(d2, float)) ** 2)


def load_pop(path, det, screen):
    L = pd.read_csv(path)
    if det and os.path.exists(det):
        D = pd.read_parquet(det)
        L = L.merge(D.groupby("grp").size().rename("n_det"), on="grp", how="left")
        L = L.merge(D.groupby("grp").date.nunique().rename("n_days"), on="grp", how="left")
        if "cc" in D.columns:
            L = L.merge(D.groupby("grp").cc.median().rename("med_cc"), on="grp", how="left")
    if screen and os.path.exists(screen):
        S = pd.read_csv(screen)
        cols = [c for c in S.columns if c != "grp" and c not in L.columns]
        L = L.merge(S[["grp"] + cols], on="grp", how="left")
    if "dz" in L.columns:
        L["abs_dz"] = L.dz.abs()
    return L


def chest(L, C, tag, lat_min=None):
    ok = L[L.located == True].copy()                                        # noqa: E712
    if lat_min is not None:
        ok = ok[ok.lat >= lat_min]
    d = dict(pop=tag, lat_min=(lat_min if lat_min is not None else -90),
             n_families=len(L), n_located=int(len(ok)))
    if not len(ok):
        return d, ok
    d3 = np.array([hyp(r.lat, r.lon, r.dep, C.latitude.values, C.longitude.values,
                       C.depth.values).min() for r in ok.itertuples()])
    de = np.array([epi_km(r.lat, r.lon, C.latitude.values, C.longitude.values).min()
                   for r in ok.itertuples()])
    ok = ok.copy()
    ok["chest_km_3d"], ok["chest_km_epi"] = d3, de
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
    Cs = C[C.latitude >= lat_min] if lat_min is not None else C
    for k in NEAR:
        hit = 0
        for r in Cs.itertuples():
            if (hyp(ok.lat.values, ok.lon.values, ok.dep.values,
                    r.latitude, r.longitude, r.depth) <= k).any():
                hit += 1
        d[f"coords_hit_{k:.0f}km"] = hit
    d["n_coords"] = len(Cs)
    return d, ok


def depth_strike(L, tag):
    """Along-strike (latitude) and depth distribution vs Slab2 and the R&B 28-45 km band."""
    ok = L[L.located == True]                                               # noqa: E712
    d = dict(pop=tag, n_located=int(len(ok)))
    if not len(ok):
        return d, pd.DataFrame(), pd.DataFrame()
    dep = ok.dep.to_numpy(float)
    d.update(dep_med=float(np.median(dep)), dep_p10=float(np.quantile(dep, .1)),
             dep_p90=float(np.quantile(dep, .9)),
             frac_in_RB_band=float(((dep >= RB_DEPTH_BAND[0])
                                    & (dep <= RB_DEPTH_BAND[1])).mean()),
             n_in_RB_band=int(((dep >= RB_DEPTH_BAND[0]) & (dep <= RB_DEPTH_BAND[1])).sum()),
             lat_med=float(np.median(ok.lat)), lat_p10=float(np.quantile(ok.lat, .1)),
             lat_p90=float(np.quantile(ok.lat, .9)))
    if "dz" in ok.columns and ok.dz.notna().any():
        d.update(dz_med=float(np.nanmedian(ok.dz)),
                 frac_abs_dz_le8=float((ok.dz.abs() <= 8).mean()),
                 frac_abs_dz_le5=float((ok.dz.abs() <= 5).mean()))
    db = np.arange(6, 80, 2.0)
    H = pd.DataFrame(dict(pop=tag, dep_lo=db[:-1], dep_hi=db[1:],
                          n=np.histogram(dep, bins=db)[0]))
    lb = np.arange(46.0, 49.2, 0.1)
    Ls = pd.DataFrame(dict(pop=tag, lat_lo=lb[:-1], lat_hi=lb[1:],
                           n=np.histogram(ok.lat.to_numpy(float), bins=lb)[0]))
    return d, H, Ls


def _in_region(L, box):
    la0, la1, lo0, lo1 = box
    return ((L.lat >= la0) & (L.lat <= la1) & (L.lon >= lo0) & (L.lon <= lo1))


def existence_row(L, tag, box):
    """THE EXISTENCE TEST, one row.  Royer & Bostock state that families exist on the interface
    at 29-46 km in Washington; this asks whether THIS population puts locatable families there.
    Nothing is cut: the row is counts and quantiles over the whole population."""
    d = dict(pop=tag, n_families=int(len(L)))
    ok = L[L.located == True]                                               # noqa: E712
    d["n_located"] = int(len(ok))
    if not len(ok):
        return d
    reg = ok[_in_region(ok, box)]
    d["n_in_region"] = int(len(reg))
    band = reg[(reg.dep >= RB_DEPTH_BAND[0]) & (reg.dep <= RB_DEPTH_BAND[1])]
    d["n_in_RB_band"] = int(len(band))
    d["frac_located_in_RB_band"] = float(len(band) / max(len(ok), 1))
    if "edge" in ok.columns:
        d["rail"] = float(ok.edge.mean())
        d["n_in_band_not_railed"] = int((~band.edge.astype(bool)).sum()) if len(band) else 0
    for col in ("boot_km", "cost"):
        if col in band.columns and len(band) and band[col].notna().any():
            d[f"band_{col}_med"] = float(np.nanmedian(band[col]))
            d[f"band_{col}_p10"] = float(np.nanquantile(band[col], 0.1))
    if "ctrl_pass" in band.columns and len(band):
        d["band_ctrl_pass_frac"] = float(band.ctrl_pass.fillna(False).astype(bool).mean())
    for t in (10.0, 20.0, 40.0):
        if "boot_km" in band.columns and len(band):
            d[f"n_in_band_boot_le{t:.0f}"] = int((band.boot_km <= t).sum())
        else:
            d[f"n_in_band_boot_le{t:.0f}"] = 0
    return d


def existence_curve(pops, box_of, out):
    """The existence test AS A CURVE, never as a count past a threshold: for every value of a
    quality criterion, how many located, in-region, in-band families each population has."""
    rows = []
    for tag, (L, box) in pops.items():
        ok = L[L.located == True]                                           # noqa: E712
        if not len(ok):
            continue
        reg = ok[_in_region(ok, box)]
        band = reg[(reg.dep >= RB_DEPTH_BAND[0]) & (reg.dep <= RB_DEPTH_BAND[1])]
        for col, sense, grid in (("boot_km", "<=", np.linspace(0, 60, 61)),
                                 ("cost", "<=", np.linspace(0, 4, 41)),
                                 ("n_sta", ">=", np.arange(3, 31, 1.0))):
            if col not in band.columns:
                continue
            v = band[col].astype(float)
            for t in grid:
                n = int((v <= t).sum()) if sense == "<=" else int((v >= t).sum())
                rows.append(dict(pop=tag, criterion=col, sense=sense, threshold=float(t),
                                 n_in_band=n, n_in_region=int(len(reg)),
                                 n_located=int(len(ok))))
    C = pd.DataFrame(rows)
    C.to_csv(out, index=False)
    return C


def separation(R, U, out):
    """Real-vs-null counts as a function of EVERY criterion, both directions.  No thresholds are
    used anywhere -- the whole curve is the evidence."""
    rows = []
    for col in HI + LO:
        if col not in R.columns or col not in U.columns:
            continue
        v = pd.concat([R[col], U[col]]).astype(float).dropna()
        if not len(v) or v.nunique() < 2:
            continue
        lo, hi = float(v.min()), float(v.max())
        for t in np.linspace(lo, hi, 41):
            if col in HI:
                r = int((R[col].astype(float) >= t).sum())
                u = int((U[col].astype(float) >= t).sum())
                sense = ">="
            else:
                r = int((R[col].astype(float) <= t).sum())
                u = int((U[col].astype(float) <= t).sum())
                sense = "<="
            rows.append(dict(criterion=col, sense=sense, threshold=float(t), real=r, null=u,
                             excess=r - u, ratio=(r / u) if u else np.nan))
    S = pd.DataFrame(rows)
    S.to_csv(out, index=False)
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--null", required=True)
    ap.add_argument("--det-real", default="")
    ap.add_argument("--det-null", default="")
    ap.add_argument("--screen-real", default="")
    ap.add_argument("--screen-null", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tag = a.tag or Y.BED
    C = pd.read_csv(Y.O.CHESTLER)
    R = load_pop(a.real, a.det_real, a.screen_real)
    U = load_pop(a.null, a.det_null, a.screen_null)
    R.insert(0, "population", "real")
    U.insert(0, "population", "null")

    crows, cfull = [], {}
    for name, L in (("real", R), ("null", U)):
        for lm in (None, OLY_LAT_MIN):
            d, ok = chest(L, C, name, lm)
            crows.append(d)
            if lm is None:
                cfull[name] = ok
    CH = pd.DataFrame(crows)
    CH.to_csv(f"{a.out}/y_chestler_{tag}.csv", index=False)

    # every family keeps its row; chestler distance is one more COLUMN
    for name, L in (("real", R), ("null", U)):
        ok = cfull.get(name)
        full = (L.merge(ok[["grp", "chest_km_3d", "chest_km_epi"]], on="grp", how="left")
                if ok is not None and len(ok) else L.assign(chest_km_3d=np.nan,
                                                            chest_km_epi=np.nan))
        full.to_csv(f"{a.out}/y_{name}_ranked_{tag}.csv", index=False)
        if name == "real":
            R = full
        else:
            U = full

    drows, HH, LL = [], [], []
    for name, L in (("real", R), ("null", U)):
        d, H, Ls = depth_strike(L, name)
        drows.append(d)
        HH.append(H)
        LL.append(Ls)
    pd.DataFrame(drows).to_csv(f"{a.out}/y_depth_strike_{tag}.csv", index=False)
    pd.concat([x for x in HH if len(x)]).to_csv(f"{a.out}/y_depth_hist_{tag}.csv", index=False)
    pd.concat([x for x in LL if len(x)]).to_csv(f"{a.out}/y_strike_hist_{tag}.csv", index=False)

    S = separation(R, U, f"{a.out}/y_separation_{tag}.csv")

    # ---------------------------------------------------------------------------------------
    # (0) THE EXISTENCE TEST, with the 2010 N-Olympic bed that WORKED in the SAME table as the
    # positive control.  Royer & Bostock's 29-46 km band and Slab2 are REFERENCE LINES on the
    # distribution; the evidence is the real-vs-null curve, not a count past a threshold.
    # ---------------------------------------------------------------------------------------
    box = REGION.get(Y.BED, REGION["2008"])
    pops = {f"{Y.BED} real": (R, box), f"{Y.BED} null (site-shift)": (U, box)}
    for name, p in CONTROL.items():
        if os.path.exists(p):
            try:
                pops[name] = (pd.read_csv(p), REGION["2010"])
            except Exception:                                               # noqa: BLE001
                pass
    E = pd.DataFrame([existence_row(L, t, b) for t, (L, b) in pops.items()])
    E.to_csv(f"{a.out}/y_existence_{tag}.csv", index=False)
    EC = existence_curve(pops, None, f"{a.out}/y_existence_curve_{tag}.csv")

    meta = dict(bed=tag, n_real=len(R), n_null=len(U),
                n_real_located=int((R.located == True).sum()),                # noqa: E712
                n_null_located=int((U.located == True).sum()),                # noqa: E712
                rb_depth_band_reference_only=list(RB_DEPTH_BAND),
                region_box=list(box), oly_overlap_lat_min=OLY_LAT_MIN,
                n_criteria_curved=int(S.criterion.nunique()) if len(S) else 0,
                existence=E.to_dict("records"))
    json.dump(meta, open(f"{a.out}/y_score_{tag}.json", "w"), indent=1, default=float)
    with pd.option_context("display.width", 260, "display.max_columns", 80):
        print("\nTHE EXISTENCE TEST (2010 N-Olympic carried as the positive control)\n"
              + E.round(3).to_string(index=False), flush=True)
        print("\nCHESTLER COMPARISON\n" + CH.round(3).to_string(index=False), flush=True)
        print("\nDEPTH / ALONG-STRIKE\n"
              + pd.DataFrame(drows).round(3).to_string(index=False), flush=True)
        print("\n" + json.dumps(meta, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
