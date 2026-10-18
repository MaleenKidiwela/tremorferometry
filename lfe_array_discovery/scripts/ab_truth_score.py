#!/usr/bin/env python
"""Score our located LFE families against the digitized Royer & Bostock (2014)
template locations (lfe_array_discovery/results/ab_royer/truth/).

Usage:
  python ab_truth_score.py OUR.csv [--truth-set Hyp2000|hypoDD|thesis_family]
                                   [--subarrays NW,SW,NVI] [--nnull 1000] [--seed 0]
OUR.csv needs columns lon, lat and optionally depth (km; column 'depth' or 'depth_km').

Per subarray:
  * candidate set = our families within REGION_KM (30 km) of any truth template
    of that subarray (so a neighbouring subarray does not steal matches);
  * one-to-one Hungarian matching on epicentral distance (km);
  * matched fraction of TRUTH templates within 3 / 5 / 10 km, median distance;
  * depth residual (ours - truth) for matched pairs within 10 km that both have depth;
  * NULL: the same candidate set randomly rotated (0-360 deg about its centroid)
    and shifted (uniform in a disc of radius = RMS radius of the truth cloud),
    re-matched; reports null mean, null 95th percentile and the one-sided
    p = P(null >= real) for each radius. Real and null are reported side by side.
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = os.path.join(HERE, "../results/ab_royer/truth/royer2014_templates_digitized.csv")
RADII = (3.0, 5.0, 10.0)
REGION_KM = 30.0
R_EARTH = 6371.0


def to_xy(lon, lat, lon0, lat0):
    x = np.radians(lon - lon0) * R_EARTH * np.cos(np.radians(lat0))
    y = np.radians(lat - lat0) * R_EARTH
    return x, y


def match(tx, ty, ox, oy):
    if len(ox) == 0:
        return np.full(len(tx), np.inf), np.full(len(tx), -1)
    D = np.hypot(tx[:, None] - ox[None, :], ty[:, None] - oy[None, :])
    i, j = linear_sum_assignment(D)
    d = np.full(len(tx), np.inf)
    jj = np.full(len(tx), -1)
    d[i] = D[i, j]
    jj[i] = j
    return d, jj


def score(truth, ours, subarrays, nnull, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for sub in subarrays:
        t = truth[truth.subarray == sub]
        if len(t) == 0:
            continue
        lon0, lat0 = t.lon.mean(), t.lat.mean()
        tx, ty = to_xy(t.lon.values, t.lat.values, lon0, lat0)
        ox, oy = to_xy(ours.lon.values, ours.lat.values, lon0, lat0)
        dmin = np.hypot(ox[:, None] - tx[None, :], oy[:, None] - ty[None, :]).min(1)
        sel = dmin <= REGION_KM
        cx, cy = ox[sel], oy[sel]
        oz = ours.depth.values[sel] if "depth" in ours else np.full(sel.sum(), np.nan)
        d, jj = match(tx, ty, cx, cy)
        row = dict(subarray=sub, n_truth=len(t), n_ours_in_region=int(sel.sum()),
                   median_match_km=float(np.median(d[np.isfinite(d)])) if np.isfinite(d).any() else np.nan)
        for r in RADII:
            row[f"frac_{r:g}km"] = float(np.mean(d <= r))
        ok = (d <= 10) & (jj >= 0)
        dz = oz[jj[ok]] - t.depth_km.values[ok]
        dz = dz[np.isfinite(dz)]
        row["n_depth_pairs"] = len(dz)
        row["dz_median_km"] = float(np.median(dz)) if len(dz) else np.nan
        row["dz_mad_km"] = float(np.median(np.abs(dz - np.median(dz)))) if len(dz) else np.nan
        # null
        R = np.sqrt(np.mean(tx ** 2 + ty ** 2))
        mx, my = (cx.mean(), cy.mean()) if len(cx) else (0, 0)
        null = np.zeros((nnull, len(RADII)))
        for k in range(nnull):
            th = rng.uniform(0, 2 * np.pi)
            rr = R * np.sqrt(rng.uniform())
            ph = rng.uniform(0, 2 * np.pi)
            x0, y0 = cx - mx, cy - my
            nx = np.cos(th) * x0 - np.sin(th) * y0 + mx + rr * np.cos(ph)
            ny = np.sin(th) * x0 + np.cos(th) * y0 + my + rr * np.sin(ph)
            dn, _ = match(tx, ty, nx, ny)
            null[k] = [np.mean(dn <= r) for r in RADII]
        for i, r in enumerate(RADII):
            row[f"null_mean_{r:g}km"] = float(null[:, i].mean())
            row[f"null_p95_{r:g}km"] = float(np.percentile(null[:, i], 95))
            row[f"p_{r:g}km"] = float((1 + np.sum(null[:, i] >= row[f"frac_{r:g}km"])) / (nnull + 1))
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ours")
    ap.add_argument("--truth", default=TRUTH)
    ap.add_argument("--truth-set", default="Hyp2000", choices=["Hyp2000", "hypoDD", "thesis_family"])
    ap.add_argument("--subarrays", default="NW,SW,NVI")
    ap.add_argument("--nnull", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None, help="optional csv for the table")
    a = ap.parse_args()
    truth = pd.read_csv(a.truth)
    truth = truth[truth.locator == a.truth_set]
    ours = pd.read_csv(a.ours)
    if "depth" not in ours and "depth_km" in ours:
        ours = ours.rename(columns={"depth_km": "depth"})
    res = score(truth, ours, a.subarrays.split(","), a.nnull, a.seed)
    pd.set_option("display.width", 250)
    print(f"truth set = {a.truth_set}; our families = {len(ours)}")
    cols_real = ["subarray", "n_truth", "n_ours_in_region", "median_match_km"] + \
        [f"frac_{r:g}km" for r in RADII] + ["n_depth_pairs", "dz_median_km", "dz_mad_km"]
    print("REAL:\n", res[cols_real].round(3).to_string(index=False))
    cols_null = ["subarray"] + [c for r in RADII for c in (f"frac_{r:g}km", f"null_mean_{r:g}km",
                                                           f"null_p95_{r:g}km", f"p_{r:g}km")]
    print("REAL vs NULL (rotated+shifted our locations):\n", res[cols_null].round(3).to_string(index=False))
    if a.out:
        res.to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
