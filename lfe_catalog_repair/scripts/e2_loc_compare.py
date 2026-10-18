#!/usr/bin/env python
"""E2 step 3b-score (added 2026-09-10, resumed run) — NonLinLoc locations vs the catalogues.

Inputs: `<events>.nll.csv` (from e2_locate_nll.py) + the event file it came from.
Comparisons
  1. TIME-MATCHED events: for every located event, the nearest Lin (2024) event within +-2 s and the
     nearest Bostock (2015) event within +-2 s (the Bostock times are shifted by the empirical
     equivalent-origin offset measured at PGC, see e2_assoc_score.bostock_offset) -> horizontal and
     depth offsets, median and IQR, split by the number of stations.
  2. NEAREST FAMILY: distance from every located epicentre to the nearest of the 122 distinct
     Bostock family hypocentres (a lower bound on the location error that needs no time match, and
     the number the "family hypocentre" comparison of the plan asks for).
  3. Internal: NonLinLoc vs PyOcto epicentre for the same event, and the formal NLL errors.
  4. Depth: how many located depths land in the 25-45 km interface band, and the S-only caveat
     (the formal ERZ).
Writes results/e2_loc_compare_<tag>.csv (per-event) and prints the summary table.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_assoc_score as AS  # noqa: E402

R_KM_DEG = 111.195
COSLAT = np.cos(np.radians(48.6))


def km(lat1, lon1, lat2, lon2):
    return np.hypot((np.asarray(lon2) - np.asarray(lon1)) * R_KM_DEG * COSLAT,
                    (np.asarray(lat2) - np.asarray(lat1)) * R_KM_DEG)


def nearest(gt_t, ev_t, tol):
    """index into gt of the nearest gt time for each ev time (or -1)."""
    o = np.argsort(gt_t); g = gt_t[o]
    i = np.searchsorted(g, ev_t)
    best = np.full(len(ev_t), -1); bd = np.full(len(ev_t), np.inf)
    for off in (-1, 0):
        j = np.clip(i + off, 0, len(g) - 1)
        d = np.abs(g[j] - ev_t)
        m = d < bd
        bd = np.where(m, d, bd); best = np.where(m, o[j], best)
    best[bd > tol] = -1
    return best, bd


def stats(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return dict(n=0)
    return dict(n=len(v), median=float(np.median(v)), q25=float(np.percentile(v, 25)),
                q75=float(np.percentile(v, 75)), mean=float(v.mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nll", required=True)
    ap.add_argument("--tol", type=float, default=2.0)
    ap.add_argument("--bostock-family-offsets", default="",
                    help="E3 step 5: e3_bostock_family_offsets.csv -> per-family Bostock times")
    ap.add_argument("--null-nll", default="",
                    help="E3 step 1 (Merlin F4): a LOCATED SHIFTED-PICK NULL catalogue "
                         "(<...>.nll.csv).  Its nearest-family distribution replaces the "
                         "uniform-over-the-bounding-box null, which is wrong because S-only "
                         "3-station solutions collapse toward the network centroid -- exactly "
                         "where the Bostock families are.  The uniform null is still printed, "
                         "labelled as the legacy number.")
    ap.add_argument("--tag-suffix", default="",
                    help="appended to the output tag so an E3 run cannot overwrite an E2 file")
    args = ap.parse_args()
    tag = os.path.basename(args.nll)
    tag = tag[tag.index("_") + 1:-len(".nll.csv")] if tag.startswith("e3_") \
        else tag[len("e2_events_"):-len(".nll.csv")]
    tag += args.tag_suffix
    loc = pd.read_csv(args.nll)
    ev = pd.read_csv(args.nll.replace(".nll.csv", ".csv"))
    loc = loc.merge(ev[["idx", "n_stations", "n_picks", "mean_prob", "latitude", "longitude",
                        "depth", "time"]].rename(columns={"latitude": "lat_oct",
                                                          "longitude": "lon_oct",
                                                          "depth": "dep_oct",
                                                          "time": "t_oct"}),
                    left_on="event_idx", right_on="idx", how="left")
    lt = pd.to_datetime(loc.time_utc).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
    print(f"[{tag}] {len(loc)} located events; NLL err_h median {loc.err_h_km.median():.1f} km, "
          f"err_z median {loc.err_z_km.median():.1f} km; rms median {loc.rms.median():.2f} s")
    print(f"[{tag}] depth: median {loc.depth.median():.1f} km, frac<5 km "
          f"{float((loc.depth<5).mean()):.3f}, frac 25-45 km {float(loc.depth.between(25,45).mean()):.3f}; "
          f"Slab2 under epicentre median {loc.slab2_km.median():.1f} km")
    dnll_oct = km(loc.lat_oct, loc.lon_oct, loc.lat, loc.lon)
    print(f"[{tag}] NonLinLoc vs PyOcto epicentre: median {np.nanmedian(dnll_oct):.1f} km "
          f"(IQR {np.nanpercentile(dnll_oct,25):.1f}-{np.nanpercentile(dnll_oct,75):.1f})")

    # ---- catalogues
    lin = pd.read_parquet(f"{RES}/e2_gt_lin.parquet")
    lin["t"] = AS.epoch(lin.OT)
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    bo["t"], blab, _ = AS.bostock_times(bo, args.bostock_family_offsets)
    fam = bo.groupby(["lat", "lon", "dep"]).size().reset_index(name="n")
    print(f"[gt] Lin {len(lin)}, Bostock {len(bo)} events / {len(fam)} distinct family hypocentres; "
          f"Bostock equivalent-origin offset: {blab}")

    rows = loc.copy()
    for name, cat, latc, lonc, depc in (("lin", lin, "lat", "lon", "depth"),
                                        ("bostock", bo, "lat", "lon", "dep")):
        j, d = nearest(cat.t.to_numpy(), lt, args.tol)
        ok = j >= 0
        rows[f"dt_{name}"] = np.where(ok, d, np.nan)
        for c, src in ((f"dh_{name}", None), (f"dz_{name}", None)):
            pass
        la = np.where(ok, cat[latc].to_numpy()[np.clip(j, 0, None)], np.nan)
        lo = np.where(ok, cat[lonc].to_numpy()[np.clip(j, 0, None)], np.nan)
        dp = np.where(ok, np.abs(cat[depc].to_numpy()[np.clip(j, 0, None)]), np.nan)
        rows[f"dh_{name}"] = km(la, lo, loc.lat.values, loc.lon.values)
        rows[f"dz_{name}"] = loc.depth.values - dp
        s = stats(rows[f"dh_{name}"]); sz = stats(rows[f"dz_{name}"])
        if s.get("n", 0):
            print(f"[vs {name}] {s['n']} of {len(loc)} matched within {args.tol} s: horizontal "
                  f"{s['median']:.1f} km (IQR {s['q25']:.1f}-{s['q75']:.1f}); depth diff "
                  f"{sz['median']:+.1f} km (IQR {sz['q25']:+.1f}..{sz['q75']:+.1f})")
    # nearest Bostock family hypocentre, no time match
    D = np.full(len(loc), np.nan)
    for k in range(len(loc)):
        D[k] = km(fam.lat.values, fam.lon.values, loc.lat.values[k], loc.lon.values[k]).min()
    rows["dh_nearest_family"] = D
    s = stats(D)
    # NULL for the "nearest family" number: the 130 families blanket the region, so a random
    # epicentre is already close to one.  Random points drawn uniformly over the bounding box of
    # the located epicentres give the distance that means nothing.
    rng = np.random.default_rng(1)
    rlat = rng.uniform(loc.lat.min(), loc.lat.max(), 4000)
    rlon = rng.uniform(loc.lon.min(), loc.lon.max(), 4000)
    Dn = np.array([km(fam.lat.values, fam.lon.values, rlat[k], rlon[k]).min()
                   for k in range(len(rlat))])
    sn = stats(Dn)
    print(f"[vs nearest family] epicentre to nearest of {len(fam)} Bostock family hypocentres: "
          f"median {s['median']:.1f} km (IQR {s['q25']:.1f}-{s['q75']:.1f}); "
          f"LEGACY uniform-random-point null median {sn['median']:.1f} km "
          f"(IQR {sn['q25']:.1f}-{sn['q75']:.1f})  [Merlin F4: this null is WRONG]")
    # E3 step 1 / F4: the DEFENSIBLE null = shifted-pick chance events located identically
    sh = None
    if args.null_nll and os.path.exists(args.null_nll):
        nl = pd.read_csv(args.null_nll)
        Ds = np.array([km(fam.lat.values, fam.lon.values, nl.lat.values[k], nl.lon.values[k]).min()
                       for k in range(len(nl))])
        sh = stats(Ds)
        print(f"[vs nearest family] SHIFTED-PICK NULL ({len(nl)} located chance events, "
              f"{os.path.basename(args.null_nll)}): median {sh['median']:.1f} km "
              f"(IQR {sh['q25']:.1f}-{sh['q75']:.1f}); observed/null median ratio "
              f"{s['median']/max(sh['median'],1e-9):.2f}")
        print(f"[null depth] shifted-pick null depth: median {nl.depth.median():.1f} km, "
              f"frac<5 km {float((nl.depth<5).mean()):.3f}, frac 25-45 km "
              f"{float(nl.depth.between(25,45).mean()):.3f}  (observed: "
              f"{loc.depth.median():.1f} km, {float((loc.depth<5).mean()):.3f}, "
              f"{float(loc.depth.between(25,45).mean()):.3f})")
        pd.DataFrame(dict(dh_nearest_family=Ds, depth=nl.depth.values,
                          err_h_km=nl.err_h_km.values, lat=nl.lat.values, lon=nl.lon.values)
                     ).to_csv(f"{RES}/e3_loc_null_{tag}.csv", index=False)
    # by station count
    print("[by n_stations] " + "; ".join(
        f"n={int(n)}: N={len(g)} dh_lin {np.nanmedian(g.dh_lin):.1f} km, "
        f"dh_fam {np.nanmedian(g.dh_nearest_family):.1f} km, errh {g.err_h_km.median():.0f} km"
        for n, g in rows.groupby("n_stations") if len(g) >= 10))
    out = (f"{RES}/e3_loc_compare_{tag}.csv" if (args.bostock_family_offsets or args.null_nll
                                                 or args.tag_suffix)
           else f"{RES}/e2_loc_compare_{tag}.csv")
    rows.to_csv(out, index=False)
    print(f"[write] {out}")


if __name__ == "__main__":
    main()
