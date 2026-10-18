#!/usr/bin/env python
"""STAGE 4 -- the CHESTLER CHECK.  Prefix `o_`.

`../codex-LFE/results/olympic_reference_discovery_v1/derived_coordinate_rows_not_a_catalog.csv`
is, as its filename says, NOT a catalogue: it is 69 coordinate rows derived from the published
Olympic LFE literature (Chestler & Creager 2017 / Sweet et al.).  It is used here as the
REPRODUCTION TARGET in exactly two descriptive ways, both with their own null:

 A. for every SCREENED family, the distance to the nearest derived coordinate, against the
    distribution obtained by drawing the same number of points from the sub-region's ACTIVE
    CELLS (the places tremor actually happens) 1,000 times -- plus the counts within 5 and 10 km;
 B. "where are their families in our output" -- for every derived coordinate, whether ANY grown
    template (screened or not) locates within 10 km of it, with the site-shift null's templates
    put through the identical question as the false-hit rate.

Nothing here is a rediscovery count in Bostock's sense (there is no event-level Washington
catalogue to score one-to-one against); it is a coordinate-proximity statement with a null.

  python scripts/o_chest_check.py --screen ... --loc ... --loc-null ... --out ...
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
import o_common as O                                                        # noqa: E402

NEAR_KM = (5.0, 10.0, 20.0)
N_DRAW = 1000


def hyp(la1, lo1, dp1, la2, lo2, dp2):
    import s_common as S
    d = S.epi_km(la1, lo1, la2, lo2)
    return np.sqrt(d ** 2 + (np.asarray(dp1) - np.asarray(dp2)) ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", required=True)
    ap.add_argument("--loc", required=True, help="locations of ALL grown real templates")
    ap.add_argument("--loc-null", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    O.activate()
    import s_common as S
    C = pd.read_csv(O.CHESTLER)
    cla = C.latitude.to_numpy(float)
    clo = C.longitude.to_numpy(float)
    cdp = C.depth.to_numpy(float)
    AC = O.active_cells()
    rng = np.random.default_rng(5)

    def nearest(la, lo, dp):
        d = np.sqrt(S.epi_km(la, lo, cla, clo) ** 2 + (dp - cdp) ** 2)
        e = S.epi_km(la, lo, cla, clo)
        return float(d.min()), float(e.min()), int(np.argmin(d))

    V = pd.read_csv(a.screen)
    col = "candidate" if "candidate" in V.columns else "pass_abcd"
    F = V[V[col].astype(bool) & np.isfinite(V.get("lat", pd.Series(dtype=float)))]
    rows = []
    for r in F.itertuples():
        d3, d2, j = nearest(float(r.lat), float(r.lon), float(r.dep))
        rows.append(dict(grp=int(r.grp), lat=float(r.lat), lon=float(r.lon), dep=float(r.dep),
                         near_km_3d=d3, near_km_epi=d2, nearest_row=int(C.iloc[j].iloc[0])))
    A_tab = pd.DataFrame(rows)

    # the null distribution: the same number of points drawn from the active cells
    null_med = []
    if len(A_tab):
        dep = AC.dep.to_numpy(float)
        for _ in range(N_DRAW):
            k = rng.integers(0, len(AC), len(A_tab))
            d = [nearest(float(AC.lat.iloc[i]), float(AC.lon.iloc[i]), float(dep[i]))[0]
                 for i in k]
            null_med.append(float(np.median(d)))
    A_sum = dict(
        n_screened_located=int(len(A_tab)),
        median_km_3d=float(A_tab.near_km_3d.median()) if len(A_tab) else None,
        min_km_3d=float(A_tab.near_km_3d.min()) if len(A_tab) else None,
        random_median_km=float(np.median(null_med)) if null_med else None,
        random_p05_km=float(np.quantile(null_med, .05)) if null_med else None,
        frac_random_at_least_as_close=(float(np.mean(np.array(null_med)
                                                     <= A_tab.near_km_3d.median()))
                                       if null_med else None),
        **{f"n_within_{int(k)}km": int((A_tab.near_km_3d <= k).sum()) if len(A_tab) else 0
           for k in NEAR_KM})

    # B: per Chestler coordinate, is ANY grown template within 10 km?
    def cover(locpath):
        if not locpath or not os.path.exists(locpath):
            return None
        L = pd.read_csv(locpath)
        L = L[L.located.astype(bool)] if "located" in L.columns else L
        if "edge" in L.columns:
            L = L[~L.edge.astype(bool)]
        if not len(L):
            return pd.DataFrame(), 0
        la, lo, dp = (L.lat.to_numpy(float), L.lon.to_numpy(float), L.dep.to_numpy(float))
        out = []
        for i in range(len(C)):
            d = np.sqrt(S.epi_km(cla[i], clo[i], la, lo) ** 2 + (cdp[i] - dp) ** 2)
            j = int(np.argmin(d))
            out.append(dict(row=int(C.iloc[i].iloc[0]), lat=cla[i], lon=clo[i], dep=cdp[i],
                            nearest_grp=int(L.grp.iloc[j]), nearest_km=float(d.min()),
                            hit_10km=bool(d.min() <= 10.0), hit_5km=bool(d.min() <= 5.0),
                            n_within_10km=int((d <= 10).sum())))
        return pd.DataFrame(out), len(L)

    Breal = cover(a.loc)
    Bnull = cover(a.loc_null)
    B_tab, n_real = Breal if Breal else (pd.DataFrame(), 0)
    Bn_tab, n_null = Bnull if Bnull else (pd.DataFrame(), 0)
    B_sum = dict(n_chestler_coords=int(len(C)),
                 n_interior_located_real=int(n_real), n_interior_located_null=int(n_null),
                 n_coords_hit_10km_real=int(B_tab.hit_10km.sum()) if len(B_tab) else 0,
                 n_coords_hit_5km_real=int(B_tab.hit_5km.sum()) if len(B_tab) else 0,
                 n_coords_hit_10km_null=int(Bn_tab.hit_10km.sum()) if len(Bn_tab) else 0,
                 n_coords_hit_5km_null=int(Bn_tab.hit_5km.sum()) if len(Bn_tab) else 0,
                 median_nearest_km_real=(float(B_tab.nearest_km.median())
                                         if len(B_tab) else None),
                 median_nearest_km_null=(float(Bn_tab.nearest_km.median())
                                         if len(Bn_tab) else None))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    A_tab.to_csv(a.out.replace(".json", "_screened.csv"), index=False)
    if len(B_tab):
        B_tab.to_csv(a.out.replace(".json", "_coverage_real.csv"), index=False)
    if len(Bn_tab):
        Bn_tab.to_csv(a.out.replace(".json", "_coverage_null.csv"), index=False)
    rep = dict(A=A_sum, B=B_sum)
    json.dump(rep, open(a.out, "w"), indent=1, default=float)
    print(json.dumps(rep, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
