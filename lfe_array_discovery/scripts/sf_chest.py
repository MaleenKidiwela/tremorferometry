#!/usr/bin/env python
"""The CHESTLER PROXIMITY test, run on the population Chestler's own method actually built.

`WA_2010_ARRAY_RESULTS.md` sec.5 ran this check on agent N's five-clause SCREEN survivors (0 of
them) and on all interior-located grown templates of arm A.  Neither is the population the
comparison is about.  This module runs it on the **deduped + located** families, real against the
matched null chain, for whichever arm it is pointed at.

The 69 rows of `derived_coordinate_rows_not_a_catalog.csv` are, as the filename says, NOT a
catalogue; they are coordinates derived from the published Olympic literature.  The only defensible
statement is a COMPARISON: does the real population sit closer to them than the null population
does, at the same sample size?  Reported as (i) the nearest-coordinate distance distribution,
(ii) counts within 5/10/20 km as a fraction of families, (iii) how many of the 69 coordinates have
a family within 5/10/20 km, real vs null, and (iv) a size-matched bootstrap so the two populations
are compared at equal n.

  python scripts/sf_chest.py --real ..._accept.csv --null ..._accept.csv --out PREFIX --tag ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import sf_common as SF                                                      # noqa: E402

CHESTLER = ("/home/jovyan/tremorferometry/codex-LFE/results/olympic_reference_discovery_v1/"
            "derived_coordinate_rows_not_a_catalog.csv")
NEAR = (5.0, 10.0, 20.0)
N_BOOT = 2000


def load_accept(p):
    A = pd.read_csv(p)
    F = A[A.ACCEPT.astype(bool)].copy()
    return F[np.isfinite(F.lat)]


def near_tab(F, C):
    cla, clo, cdp = (C.latitude.to_numpy(float), C.longitude.to_numpy(float),
                     C.depth.to_numpy(float))
    rows = []
    for r in F.itertuples():
        e = SF.epi_km(float(r.lat), float(r.lon), cla, clo)
        d3 = np.sqrt(e ** 2 + (float(r.dep) - cdp) ** 2)
        j = int(np.argmin(d3))
        rows.append(dict(grp=int(r.grp), lat=float(r.lat), lon=float(r.lon), dep=float(r.dep),
                         n_det=int(r.n_det), n_days=int(r.n_days), edge=bool(r.edge),
                         ctrl_pass=bool(r.ctrl_pass),
                         near_km_3d=float(d3.min()), near_km_epi=float(e.min()),
                         nearest_row=int(C.iloc[j].iloc[0])))
    return pd.DataFrame(rows)


def coord_hits(F, C, km):
    cla, clo, cdp = (C.latitude.to_numpy(float), C.longitude.to_numpy(float),
                     C.depth.to_numpy(float))
    hit = np.zeros(len(C), bool)
    for r in F.itertuples():
        e = SF.epi_km(float(r.lat), float(r.lon), cla, clo)
        d3 = np.sqrt(e ** 2 + (float(r.dep) - cdp) ** 2)
        hit |= d3 <= km
    return int(hit.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--null", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    C = pd.read_csv(CHESTLER)
    FR, FN = load_accept(a.real), load_accept(a.null)
    TR, TN = near_tab(FR, C), near_tab(FN, C)
    TR.to_csv(f"{a.out}_real.csv", index=False)
    TN.to_csv(f"{a.out}_null.csv", index=False)
    rng = np.random.default_rng(7)

    def stats(T, F):
        d = dict(n=len(T))
        for km in NEAR:
            d[f"n_within_{km:.0f}km"] = int((T.near_km_3d <= km).sum())
            d[f"frac_within_{km:.0f}km"] = float((T.near_km_3d <= km).mean()) if len(T) else 0.0
            d[f"n_coords_hit_{km:.0f}km"] = coord_hits(F, C, km)
        d["median_near_km_3d"] = float(T.near_km_3d.median()) if len(T) else float("nan")
        d["median_near_km_epi"] = float(T.near_km_epi.median()) if len(T) else float("nan")
        d["min_near_km_3d"] = float(T.near_km_3d.min()) if len(T) else float("nan")
        return d

    R = dict(tag=a.tag, n_chestler_coords=int(len(C)),
             real=stats(TR, FR), null=stats(TN, FN))
    # size-matched bootstrap: the larger population subsampled to the smaller one's size
    if len(TR) and len(TN):
        big, small = (TR, TN) if len(TR) >= len(TN) else (TN, TR)
        k = len(small)
        boot = {f"{km:.0f}": [] for km in NEAR}
        med = []
        for _ in range(N_BOOT):
            s = big.iloc[rng.choice(len(big), size=k, replace=False)]
            for km in NEAR:
                boot[f"{km:.0f}"].append(int((s.near_km_3d <= km).sum()))
            med.append(float(s.near_km_3d.median()))
        R["bootstrap"] = dict(
            resampled=("real" if len(TR) >= len(TN) else "null"), n=k,
            median_km_mean=float(np.mean(med)),
            **{f"n_within_{km:.0f}km_mean": float(np.mean(boot[f'{km:.0f}'])) for km in NEAR},
            **{f"p_ge_other_{km:.0f}km": float(np.mean(
                np.array(boot[f"{km:.0f}"]) >= int((small.near_km_3d <= km).sum())))
                for km in NEAR})
    json.dump(R, open(f"{a.out}.json", "w"), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
