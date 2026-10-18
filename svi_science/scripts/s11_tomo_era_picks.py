#!/usr/bin/env python
"""S11 tomo (Amendment 8): picks on the SVI 2003-06 era stacks.
  python s11_tomo_era_picks.py -> tables/picks_EallV3.csv, picks_E60V3.csv, picks_E40V3.csv (hypocentral km limits)
Frame: t_det + offset + i/fs (offset = ttS_model - 16 s), common to all stations of a family."""
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import resample_poly

ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s11_tomo_svi_picks import picks_for  # noqa: E402

T = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"


def main():
    I = pd.read_csv(f"{T}/era_stack_index.csv")
    off = {(int(r.fam), r.netsta): float(r.offset) for r in I.itertuples()}
    nrow = {(int(r.fam), r.netsta): int(r.n) for r in I.itertuples()}
    ST = {}
    for p in sorted(glob.glob(f"{T}/era_stacks/*.npz")):
        ns = os.path.basename(p)[:-4]
        z = np.load(p)
        for f in z.files:
            ST.setdefault(int(f), {})[ns] = resample_poly(z[f].astype(np.float64), 2, 1, axis=-1)
    P = picks_for(ST, None)
    o = np.array([off[(g, s)] for g, s in zip(P.grp, P.sta)])
    for c in ("tP", "tS", "tP_v0"):
        P[c] = P[c] + o
    P["frame_offset"] = o
    P["n_win"] = [nrow[(g, s)] for g, s in zip(P.grp, P.sta)]
    F = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv").set_index("bkey")
    S = pd.read_csv(f"{ROOT}/svi_science/results/backward_scope/inventory_stations.csv")
    S["ns"] = S.net + "." + S.sta
    S = S.drop_duplicates("ns").set_index("ns")
    P["fam"] = P.grp.astype(str)
    P["lat"], P["lon"], P["depth_km"] = P.grp.map(F.lat), P.grp.map(F.lon), P.grp.map(F.dep)
    P["n_detections"] = P.grp.map(I.groupby("fam").n.max())
    slat, slon = P.sta.map(S.lat), P.sta.map(S.lon)
    P["epi_km"] = np.hypot((slon - P.lon) * 111.195 * np.cos(np.radians(P.lat)), (slat - P.lat) * 111.195)
    P["hyp_km"] = np.hypot(P.epi_km, P.depth_km)
    P["src"] = "era stacks 2003-06"
    P.to_csv(f"{T}/picks_EallV3.csv", index=False)
    for lim in (60, 40):
        Q = P[P.hyp_km <= lim]
        Q.to_csv(f"{T}/picks_E{lim}V3.csv", index=False)
        ok = (Q.rP >= 1.2) & (Q.rS >= 1.2)
        n = Q[ok].groupby("fam").size()
        print(f"[E{lim}V3] rows {len(Q)}, ok {ok.sum()}, stations {Q.sta.nunique()}; families with >=3/4/6 ok rows: "
              f"{(n >= 3).sum()}/{(n >= 4).sum()}/{(n >= 6).sum()}", flush=True)
    ok = (P.rP >= 1.2) & (P.rS >= 1.2)
    print(f"[EallV3] {P.fam.nunique()} families, {len(P)} rows, ok {ok.sum()}", flush=True)


if __name__ == "__main__":
    main()
