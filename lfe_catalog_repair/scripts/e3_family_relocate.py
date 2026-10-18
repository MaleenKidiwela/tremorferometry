#!/usr/bin/env python
"""E3 step 2, configs A and B — relocate the 129 Bostock families from the SeisBench stack picks.

WHY THIS IS NOT A NonLinLoc RUN (data fact found 2026-09-11, contradicting the audit's premise).
`seisbench_lfestacks_cascadia/metadata.csv` has `trace_P_predicted_arrival_sample` == 800.0 for all
750 traces that carry it and `trace_P_arrival_sample` == 800.0 for all 1,067 that carry it: EVERY
trace is individually aligned on its own P arrival at sample 800, there is no trace start time
anywhere in the dataset, and cross-station moveout has been destroyed by that alignment.  The only
timing information that survives per (family, station) is the **S-P differential time**
(S_sample - 800) / 20 s.  Writing NLL obs files with "P at a common origin + 40 s, S at
40 + (S-P)" would impose the false constraint that every P arrival is simultaneous, i.e. that the
family is equidistant from all stations, and would drag every solution to the station centroid.
So configs A and B are an S-P differential-time grid search on the SAME Crosson (1976) 1-D model
that `e2_locate_nll.py` uses; the origin time cancels identically.  Config C (the geometry /
S-only degeneracy term) DOES go through the exact NLL control, from synthetic times -- see
`e3_geom_synth.py`.

Config A: all stations with a BOTH-PICKED trace (P hand-picked at 800 AND S hand-picked), 20
  stations -- SHVB excluded (its metadata coordinates disagree with `path_ep_distance_km` by up to
  19 km where all 20 others agree to < 2 km).  Traces with only a PREDICTED P are excluded: their
  P sample is Bostock's own model prediction at the published hypocentre, so their S-P would be
  partly circular.
Config B: the same search masked to the six vi2010-overlap stations (PGC, LZB, VGZ, YOUB, SNB,
  PFB) -- reported, not gated on, because S-P at 6 stations is not E2's observable.

Elevation: `e2_locate_nll.py` places receivers at their true elevation (GTSRCE, grid top -3 km)
while the `e2_ttable_crosson.npz` table puts them at z = 0, so a vertical-incidence correction
elev * (1/Vs0 - 1/Vp0) is added to the predicted S-P.

Outputs results/e3_family_reloc_{A,B}.csv and results/e3_family_reloc_summary.csv.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
MD = f"{RES}/external/other/seisbench_lfestacks_cascadia/metadata.csv"
TT = f"{RES}/e2_ttable_crosson.npz"
R_KM_DEG = 111.195
OVERLAP = ["PGC", "LZB", "VGZ", "YOUB", "SNB", "PFB"]
BAD_STATIONS = ["SHVB"]          # metadata coordinates inconsistent with path_ep_distance_km
VP0, VPVS = 5.40, 1.78           # Crosson surface layer, for the elevation correction


def load_tt():
    z = np.load(TT)
    return z["x"], z["z"], z["P"], z["S"]


def interp2(x, zg, G, xq, zq):
    """bilinear interpolation of G(x, z) at (xq, zq), clipped to the table."""
    xq = np.clip(xq, x[0], x[-1]); zq = np.clip(zq, zg[0], zg[-1])
    i = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    j = np.clip(np.searchsorted(zg, zq) - 1, 0, len(zg) - 2)
    tx = (xq - x[i]) / (x[i + 1] - x[i])
    tz = (zq - zg[j]) / (zg[j + 1] - zg[j])
    return ((1 - tx) * (1 - tz) * G[i, j] + tx * (1 - tz) * G[i + 1, j]
            + (1 - tx) * tz * G[i, j + 1] + tx * tz * G[i + 1, j + 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vpvs", type=float, default=VPVS,
                    help="only used for the elevation correction; the table fixes the ratio")
    ap.add_argument("--half-km", type=float, default=60.0)
    ap.add_argument("--dh", type=float, default=1.0)
    ap.add_argument("--zmin", type=float, default=5.0)
    ap.add_argument("--zmax", type=float, default=70.0)
    ap.add_argument("--dz", type=float, default=1.0)
    ap.add_argument("--min-stations", type=int, default=3)
    args = ap.parse_args()

    x, zg, TP, TS = load_tt()
    print(f"[tt] Crosson table x {x[0]}..{x[-1]} km ({len(x)}), z {zg[0]}..{zg[-1]} km ({len(zg)});"
          f" finite P {np.isfinite(TP).mean():.3f}, S {np.isfinite(TS).mean():.3f}")

    m = pd.read_csv(MD)
    both = m[(m.trace_P_arrival_sample == 800) & (m.trace_S_arrival_sample.notna())].copy()
    both = both[~both.station_code.isin(BAD_STATIONS)]
    both["sp_obs"] = (both.trace_S_arrival_sample - 800.0) / both.trace_sampling_rate_hz
    # one row per (family, station): identical stacks are not duplicated in this dataset, but guard
    both = both.drop_duplicates(["source_family_id", "station_code"])
    print(f"[stacks] {len(both)} both-picked traces, {both.source_family_id.nunique()} families, "
          f"{both.station_code.nunique()} stations (SHVB excluded); "
          f"S-P {both.sp_obs.min():.2f}..{both.sp_obs.max():.2f} s")

    # elevation correction: the table has receivers at z=0, NLL/GTSRCE has them at elev
    VS0 = VP0 / args.vpvs
    both["elev_corr"] = (both.station_elevation_m / 1000.0) * (1.0 / VS0 - 1.0 / VP0)

    cos0 = np.cos(np.radians(48.6))
    rows = []
    for cfg, mask in (("A", np.ones(len(both), bool)),
                      ("B", both.station_code.isin(OVERLAP).to_numpy())):
        sub = both[mask]
        for fid, g in sub.groupby("source_family_id"):
            if len(g) < args.min_stations:
                continue
            lat0 = float(g.source_latitude_deg.iloc[0])
            lon0 = float(g.source_longitude_deg.iloc[0])
            dep0 = float(g.source_depth_km.iloc[0])
            # grid centred on the PUBLISHED epicentre (a search window, not a prior: +-60 km)
            n = int(2 * args.half_km / args.dh) + 1
            gx = np.linspace(-args.half_km, args.half_km, n)
            GX, GY = np.meshgrid(gx, gx, indexing="ij")
            glat = lat0 + GY / R_KM_DEG
            glon = lon0 + GX / (R_KM_DEG * cos0)
            gz = np.arange(args.zmin, args.zmax + args.dz / 2, args.dz)
            slat = g.station_latitude_deg.to_numpy()
            slon = g.station_longitude_deg.to_numpy()
            sp = g.sp_obs.to_numpy() - g.elev_corr.to_numpy()
            # horizontal distance station -> every node, (nsta, n, n)
            dx = (slon[:, None, None] - glon[None]) * R_KM_DEG * cos0
            dy = (slat[:, None, None] - glat[None]) * R_KM_DEG
            dh = np.hypot(dx, dy)
            best = (np.inf, None)
            for zq in gz:
                pred = interp2(x, zg, TS, dh, zq) - interp2(x, zg, TP, dh, zq)
                # L1 (Huber-like robustness without a scale parameter)
                mis = np.abs(pred - sp[:, None, None]).sum(axis=0) / len(sp)
                k = int(np.argmin(mis))
                if mis.flat[k] < best[0]:
                    ii, jj = np.unravel_index(k, mis.shape)
                    best = (float(mis.flat[k]), (float(glat[ii, jj]), float(glon[ii, jj]),
                                                 float(zq)))
            mis, (blat, blon, bz) = best
            ddh = np.hypot((blon - lon0) * R_KM_DEG * cos0, (blat - lat0) * R_KM_DEG)
            rows.append(dict(config=cfg, fam_id=float(fid), n_sta=len(g),
                             lat=blat, lon=blon, dep=bz,
                             lat_pub=lat0, lon_pub=lon0, dep_pub=dep0,
                             dh_km=ddh, dz_km=bz - dep0, l1_resid_s=mis,
                             gap_deg=float(g.source_azimuthal_gap_deg.iloc[0]),
                             h_unc_pub=float(g.source_horizontal_uncertainty_km.iloc[0]),
                             z_unc_pub=float(g.source_vertical_uncertainty_km.iloc[0]),
                             at_edge=bool(ddh > args.half_km - 2 * args.dh
                                          or bz <= args.zmin + args.dz
                                          or bz >= args.zmax - args.dz)))
        print(f"[config {cfg}] {sum(r['config']==cfg for r in rows)} families located")

    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/e3_family_reloc.csv", index=False)

    srows = []
    for cfg, g in out.groupby("config"):
        for lab, sel in (("all", g), ("n>=6", g[g.n_sta >= 6]),
                         ("n>=6, gap<200", g[(g.n_sta >= 6) & (g.gap_deg < 200)])):
            if len(sel) == 0:
                continue
            srows.append(dict(config=cfg, subset=lab, n_fam=len(sel),
                              dh_median=sel.dh_km.median(),
                              dh_q25=sel.dh_km.quantile(.25), dh_q75=sel.dh_km.quantile(.75),
                              dz_median=sel.dz_km.median(),
                              dz_q25=sel.dz_km.quantile(.25), dz_q75=sel.dz_km.quantile(.75),
                              abs_dz_median=sel.dz_km.abs().median(),
                              l1_median=sel.l1_resid_s.median(),
                              frac_at_edge=float(sel.at_edge.mean())))
    s = pd.DataFrame(srows)
    s.to_csv(f"{RES}/e3_family_reloc_summary.csv", index=False)
    print()
    print(s.round(2).to_string(index=False))
    print(f"\n[write] {RES}/e3_family_reloc.csv, e3_family_reloc_summary.csv")


if __name__ == "__main__":
    main()
