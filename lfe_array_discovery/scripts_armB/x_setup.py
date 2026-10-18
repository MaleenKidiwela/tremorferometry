#!/usr/bin/env python
"""Step 0 (a) -- build the WASHINGTON N-Olympic 2010 BED and write `results/step5_bed.json`.

Everything the Step-5 chain needs that is DECIDED rather than measured is decided here, once,
before the freeze is written and before a single WA byte is scored:

  * the admitted site set (the Step-4 site probe's admitted stations, minus the CL array, minus
    PFB, minus the withheld set), merged at 1 km into sites BY MEMBER STATION;
  * the withheld set (P_STEP2_SITES.md sec.5 Option B: B001, the B005 cluster B005/B006/B007/LRIV,
    and VGZ) -- never in any autocorrelation sum, growth sum or stack;
  * the moveout spread of the admitted sites at the sub-region centre, against the frozen 10-s
    window.  If it exceeds AC_WIN_S the SAME written rule used for Step 1b/1c is applied: sites are
    added in order of increasing median epicentral distance to the active cells while the MEDIAN
    S-moveout spread over the active cells stays <= AC_WIN_S;
  * the active cells, the candidate (tremor) hours, the quiet days and quiet hours;
  * the HF-veto station (a WA 100-Hz borehole that the probe admits).

The HF cut itself is measured by `x_hf_cut.py`, which updates this file before the freeze.

    python scripts/x_setup.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x_common as X                                                        # noqa: E402

RES, CACHE = X.RES, X.CACHE
CL_ARRAY_PREFIX = "CL"          # XG.CL: 0/16 admitted, median contrast 1.9 dB (P_STEP2 sec.3.1)
DROP_STATIONS = ["PFB"]         # CN.PFB fails local-far (-10.6 dB); the S-WA "far" days are not
                                # far from Vancouver Island (P_STEP2 sec.3.1)


def cached_stations(year=X.YEAR):
    import glob
    out = {}
    for p in glob.glob(f"{CACHE}/*/{year}/*.npz"):
        ns = p.split(os.sep)[-3]
        out.setdefault(ns.split(".")[1], set()).add(int(os.path.basename(p)[:3]))
    return out


def single_linkage(stas, ST, km=1.0):
    import s_common as S
    stas = [s for s in stas if s in ST.index]
    n = len(stas)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            d = float(S.epi_km(ST.lat[stas[i]], ST.lon[stas[i]],
                               ST.lat[stas[j]], ST.lon[stas[j]]))
            if d <= km:
                parent[find(i)] = find(j)
    grp = {}
    for i in range(n):
        grp.setdefault(find(i), []).append(stas[i])
    return [sorted(v) for v in grp.values()]


def main():
    os.makedirs(X.OUT, exist_ok=True)
    ST = X.wa_station_table()
    import s_common as S
    import g_common as G

    # ---------------------------------------------------------------- admitted stations
    P = pd.read_csv(f"{RES}/site_probe_2010.csv")
    have = cached_stations()
    window = pd.date_range(X.DAY0, X.DAY1)
    win_doy = {int(d.dayofyear) for d in window}
    adm = []
    for r in P.itertuples():
        s = str(r.sta)
        rec = dict(sta=s, net=str(r.net), admit=bool(r.admit),
                   contrast_db=float(r.contrast_db), local_minus_far_db=float(r.local_minus_far_db),
                   n_days_window=len(have.get(s, set()) & win_doy))
        drop = None
        if s in X.WITHHELD_STATIONS:
            drop = "withheld"
        elif s in DROP_STATIONS:
            drop = "probe_far_clause"
        elif s.startswith(CL_ARRAY_PREFIX) and str(r.net) == "XG":
            drop = "XG.CL array flagged (<5 admitted members)"
        elif not bool(r.admit):
            drop = "probe_not_admitted"
        elif rec["n_days_window"] == 0:
            drop = "no cached window day"
        elif s not in ST.index:
            drop = "no station metadata"
        rec["drop_reason"] = drop
        adm.append(rec)
    A = pd.DataFrame(adm)
    A.to_csv(f"{X.OUT}/admitted_stations.csv", index=False)
    admitted = sorted(A[A.drop_reason.isna()].sta)
    X.log(f"[setup] {len(admitted)} admitted stations of {len(A)} probed "
          f"({A.drop_reason.value_counts().to_dict()})")

    # ---------------------------------------------------------------- sites (1-km merge)
    groups = single_linkage(admitted, ST, 1.0)
    cells = X.active_cells()
    X.log(f"[setup] {len(cells)} active cells in {X.SUBREGION} "
          f"(Slab2 {X.CELL_DEP[0]:.0f}-{X.CELL_DEP[1]:.0f} km, >= {X.CELL_MIN_EPI} epicentres "
          f"within {X.CELL_KM:.0f} km, {X.CELL_YEARS[0]}-{X.CELL_YEARS[1]})")
    xs, zs, TP, TS = G.tt_table()
    sites = []
    for g in groups:
        key = sorted(g, key=lambda s: (-len(have.get(s, set()) & win_doy), s))[0]
        lat = float(np.median([ST.lat[s] for s in g]))
        lon = float(np.median([ST.lon[s] for s in g]))
        d = S.epi_km(cells.lat.to_numpy(), cells.lon.to_numpy(), lat, lon)
        tt = G._interp2(xs, zs, TS, d, cells.dep.to_numpy())
        sites.append(dict(site=key, members=g, n_members=len(g), lat=lat, lon=lon,
                          epi_med=float(np.median(d)),
                          n_cells_60km=int((np.hypot(d, cells.dep.to_numpy()) <= 60.0).sum()),
                          tt=tt))
    sites.sort(key=lambda r: r["epi_med"])

    # ---------------------------------------------------------------- moveout spread
    clat, clon = 0.5 * (X.LAT0 + X.LAT1), 0.5 * (X.LON0 + X.LON1)
    cdep = float(G.slab2_depth([clat], [clon])[0])
    dc = np.array([float(S.epi_km(clat, clon, r["lat"], r["lon"])) for r in sites])
    ttc = G._interp2(xs, zs, TS, dc, np.full(len(dc), cdep))
    spread_centre = float(ttc.max() - ttc.min())
    TTC = np.vstack([r["tt"] for r in sites])                 # (n_site, n_cell)
    spread_cells_all = float(np.median(TTC.max(0) - TTC.min(0)))
    X.log(f"[setup] centre {clat:.3f}/{clon:.3f} Slab2 {cdep:.1f} km -> admitted-site moveout "
          f"spread {spread_centre:.2f} s (median over cells {spread_cells_all:.2f} s); "
          f"AC_WIN_S = 10.0 s")

    rule = "all admitted sites (spread at the sub-region centre <= AC_WIN_S)"
    keep = list(range(len(sites)))
    if spread_centre > 10.0:
        # The Step-1b/1c rule, verbatim: add sites in order of increasing median epicentral
        # distance while the MEDIAN S-moveout spread over the active cells stays <= AC_WIN_S.
        rule = ("spread rule (Step 1b/1c): sites added by increasing median epicentral distance "
                "while the median S-moveout spread over the active cells stays <= AC_WIN_S")
        keep = []
        for i in range(len(sites)):
            trial = keep + [i]
            if len(trial) < 2:
                keep = trial
                continue
            M = TTC[trial]
            if float(np.median(M.max(0) - M.min(0))) <= 10.0:
                keep = trial
        X.log(f"[setup] spread rule applied: {len(keep)}/{len(sites)} sites kept")
    dropped = [dict(site=sites[i]["site"], members=sites[i]["members"],
                    epi_med_km=sites[i]["epi_med"])
               for i in range(len(sites)) if i not in set(keep)]
    if dropped:
        X.log("[setup] dropped by the spread rule: "
              + ", ".join(f"{d['site']}({d['epi_med_km']:.0f} km)" for d in dropped))
    sites = [sites[i] for i in keep]
    TTC = np.vstack([r["tt"] for r in sites])
    dc = np.array([float(S.epi_km(clat, clon, r["lat"], r["lon"])) for r in sites])
    ttc = G._interp2(xs, zs, TS, dc, np.full(len(dc), cdep))
    spread_final = float(ttc.max() - ttc.min())
    spread_cells = float(np.median(TTC.max(0) - TTC.min(0)))

    # ---------------------------------------------------------------- coverage
    hyp = np.zeros((len(sites), len(cells)))
    for i, r in enumerate(sites):
        d = S.epi_km(cells.lat.to_numpy(), cells.lon.to_numpy(), r["lat"], r["lon"])
        hyp[i] = np.hypot(d, cells.dep.to_numpy())
    n60 = (hyp <= 60.0).sum(0)
    cov = dict(n_cells=int(len(cells)),
               pct_cells_ge5_sites_60km=float(100.0 * (n60 >= 5).mean()),
               median_sites_60km=float(np.median(n60)))
    # withheld-site visibility (fraction of cells within 60 km hypocentral)
    vis = {}
    for q, mm in X.WITHHELD_SITES.items():
        m0 = [m for m in mm if m in ST.index]
        if not m0:
            continue
        d = S.epi_km(cells.lat.to_numpy(), cells.lon.to_numpy(),
                     float(ST.lat[m0[0]]), float(ST.lon[m0[0]]))
        vis[q] = float((np.hypot(d, cells.dep.to_numpy()) <= 60.0).mean())

    # ---------------------------------------------------------------- HF-veto station
    hf_sta = None
    for s in X.HF_STATION_CANDIDATES:
        row = P[P.sta == s]
        if len(row) and bool(row.admit.iloc[0]) and float(row.fs_native.iloc[0]) >= 100.0:
            hf_sta = s
            break
    if hf_sta is None:
        raise SystemExit("no probe-admitted 100-Hz WA borehole for the HF veto")
    X.log(f"[setup] HF-veto station = PB.{hf_sta} "
          f"(probe contrast {float(P[P.sta == hf_sta].contrast_db.iloc[0]):.1f} dB, "
          f"{float(P[P.sta == hf_sta].fs_native.iloc[0]):.0f} Hz native)")

    # ---------------------------------------------------------------- quiet days
    FW = json.load(open(f"{RES}/freeze_windows.json"))
    quiet = [d for d in FW["quiet_days"]["y2010"]["days"]]
    cfg = dict(
        study="lfe_array_discovery / Step 5 -- N-Olympic 2010 discovery bed",
        subregion=X.SUBREGION, lat0=X.LAT0, lat1=X.LAT1, lon0=X.LON0, lon1=X.LON1,
        year=X.YEAR, day0=X.DAY0, day1=X.DAY1,
        centre=dict(lat=clat, lon=clon, slab2_dep_km=cdep),
        grid=dict(lat0=X.GRID_LAT0, lon0=X.GRID_LON0, halfwidth_km=X.GRID_HALFWIDTH),
        active_cell_rule=(f"Slab2 interface points in the sub-region box, depth "
                          f"{X.CELL_DEP[0]:.0f}-{X.CELL_DEP[1]:.0f} km, with >= {X.CELL_MIN_EPI} "
                          f"PNSN {X.CELL_YEARS[0]}-{X.CELL_YEARS[1]} epicentres within "
                          f"{X.CELL_KM:.0f} km"),
        n_active_cells=int(len(cells)),
        admitted_rule=("Step-4 site probe (>= 4/6 days >= +3 dB contrast AND local-far >= +3 dB), "
                       "minus the XG.CL array (0/16 admitted), minus CN.PFB (far clause), minus "
                       "the withheld set, and requiring >= 1 cached window day"),
        n_admitted_stations=len(admitted),
        site_merge_km=1.0,
        site_rule=rule,
        sites_dropped_by_spread_rule=dropped,
        site_order=[r["site"] for r in sites],
        site_members={r["site"]: r["members"] for r in sites},
        site_table=[dict(site=r["site"], members=r["members"], n_members=r["n_members"],
                         lat=r["lat"], lon=r["lon"], epi_med_km=r["epi_med"],
                         n_cells_60km=r["n_cells_60km"]) for r in sites],
        n_sites=len(sites),
        n_detection_stations=int(sum(len(r["members"]) for r in sites)),
        moveout_spread_s_centre_all_admitted=spread_centre if rule.startswith("all") else None,
        moveout_spread_s_centre=spread_final,
        moveout_spread_s_median_over_cells=spread_cells,
        ac_win_s=10.0,
        withheld_stations=list(X.WITHHELD_STATIONS),
        withheld_sites=X.WITHHELD_SITES,
        withheld_rule=("P_STEP2_SITES.md sec.5 Option B, adopted in PLAN sec.10: B001 + the B005 "
                       "cluster (B005/B006/B007/LRIV) + VGZ are barred from EVERY autocorrelation "
                       "sum, growth sum and stack for the whole run; they are registered post hoc "
                       "for the extension"),
        withheld_cell_visibility_frac=vis,
        coverage=cov,
        quiet_days=quiet, n_quiet_days=len(quiet),
        quiet_rule=FW["quiet_days"]["y2010"]["rule"],
        hf_station=hf_sta, hf_veto=None,
        hf_veto_note="set by x_hf_cut.py before the freeze",
    )
    json.dump(cfg, open(X.BED_JSON, "w"), indent=1)
    cells.to_csv(f"{X.OUT}/active_cells.csv", index=False)
    X.log(f"[setup] {len(sites)} sites / {cfg['n_detection_stations']} detection stations; "
          f"spread {spread_final:.2f} s; coverage {cov['pct_cells_ge5_sites_60km']:.1f} % of "
          f"{cov['n_cells']} cells with >= 5 sites <= 60 km (median {cov['median_sites_60km']:.0f})")
    X.log(f"[setup] withheld visibility {vis}")
    X.log(f"[setup] wrote {X.BED_JSON}")


if __name__ == "__main__":
    main()
