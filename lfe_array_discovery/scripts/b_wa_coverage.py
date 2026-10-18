#!/usr/bin/env python
"""WA preparation (1b) -- COVERAGE RECOMPUTED without the unhealthy stations.  Prefix `b_`.

Merlin: "if any XG array has > 20 % toned/degenerate station-hours, the site-probe admission was
blind to it and the WA coverage numbers in `results/P_STEP2_SITES.md` are overstated" -- then
recompute coverage excluding stations with > 30 % vetoed hours.

Two numbers are recomputed, both with the SAME code the originals used
(`d_coverage.site_frame` / `d_coverage.cover`, and the freeze's own site table):

  (a) the P_STEP2_SITES sub-region table (admitted sites, >= 5 sites <= 60 km hypocentral),
  (b) the freeze_wa_v1 N-Olympic 11-site coverage over its 145 active cells.

-> results/wa_coverage_health_adjusted_<year>.csv + .json
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
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
import d_common as D                                                        # noqa: E402
import d_coverage as C                                                      # noqa: E402

VETO_PCT_DROP = 30.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True, choices=[2008, 2010])
    a = ap.parse_args()
    S = pd.read_csv(f"{RES}/wa_station_veto_summary_{a.year}.csv")
    drop = set(S[S.pct_veto > VETO_PCT_DROP].netsta)
    print(f"[drop] {len(drop)} stations with > {VETO_PCT_DROP:.0f} % vetoed hours: "
          f"{sorted(drop)}")
    G = D.active_cells()
    A = G[G.tremor_active].reset_index(drop=True)
    ST = pd.read_csv(f"{RES}/sites_array_era_final.csv", keep_default_na=False)
    probe = pd.read_csv(f"{RES}/site_probe_{a.year}.csv")
    admitted = set(probe[probe.admit].netsta)
    wcol = "w2010_08_days_3c" if a.year == 2010 else "w2008_05_days_3c"
    subs = (["N-Olympic", "S-Olympic"] if a.year == 2010 else ["Puget-SWA", "S-Olympic"])
    rows = []
    for tag, adm in (("as_published", admitted), ("health_adjusted", admitted - drop)):
        SF = C.site_frame(ST, adm, wcol)
        for sub in subs:
            lo, hi = C.SUBREG[sub]
            cells = A[(A.lat >= lo) & (A.lat <= hi)].reset_index(drop=True)
            adm_s = SF[SF.admitted]
            fh, fe, nh, _ = C.cover(cells, adm_s)
            rows.append(dict(variant=tag, year=a.year, subregion=sub, n_cells=len(cells),
                             n_sites_admitted=len(adm_s),
                             pct_cells_ge5_sites_60km=round(100 * fh, 1),
                             pct_cells_ge5_sites_30km=round(100 * fe, 1),
                             median_sites_60km=float(np.median(nh)) if len(nh) else 0.0))
    R = pd.DataFrame(rows)
    out = f"{RES}/wa_coverage_health_adjusted_{a.year}.csv"
    R.to_csv(out, index=False)
    print(R.to_string(index=False))

    meta = dict(year=a.year, dropped_stations=sorted(drop), veto_pct_drop=VETO_PCT_DROP)
    # ---- (b) the freeze_wa_v1 11-site N-Olympic coverage over its own 145 cells
    if a.year == 2010 and os.path.exists(f"{RES}/freeze_wa_v1.json"):
        F = json.load(open(f"{RES}/freeze_wa_v1.json"))
        cells = pd.read_csv(f"{RES}/step5/active_cells.csv")
        cells = cells.rename(columns={"dep": "depth_km"})
        keep_codes = {s.split(".")[-1] for s in drop}
        surv, died = [], []
        for s in F["sites"]["site_table"]:
            mem = [m for m in s["members"] if m not in keep_codes]
            (surv if mem else died).append(dict(site=s["site"], lat=s["lat"], lon=s["lon"],
                                                n_members=len(mem)))
        SV = pd.DataFrame(surv)
        AL = pd.DataFrame([dict(site=s["site"], lat=s["lat"], lon=s["lon"])
                           for s in F["sites"]["site_table"]])
        f0, e0, n0, _ = C.cover(cells, AL)
        f1, e1, n1, _ = C.cover(cells, SV) if len(SV) else (0.0, 0.0, np.zeros(len(cells)), None)
        meta["freeze_wa_v1_11site"] = dict(
            n_cells=int(len(cells)), n_sites=len(AL), n_sites_surviving=len(SV),
            sites_lost=[d["site"] for d in died],
            pct_cells_ge5_sites_60km_published=round(100 * f0, 1),
            pct_cells_ge5_sites_60km_health_adjusted=round(100 * f1, 1),
            median_sites_60km_published=float(np.median(n0)),
            median_sites_60km_health_adjusted=float(np.median(n1)) if len(SV) else 0.0)
        print("[freeze_wa_v1 11-site] ", json.dumps(meta["freeze_wa_v1_11site"], indent=1))
    json.dump(meta, open(out.replace(".csv", ".json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
