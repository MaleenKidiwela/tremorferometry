#!/usr/bin/env python
"""STAGE 4 -- CONFIRMATION of the 2008 survivors on the 2007-01 window.  Prefix `n_`.

The rule is fixed here, before the 2007 numbers exist, and is the 2011 re-fire's rule with the
bed's own quiet days substituted (results/WA_2011_REFIRE_RESULTS.md sec.5):

    >= 3 distinct 2007 tremor days with detections
    AND quiet-day false alarms <= 0.1 per template-day on the 2007 literature-quiet days
        (the 15 out-of-episode days of tranche T5: 2006-11-xx and 2007-04-xx)
    AND a 2007 location that is INTERIOR to the wide search box and within 10 km (3-D) of the
        2008 location.

The site-shift-null templates that survived the 2008 screen go through byte-identical code on the
identical 2007 data: their confirmation count IS the false-confirmation rate.

  python scripts/n_confirm.py --screen ... --det-tremor ... --det-fa ... --loc2007 ... --out ...
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
import n_screen as NS                                                       # noqa: E402

MIN_DAYS = 3
FA_MAX = 0.1
MOVE_KM = 10.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", required=True, help="the 2008 screen table (2008 locations)")
    ap.add_argument("--grps", default="", help="which grps were re-fired (default: candidates)")
    ap.add_argument("--det-tremor", required=True)
    ap.add_argument("--det-fa", required=True)
    ap.add_argument("--loc2007", required=True)
    ap.add_argument("--tremor-dates", required=True)
    ap.add_argument("--quiet-dates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="real")
    a = ap.parse_args()

    V = pd.read_csv(a.screen)
    if a.grps and os.path.exists(a.grps):
        G = pd.read_csv(a.grps)
        keep = set(int(x) for x in G.grp)
    else:
        col = "candidate" if "candidate" in V.columns else "pass_abcd"
        keep = set(int(x) for x in V[V[col].astype(bool)].grp)
    V = V[V.grp.isin(keep)].set_index("grp")

    def rd(p):
        return (pd.read_parquet(p) if os.path.exists(p)
                else pd.DataFrame(columns=["grp", "t", "cc", "nmad", "date"]))

    DT, DF = rd(a.det_tremor), rd(a.det_fa)
    td = [x for x in open(a.tremor_dates).read().split(",") if x]
    qd = [x for x in open(a.quiet_dates).read().split(",") if x]
    L = pd.read_csv(a.loc2007) if os.path.exists(a.loc2007) else pd.DataFrame(columns=["grp"])
    L = L.set_index("grp") if len(L) else L

    nq = max(len(qd), 1)
    rows = []
    for g in sorted(keep):
        d = DT[DT.grp == g]
        f = DF[DF.grp == g]
        rec = dict(grp=int(g),
                   lat_2008=float(V.lat.get(g, np.nan)), lon_2008=float(V.lon.get(g, np.nan)),
                   dep_2008=float(V.dep.get(g, np.nan)),
                   n_det_2007=len(d), n_days_2007=int(d.date.nunique()) if len(d) else 0,
                   n_det_quiet_2007=len(f), fa_2007_per_tpl_day=len(f) / nq)
        rec["c_days"] = bool(rec["n_days_2007"] >= MIN_DAYS)
        rec["c_fa"] = bool(rec["fa_2007_per_tpl_day"] <= FA_MAX)
        loc = L.loc[g] if (len(L) and g in L.index) else None
        if loc is not None and bool(loc.get("located", False)):
            rec.update(lat_2007=float(loc.lat), lon_2007=float(loc.lon), dep_2007=float(loc.dep),
                       edge_2007=bool(loc.edge) if loc.edge == loc.edge else None,
                       ctrl_pass_2007=bool(loc.ctrl_pass))
            rec["move_km"] = (NS.hyp_km(rec["lat_2008"], rec["lon_2008"], rec["dep_2008"],
                                        rec["lat_2007"], rec["lon_2007"], rec["dep_2007"])
                              if np.isfinite(rec["lat_2008"]) else np.nan)
            rec["c_loc"] = bool((rec["edge_2007"] is False)
                                and np.isfinite(rec["move_km"]) and rec["move_km"] <= MOVE_KM)
        else:
            rec.update(lat_2007=np.nan, lon_2007=np.nan, dep_2007=np.nan, edge_2007=None,
                       ctrl_pass_2007=False, move_km=np.nan, c_loc=False)
        rec["CONFIRMED"] = bool(rec["c_days"] and rec["c_fa"] and rec["c_loc"])
        rows.append(rec)
    C = pd.DataFrame(rows)
    C.to_csv(a.out, index=False)
    rep = dict(tag=a.tag, n_refired=int(len(C)), n_tremor_days_2007=len(td),
               n_quiet_days_2007=len(qd),
               n_det_2007=int(len(DT)), n_det_quiet_2007=int(len(DF)),
               rate_tremor_per_tpl_day=float(len(DT) / max(len(C) * max(len(td), 1), 1)),
               rate_quiet_per_tpl_day=float(len(DF) / max(len(C) * nq, 1)),
               n_ge3_days=int(C.c_days.sum()), n_fa_ok=int(C.c_fa.sum()),
               n_loc_ok=int(C.c_loc.sum()), n_CONFIRMED=int(C.CONFIRMED.sum()))
    json.dump(rep, open(a.out.replace(".csv", "_summary.json"), "w"), indent=1, default=float)
    print(json.dumps(rep, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
