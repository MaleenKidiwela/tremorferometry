#!/usr/bin/env python
"""Apply the STAGE-1 CALIBRATED confirmation rule to a sealed-window re-fire.

The rule comes from `results/sf_final/calib.json`, which was fixed by measuring Bostock's own
rediscovered families on the north half-split BEFORE any Washington number was looked at
(`results/CONFIRMATION_CALIBRATION_NORTH.md`).  Nothing here is tuned per bed: the two (or three)
thresholds are read from that file and applied mechanically to the real chain and to the null
chain through identical code, so the null's confirmation count IS the false-confirmation rate.

  python scripts/sf_confirm.py --grps ... --det-tremor ... --det-quiet ... \
      --tremor-dates F --quiet-dates F --loc ... [--loc-sealed ...] --out CSV --tag real
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


def read_dates(p):
    if not p or not os.path.exists(p):
        return []
    return [x.strip() for x in open(p).read().replace("\n", ",").split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grps", required=True)
    ap.add_argument("--det-tremor", required=True)
    ap.add_argument("--det-quiet", default="")
    ap.add_argument("--tremor-dates", required=True)
    ap.add_argument("--quiet-dates", default="")
    ap.add_argument("--loc", required=True, help="discovery-window locations")
    ap.add_argument("--loc-sealed", default="", help="sealed-window locations (optional)")
    ap.add_argument("--calib", default=f"{SF.RES}/sf_final/calib.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="real")
    a = ap.parse_args()
    R = json.load(open(a.calib))
    min_days = R.get("min_days_h2")
    fa_max = R.get("fa_max_per_tpl_day")
    use_move = bool(R.get("use_move"))
    mv = R.get("move_clause") or {}

    G = pd.read_csv(a.grps)
    keep = sorted(int(x) for x in G.grp)
    def rd(p):
        return (pd.read_parquet(p) if p and os.path.exists(p)
                else pd.DataFrame(columns=["grp", "t", "cc", "nmad", "date"]))
    DT, DQ = rd(a.det_tremor), rd(a.det_quiet)
    td, qd = read_dates(a.tremor_dates), read_dates(a.quiet_dates)
    L = pd.read_csv(a.loc).set_index("grp")
    LS = (pd.read_csv(a.loc_sealed).set_index("grp")
          if a.loc_sealed and os.path.exists(a.loc_sealed) else None)
    nq, nt = max(len(qd), 1), max(len(td), 1)
    ndet = DT.groupby("grp").size()
    ndays = DT.groupby("grp").date.nunique()
    nq_g = DQ.groupby("grp").size()
    rows = []
    for g in keep:
        rec = dict(grp=g,
                   lat=float(L.lat.get(g, np.nan)) if g in L.index else np.nan,
                   lon=float(L.lon.get(g, np.nan)) if g in L.index else np.nan,
                   dep=float(L.dep.get(g, np.nan)) if g in L.index else np.nan,
                   n_det_disc=int(L.n_det.get(g, 0)) if g in L.index else 0,
                   n_days_disc=int(L.n_days.get(g, 0)) if g in L.index else 0,
                   n_det_sealed=int(ndet.get(g, 0)),
                   n_days_sealed=int(ndays.get(g, 0)),
                   n_det_quiet_sealed=int(nq_g.get(g, 0)),
                   fa_sealed=float(nq_g.get(g, 0)) / nq,
                   rate_sealed=float(ndet.get(g, 0)) / nt)
        rec["c_days"] = bool(min_days is None or rec["n_days_sealed"] >= min_days)
        rec["c_fa"] = bool(fa_max is None or rec["fa_sealed"] <= fa_max)
        rec["move_km"] = np.nan
        rec["c_move"] = True
        if LS is not None and g in LS.index and bool(LS.located.get(g, False)) \
                and np.isfinite(rec["lat"]):
            rec.update(lat_sealed=float(LS.lat[g]), lon_sealed=float(LS.lon[g]),
                       dep_sealed=float(LS.dep[g]))
            rec["move_km"] = float(SF.hyp_km(rec["lat"], rec["lon"], rec["dep"],
                                             rec["lat_sealed"], rec["lon_sealed"],
                                             rec["dep_sealed"]))
            rec["move_epi_km"] = float(SF.epi_km(rec["lat"], rec["lon"],
                                                 rec["lat_sealed"], rec["lon_sealed"]))
        if use_move and mv and LS is not None:
            col = "move_km" if mv.get("col") == "move_3d_km" else "move_epi_km"
            v = rec.get(col, np.nan)
            rec["c_move"] = bool(np.isfinite(v) and v <= float(mv["thr"]))
        rec["CONFIRMED"] = bool(rec["c_days"] and rec["c_fa"] and rec["c_move"])
        rows.append(rec)
    if use_move and LS is None:
        print("NOTE: the calibrated rule keeps a location-move clause but no sealed-window "
              "locations were supplied for this chain -- the clause is reported as not "
              "evaluated, not as a failure.", flush=True)
    C = pd.DataFrame(rows)
    C.to_csv(a.out, index=False)
    summ = dict(tag=a.tag, rule=dict(min_days=min_days, fa_max=fa_max, use_move=use_move,
                                     move=mv),
                n_refired=int(len(C)), n_tremor_days=len(td), n_quiet_days=len(qd),
                n_det_sealed=int(len(DT)), n_det_quiet_sealed=int(len(DQ)),
                rate_tremor_per_tpl_day=float(len(DT) / max(len(C) * nt, 1)),
                rate_quiet_per_tpl_day=float(len(DQ) / max(len(C) * nq, 1)),
                n_days_ok=int(C.c_days.sum()) if len(C) else 0,
                n_fa_ok=int(C.c_fa.sum()) if len(C) else 0,
                move_evaluated=bool(use_move and LS is not None),
                n_move_ok=int(C.c_move.sum()) if len(C) else 0,
                n_CONFIRMED=int(C.CONFIRMED.sum()) if len(C) else 0)
    json.dump(summ, open(a.out.replace(".csv", "_summary.json"), "w"), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
