#!/usr/bin/env python
"""ONE ROW PER GROWN TEMPLATE, nothing filtered.  Prefix `sf_`.

Binding rule (user, 2026-09-23, memory `rank-dont-discard`): calibrating acceptance criteria
against the 58 known Bostock families can only detect criteria that are too strict FOR THAT
POPULATION; it gives no protection to a real family that is smaller, deeper, noisier or seen by
fewer stations -- which is exactly what an undiscovered Washington family may be.  So:

  * every grown family keeps its row, permanently.  There is no "candidates after the screen"
    number and no template is deleted for failing a clause;
  * every criterion is a COLUMN carrying its measured value, never a gate;
  * the whole-chain site-shift NULL gets the byte-identical table -- the only honest comparator;
  * the two steps the published chains actually perform (the dedupe and the location) are recorded
    as FLAGS with counts before and after; merged and unlocatable rows stay in the table.

This module builds that table for one chain.  It reads only saved outputs, so it re-scores an
already-run chain without touching the discovery recipe.

  python scripts/sf_table.py --key oly2010 --chain real --det ... --out ...
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
WA_RAILS = (15.0, 59.0)
N_BASE = 200


def dates_of(spec):
    if not spec:
        return []
    if os.path.exists(spec):
        spec = open(spec).read()
    return [x.strip() for x in spec.replace("\n", ",").split(",") if x.strip()]


def rd(p):
    return (pd.read_parquet(p) if p and os.path.exists(p)
            else pd.DataFrame(columns=["grp", "t", "cc", "nmad", "date"]))


def co_timing(D, hours_csv, rng):
    """fraction of a template's detections inside an audited tremor hour, and the same-day
    uniform-redraw base rate -- n_screen.py's statistic, recomputed here so every chain (Arm B
    included) has the column."""
    if not hours_csv or not os.path.exists(hours_csv) or not len(D):
        return {}
    H = pd.read_csv(hours_csv)
    THR = {int(pd.Timestamp(f"{r.date} {int(r.hour):02d}:00:00").timestamp()) // 3600
           for r in H.itertuples()}
    out = {}
    for g, d in D.groupby("grp"):
        tv = d.t.to_numpy(np.float64)

        def frac(ts):
            ep = np.floor_divide(ts.astype(np.int64), 3600)
            return float(np.mean([int(x) in THR for x in ep]))

        obs = frac(tv)
        dsec = np.array([pd.Timestamp(x).timestamp() for x in d.date], dtype=np.float64)
        bs = np.array([frac(dsec + rng.random(len(dsec)) * 86400.0) for _ in range(N_BASE)])
        out[int(g)] = (obs, float(bs.mean()), float((bs >= obs).mean()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--key", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--det", required=True, help="discovery-window detections")
    ap.add_argument("--det-fa", default="", help="discovery-era quiet-day detections")
    ap.add_argument("--loc", required=True)
    ap.add_argument("--dedupe", default="")
    ap.add_argument("--screen", default="", help="an existing screen_final.csv, merged if present")
    ap.add_argument("--hours", default="", help="audited tremor hours CSV (date,hour)")
    ap.add_argument("--tremor-dates", default="")
    ap.add_argument("--quiet-dates", default="")
    ap.add_argument("--sealed-det", default="")
    ap.add_argument("--sealed-fa", default="")
    ap.add_argument("--sealed-loc", default="")
    ap.add_argument("--sealed-tremor-dates", default="")
    ap.add_argument("--sealed-quiet-dates", default="")
    ap.add_argument("--chestler", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    D = rd(a.det)
    Q = rd(a.det_fa)
    L = pd.read_csv(a.loc)
    DD = pd.read_csv(a.dedupe).set_index("grp") if a.dedupe and os.path.exists(a.dedupe) else None
    SC = (pd.read_csv(a.screen).set_index("grp")
          if a.screen and os.path.exists(a.screen) else None)
    SD = rd(a.sealed_det)
    SQ = rd(a.sealed_fa)
    SL = (pd.read_csv(a.sealed_loc).set_index("grp")
          if a.sealed_loc and os.path.exists(a.sealed_loc) else None)
    nt = max(len(dates_of(a.tremor_dates)), 1)
    nq = max(len(dates_of(a.quiet_dates)), 1)
    snt = max(len(dates_of(a.sealed_tremor_dates)), 1)
    snq = max(len(dates_of(a.sealed_quiet_dates)), 1)

    nd = D.groupby("grp").size()
    ndy = D.groupby("grp").date.nunique()
    mcc = D.groupby("grp").cc.median()
    mnm = D.groupby("grp").nmad.max()
    nq_g = Q.groupby("grp").size()
    snd = SD.groupby("grp").size()
    sndy = SD.groupby("grp").date.nunique()
    snq_g = SQ.groupby("grp").size()
    rng = np.random.default_rng(17)
    CO = co_timing(D, a.hours, rng)

    C = pd.read_csv(CHESTLER) if a.chestler else None
    rows = []
    for r in L.itertuples():
        g = int(r.grp)
        rec = dict(key=a.key, chain=a.chain, grp=g,
                   n_det=int(nd.get(g, 0)), n_days=int(ndy.get(g, 0)),
                   med_cc=float(mcc.get(g, np.nan)), max_nmad=float(mnm.get(g, np.nan)),
                   rate_tremor=float(nd.get(g, 0)) / nt,
                   n_det_quiet=int(nq_g.get(g, 0)),
                   fa_quiet=float(nq_g.get(g, 0)) / nq)
        rec["tremor_quiet_ratio"] = (rec["rate_tremor"] / rec["fa_quiet"]
                                     if rec["fa_quiet"] > 0 else np.nan)
        co = CO.get(g)
        rec.update(co_timing=co[0] if co else np.nan, co_base=co[1] if co else np.nan,
                   co_p=co[2] if co else np.nan,
                   co_excess=(co[0] - co[1]) if co else np.nan)
        # ---- the dedupe: a FLAG plus what it merged on; merged rows are kept
        if DD is not None and g in DD.index:
            rec.update(dedupe_kept=bool(DD.kept[g]), dedupe_merged_into=int(DD.merged_into[g]),
                       dedupe_max_share=float(DD.max_share[g]),
                       dedupe_max_cc=float(DD.max_cc[g]),
                       dedupe_n_partners=int(DD.n_partners[g]))
        else:
            rec.update(dedupe_kept=True, dedupe_merged_into=-1, dedupe_max_share=np.nan,
                       dedupe_max_cc=np.nan, dedupe_n_partners=0)
        # ---- location: a FLAG plus every statistic; unlocatable rows are kept
        rec["located"] = bool(r.located)
        rec["n_sta"] = int(r.n_sta)
        if bool(r.located):
            rec.update(lat=float(r.lat), lon=float(r.lon), dep=float(r.dep), rms=float(r.rms),
                       rms_shuf_p05=float(r.rms_shuf_p05), rms_shuf_med=float(r.rms_shuf_med),
                       ctrl_pass=bool(r.ctrl_pass),
                       ctrl_margin=float(r.rms_shuf_p05) / max(float(r.rms), 1e-9),
                       edge=bool(r.edge), dep_rail=bool(float(r.dep) in WA_RAILS),
                       slab2=float(r.slab2), dz=float(r.dz), abs_dz=abs(float(r.dz)),
                       depth_window_15_58=bool(15.0 <= float(r.dep) <= 58.0))
        else:
            for k in ("lat", "lon", "dep", "rms", "rms_shuf_p05", "rms_shuf_med", "ctrl_margin",
                      "slab2", "dz", "abs_dz"):
                rec[k] = np.nan
            rec.update(ctrl_pass=False, edge=False, dep_rail=False, depth_window_15_58=False)
        # ---- distance to the nearest published Olympic coordinate
        if C is not None and rec["located"]:
            e = SF.epi_km(rec["lat"], rec["lon"], C.latitude.to_numpy(float),
                          C.longitude.to_numpy(float))
            d3 = np.sqrt(e ** 2 + (rec["dep"] - C.depth.to_numpy(float)) ** 2)
            rec.update(near_chestler_km=float(d3.min()), near_chestler_epi_km=float(e.min()),
                       near_chestler_row=int(C.iloc[int(np.argmin(d3))].iloc[0]))
        else:
            rec.update(near_chestler_km=np.nan, near_chestler_epi_km=np.nan,
                       near_chestler_row=-1)
        # ---- the split-half re-stack move, if the earlier screen computed it
        if SC is not None and g in SC.index:
            for src, dst in (("h1_h2_km", "restack_move_km"), ("h2_km", "restack_h2_km")):
                if src in SC.columns:
                    v = SC[src].get(g, np.nan)
                    rec[dst] = float(v) if pd.notna(v) else np.nan
            for c in ("c_a", "c_b", "c_c1", "c_c2", "c_c3", "c_c4", "c_d", "c_e"):
                if c in SC.columns:
                    rec[f"old_{c}"] = bool(SC[c].get(g, False))
        # ---- the sealed-window re-fire
        rec.update(sealed_n_det=int(snd.get(g, 0)), sealed_n_days=int(sndy.get(g, 0)),
                   sealed_rate=float(snd.get(g, 0)) / snt,
                   sealed_n_det_quiet=int(snq_g.get(g, 0)),
                   sealed_fa=float(snq_g.get(g, 0)) / snq)
        rec["sealed_ratio"] = (rec["sealed_rate"] / rec["sealed_fa"]
                               if rec["sealed_fa"] > 0 else np.nan)
        if SL is not None and g in SL.index and bool(SL.located.get(g, False)) \
                and rec["located"]:
            rec.update(sealed_lat=float(SL.lat[g]), sealed_lon=float(SL.lon[g]),
                       sealed_dep=float(SL.dep[g]))
            rec["sealed_move_km"] = float(SF.hyp_km(rec["lat"], rec["lon"], rec["dep"],
                                                    rec["sealed_lat"], rec["sealed_lon"],
                                                    rec["sealed_dep"]))
        else:
            rec.update(sealed_lat=np.nan, sealed_lon=np.nan, sealed_dep=np.nan,
                       sealed_move_km=np.nan)
        rows.append(rec)
    T = pd.DataFrame(rows)
    # tier A/B/C -- SUPERSEDED comparison columns only, never a filter
    T["tierA"] = (T.located & (T.abs_dz <= 19) & (T.n_days >= 3)
                  & (T.fa_quiet <= 0.25) & (T.rate_tremor >= 3 * T.fa_quiet) & (T.n_det > 0))
    T["tierB"] = T.tierA & T.ctrl_pass
    T["tierC"] = T.tierB & (T.co_timing > T.co_base)
    T.to_csv(a.out, index=False)
    summ = dict(key=a.key, chain=a.chain, n_rows=int(len(T)),
                n_dedupe_kept=int(T.dedupe_kept.sum()),
                n_located=int(T.located.sum()),
                n_dedupe_and_located=int((T.dedupe_kept & T.located).sum()),
                n_on_grid_edge=int(T.edge.sum()), n_on_depth_rail=int(T.dep_rail.sum()),
                n_depth_window=int(T.depth_window_15_58.sum()),
                n_tierA=int(T.tierA.sum()), n_tierB=int(T.tierB.sum()),
                n_tierC=int(T.tierC.sum()),
                n_sealed_refired=int((T.sealed_n_det > 0).sum()),
                n_tremor_days=nt, n_quiet_days=nq,
                n_sealed_tremor_days=snt, n_sealed_quiet_days=snq)
    json.dump(summ, open(a.out.replace(".csv", "_summary.json"), "w"), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
