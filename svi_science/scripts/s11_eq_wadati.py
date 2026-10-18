#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6.4): earthquake Wadati κ per region / 30-km zone / class, same estimators and gates
as the LFE κ (s11_tomo_pooled.zone with E1 tP-on-tS + clipping, OLS column; s11_tomo_profile.zone inlier profile).

  python s11_eq_wadati.py [--picks data/eq_picks/picks.parquet] [--out eq_kappa] [--workers 8] [--quick]
writes tomo/eq/tables/<out>_rows.csv (event-station P/S rows), <out>_zones.csv (one row per cell x estimator x variant)
"""
import argparse
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
EQ = f"{ROOT}/svi_science/results/s11_depth/tomo/eq"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s11_tomo_pooled as PO  # noqa: E402
import s11_tomo_profile as PR  # noqa: E402

MIN_ROWS, MIN_EV = 4, 5


def build_rows(picks):
    P = pd.read_parquet(picks)
    P = P[P.phase1.isin(["P", "S"]) & P.t_rel.notna() & P.event.notna()].copy()
    # one pick per (event, station, phase): manual first, then the earliest
    P["man"] = (P["mode"] != "manual").astype(int)
    P = P.sort_values(["man", "t_rel"]).drop_duplicates(["cc_id", "stakey", "phase1"])
    key = ["cc_id", "stakey"]
    a = P[P.phase1 == "P"].set_index(key)
    b = P[P.phase1 == "S"].set_index(key)
    R = a[["event", "reg", "zone", "kind", "tier", "mag", "o_depth_km", "cat_depth_km", "along", "o_lat", "o_lon",
           "dist_km", "t_rel", "unc", "time_weight", "time_residual", "cha"]].join(
        b[["t_rel", "unc", "time_weight", "time_residual", "cha", "dist_km"]], rsuffix="_S", how="inner")
    R = R.rename(columns=dict(t_rel="tP", t_rel_S="tS", unc="uncP", unc_S="uncS", time_weight="wP",
                              time_weight_S="wS", time_residual="resP", time_residual_S="resS"))
    R = R.reset_index()
    R["dist_km"] = R.dist_km.fillna(R.dist_km_S)
    R = R[(R.tS > R.tP) & (R.tP > 0)]
    R["fam"] = R.cc_id.astype(str)
    R["sta"] = R.stakey
    R["zone"] = R.zone.astype(int)
    return R


def cells(R, quick=False):
    """(label, variant, frame) list. Headline rows: dist <= 100 km, all weights."""
    out = []
    var = {"d100": R.dist_km <= 100, "d60": R.dist_km <= 60, "d150": R.dist_km <= 150,
           "d100_w": (R.dist_km <= 100) & (R.wP > 0) & (R.wS > 0)}
    for v, m in var.items():
        Q = R[m]
        for kind in ("crust", "slab", "layer"):
            for reg in ("SVI", "WA_NW", "WA_SW"):
                q = Q[(Q.kind == kind) & (Q.reg == reg)]
                out.append((f"{kind}|{reg}|all", v, q))
                if v != "d100" or quick:
                    continue
                for z in sorted(q.zone.unique()):
                    out.append((f"{kind}|{reg}|{z:+d}", v, q[q.zone == z]))
    return out


def enough(q):
    n = q.groupby("fam").size()
    return int((n >= MIN_ROWS).sum())


def run_cell(arg):
    lab, var, q, seed, nb, ns = arg
    kind, reg, zone = lab.split("|")
    base = dict(kind=kind, reg=reg, zone=zone, variant=var, n_ev_in=q.fam.nunique(), n_rows_in=len(q),
                n_ev_ge4=enough(q), dist_med=float(q.dist_km.median()) if len(q) else np.nan,
                depth_med=float(q.o_depth_km.median()) if len(q) else np.nan,
                along_med=float(q.along.median()) if len(q) else np.nan)
    if base["n_ev_ge4"] < MIN_EV:
        return [dict(base, est="none", note=f"< {MIN_EV} events with >= {MIN_ROWS} P+S stations")]
    rng = np.random.default_rng(seed)
    q = q.groupby("fam").filter(lambda g: len(g) >= MIN_ROWS)
    rows = []
    PO.NB, PO.NS, PO.NSTA_MIN = nb, ns, MIN_ROWS
    for est in ("inv", "ols"):
        PO.EST = est
        o, _ = PO.zone(q, rng, lab, True)
        rows.append({**base, **{k: v for k, v in o.items() if k not in base and k not in ("zone", "est")},
                     "est": f"clip_{est}"})
    sigP = float(np.nanmedian(q.uncP)) if q.uncP.notna().any() else 0.05
    sigS = float(np.nanmedian(q.uncS)) if q.uncS.notna().any() else 0.1
    PR.NB, PR.NSC, PR.MINROWS = nb, max(ns // 5, 100), MIN_ROWS
    for thr in (0.25, 0.5):
        PR.THR = thr
        o, _ = PR.zone(q, rng, max(sigP, 0.025), max(sigS, 0.025))
        rows.append({**base, **{k: v for k, v in o.items() if k not in base and k not in ("zone", "est")},
                     "est": f"profile_{thr}"})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--picks", default=f"{ROOT}/data/eq_picks/picks.parquet")
    ap.add_argument("--out", default="eq_kappa")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--nb", type=int, default=1000)
    ap.add_argument("--ns", type=int, default=1000)
    ap.add_argument("--quick", action="store_true", help="region-level cells only")
    a = ap.parse_args()
    os.makedirs(f"{EQ}/tables", exist_ok=True)
    R = build_rows(a.picks)
    R.to_csv(f"{EQ}/tables/{a.out}_rows.csv", index=False)
    print(f"rows {len(R)} events {R.fam.nunique()} stations {R.sta.nunique()}", flush=True)
    jobs = [(lab, v, q, 1000 + i, a.nb, a.ns) for i, (lab, v, q) in enumerate(cells(R, a.quick))]
    with Pool(min(a.workers, 8)) as pool:
        res = pool.map(run_cell, jobs, chunksize=1)
    Z = pd.DataFrame([r for rr in res for r in rr])
    Z.to_csv(f"{EQ}/tables/{a.out}_zones.csv", index=False)
    show = ["kind", "reg", "zone", "variant", "est", "n_fam", "n_rows", "kappa", "sig", "rms", "z_vs_scr", "loso_max",
            "d_site", "ctrl_pass", "pass", "note"]
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(Z[[c for c in show if c in Z]].round(3).to_string())


if __name__ == "__main__":
    main()
