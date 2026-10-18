#!/usr/bin/env python
"""P4-2b step 5: completeness factors k from the ladder (old station set run 60 d after each boundary) and from the
era/frozen overlap.  WAITS FOR P6: needs a confidence P per candidate.  Rules: results/backward/p4b/P4B.md.

  excess rate  = (sum P_data - sum P_control) / observable data family-days, over a day set
  k(new/old)   = excess rate of the new detector / excess rate of the old one on the SAME days (paired day bootstrap,
                 2000 draws, 68 % and 95 % CI)
  day sets     = env_quiet (quiet_env3_CN4), ets (inside a Bostock window, dates only), between (> 3 d from every window),
                 all;  rows: all candidates, and multi_site_z2 only (n_z2 >= 2; B2/B3 bracket [all, n_z2 >= 2])
  agreement    = k estimates of the same pair within x2 of each other, else bracket [min, max]

P input (from P6):  --p-dir DIR with DIR/<run>_P.parquet, columns run, date, row, P  (row = 0-based row of the feature
file <date>.parquet, the same key as the p4b/sidecar tables).  --smoke: P = 1 for every candidate (mechanics test only,
writes k_smoke.csv, NOT a result).  Uses the sidecar tables (bw_p4b_flags.py) for the row keys and family-day counts.

  python bw_p4b_kladder.py --p-dir <dir> | --smoke   -> p4b/tables/k_ladder.csv (or k_smoke.csv)
"""
import os, sys, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
SC = f"{BW}/p4b/sidecar"; OUT = f"{BW}/p4b/tables"
LADDER = [("CN4", "B1", "lad_CN4_20050701", "2005-07-01"), ("B1", "B2", "lad_B1_20051102", "2005-11-02"),
          ("B2", "B3", "lad_B2_20071113", "2007-11-13"), ("B3", "S15", "lad_B3_20090827", "2009-08-27")]
OVERLAP = [("CN4", "2003-02-01", "2005-06-30"), ("B1", "2005-07-01", "2005-11-01"), ("B2", "2005-11-02", "2006-06-30")]
ERA_RUN = {"2003": "p5e_2003_kelb", "2004": "p5e_2004", "2005": "p5e_2005", "2006": "p5e_2006"}


def day_classes():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import bw_p4_common as P4
    F = pd.read_csv(f"{BW}/p4/envflag/tables/quiet_flag_env.csv", dtype={"date": str})
    F["ets"] = P4.bostock_near(F.date, pad=0); F["near3"] = P4.bostock_near(F.date, pad=3)
    return dict(env_quiet=set(F[F.quiet_env3_CN4 == True].date), ets=set(F[F.ets].date),       # noqa: E712
                between=set(F[~F.near3].date), all=set(F.date))


def per_day(runs, days, pdir, smoke):
    """-> DataFrame date, spd_all, spc_all, spd_ms, spc_ms, m_data (summed over the runs that hold each day)."""
    out = []
    for run in runs:
        fp = f"{SC}/{run}_flags.parquet"
        if not os.path.exists(fp):
            raise SystemExit(f"[k] missing sidecar {fp} (run bw_p4b_flags.py after the run is finished)")
        A = pd.read_parquet(fp, columns=["date", "row", "series", "multi_site_z2"])
        A = A[A.date.isin(days)]
        if smoke:
            A["P"] = 1.0
        else:
            Pp = pd.read_parquet(f"{pdir}/{run}_P.parquet", columns=["date", "row", "P"])
            A = A.merge(Pp, on=["date", "row"], how="left")
            if A.P.isna().any():
                raise SystemExit(f"[k] {run}: {int(A.P.isna().sum())} candidates without P")
        Dd = pd.read_csv(f"{SC}/{run}_days.csv", dtype={"date": str})
        Dd = Dd[Dd.date.isin(days)][["date", "m_data"]]
        g = lambda m: A[m].groupby("date").P.sum()
        d = pd.DataFrame(dict(spd_all=g(A.series == "data"), spc_all=g(A.series == "c1"),
                              spd_ms=g((A.series == "data") & A.multi_site_z2), spc_ms=g((A.series == "c1") & A.multi_site_z2)))
        out.append(Dd.set_index("date").join(d).fillna(0.0))
    return pd.concat(out).groupby(level=0).sum() if out else pd.DataFrame()


def k_boot(N, O, rows, nb=2000, seed=0):
    days = sorted(set(N.index) & set(O.index))
    days = [d for d in days if N.loc[d, "m_data"] > 0 and O.loc[d, "m_data"] > 0]
    if len(days) < 5:
        return dict(n_days=len(days))
    n = N.loc[days]; o = O.loc[days]
    en = (n[f"spd_{rows}"] - n[f"spc_{rows}"]).to_numpy(); mn = n.m_data.to_numpy(float)
    eo = (o[f"spd_{rows}"] - o[f"spc_{rows}"]).to_numpy(); mo = o.m_data.to_numpy(float)
    rn, ro = en.sum() / mn.sum(), eo.sum() / mo.sum()
    rng = np.random.default_rng(seed); ks = []
    for _ in range(nb):
        i = rng.integers(0, len(days), len(days))
        a, b = en[i].sum() / mn[i].sum(), eo[i].sum() / mo[i].sum()
        ks.append(a / b if b > 0 else np.nan)
    ks = np.array(ks); f = np.isfinite(ks)
    q = np.nanpercentile(ks[f], [2.5, 16, 50, 84, 97.5]) if f.sum() > 10 else [np.nan] * 5
    return dict(n_days=len(days), excess_new=rn, excess_old=ro, k=rn / ro if ro > 0 else np.nan,
                k_lo95=q[0], k_lo68=q[1], k_med=q[2], k_hi68=q[3], k_hi95=q[4], frac_boot_finite=float(f.mean()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--p-dir", default=""); ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    if not a.smoke and not a.p_dir:
        raise SystemExit("[k] needs --p-dir (P6 output) or --smoke")
    DC = day_classes(); rows = []
    for old, new, lad, b in LADDER:
        days = [str(x.date()) for x in pd.date_range(b, periods=60)]
        if not os.path.exists(f"{SC}/{lad}_flags.parquet"):
            print(f"[k] {lad}: ladder sidecar not there yet, skipped"); continue
        newruns = sorted({f"p5f_{d[:4]}" for d in days})
        if not all(os.path.exists(f"{SC}/{r}_flags.parquet") for r in newruns):
            print(f"[k] {lad}: new-set sidecars {newruns} not all there yet, skipped"); continue
        O = per_day([lad], days, a.p_dir, a.smoke); N = per_day(newruns, days, a.p_dir, a.smoke)
        for cl, S in DC.items():
            for rw in ("all", "ms"):
                r = k_boot(N.loc[[d for d in N.index if d in S]], O.loc[[d for d in O.index if d in S]], rw)
                rows.append(dict(kind="ladder", pair=f"{new}/{old}", boundary=b, days=cl, rows=rw, **r))
    for fset, d0, d1 in OVERLAP:
        days = [str(x.date()) for x in pd.date_range(d0, d1)]
        yrs = sorted({d[:4] for d in days})
        er = [ERA_RUN[y] for y in yrs]; fr = [f"p5f_{y}" for y in yrs]
        if not all(os.path.exists(f"{SC}/{r}_flags.parquet") for r in er + fr):
            continue
        E = per_day(er, days, a.p_dir, a.smoke); Fz = per_day(fr, days, a.p_dir, a.smoke)
        for cl, S in DC.items():
            for rw in ("all", "ms"):
                r = k_boot(E.loc[[d for d in E.index if d in S]], Fz.loc[[d for d in Fz.index if d in S]], rw)
                rows.append(dict(kind="overlap", pair=f"era/{fset}", boundary=f"{d0}..{d1}", days=cl, rows=rw, **r))
    R = pd.DataFrame(rows)
    # agreement within x2 between the day classes of the same pair / rows (env_quiet is the primary estimate)
    ag = []
    for (kind, pair, rw), g in R.groupby(["kind", "pair", "rows"]):
        ks = g.set_index("days").k.dropna(); ks = ks[ks > 0]
        if len(ks) >= 2:
            ag.append(dict(kind=kind, pair=pair, rows=rw, k_min=ks.min(), k_max=ks.max(),
                           agree_x2=bool(ks.max() / ks.min() <= 2.0), bracket=[ks.min(), ks.max()]))
    fn = "k_smoke" if a.smoke else "k_ladder"
    R.to_csv(f"{OUT}/{fn}.csv", index=False); pd.DataFrame(ag).to_csv(f"{OUT}/{fn}_agreement.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.round(3).to_string(index=False)); print(pd.DataFrame(ag).round(3).to_string(index=False))
    if a.smoke:
        print("[k] SMOKE (P = 1): mechanics only, not a result")


if __name__ == "__main__":
    main()
