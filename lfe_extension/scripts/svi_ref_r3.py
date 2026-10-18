#!/usr/bin/env python
"""SVI_REF R3 -- A/B test of the southern cleaning (svi_ref_prep.py: flat / dead masks, burst mute, line notch) against the
frozen v1.1 data path, on the top-flagged production months and on the truth years.  Decides whether the production
feature pass applies the cleaning (that pass IS the 'R4 re-run' when the gate passes).

Runs compared (svi_ref_mf.py, S15 / all stations, the same scrambled control): <run>_nc (cleaning flags only) vs <run>_cl
(cleaning applied).  Truth days scored = 2011 inside Bostock's window + 2012 TEST blocks (build blocks are in-sample for
the templates, block edges excluded); 2010 (9 stations) reported separately; 2013 is NOT used (sealed).

Gate (Merlin): adopt the cleaning only if
  A  truth recall (one-to-one, own family, +-2 s, all >= T candidates) at a MATCHED scrambled-control rate is >= v1.1's,
     at the two v1.1 operating rates (control rate of the nc run at 8.25 and at 8.75 MAD per observable family-day), AND
  V  the known-family veto holds: >= 85 % of v1.1's matched events and >= 85 % of v1.1's rediscovered families (>= 5
     matches) are kept at the matched rate (a veto floor, not a target), AND
  B  on the flagged months the scrambled-control rate falls MORE than the data rate (ratio cl / nc at the month's v1.1
     T_ext, pooled over the months).
Outputs OUT/r3/: r3_recall.csv, r3_families.csv, r3_months.csv, r3_decision.json.

  python scripts/svi_ref_r3.py months          -> OUT/r3/flagged_months.json (top-flagged months from qc_summary.csv)
  python scripts/svi_ref_r3.py evaluate
"""
from __future__ import annotations
import os, sys, json, glob, argparse, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C

O3 = f"{R.OUT}/r3"
NMONTHS = 6


def flagged_months(n=NMONTHS):
    q = pd.read_csv(f"{C.TRC}/qc_summary.csv")
    s = q[q.kind == "site"].groupby("ym").flag_2lambda.apply(lambda v: int(pd.Series(v).fillna(False).astype(bool).sum()))
    h = q[q.kind == "network"].set_index("ym").flag_hfr_shift.fillna(False).astype(bool).astype(int)
    T = pd.DataFrame(dict(n_2lambda=s, hfr_shift=h)).fillna(0)
    T["score"] = T.n_2lambda + T.hfr_shift
    T = T[T.index >= "2014-01"].sort_values(["score", "n_2lambda"], ascending=False)
    return T.head(n)


def load_run(run):
    fs = sorted(p for p in glob.glob(f"{R.OUT}/feat/{run}/*.parquet") if not p.endswith("_fd.parquet"))
    import pyarrow.parquet as _pq   # 2026-10-02 fix: no-data days are written as empty, column-less parquet files -> skip them
    fs = [p for p in fs if os.path.getsize(p) > 0 and _pq.read_metadata(p).num_rows > 0]
    A = pd.concat([pd.read_parquet(p, columns=["series", "fam", "date", "t", "net_mad"]) for p in fs],
                  ignore_index=True)
    fd = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{R.OUT}/feat/{run}/*_fd.parquet")) if os.path.getsize(p) > 0],
                   ignore_index=True)
    return A, fd


def scored_rows(A, year):
    DT = R.load_day_table()
    roles = DT[DT.kind == f"truth{year}"].set_index("date").role.to_dict()
    A = A[A.date.map(roles).isin(["fit", "s9"])].copy()
    if year == 2012:
        A = A[R.block_role(A.t.to_numpy()) == "fit"]
    return A


def ctl_rate(A, fd, T, days):
    n = int(((A.series == "c1") & (A.net_mad >= T) & A.date.isin(days)).sum())
    m = int(((fd.series == "c1") & fd.observable.astype(bool) & fd.date.isin(days)).sum())
    return n / max(m, 1)


def matched_T(A, fd, target, days):
    grid = np.round(np.arange(8.0, 14.0, 0.01), 2)
    c = A[(A.series == "c1") & A.date.isin(days)].net_mad.to_numpy()
    m = int(((fd.series == "c1") & fd.observable.astype(bool) & fd.date.isin(days)).sum())
    for T in grid:
        if (c >= T).sum() / max(m, 1) <= target + 1e-12:
            return float(T)
    return float(grid[-1])


def recall(A, year, T, B):
    D = A[(A.series == "data") & (A.net_mad >= T)]
    hit = R.one2one_recall(D.t.to_numpy(), D.fam.to_numpy(), B.t_origin.to_numpy(), B.fam.to_numpy())
    return hit


def evaluate(a):
    R.ensure(O3)
    rows, frows = [], []
    dec = dict()
    for year in (2011, 2012, 2010):
        rn, rc = f"r3t{year}_nc", f"r3t{year}_cl"
        if not os.path.exists(f"{R.OUT}/feat/{rn}") or not os.path.exists(f"{R.OUT}/feat/{rc}"):
            continue
        An, fdn = load_run(rn); Ac, fdc = load_run(rc)
        An = scored_rows(An, year); Ac = scored_rows(Ac, year)
        days = sorted(set(An.date))
        B = R.bostock([year]); B = B[B.date.isin(days)]
        if year == 2012:
            B = B[R.block_role(B.t_origin.to_numpy()) == "fit"]
        for Tref in (8.25, 8.75):
            target = ctl_rate(An, fdn, Tref, days)
            Tc = matched_T(Ac, fdc, target, days)
            hn = recall(An, year, Tref, B); hc = recall(Ac, year, Tc, B)
            r = dict(year=year, T_v11=Tref, ctl_rate_target=target, T_clean=Tc, ctl_rate_clean=ctl_rate(Ac, fdc, Tc, days),
                     n_events=len(B), recall_v11=hn.mean() if len(B) else np.nan, recall_clean=hc.mean() if len(B) else np.nan,
                     matched_v11=int(hn.sum()), matched_clean=int(hc.sum()),
                     kept_of_v11=float((hn & hc).sum() / max(hn.sum(), 1)),
                     n_data_v11=int(((An.series == "data") & (An.net_mad >= Tref)).sum()),
                     n_data_clean=int(((Ac.series == "data") & (Ac.net_mad >= Tc)).sum()))
            fam = B.fam.to_numpy()
            fn = pd.Series(hn).groupby(fam).sum(); fc = pd.Series(hc).groupby(fam).sum()
            red = fn[fn >= 5].index
            r["n_fam_rediscovered_v11"] = len(red)
            r["frac_fam_kept"] = float((fc.reindex(red).fillna(0) >= 5).mean()) if len(red) else np.nan
            rows.append(r)
            for f_ in fn.index:
                frows.append(dict(year=year, T_v11=Tref, fam=int(f_), n_events=int((fam == f_).sum()), matched_v11=int(fn[f_]),
                                  matched_clean=int(fc.get(f_, 0))))
    Rc = pd.DataFrame(rows); Rc.to_csv(f"{O3}/r3_recall.csv", index=False)
    pd.DataFrame(frows).to_csv(f"{O3}/r3_families.csv", index=False)
    # flagged months
    mrows = []
    TB = C.tiers_by_epoch()
    if os.path.exists(f"{R.OUT}/feat/r3m_nc") and os.path.exists(f"{R.OUT}/feat/r3m_cl"):
        An, fdn = load_run("r3m_nc"); Ac, fdc = load_run("r3m_cl")
        for m in sorted(set(An.date.str[:7])):
            T = TB.get(C.epoch_of_month(m), (8.75, 8.25))[1]
            r = dict(month=m, T_ext=T)
            for ser in ("data", "c1"):
                nn = int(((An.series == ser) & (An.net_mad >= T) & An.date.str.startswith(m)).sum())
                nc = int(((Ac.series == ser) & (Ac.net_mad >= T) & Ac.date.str.startswith(m)).sum())
                r[f"{ser}_v11"] = nn; r[f"{ser}_clean"] = nc
            mrows.append(r)
    Mo = pd.DataFrame(mrows); Mo.to_csv(f"{O3}/r3_months.csv", index=False)
    fit = Rc[Rc.year.isin([2011, 2012])]
    if len(fit):
        A_ok = bool(all(fit.groupby("T_v11").apply(lambda g: g.matched_clean.sum() / g.n_events.sum() >= g.matched_v11.sum() / g.n_events.sum())))
        V_ok = bool(all(fit.groupby("T_v11").apply(lambda g: (g.kept_of_v11 * g.matched_v11).sum() / max(g.matched_v11.sum(), 1) >= 0.85)) and
                    (fit.frac_fam_kept.fillna(1) >= 0.85).all())
    else:
        A_ok = V_ok = False
    if len(Mo):
        dr = Mo.data_clean.sum() / max(Mo.data_v11.sum(), 1); cr = Mo.c1_clean.sum() / max(Mo.c1_v11.sum(), 1)
        B_ok = bool(cr < dr)
    else:
        dr = cr = np.nan; B_ok = False
    dec = dict(adopt_cleaning=bool(A_ok and V_ok and B_ok), A_recall_at_matched_control=A_ok, V_known_family_veto=V_ok,
               B_control_falls_more=B_ok, months_data_ratio=dr, months_control_ratio=cr,
               recall_table=Rc.round(4).to_dict("records"), months=Mo.to_dict("records"))
    json.dump(dec, open(f"{O3}/r3_decision.json", "w"), indent=1, default=float)
    R.log(f"[r3] adopt_cleaning={dec['adopt_cleaning']} (A {A_ok}, V {V_ok}, B {B_ok}; months data x{dr:.3f} control x{cr:.3f})\n"
          + Rc.round(4).to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["months", "evaluate"])
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    R.ensure(O3)
    if a.what == "months":
        T = flagged_months(); T.to_csv(f"{O3}/flagged_months.csv")
        json.dump(list(T.index), open(f"{O3}/flagged_months.json", "w"))
        R.log(f"[r3] flagged months:\n{T}")
    else:
        evaluate(a)
