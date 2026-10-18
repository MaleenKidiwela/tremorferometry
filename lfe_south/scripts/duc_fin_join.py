#!/usr/bin/env python
"""DUC_FIN join stage (CPU): check the feature re-run against duc_fix pass 1, assemble the final candidate list
2008-01..2012-07 and re-join with his catalogue.

  1. mfcheck: on every day where no re-templated station took part (n_new_tpl == 0 for all families), the feature
     run's data candidates must be IDENTICAL to results/duc_fix/mf1 (same (family, sample) keys, |cc diff| < 1e-9);
     exit 7 when > 3 % of those days differ or a day file is missing.
  2. OUT/catalogue_final.parquet: every data candidate (all >= 8 x MAD, nothing cut) with its detector statistics;
     src = 'pass1_identical' (no new template that day) or 'rescan_new_templates' (2011-10-27 .. 2012-07 with the
     DUC_FIN re-templated epochs = the pass-2 re-scan).
  3. joins: monthly / yearly ours vs his (his rule both sides, raw, 60 kept families; one-to-one +-2 s same family),
     templates per family-day by month (before = mf1 diag, after = this run), 2011-10..2012-07 rule counts before / after.
"""
import os, sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
from t_common import log


def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--days", default="", help="smoke: only these days"); a = ap.parse_args()
    DT = F.day_table()
    if a.days: DT = DT[DT.date.isin(a.days.split(","))]
    rows = []; parts = []; nt = []
    mf1 = {}
    for m in F.MONTHS:
        p = f"{F.FIX}/mf1/{m}.parquet"
        if os.path.exists(p):
            M = pd.read_parquet(p); M["fam"] = M.fam.astype(str); M["date"] = M.t.dt.strftime("%Y-%m-%d"); mf1[m] = M
    miss = []
    for ds in DT.date:
        pj = f"{F.OUT}/feat/{ds}.json"; pq = f"{F.OUT}/feat/{ds}.parquet"
        if not os.path.exists(pj): miss.append(ds); continue
        meta = json.load(open(pj))
        for f, (n, nn) in meta.get("ntpl", {}).items(): nt.append((ds, f, n, nn))
        A = pd.read_parquet(pq) if os.path.getsize(pq) > 0 else pd.DataFrame()
        a = A[A.series == "data"].copy() if len(A) else pd.DataFrame(columns=["fam", "t_epoch", "cc", "n_new_tpl"])
        newday = any(v[1] > 0 for v in meta.get("ntpl", {}).values())
        r = dict(date=ds, n_feat=len(a), new_tpl_day=newday)
        M = mf1.get(ds[:7]); m = M[M.date == ds] if M is not None else None
        if m is not None:
            r["n_mf1"] = len(m)
            if not newday:
                ka = dict(zip(zip(a.fam.astype(str), np.round(a.t_epoch.to_numpy() * 20).astype(np.int64)), a.cc))
                km = dict(zip(zip(m.fam, np.round(m.t_epoch.to_numpy() * 20).astype(np.int64)), m.cc))
                same = set(ka) == set(km)
                dcc = max([abs(ka[k] - km[k]) for k in set(ka) & set(km)] or [0.0])
                r.update(identical=bool(same and dcc < 1e-9), max_dcc=dcc, n_only_feat=len(set(ka) - set(km)), n_only_mf1=len(set(km) - set(ka)))
        if len(a):
            a = a.drop(columns=[c for c in a.columns if c.startswith(("cc_", "nc_", "z_"))] + ["series", "k"], errors="ignore")
            a["date"] = ds; a["src"] = "rescan_new_templates" if newday else "pass1_identical"
            parts.append(a)
        rows.append(r)
    C = pd.DataFrame(rows); C.to_csv(f"{F.OUT}/mfcheck_by_day.csv", index=False)
    cmp_ = C[C.identical.notna()] if "identical" in C else C.iloc[:0]
    nbad = int((cmp_.identical == False).sum()) if len(cmp_) else 0
    msg = (f"[mfcheck] days {len(DT)}, feature files missing {len(miss)}; days without new templates compared with mf1 "
           f"{len(cmp_)}: identical {int((cmp_.identical == True).sum())}, different {nbad}; "
           f"max |cc diff| {cmp_.max_dcc.max() if len(cmp_) else np.nan:.2e}")
    log(msg); os.makedirs(f"{F.OUT}/status", exist_ok=True); open(f"{F.OUT}/status/MFCHECK.txt", "w").write(msg + "\n")
    if len(miss) + nbad > F.ERR_MAX * len(DT):
        log(f"[mfcheck] FAIL: missing {miss[:5]} different {cmp_[cmp_.identical == False].date.tolist()[:5]}"); sys.exit(7)
    D = pd.concat(parts, ignore_index=True)
    D["fam"] = D.fam.astype(str); D["t"] = pd.to_datetime(D.t_epoch, unit="s").astype("datetime64[ns]")
    D = D.sort_values(["t", "fam"]).reset_index(drop=True)
    D.to_parquet(f"{F.OUT}/catalogue_final.parquet", index=False)
    NT = pd.DataFrame(nt, columns=["date", "fam", "n_tpl", "n_new_tpl"]); NT.to_csv(f"{F.OUT}/templates_per_famday.csv", index=False)
    # ---- joins
    H = F.his_catalogue(); H = H[(H.t >= "2008-01-01") & (H.t < "2012-01-01")].reset_index(drop=True)
    H["m"] = H.t.dt.strftime("%Y-%m"); D["m"] = D.t.dt.strftime("%Y-%m")
    O = D[(D.pass_his == 1) & (D.t < "2012-01-01")].reset_index(drop=True); Hh = H[H.hp].reset_index(drop=True)
    ia, ib, _ = F.match_one_to_one(O.t_epoch.to_numpy(), O.fam.to_numpy(), Hh.t_epoch.to_numpy(), Hh.fam.to_numpy())
    O["mt"] = False; O.loc[ia, "mt"] = True; Hh["mt"] = False; Hh.loc[ib, "mt"] = True
    jr = []
    for m in sorted(set(D.m)):
        d = D[D.m == m]; h = H[H.m == m]; o = O[O.m == m]; hh = Hh[Hh.m == m]
        before = mf1.get(m)
        jr.append(dict(month=m, ours_rule=int(d.pass_his.sum()), his_rule=int(h.hp.sum()) if m < "2012" else np.nan,
                       ratio_rule=d.pass_his.sum() / max(h.hp.sum(), 1) if m < "2012" else np.nan,
                       recall_same=hh.mt.mean() if len(hh) else np.nan, precision_same=o.mt.mean() if len(o) else np.nan,
                       ours_raw=len(d), his_raw=len(h) if m < "2012" else np.nan,
                       ours_rule_mf1=int(before.pass_his.sum()) if before is not None else np.nan,
                       ours_raw_mf1=len(before) if before is not None else np.nan,
                       ntpl_famday=NT[NT.date.str[:7] == m].n_tpl.mean(),
                       ntpl_famday_mf1=pd.read_csv(f"{F.FIX}/mf1/{m}_diag.csv").n_tpl.mean() if os.path.exists(f"{F.FIX}/mf1/{m}_diag.csv") else np.nan,
                       new_tpl_famday=NT[NT.date.str[:7] == m].n_new_tpl.mean()))
    J = pd.DataFrame(jr); J.to_csv(f"{F.OUT}/joins_by_month.csv", index=False)
    J["year"] = J.month.str[:4]
    Y = J.groupby("year").agg(ours_rule=("ours_rule", "sum"), his_rule=("his_rule", "sum"), ours_raw=("ours_raw", "sum"),
                              his_raw=("his_raw", "sum"), ours_rule_mf1=("ours_rule_mf1", "sum")).reset_index()
    Y["ratio_rule"] = Y.ours_rule / Y.his_rule.replace(0, np.nan)
    yy = O.assign(y=O.t.dt.year).groupby("y").mt.mean(); yh = Hh.assign(y=Hh.t.dt.year).groupby("y").mt.mean()
    Y["precision_same"] = Y.year.astype(int).map(yy); Y["recall_same"] = Y.year.astype(int).map(yh)
    Y.to_csv(f"{F.OUT}/joins_by_year.csv", index=False)
    log(f"[join] catalogue_final {len(D)} candidates ({int(D.pass_his.sum())} pass his rule; src {D.src.value_counts().to_dict()}); "
        f"years:\n{Y.round(3).to_string(index=False)}")


if __name__ == "__main__":
    main()
