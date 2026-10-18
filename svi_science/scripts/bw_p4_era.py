#!/usr/bin/env python
"""P4-3 / P4-4: thresholds of the era detector (2003-02-01..2006-06-30) by the forward rule (quiet bar) on its own
quiet days, and the cross-check of the frozen network on the same days against its P4-2 (E0) value.
Quiet days: option 3 (CN4 energy + Bostock +-3 d) and, beside it, Bostock-only flags (> 3 d, > 30 d from every
window), from p4/tables/quiet_flag_cn4.csv.  Also a diagnostic on B4 2010-13: the same flags against PNSN.
Rules: p4/P4.md.  Works on whatever runs exist (n_days reported); the queue calls it when all are in.

  python bw_p4_era.py   -> p4/tables/p4_3_era.csv, p4_4_crosscheck.csv, p4_flagdiag_B4.csv, p4_3_decision.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import bw_p4_common as P

F = pd.read_csv(f"{P.TAB}/quiet_flag_cn4.csv", dtype={"date": str})
FLAGS = ["quiet_opt3", "quiet_bostock3", "quiet_bostock30"]
QS = {f: set(F[F[f] == True].date) for f in FLAGS}                       # noqa: E712
ok_val = False
if os.path.exists(f"{P.TAB}/qflag_validation.json"):
    ok_val = bool(json.load(open(f"{P.TAB}/qflag_validation.json")).get("pass_", False))


def load_many(dirs, d0, d1):
    As, FDs, dss = [], [], []
    for d in dirs:
        if os.path.isdir(d):
            A, FD, ds = P.load_run(d)
            keep = [x for x in ds if d0 <= x <= d1]
            As.append(A[A.date.isin(keep)]); FDs.append(FD[FD.date.isin(keep)]); dss += keep
    if not As:
        return None, None, []
    return pd.concat(As, ignore_index=True), pd.concat(FDs, ignore_index=True), sorted(set(dss))


rows = []
# ---- diagnostic: B4 (p5f 2010-13) with each flag vs PNSN
Q = P.pnsn_quiet("2010-01-01", "2013-12-31"); pq = set(Q[Q.quiet == True].date)        # noqa: E712
yrs = [y for y in (2010, 2011, 2012, 2013) if os.path.exists(f"{P.BW}/P5F_{y}_DONE")]
A, FD, ds = load_many([f"{P.BW}/mf/feat/p5f_{y}" for y in yrs], "2010-01-10", "2013-12-31")
diag = []
if ds:
    known = set(F[F.Z.notna()].date)
    for lab, qs in [("pnsn", pq)] + [(f, QS[f]) for f in FLAGS]:
        dd = sorted(set(ds) & qs & (known if lab == "quiet_opt3" else set(ds)))
        r, _, _ = P.forward_rule(A, FD, dd)
        r.update(flag=lab, years=",".join(map(str, yrs)), n_quiet_days=len(dd),
                 n_quiet_days_pnsn_active=len(set(dd) - pq))
        diag.append(r)
    pd.DataFrame(diag).to_csv(f"{P.TAB}/p4_flagdiag_B4.csv", index=False)
# ---- P4-3: era detector, and the frozen network on the same days
E0, E1 = "2003-02-01", "2006-06-30"
era_dirs = [f"{P.BW}/mf/feat/{r}_{y}" for y in (2003, 2004, 2005, 2006) for r in ("p5e", "p5e2")]
frz_dirs = [f"{P.BW}/mf/feat/p5f_{y}" for y in (2003, 2004, 2005, 2006)]
SUB = [("CN4", "2003-02-01", "2005-06-30"), ("B1", "2005-07-01", "2005-10-31"), ("B2", "2005-11-01", "2006-06-30")]
S2 = pd.read_csv(f"{P.TAB}/p4_2_sets.csv").set_index("set") if os.path.exists(f"{P.TAB}/p4_2_sets.csv") else None
cross = []
Ae, FDe, dse = load_many(era_dirs, E0, E1)
Af, FDf, dsf = load_many(frz_dirs, E0, E1)
for fl in FLAGS:
    if dse:
        dd = sorted(set(dse) & QS[fl])
        r, _, _ = P.forward_rule(Ae, FDe, dd)
        r.update(detector="era (S15ERA)", span=f"{E0}..{E1}", flag=fl, n_days=len(dse), n_quiet_days=len(dd))
        rows.append(r)
        for y in (2003, 2004, 2005, 2006):
            ddy = [x for x in dd if x[:4] == str(y)]
            if ddy:
                ry, _, _ = P.forward_rule(Ae, FDe, ddy)
                ry.update(detector="era (S15ERA)", span=str(y), flag=fl, n_days=len([x for x in dse if x[:4] == str(y)]),
                          n_quiet_days=len(ddy))
                rows.append(ry)
    if dsf:
        for name, a, b in SUB:
            dd = sorted(x for x in set(dsf) & QS[fl] if a <= x <= b)
            if not dd:
                continue
            r, _, _ = P.forward_rule(Af, FDf, dd)
            r.update(detector=f"frozen ({name})", span=f"{a}..{b}", flag=fl,
                     n_days=len([x for x in dsf if a <= x <= b]), n_quiet_days=len(dd))
            rows.append(r)
            if S2 is not None and name in S2.index:
                for k in ("T_core", "T_ext"):
                    t42 = float(S2.loc[name, k]); tq = r[k]
                    cross.append(dict(flag=fl, set=name, tier=k, T_option_days=tq, T_p4_2=t42, diff=tq - t42,
                                      within_one_step=bool(np.isfinite(tq) and abs(tq - t42) <= 0.25 + 1e-9),
                                      n_quiet_family_days=r["n_quiet_family_days"], status=r["status"]))
R = pd.DataFrame(rows); R.to_csv(f"{P.TAB}/p4_3_era.csv", index=False)
X = pd.DataFrame(cross); X.to_csv(f"{P.TAB}/p4_4_crosscheck.csv", index=False)
# ---- decision: era thresholds per flag; bracket when the frozen cross-check is off by more than one step
dec = dict(flag_validation_pass=ok_val, primary_flag="quiet_opt3" if ok_val else "none (flag failed: all reported)", flags={})
for fl in FLAGS:
    e = R[(R.detector == "era (S15ERA)") & (R.span == f"{E0}..{E1}") & (R.flag == fl)] if len(R) else R
    if not len(e):
        continue
    e = e.iloc[0]
    out = dict(T_core=e.T_core, T_ext=e.T_ext, n_quiet_family_days=int(e.n_quiet_family_days), status=e.status)
    if len(X):
        for k in ("T_core", "T_ext"):
            x = X[(X.flag == fl) & (X.set == "CN4") & (X.tier == k)]      # the frozen set of most of the era
            if len(x):
                off = float(-x["diff"].iloc[0])
                out[f"{k}_crosscheck_CN4_diff"] = float(x["diff"].iloc[0])
                out[f"{k}_bracket"] = sorted([float(e[k]), float(e[k]) + off]) if not bool(x.within_one_step.iloc[0]) \
                    else [float(e[k]), float(e[k])]
    dec["flags"][fl] = out
json.dump(dec, open(f"{P.TAB}/p4_3_decision.json", "w"), indent=1, default=float)
pd.set_option("display.width", 250)
if diag:
    print(pd.DataFrame(diag)[["flag", "years", "n_quiet_days", "n_quiet_days_pnsn_active", "n_quiet_family_days", "T_core",
                              "T_ext", "fa_quiet_at_8p75", "fa_quiet_at_9p00"]].round(4).to_string(index=False))
if len(R):
    print(R[["detector", "span", "flag", "n_days", "n_quiet_days", "n_quiet_family_days", "T_core", "T_ext",
             "fa_quiet_at_8p75", "status"]].round(4).to_string(index=False))
if len(X):
    print(X.to_string(index=False))
print(json.dumps(dec, indent=1, default=float))
