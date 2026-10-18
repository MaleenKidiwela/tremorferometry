#!/usr/bin/env python
"""P4-1 anchor: the forward false-alarm rule (c_tiers_by_epoch.py) on B4 = the P5 frozen-network runs p5f_2010..2013
(PNSN quiet days), plus the self-check on the forward production E0 (refine/feat/prod, 2014-01..2017-02).
Pass (P4.md): |T_core - 8.75| <= 0.25 and |T_ext - 8.25| <= 0.25.  Years used = those with P5F_<Y>_DONE.

  python bw_p4_anchor.py            -> p4/tables/p4_1_anchor.csv, p4_1_curves.csv, p4_1_decision.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import bw_p4_common as P

Q = P.pnsn_quiet("2010-01-01", "2017-02-28")
quiet = set(Q[Q.quiet == True].date); post = set(Q[Q.post == True].date)            # noqa: E712
rows, curves = [], []


def do(label, featdir, dates, years):
    A, FD, ds = P.load_run(featdir, dates)
    qd = sorted(set(ds) & quiet); pdd = sorted(set(ds) & post)
    r, q, p_ = P.forward_rule(A, FD, qd, pdd)
    r.update(run=label, years=years, n_days=len(ds), n_quiet_days=len(qd), n_post_days=len(pdd))
    rows.append(r)
    curves.append(pd.DataFrame(dict(run=label, T=P.GRID, fa_quiet=q, fa_post=p_)))
    P.log(f"[anchor] {label}: T_core {r['T_core']} T_ext {r['T_ext']} (quiet fd {r['n_quiet_family_days']}, "
          f"FA@8.75 {r['fa_quiet_at_8p75']:.4f}, FA@9.00 {r['fa_quiet_at_9p00']:.4f})")
    return r


# self-check: forward production E0
e0 = [str(d.date()) for d in pd.date_range("2014-01-01", "2017-02-28")]
do("selfcheck_forward_prod_E0", f"{P.REFINE}/feat/prod", e0, "2014-01..2017-02")
# B4 = p5f per year, then pooled
yrs = [y for y in (2010, 2011, 2012, 2013) if os.path.exists(f"{P.BW}/P5F_{y}_DONE")]
for y in yrs:
    do(f"B4_p5f_{y}", f"{P.BW}/mf/feat/p5f_{y}", None, str(y))
if yrs:
    # pooled: concatenate the year runs
    As, FDs, dss = [], [], []
    for y in yrs:
        A, FD, ds = P.load_run(f"{P.BW}/mf/feat/p5f_{y}")
        As.append(A); FDs.append(FD); dss += ds
    A = pd.concat(As, ignore_index=True); FD = pd.concat(FDs, ignore_index=True)
    qd = sorted(set(dss) & quiet); pdd = sorted(set(dss) & post)
    r, q, p_ = P.forward_rule(A, FD, qd, pdd)
    lab = f"B4_p5f_{yrs[0]}-{yrs[-1]}"
    r.update(run=lab, years=",".join(map(str, yrs)), n_days=len(dss), n_quiet_days=len(qd), n_post_days=len(pdd))
    rows.append(r); curves.append(pd.DataFrame(dict(run=lab, T=P.GRID, fa_quiet=q, fa_post=p_)))
    ok = (abs(r["T_core"] - 8.75) <= 0.25 + 1e-9) and (abs(r["T_ext"] - 8.25) <= 0.25 + 1e-9) and r["status"] == "derived"
    final = yrs == [2010, 2011, 2012, 2013]
    dec = dict(years=yrs, final=final, T_core=r["T_core"], T_ext=r["T_ext"], n_quiet_family_days=r["n_quiet_family_days"],
               n_post_family_days=r["n_post_family_days"], fa_quiet_at_8p75=r["fa_quiet_at_8p75"],
               fa_quiet_at_8p25=r["fa_quiet_at_8p25"], fa_post_at_8p25=r["fa_post_at_8p25"],
               rule="|T_core-8.75|<=0.25 and |T_ext-8.25|<=0.25", pass_=bool(ok),
               verdict=("PASS" if ok else "FAIL") + ("" if final else " (provisional: not all of 2010-13 yet)"))
    json.dump(dec, open(f"{P.TAB}/p4_1_decision.json", "w"), indent=1)
    P.log(f"[anchor] {lab}: {dec['verdict']}  T_core {r['T_core']} T_ext {r['T_ext']}")
R = pd.DataFrame(rows)
R.to_csv(f"{P.TAB}/p4_1_anchor.csv", index=False)
pd.concat(curves).to_csv(f"{P.TAB}/p4_1_curves.csv", index=False)
pd.set_option("display.width", 250)
print(R[["run", "n_days", "n_quiet_days", "n_quiet_family_days", "n_post_family_days", "T_core", "T_ext",
         "fa_quiet_at_8p25", "fa_quiet_at_8p75", "fa_quiet_at_9p00", "fa_post_at_8p25", "status"]].round(4).to_string(index=False))
