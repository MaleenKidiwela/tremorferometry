#!/usr/bin/env python
"""P4-2 scoring: per reduced station set (CN4 / B1 / B2 / B3; extra S15 fixed set) on the E0 PNSN-quiet days:
the forward-rule T_core / T_ext (quiet FA <= 0.1 / <= 1.0 per template-day, hf_pass) and, beside it, the
control-matched T (G1 procedure: scrambled-control rate per observable control family-day = the reference's rate at
8.75 / 8.25 on the same days).  Reference = forward production refine/feat/prod (S15 + R2 exclusions); the S15 fixed
set (no exclusions) is reported as a second reference.  Rules: p4/P4.md.  Safe to run on partial runs (n_days shown).

  python bw_p4_sets.py   -> p4/tables/p4_2_sets.csv, p4_2_curves.csv, p4_2_decision.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import bw_p4_common as P

days = json.load(open(f"{P.TAB}/p4_e0_quiet_days.json"))
REF, rref = {}, {}
for lab, d in (("prod", f"{P.REFINE}/feat/prod"), ("S15fix", f"{P.P4}/mf/feat/p4e0_S15fix")):
    if os.path.isdir(d):
        REF[lab] = P.load_run(d, days)
rows, curves = [], []
ERA_USE = dict(CN4="2003-01..2005-06 (frozen network)", B1="2005-07..2005-10", B2="2005-11..2007-11",
               B3="2007-11..2009-08", S15fix="2009-08..2013 (= B4, no exclusions)")
for SET in ("CN4", "B1", "B2", "B3", "S15fix"):
    d = f"{P.P4}/mf/feat/p4e0_{SET}"
    if not os.path.isdir(d):
        continue
    A, FD, ds = P.load_run(d, days)
    if not ds:
        continue
    r, q, _ = P.forward_rule(A, FD, ds)
    r.update(set=SET, used_for=ERA_USE[SET], n_days=len(ds), complete=len(ds) == len(days))
    curves.append(pd.DataFrame(dict(set=SET, T=P.GRID, fa_quiet=q)))
    for rl, (RA, RFD, rds) in REF.items():
        if rl == SET:
            continue
        common = sorted(set(ds) & set(rds))
        for tier, Tref in (("core", 8.75), ("ext", 8.25)):
            tgt = P.ctl_rate(RA, RFD, Tref, common)
            tm = P.matched_T(A, FD, tgt, common)
            r[f"ctl_rate_ref_{rl}_{tier}"] = tgt
            r[f"T_matched_{rl}_{tier}"] = tm
            r[f"T_matched_{rl}_{tier}_grid"] = P.up_to_grid(tm)
            r[f"T_matched_{rl}_{tier}_floor"] = bool(tm <= 8.0 + 1e-9)
            r[f"n_common_days_{rl}"] = len(common)
    for tier, Tref in (("core", 8.75), ("ext", 8.25)):
        r[f"ctl_rate_self_at_{Tref}"] = P.ctl_rate(A, FD, Tref, ds)
    # one-step rule (primary reference = forward production)
    for tier, k in (("core", "T_core"), ("ext", "T_ext")):
        tg = r.get(f"T_matched_prod_{tier}_grid", np.nan)
        r[f"diff_rule_vs_matched_{tier}"] = r[k] - tg if np.isfinite(tg) else np.nan
        r[f"report_both_{tier}"] = bool(np.isfinite(tg) and abs(r[k] - tg) > 0.25 + 1e-9)
    rows.append(r)
    P.log(f"[sets] {SET}: {len(ds)}/{len(days)} days  rule {r['T_core']:.2f}/{r['T_ext']:.2f}  matched(prod) "
          f"{r.get('T_matched_prod_core', np.nan):.2f}/{r.get('T_matched_prod_ext', np.nan):.2f}")
# the references themselves (forward rule on the same days)
for rl, (RA, RFD, rds) in REF.items():
    if rl == "S15fix":
        continue
    r, q, _ = P.forward_rule(RA, RFD, rds)
    r.update(set=f"ref_{rl}", used_for="forward E0 reference", n_days=len(rds), complete=len(rds) == len(days))
    for tier, Tref in (("core", 8.75), ("ext", 8.25)):
        r[f"ctl_rate_self_at_{Tref}"] = P.ctl_rate(RA, RFD, Tref, rds)
    rows.append(r); curves.append(pd.DataFrame(dict(set=f"ref_{rl}", T=P.GRID, fa_quiet=q)))
R = pd.DataFrame(rows)
R.to_csv(f"{P.TAB}/p4_2_sets.csv", index=False)
pd.concat(curves).to_csv(f"{P.TAB}/p4_2_curves.csv", index=False)
cols = ["set", "n_days", "n_quiet_family_days", "T_core", "T_ext", "fa_quiet_at_8p75", "fa_quiet_at_9p00",
        "T_matched_prod_core", "T_matched_prod_core_grid", "T_matched_prod_ext", "T_matched_prod_ext_grid",
        "T_matched_S15fix_core", "T_matched_S15fix_ext", "report_both_core", "report_both_ext"]
pd.set_option("display.width", 300)
print(R[[c for c in cols if c in R]].round(4).to_string(index=False))
json.dump(dict(sets={r["set"]: {k: (v if not isinstance(v, (np.floating, np.bool_)) else v.item()) for k, v in r.items()}
                     for r in rows}, n_quiet_days=len(days)),
          open(f"{P.TAB}/p4_2_decision.json", "w"), indent=1, default=str)
