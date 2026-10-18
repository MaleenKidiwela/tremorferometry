#!/usr/bin/env python
"""P4-2b step 1: burst-robust control-matched T_core for CN4 / B1 / B2 / B3 (and S15fix itself) on the 509 E0 PNSN-quiet
days, statistics MED and TRIM (top-5 % S15fix control days dropped), POOLED (P4-2) beside them; pass rule of P4B.md.
Also verifies the "25 PGC burst days carry 85 % of the CN4 control" claim and the 2003-05 CN4 control rate at 8.75.
  python bw_p4b_e0.py -> p4b/tables/p4b_e0_thresholds.csv, p4b_e0_curves.csv, p4b_e0_burstcheck.json, p4b_e0_decision.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import bw_p4b_common as B

P4 = f"{B.BW}/p4"; OUT = f"{B.BW}/p4b/tables"
days = json.load(open(f"{P4}/tables/p4_e0_quiet_days.json"))
ZC = {"CN4": ["z_PGC", "z_VGZ", "z_NLLB", "z_SNB"]}
L = {}
for s in ("S15fix", "CN4", "B1", "B2", "B3"):
    L[s] = B.load_ctl(f"{P4}/mf/feat/p4e0_{s}", days, zcols=ZC.get(s))
    print(s, len(L[s][0]), "days", int(L[s][0].m_c1.sum()), "c1 family-days", flush=True)
D0, C0 = L["S15fix"][0], L["S15fix"][1]
rate_ref = B.day_rate(D0, C0, 8.75)
drop = B.top_days(rate_ref, 0.05)
keep = [d for d in D0.date if d not in set(drop)]
REF = B.stats(D0, C0); REFT = B.stats(D0, C0, keep)
tgt = dict(MED_core=B.at(REF["MED"], 8.75), MED_ext=B.at(REF["MED"], 8.25), TRIM_core=B.at(REFT["POOLED"], 8.75),
           TRIM_ext=B.at(REFT["POOLED"], 8.25), POOLED_core=B.at(REF["POOLED"], 8.75), POOLED_ext=B.at(REF["POOLED"], 8.25))
print("targets", {k: round(v, 4) for k, v in tgt.items()}, "dropped", len(drop))
P42 = json.load(open(f"{P4}/tables/p4_2_decision.json"))["sets"]
rows, curves = [], []
for s in ("CN4", "B1", "B2", "B3", "S15fix"):
    D, C = L[s][0], L[s][1]
    a = B.stats(D, C); t = B.stats(D, C, keep)
    r = dict(set=s, n_days=a["n_days"], n_c1_fd=a["n_fd"], n_days_trim=t["n_days"])
    for st, cur, key in (("MED", a["MED"], "MED"), ("TRIM", t["POOLED"], "TRIM"), ("POOLED", a["POOLED"], "POOLED")):
        for tier in ("core", "ext"):
            T = B.matched(cur, tgt[f"{key}_{tier}"])
            r[f"T_{st}_{tier}"] = T; r[f"T_{st}_{tier}_grid"] = B.up25(T); r[f"T_{st}_{tier}_floor"] = bool(T <= 8.0 + 1e-9)
        r[f"{st}_at_8p75"] = B.at(cur, 8.75); r[f"{st}_at_8p25"] = B.at(cur, 8.25)
    r["T_P42_pooled_S15fix_core"] = P42[s].get("T_matched_S15fix_core", 8.75 if s == "S15fix" else np.nan)
    r["T_P42_rule_core"] = P42[s]["T_core"]
    r["in_window_MED"] = bool(8.50 - 1e-9 <= r["T_MED_core"] <= 9.00 + 1e-9)
    r["in_window_TRIM"] = bool(8.50 - 1e-9 <= r["T_TRIM_core"] <= 9.00 + 1e-9)
    r["MED_vs_TRIM"] = r["T_MED_core"] - r["T_TRIM_core"]
    r["bracket"] = bool(abs(r["MED_vs_TRIM"]) > 0.5 + 1e-9)
    rows.append(r)
    curves.append(pd.DataFrame(dict(set=s, T=B.TG, MED=a["MED"], TRIM=t["POOLED"], POOLED=a["POOLED"])))
# diagnostic only (added after the first computation, not a pass rule): TRIM with each set's OWN top-5 % control days
for r in rows:
    D, C = L[r["set"]][0], L[r["set"]][1]
    own = B.top_days(B.day_rate(D, C, 8.75), 0.05)
    t2 = B.stats(D, C, [d for d in D.date if d not in set(own)])
    r["T_TRIMown_core_diag"] = B.matched(t2["POOLED"], tgt["TRIM_core"])
    r["TRIMown_overlap_S15fix_days"] = len(set(own) & set(drop))
R = pd.DataFrame(rows)
R.to_csv(f"{OUT}/p4b_e0_thresholds.csv", index=False)
pd.concat(curves).round(6).to_csv(f"{OUT}/p4b_e0_curves.csv", index=False)
fr = R[R.set.isin(["CN4", "B1", "B2", "B3"])]
passed = bool((fr.in_window_MED & fr.in_window_TRIM).all())
dec = dict(rule="P4B.md step 1", pass_=passed, targets=tgt, n_days=len(days), dropped_days=drop,
           adopted=(dict(T_core={k: 8.75 for k in ("CN4", "B1", "B2", "B3", "S15")}, T_ext=8.00) if passed else None),
           brackets={r.set: [min(r.T_MED_core, r.T_TRIM_core), max(r.T_MED_core, r.T_TRIM_core)] for r in fr.itertuples() if r.bracket},
           not_adopted=[r.set for r in fr.itertuples() if not (r.in_window_MED and r.in_window_TRIM) and not r.bracket])
json.dump(dec, open(f"{OUT}/p4b_e0_decision.json", "w"), indent=1, default=float)
pd.set_option("display.width", 250)
print(R[["set", "T_MED_core", "T_TRIM_core", "T_POOLED_core", "T_P42_pooled_S15fix_core", "T_MED_ext", "T_TRIM_ext",
         "MED_at_8p75", "TRIM_at_8p75", "POOLED_at_8p75", "T_TRIMown_core_diag", "TRIMown_overlap_S15fix_days", "in_window_MED", "in_window_TRIM", "bracket"]].round(4).to_string(index=False))
print("PASS" if passed else "FAIL", json.dumps({k: dec[k] for k in ("brackets", "not_adopted")}))

# ---- verification of the 25-burst-day claim
D, C, Zs = L["CN4"][0], L["CN4"][1], L["CN4"][2]
bc = {}
for T in (8.0, 8.75):
    n = pd.Series({d: int(len(C[d]) - np.searchsorted(C[d], T, side="left")) for d in D.date})
    top = list(n.sort_values(ascending=False, kind="mergesort").index[:25])
    bc[f"T{T}"] = dict(total=int(n.sum()), top25=int(n[top].sum()), share_top25=float(n[top].sum() / max(n.sum(), 1)),
                       top25_days=top, n_feb2015=sum(d[:7] == "2015-02" for d in top),
                       n_after_2015_01_12=sum(d >= "2015-01-12" for d in top),
                       months=pd.Series([d[:7] for d in top]).value_counts().to_dict(),
                       overlap_with_S15fix_trim=len(set(top) & set(drop)))
    zz = pd.concat([Zs[d].assign(date=d) for d in top if Zs[d] is not None and len(Zs[d])], ignore_index=True)
    zz = zz[zz.index.isin(zz.index)]
    # keep only candidates >= T: Zs rows are aligned with all c1 rows of the day (unsorted); recompute from files
    bc[f"T{T}"]["note_z"] = "PGC-top fraction over all stored c1 candidates (>= 8.00) on the top-25 days"
    zc = ["z_PGC", "z_VGZ", "z_NLLB", "z_SNB"]
    top_site = zz[zc].fillna(-99).to_numpy().argmax(1)
    bc[f"T{T}"]["frac_PGC_top_site"] = float((top_site == 0).mean())
    bc[f"T{T}"]["median_dom_top"] = float(zz.dom_top.median()) if "dom_top" in zz else None
    rest = [d for d in D.date if d not in set(top) and Zs[d] is not None and len(Zs[d])]
    zr = pd.concat([Zs[d] for d in rest], ignore_index=True)
    bc[f"T{T}"]["frac_PGC_top_site_other_days"] = float((zr[zc].fillna(-99).to_numpy().argmax(1) == 0).mean())
    bc[f"T{T}"]["median_dom_top_other_days"] = float(zr.dom_top.median())
# 2003-05 CN4 control rate at 8.75 (p5f, CN4-only days)
real = {}
for y in (2003, 2004, 2005):
    Dy, Cy, _, _ = B.load_ctl(f"{B.BW}/mf/feat/p5f_{y}")
    Dy = Dy[Dy.stations.apply(lambda s: bool(s) and set(s.split(",")) <= {"PGC", "VGZ", "NLLB", "SNB"})]
    sy = B.stats(Dy, Cy)
    real[y] = dict(n_days=sy["n_days"], pooled_8p75=B.at(sy["POOLED"], 8.75), med_8p75=B.at(sy["MED"], 8.75))
bc["p5f_CN4_2003_2005_ctl_rate_8p75"] = real
bc["E0_CN4_ctl_rate_8p75"] = dict(pooled=float(R.set_index("set").loc["CN4", "POOLED_at_8p75"]),
                                  med=float(R.set_index("set").loc["CN4", "MED_at_8p75"]),
                                  trim=float(R.set_index("set").loc["CN4", "TRIM_at_8p75"]))
bc["E0_S15fix_ctl_rate_8p75"] = dict(pooled=tgt["POOLED_core"], med=tgt["MED_core"], trim=tgt["TRIM_core"])
json.dump(bc, open(f"{OUT}/p4b_e0_burstcheck.json", "w"), indent=1, default=float)
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("top25_days",)} for k, v in bc.items()}, indent=1, default=float))
