#!/usr/bin/env python
"""P4-2b step 2: era-detector thresholds on its env-quiet days (quiet_env3_CN4, 2003-02-01..2006-06-30), control-matched
to the S15fix E0 levels (MED, TRIM; POOLED beside), and the frozen-CN4 same-day cross-check (8.75 +- 0.25, else bracket).
Scrambled control only; data FA never used.  Rules: results/backward/p4b/P4B.md.
  python bw_p4b_era.py -> p4b/tables/p4b_era_thresholds.csv, p4b_era_decision.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import bw_p4b_common as B

OUT = f"{B.BW}/p4b/tables"
tgt = json.load(open(f"{OUT}/p4b_e0_decision.json"))["targets"]
F = pd.read_csv(f"{B.BW}/p4/envflag/tables/quiet_flag_env.csv", dtype={"date": str})
Q = set(F[(F.quiet_env3_CN4 == True) & F.date.between("2003-02-01", "2006-06-30")].date)       # noqa: E712
ERA = [f"{B.BW}/mf/feat/{r}" for r in ("p5e_2003_kelb", "p5e_2004", "p5e_2005", "p5e_2006")]
FRZ = [f"{B.BW}/mf/feat/p5f_{y}" for y in (2003, 2004, 2005, 2006)]


def load(dirs):
    Ds, Cs = [], {}
    for d in dirs:
        D, C, _, _ = B.load_ctl(d, sorted(Q)); Ds.append(D); Cs.update(C)
    return pd.concat(Ds, ignore_index=True), Cs


De, Ce = load(ERA); Df, Cf = load(FRZ)
common = sorted(set(De[De.m_c1 > 0].date) & set(Df[Df.m_c1 > 0].date))
De = De[De.date.isin(common)].reset_index(drop=True); Df = Df[Df.date.isin(common)].reset_index(drop=True)
print("env-quiet days", len(Q), "with both runs", len(common), flush=True)
cn4 = set(Df[(Df.date <= "2005-06-30") & Df.stations.apply(lambda s: bool(s) and set(s.split(",")) <= {"PGC", "VGZ", "NLLB", "SNB"})].date)


def derive(D, C, days, rank_D, rank_C, label, detector):
    E = D[D.date.isin(set(days))]
    rr = B.day_rate(rank_D[rank_D.date.isin(set(days))], rank_C, 8.75)
    drop = set(B.top_days(rr, 0.05))
    a = B.stats(E, C); t = B.stats(E, C, [d for d in E.date if d not in drop])
    r = dict(detector=detector, span=label, n_days=a["n_days"], n_c1_fd=a["n_fd"], n_dropped=len(drop))
    for st, cur, key in (("MED", a["MED"], "MED"), ("TRIM", t["POOLED"], "TRIM"), ("POOLED", a["POOLED"], "POOLED")):
        for tier in ("core", "ext"):
            T = B.matched(cur, tgt[f"{key}_{tier}"])
            r[f"T_{st}_{tier}"] = T; r[f"T_{st}_{tier}_grid"] = B.up25(T)
        r[f"{st}_at_8p75"] = B.at(cur, 8.75)
    return r, drop


rows = []
full = sorted(common)
r_era, drop_full = derive(De, Ce, full, Df, Cf, "2003-02..2006-06", "era (S15ERA+KELB)"); rows.append(r_era)
r_frz_full, _ = derive(Df, Cf, full, Df, Cf, "2003-02..2006-06", "frozen (as alive)"); rows.append(r_frz_full)
# diagnostic: era trimmed by its own control ranking
r_d, _ = derive(De, Ce, full, De, Ce, "2003-02..2006-06 (diag: trim by era's own ranking)", "era (S15ERA+KELB)"); rows.append(r_d)
sub = sorted(cn4)
r_xf, _ = derive(Df, Cf, sub, Df, Cf, "CN4 days 2003-02..2005-06", "frozen CN4"); rows.append(r_xf)
r_xe, _ = derive(De, Ce, sub, Df, Cf, "CN4 days 2003-02..2005-06", "era (S15ERA+KELB)"); rows.append(r_xe)
for y in ("2003", "2004", "2005", "2006"):
    dy = [d for d in full if d[:4] == y]
    ry, _ = derive(De, Ce, dy, Df, Cf, y, "era (S15ERA+KELB)"); rows.append(ry)
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/p4b_era_thresholds.csv", index=False)
dec = dict(rule="P4B.md step 2", n_env_quiet_days=len(Q), n_days_both_runs=len(full), n_days_cn4_crosscheck=len(sub), stats={})
for st in ("MED", "TRIM"):
    te, tf = r_era[f"T_{st}_core"], r_xf[f"T_{st}_core"]
    ok = bool(abs(tf - 8.75) <= 0.25 + 1e-9)
    dec["stats"][st] = dict(era_T_core=te, frozen_CN4_T_core=tf, crosscheck_pass=ok,
                            era_T_core_final=[te, te] if ok else sorted([te, te + (8.75 - tf)]),
                            era_T_ext=r_era[f"T_{st}_ext"])
dec["MED_vs_TRIM_era"] = r_era["T_MED_core"] - r_era["T_TRIM_core"]
dec["bracket_MED_TRIM"] = bool(abs(dec["MED_vs_TRIM_era"]) > 0.5 + 1e-9)
lo = min(min(v["era_T_core_final"]) for v in dec["stats"].values()); hi = max(max(v["era_T_core_final"]) for v in dec["stats"].values())
dec["era_T_core_adopted"] = [lo, hi] if (dec["bracket_MED_TRIM"] or not all(v["crosscheck_pass"] for v in dec["stats"].values())) \
    else [r_era["T_MED_core"], r_era["T_TRIM_core"]]
dec["targets"] = tgt
json.dump(dec, open(f"{OUT}/p4b_era_decision.json", "w"), indent=1, default=float)
pd.set_option("display.width", 250)
print(R[["detector", "span", "n_days", "n_c1_fd", "n_dropped", "T_MED_core", "T_TRIM_core", "T_POOLED_core", "T_MED_ext",
         "T_TRIM_ext", "MED_at_8p75", "TRIM_at_8p75", "POOLED_at_8p75"]].round(4).to_string(index=False))
print(json.dumps(dec, indent=1, default=float))
