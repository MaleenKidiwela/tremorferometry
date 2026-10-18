#!/usr/bin/env python
"""Assemble the gate validation tables -> results/gates/g_tables.md (+ the CSVs behind them)."""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402

SETS = ["pos_bost_N20", "pos_bost_N50", "pos_bost_N100", "pos_bost_N200",
        "pos_tpl_N20", "pos_tpl_N50", "pos_tpl_N100",
        "neg_moveout_N50", "neg_timeshift_N50", "neg_dvv_N50", "neg_hf_N50", "alias_N50"]
G3 = ["g3_pos", "g3_neg_random", "g3_neg_activity"]


def load(s):
    p = f"{G.GATES}/g_val_{s}.csv"
    if not os.path.exists(p):
        return None
    D = pd.read_csv(p)
    D["set"] = s
    return D


def frac(d, col, val):
    if col not in d or d[col].notna().sum() == 0:
        return np.nan
    return float((d[col] == val).mean())


def md(df):
    """Minimal markdown table (no tabulate dependency)."""
    def f(v):
        if isinstance(v, float):
            return "" if not np.isfinite(v) else (f"{v:.3f}" if abs(v) < 1e4 else f"{v:.3g}")
        return str(v)
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(f(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def main():
    out = []
    A = pd.concat([x for x in (load(s) for s in SETS) if x is not None], ignore_index=True)
    A.to_csv(f"{G.GATES}/g_val_all.csv", index=False)

    # ---------------- table 1: verdicts per set
    rows = []
    for s in SETS:
        d = A[A.set == s]
        if not len(d):
            continue
        r = dict(set=s, n=len(d), N=int(d.n_members.median()))
        for g in ("g0", "g1", "g2"):
            for v in ("PASS", "FAIL", "UNTESTABLE"):
                r[f"{g}_{v[:4]}"] = int((d[g] == v).sum()) if g in d else 0
        r["overall_PASS"] = int(d.verdict.isin(["PASS", "PASS_inter_ETS_candidate"]).sum())
        r["runtime_s_med"] = float(d.runtime_s.median())
        rows.append(r)
    T1 = pd.DataFrame(rows)
    T1.to_csv(f"{G.GATES}/g_table_verdicts.csv", index=False)
    out += ["## Table 1 — verdicts per proposal set", "", md(T1), ""]

    # ---------------- table 2: gate 1 power / false pass
    rows = []
    for s in SETS:
        d = A[A.set == s]
        if not len(d) or "g1" not in d:
            continue
        t = d[d.g1.isin(["PASS", "FAIL"])]
        rows.append(dict(set=s, N=int(d.n_members.median()), n=len(d),
                         n_testable=len(t), untestable=int((d.g1 == "UNTESTABLE").sum()),
                         g1_pass=int((d.g1 == "PASS").sum()),
                         rate_of_testable=float((t.g1 == "PASS").mean()) if len(t) else np.nan,
                         rate_of_all=float((d.g1 == "PASS").mean()),
                         med_excess_db=float(d.g1_best_excess_db.median()),
                         med_best_dist_km=float(d.g1_best_dist_km.median())))
    T2 = pd.DataFrame(rows)
    T2.to_csv(f"{G.GATES}/g_table_gate1.csv", index=False)
    out += ["## Table 2 — Gate 1 (withheld-station stack): power and false-pass rate", "", md(T2), ""]

    # ---------------- table 3: gate 0 power by nearest-neighbour distance
    rows = []
    for s in ("pos_bost_N50", "pos_bost_N200", "pos_tpl_N50", "alias_N50"):
        d = A[A.set == s]
        if not len(d) or "g0" not in d:
            continue
        bins = [(0, 2.5), (2.5, 5.0), (5.0, 999)]
        for lo, hi in bins:
            k = d[(d.g0_min_dist_km >= lo) & (d.g0_min_dist_km < hi)]
            if not len(k):
                continue
            rows.append(dict(set=s, nn_km=f"{lo}-{hi if hi < 999 else 'inf'}", n=len(k),
                             g0_pass=int((k.g0 == "PASS").sum()),
                             g0_pass_rate=float((k.g0 == "PASS").mean()),
                             literal_rule_pass_rate=float(
                                 ((k.g0_min_dist_km >= 5.0) | k.g0_literal.astype(bool)).mean()),
                             med_frac_above=float(k.g0_frac_above.median()),
                             med_hf=float(k.g0_hf_median.median())))
    T3 = pd.DataFrame(rows)
    T3.to_csv(f"{G.GATES}/g_table_gate0.csv", index=False)
    out += ["## Table 3 — Gate 0 (identity/class) by distance to the nearest ACCEPTED family "
            "(leave-one-out for the positives)", "", md(T3), ""]

    # ---------------- table 4: gate 2
    rows = []
    for s in SETS:
        d = A[A.set == s]
        if not len(d) or "g2" not in d:
            continue
        t = d[d.g2.isin(["PASS", "FAIL"])]
        rows.append(dict(set=s, N=int(d.n_members.median()), n=len(d),
                         testable=len(t), untestable=int((d.g2 == "UNTESTABLE").sum()),
                         g2_pass=int((d.g2 == "PASS").sum()),
                         med_sites_60km=float(d.g2_n_sites_60km.median()),
                         med_dist_to_proposal_km=float(d.g2_dist_km.median()),
                         med_null_ratio=float(d.g2_null_ratio.median()),
                         med_synth_km=float(d.g2_synth_km.median()),
                         med_loso_km=float(d.g2_loso_km.median())))
    T4 = pd.DataFrame(rows)
    T4.to_csv(f"{G.GATES}/g_table_gate2.csv", index=False)
    out += ["## Table 4 — Gate 2 (location with controls)", "", md(T4), ""]

    # ---------------- table 5: gate 3
    rows = []
    for s in G3:
        d = load(s)
        if d is None:
            continue
        rows.append(dict(set=s, n=len(d), N=int(d.n_members.median()),
                         g3_pass=int((d.g3 == "PASS").sum()),
                         rate=float((d.g3 == "PASS").mean()),
                         med_E_local=float(d.g3_E_local.median()),
                         med_E_far=float(d.g3_E_far.median()),
                         med_ratio=float(d.g3_ratio.median()),
                         frac_p_lt_0p01=float((d.g3_p < 0.01).mean())))
    T5 = pd.DataFrame(rows)
    T5.to_csv(f"{G.GATES}/g_table_gate3.csv", index=False)
    out += ["## Table 5 — Gate 3 (tremor timing, long-span validation on catalog v1 2014-2026)", "",
            md(T5), ""]

    # ---------------- table 6: gate 1 controls
    C = pd.read_csv(f"{G.GATES}/g_gate1_controls.csv")
    T6 = (C.groupby(["sta", "N"]).apply(
        lambda g: pd.Series(dict(n=len(g), pass_=int((g.verdict == "PASS").sum()),
                                 rate=float((g.verdict == "PASS").mean()),
                                 med_excess_db=float(g.excess_db.median()),
                                 med_dist_km=float(g.dist_km.median()))), include_groups=False)
          .reset_index())
    T6.to_csv(f"{G.GATES}/g_table_controls.csv", index=False)
    out += ["## Table 6 — Gate 1 POSITIVE-CONTROL table (known Bostock families, their own event "
            "times)", "", md(T6), ""]

    # ---------------- runtime
    out += ["## Runtime", "",
            f"median per proposal (gates 0+1+2): **{A.runtime_s.median():.1f} s**, "
            f"p90 {A.runtime_s.quantile(.9):.1f} s, max {A.runtime_s.max():.1f} s "
            f"(n = {len(A)} proposals)", ""]
    open(f"{G.GATES}/g_tables.md", "w").write("\n".join(out))
    print("\n".join(out))


if __name__ == "__main__":
    main()
