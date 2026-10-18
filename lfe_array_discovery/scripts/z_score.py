#!/usr/bin/env python
"""STEP 1c' STAGE 10 -- SCORING and the PRE-REGISTERED RULE.

Same relative-clause scorer as Step 1c (`q_rediscover.score`, 3-km cluster truth reported
alongside per the `scoring-truth-colocated-families` memory), with Merlin's two corrections:

  * the relative clauses use the **2005 ORACLE REFERENCE** (results/step1c_prime/
    oracle_reference_2005.csv) instead of the transferred 2012 truth-template table;
  * the site-shift null is judged on SURVIVORS as well as gated families.

RULE (pre-registered, written before any 1c' discovery output was read):
  PASS       >= 3/10 isolated (centroid <= 5 km, recall >= 0.5 x and purity >= 0.8 x the oracle)
             -- or, when the oracle's C_iso < 4, the fallback bar >= max(ceil(0.5 C_30), 5) of the
             30 hidden families -- AND spurious split <= 0.10 AND null gated <= 1 per realisation
             AND null survivors <= 10 % of real survivors.
  FAIL       < 2/10 (or < 0.4 x C_iso) with clean nulls  -> STOP Arm A (north AND WA).
  FA_BROKEN  null survivors >= 50 % of real survivors    -> STOP Arm A.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import w_score_step1 as SC                                                  # noqa: E402

NULL_SURV_PASS_FRAC = 0.10
NULL_SURV_BROKEN_FRAC = 0.50
FAIL_FRAC_OF_CEILING = 0.40


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=f"{W.RES}/step1c_prime")
    ap.add_argument("--oracle", default=f"{W.RES}/step1c_prime/oracle_reference_2005.csv")
    ap.add_argument("--out", default=f"{W.RES}/step1c_prime_score")
    a = ap.parse_args()
    import q_rediscover as Q
    import g_rediscover as GR
    import p_common as P
    import s_common as S

    F = P.families()
    fam = {int(b): (float(F.lat[b]), float(F.lon[b])) for b in F.index}
    B = P.bostock()
    B = B[B.year == W.STEP1_YEAR]
    events = {int(b): np.sort(g.t_origin.to_numpy(float)) for b, g in B.groupby("bkey")}
    hid, iso, _ = S.holdouts()
    hidden, isolated = sorted(hid), sorted(iso)
    assert isolated == W.ISOLATED_HIDDEN, (isolated, W.ISOLATED_HIDDEN)
    windows = [(pd.Timestamp(W.STEP1_DAY0).timestamp(),
                pd.Timestamp(W.STEP1_DAY1).timestamp() + 86400)]

    ORA = pd.read_csv(a.oracle).set_index("bkey")
    oj = json.load(open(a.oracle.replace(".csv", ".json")))
    ref = ORA[["recall_ref", "purity_ref"]].copy()

    C = (pd.read_parquet(f"{a.run}/lg_clusters.parquet")
         if os.path.exists(f"{a.run}/lg_clusters.parquet")
         else pd.DataFrame(columns=["cid", "t", "lat", "lon"]))
    clusters = {c: dict(times=g.t.to_numpy(float), lat=float(g.lat.iloc[0]),
                        lon=float(g.lon.iloc[0])) for c, g in C.groupby("cid")} if len(C) else {}
    R, SP, summ = Q.score(clusters, events, fam, hidden, ref, windows,
                          tol=W.SCORE_TOL_S, r_km=W.SCORE_R_KM,
                          recall_frac=W.SCORE_RECALL_FRAC, purity_frac=W.SCORE_PURITY_FRAC,
                          max_explained=W.SCORE_MAX_EXPLAINED, split_frac=W.SCORE_SPLIT_FRAC)
    RI = R[R.hidden.isin(isolated)] if len(R) else R
    n_iso = int(RI.rediscovered.sum()) if len(RI) else 0
    n_all = int(R.rediscovered.sum()) if len(R) else 0

    CT = SC.cluster_truth(fam, events, hidden)
    rows = []
    for h in isolated:
        best = None
        for c, v in clusters.items():
            t = v["times"]
            if not len(t):
                continue
            d = float(S.epi_km(v["lat"], v["lon"], *fam[h]))
            mr, _ = GR.one2one(events[h], t, W.SCORE_TOL_S)
            mp, _ = GR.one2one(t, events[h], W.SCORE_TOL_S)
            mc, _ = GR.one2one(t, CT[h], W.SCORE_TOL_S)
            r = dict(hidden=h, cluster=c, dist_km=d, n_cluster=len(t),
                     recall=float(mr.mean()) if len(events[h]) else 0.0,
                     purity=float(mp.mean()), purity_cluster3km=float(mc.mean()),
                     oracle_recall=float(ORA.recall_ref.get(h, np.nan)),
                     oracle_purity=float(ORA.purity_ref.get(h, np.nan)))
            if best is None or r["recall"] > best["recall"]:
                best = r
        rows.append(best or dict(hidden=h, cluster=None, recall=0.0, purity=np.nan,
                                 purity_cluster3km=np.nan, dist_km=np.nan, n_cluster=0,
                                 oracle_recall=float(ORA.recall_ref.get(h, np.nan)),
                                 oracle_purity=float(ORA.purity_ref.get(h, np.nan))))
    CTR = pd.DataFrame(rows)
    CTR.to_csv(a.out + "_isolated_detail.csv", index=False)

    # ---------------------------------------------------------------- nulls and survivors
    def n_surv(p):
        S_ = W.read_csv_safe(p, W.SURV_COLS)
        return int(S_.survives.astype(bool).sum()) if len(S_) else 0

    real_p0 = n_surv(f"{a.run}/p0_survival.csv")
    real_fin = n_surv(f"{a.run}/p{W.N_ITERATIONS}_survival.csv")
    N = W.read_csv_safe(f"{a.run}/null_summary.csv",
                        ["realisation", "n_gate_pass", "n_survivors_p0", "n_survivors_final"])
    def col(name):
        if not len(N) or name not in N.columns:
            return []
        return [int(x) for x in pd.to_numeric(N[name], errors="coerce").fillna(-1)]

    null_gated = col("n_gate_pass")
    null_p0 = [x for x in col("n_survivors_p0") if x >= 0]
    null_fin = [x for x in col("n_survivors_final") if x >= 0]
    if null_gated and min(null_gated) < 0:
        raise SystemExit("null_summary.csv has a missing n_gate_pass -- re-run stage 9")
    if len(null_p0) != len(null_gated):
        raise SystemExit("null_summary.csv carries no survivor counts -- it was written by the "
                         "1c script, not z_null.py; re-run stage 9")
    max_null_p0 = max(null_p0) if null_p0 else 0
    surv_frac = (max_null_p0 / real_p0) if real_p0 > 0 else (np.inf if max_null_p0 else 0.0)

    c_iso, c_30 = int(oj["C_iso"]), int(oj["C_30"])
    fb = bool(oj["fallback_active"])
    bar = int(oj["pass_bar"]["value"])
    got = n_all if fb else n_iso
    fail_below = (max(2, int(np.ceil(FAIL_FRAC_OF_CEILING * c_30)))
                  if fb else max(2, int(np.ceil(FAIL_FRAC_OF_CEILING * c_iso))))
    split = summ.get("spurious_split_rate")
    split_ok = (split is None or not np.isfinite(split) or split <= W.ACCEPT_SPLIT_RATE)
    null_gated_ok = (not null_gated) or max(null_gated) <= W.ACCEPT_NULL_GATED
    null_surv_ok = surv_frac <= NULL_SURV_PASS_FRAC
    fa_broken = surv_frac >= NULL_SURV_BROKEN_FRAC

    if fa_broken:
        verdict = "FA_BROKEN"
    elif got >= bar and split_ok and null_gated_ok and null_surv_ok:
        verdict = "PASS"
    elif got < fail_below and null_gated_ok and null_surv_ok:
        verdict = "FAIL"
    else:
        verdict = "BELOW_BAR"

    out = dict(
        year=int(W.STEP1_YEAR), run=a.run, bed_mode=W.BED_MODE,
        repair_flags=W.repair_flags(),
        reference="results/step1c_prime/oracle_reference_2005.csv (2005 ORACLE, not the "
                  "transferred 2012 truth-template table)",
        oracle=dict(C_iso=c_iso, C_30=c_30, fallback_active=fb, bar=bar,
                    mode=oj["pass_bar"]["mode"]),
        n_clusters=len(clusters),
        n_isolated_rediscovered=n_iso, n_hidden_rediscovered_all30=n_all,
        scored_against_bar=got, bar=bar, fail_below=fail_below,
        isolated_no_reference=int(RI.no_reference.sum())
        if len(RI) and "no_reference" in RI else 0,
        spurious_split_rate=split, n_kept_split=summ.get("n_kept_split"),
        median_centroid_km=float(np.nanmedian(CTR.dist_km)) if len(CTR) else None,
        median_recall=float(np.nanmedian(CTR.recall)) if len(CTR) else None,
        median_purity=float(np.nanmedian(CTR.purity)) if len(CTR) else None,
        median_purity_cluster3km=float(np.nanmedian(CTR.purity_cluster3km)) if len(CTR) else None,
        real_survivors_p0=real_p0, real_survivors_final=real_fin,
        null_gated_per_realisation=null_gated,
        null_survivors_p0=null_p0, null_survivors_final=null_fin,
        null_survivor_fraction=float(surv_frac) if np.isfinite(surv_frac) else None,
        accept=dict(bar=bar, split_le=W.ACCEPT_SPLIT_RATE,
                    null_gated_le=W.ACCEPT_NULL_GATED,
                    null_survivors_le_frac=NULL_SURV_PASS_FRAC,
                    fa_broken_ge_frac=NULL_SURV_BROKEN_FRAC,
                    fail_below=fail_below),
        verdict=verdict,
        verdict_line=(f"STEP 1c' {verdict}: {got} of "
                      f"{'30 hidden (fallback bar %d)' % bar if fb else '10 isolated (bar %d)' % bar}"
                      f"; split {split}; null gated {null_gated}; null survivors {null_p0} vs real "
                      f"{real_p0}"),
        consequence=dict(
            PASS="WA_GO may be created once freeze_wa_v2 and check_constants pass (Merlin Q3)",
            FAIL="STOP Arm A for north AND WA; no further beds, no further seed rules",
            FA_BROKEN="STOP Arm A: the family-level false-alarm control is broken",
            BELOW_BAR="neither PASS nor the STOP clause; report and hand to Merlin"
        )[verdict])
    if len(R):
        R.to_csv(a.out + "_rediscovery.csv", index=False)
    if len(SP):
        SP.to_csv(a.out + "_splits.csv", index=False)
    json.dump(out, open(a.out + ".json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
