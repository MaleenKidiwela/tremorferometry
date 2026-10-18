#!/usr/bin/env python
"""P4-0b POST-HOC RE-SCORE with a CLUSTER truth -- NOT THE FROZEN VERDICT.

Merlin's close-out F1 (results/merlin_review_p4_closeout_2026-09-20.md): arm (i)'s purity counted a
hit correct only if it one-to-one matched a p_real >= 0.8 detection labelled to ONE family, but
97.5 % of hp events are co-detected by another family's template and the 18 eligible families sit
0.12-1.6 km from their neighbours.  A per-family single-label purity is therefore not a measure of
noise admission.  This re-score replaces the truth with a CLUSTER truth and keeps everything else
(the templates, the hits, the 0.5x / 0.8x relative clauses, the reference) identical.

  truth(X) = frozen-catalogue detections (realness_v1_1) on the usable episode days of every
             NON-HIDDEN family within 3.0 km hypocentral of X, deduplicated at 2 s
  primary  = p_real >= 0.2      secondary = any p_real
  scored   = the frozen it2 grown templates, the A2 ablation's it2 templates, and REF-X, all
             against the same truth, with recall >= 0.5 x recall_ref and purity >= 0.8 x purity_ref

Accept (Merlin): <= 6/18 confirms closure.  >= 12/18 is a FAILURE SIGNATURE and reverses the
CONSTRUCTION verdict (not the Puget decision).

  python scripts/v_rescore.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

FRZ = V.OUT
ABL = f"{V.RES}/p4/p4_0b_abl"
OUT = f"{V.RES}/p4/p4_0b_rescore"
CLUSTER_KM = 3.0
DEDUP_S = 2.0
LABEL = "POST-HOC RE-SCORE -- not the frozen verdict"


def dedup(t, sep=DEDUP_S):
    t = np.sort(np.asarray(t, float))
    if len(t) == 0:
        return t
    keep, last = [], -1e18
    for x in t:
        if x - last >= sep:
            keep.append(x)
            last = x
    return np.asarray(keep)


def build_truth():
    hid, _, _ = S.holdouts()
    F = S.families()
    keep, _, _ = V.eligible_families()
    A = pd.read_csv(f"{FRZ}/availability_episode_days.csv")
    usable = {str(pd.Timestamp(d).date()) for d, u in zip(A.date, A.usable) if u}
    D = pd.read_parquet("/home/jovyan/tremorferometry/lfe_extension/results/trackE/"
                        "realness_v1_1.parquet",
                        columns=["family_id", "t", "date", "p_real"])
    D = D.rename(columns={"family_id": "bkey"})
    D["bkey"] = D.bkey.astype(int)
    D = D[~D.bkey.isin(hid)]                     # hidden families never enter the truth
    D = D[D.date.astype(str).isin(usable)]
    cl, truth_p, truth_any = {}, {}, {}
    for X in keep:
        d = S.hyp_km(F.lat.values, F.lon.values, F.dep.values, F.lat[X], F.lon[X])
        # hypocentral distance between two family hypocentres
        dk = np.sqrt(V_epi(F, X) ** 2 + (F.dep.values - float(F.dep[X])) ** 2)
        mem = [int(b) for b, dd in zip(F.index, dk) if dd <= CLUSTER_KM and int(b) not in hid]
        cl[X] = mem
        g = D[D.bkey.isin(mem)]
        truth_p[X] = dedup(g[g.p_real >= 0.2].t.to_numpy(float))
        truth_any[X] = dedup(g.t.to_numpy(float))
    return keep, cl, truth_p, truth_any


def V_epi(F, X):
    return S.epi_km(F.lat.values, F.lon.values, float(F.lat[X]), float(F.lon[X]))


def rp(hits, truth):
    if len(truth) == 0:
        return np.nan, np.nan
    if len(hits) == 0:
        return 0.0, np.nan
    mr, _ = V.one2one(truth, hits, V.MATCH_TOL_S)
    mp, _ = V.one2one(hits, truth, V.MATCH_TOL_S)
    return float(mr.mean()), float(mp.mean())


def main():
    os.makedirs(OUT, exist_ok=True)
    keep, cl, TP, TA = build_truth()
    REF = pd.read_parquet(f"{FRZ}/run_ref_indep.parquet")
    FRO = pd.read_parquet(f"{FRZ}/run_given_it2_indep.parquet")
    A2 = pd.read_parquet(f"{ABL}/A2_run_it2_indep.parquet")
    rows = []
    for tag, TR in (("primary_p0.2", TP), ("secondary_anyp", TA)):
        for X in keep:
            tr = TR[X]
            rt = np.sort(REF[REF.tid == f"REF{X}"].t.to_numpy(float))
            rrec, rpur = rp(rt, tr)
            for src, G in (("frozen", FRO), ("A2", A2)):
                cand = [t for t in sorted(set(G.tid)) if str(t).startswith(f"G{X}_")]
                best = None
                for tid in cand:
                    ht = np.sort(G[G.tid == tid].t.to_numpy(float))
                    rec, pur = rp(ht, tr)
                    ok = bool(np.isfinite(rrec) and rrec > 0 and rec >= V.REL_RECALL * rrec
                              and np.isfinite(pur) and pur >= V.REL_PURITY * rpur)
                    r = dict(tid=tid, n_hits=len(ht), recall=rec, purity=pur,
                             rel_recall=rec / rrec if rrec else np.nan,
                             rel_purity=pur / rpur if rpur else np.nan, assembled=ok)
                    if best is None or (ok, r["recall"]) > (best["assembled"], best["recall"]):
                        best = r
                row = dict(truth=tag, source=src, bkey=int(X), n_cluster=len(cl[X]),
                           cluster=",".join(map(str, cl[X])), n_truth=len(tr),
                           n_ref=len(rt), recall_ref=rrec, purity_ref=rpur,
                           need_recall=V.REL_RECALL * rrec, need_purity=V.REL_PURITY * rpur)
                row.update(best or dict(tid=None, n_hits=0, recall=0.0, purity=np.nan,
                                        rel_recall=np.nan, rel_purity=np.nan, assembled=False))
                rows.append(row)
    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/rescore_cluster_truth.csv", index=False)
    # quiet-day FA at the 8.0 growth floor (Merlin F4)
    Q = pd.read_parquet(f"{FRZ}/run_given0_quiethits.parquet")
    nq = len(pd.read_csv(f"{FRZ}/quiet_days.csv"))
    ntpl = len(pd.read_csv(f"{FRZ}/run_given0_thresholds.csv"))
    fa = float((Q.net_mad >= 8.0).sum()) / (ntpl * nq)
    summ = dict(label=LABEL, cluster_km=CLUSTER_KM, dedup_s=DEDUP_S,
                frozen_verdict="arm (i) 0/18 assembled -- unchanged by this re-score",
                counts={f"{t}/{s}": int(g.assembled.sum())
                        for (t, s), g in R.groupby(["truth", "source"])},
                n_eligible=len(keep),
                cluster_sizes={int(X): len(cl[X]) for X in keep},
                quiet_day_FA_at_T8=dict(n_hits=int((Q.net_mad >= 8.0).sum()),
                                        template_days=int(ntpl * nq), fa_per_template_day=fa,
                                        allowance=V.FA_PER_TPL_DAY,
                                        note="Merlin F4: the 8.0 growth floor was never exercised "
                                             "by the FA-only rule -- the measured rate is ~8x "
                                             "below the allowance, so T was the recording floor, "
                                             "not a measured noise limit"),
                accept="<= 6/18 confirms closure; >= 12/18 is a failure signature",
                failure_signature_fired=bool(max(int(g.assembled.sum())
                                                 for _, g in R.groupby(["truth", "source"])) >= 12))
    json.dump(summ, open(f"{OUT}/rescore_summary.json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))
    sv = set(pd.read_csv(f"{FRZ}/run_given0_survival.csv").query("survives").tid)
    fam5 = sorted({int(str(t)[1:].split("_")[0]) for t in sv})
    cols = ["truth", "source", "bkey", "n_cluster", "n_truth", "tid", "n_hits", "recall",
            "need_recall", "purity", "need_purity", "rel_purity", "assembled"]
    print("\n--- the 5 frozen survivors' families ---")
    print(R[(R.bkey.isin(fam5)) & (R.source == "frozen")][cols].to_string(index=False))


if __name__ == "__main__":
    main()
