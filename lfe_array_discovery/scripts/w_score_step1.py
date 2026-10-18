#!/usr/bin/env python
"""Stage 8 -- STEP-1 SCORING against the 10 ISOLATED hidden Bostock families.

RELATIVE clauses (`q_rediscover`, Merlin flaw #5): a hidden family counts as REDISCOVERED when a
cluster is within 5 km of it AND its one-to-one recall (+- 2 s) is >= 0.5 x the FROZEN TEMPLATE
detector's own recall for that family AND its purity >= 0.8 x the template's own purity AND it is
not > 50 % explained by a kept family.  The reference is
`lfe_picker/results/p3/truth_template_reference.csv` (frozen template, 2012, same window, SNB and
B926 withheld).

PURITY IS ALSO REPORTED AGAINST A <= 3 km CLUSTER TRUTH (memory `scoring-truth-colocated-families`:
per-family one-to-one purity is not a noise measure when families sit 0.1-1.6 km apart).  The 10
isolated families have no neighbour within 5 km, so for them the two truths coincide -- which is
exactly why the isolated stratum is the stratum Step 1 is scored on.

The hidden families are NEVER seeded, NEVER stacked and NEVER matched anywhere else in the chain.

Accept (pre-registered): >= 3/10 isolated rediscovered, spurious-split rate <= 0.1, whole-chain
site-shift null <= 1 gated family per realisation.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402

PICKREF = "/home/jovyan/tremorferometry/lfe_picker/results/p3/truth_template_reference.csv"


def cluster_truth(families, events, hidden, km=W.CLUSTER_TRUTH_KM, tol=2.0):
    """{fid: merged event times of every family within `km` of fid} (2-s dedup)."""
    import s_common as S
    out = {}
    for f, (la, lo) in families.items():
        t = [events.get(f, np.zeros(0))]
        for g, (lb, ob) in families.items():
            if g == f:
                continue
            if float(S.epi_km(la, lo, lb, ob)) <= km:
                t.append(events.get(g, np.zeros(0)))
        v = np.sort(np.concatenate([x for x in t if len(x)])) if any(len(x) for x in t) \
            else np.zeros(0)
        if len(v):
            keep = np.concatenate([[True], np.diff(v) > tol])
            v = v[keep]
        out[f] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--gates", default="")
    ap.add_argument("--nulls", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--year", type=int, default=0)
    ap.add_argument("--windows", default="")
    a = ap.parse_args()
    year = a.year or W.STEP1_YEAR
    sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
    import q_rediscover as Q
    import g_rediscover as GR
    import p_common as P
    import s_common as S

    F = P.families()
    families = {int(b): (float(F.lat[b]), float(F.lon[b])) for b in F.index}
    B = P.bostock()
    B = B[B.year == year]
    events = {int(b): np.sort(g.t_origin.to_numpy(float)) for b, g in B.groupby("bkey")}
    hid, iso, _ = S.holdouts()
    hidden = sorted(hid)
    isolated = sorted(iso)
    assert isolated == W.ISOLATED_HIDDEN, (isolated, W.ISOLATED_HIDDEN)
    ref = pd.read_csv(PICKREF).set_index("bkey")
    wtxt = a.windows or f"{W.STEP1_DAY0},{W.STEP1_DAY1}"
    windows = [(pd.Timestamp(wtxt.split(",")[0]).timestamp(),
                pd.Timestamp(wtxt.split(",")[1]).timestamp() + 86400)]

    C = (pd.read_parquet(a.clusters) if os.path.exists(a.clusters)
         else pd.DataFrame(columns=["cid", "t", "lat", "lon"]))
    clusters = {c: dict(times=g.t.to_numpy(float), lat=float(g.lat.iloc[0]),
                        lon=float(g.lon.iloc[0])) for c, g in C.groupby("cid")} if len(C) else {}
    R, SP, summ = Q.score(clusters, events, families, hidden, ref, windows,
                          tol=W.SCORE_TOL_S, r_km=W.SCORE_R_KM,
                          recall_frac=W.SCORE_RECALL_FRAC, purity_frac=W.SCORE_PURITY_FRAC,
                          max_explained=W.SCORE_MAX_EXPLAINED, split_frac=W.SCORE_SPLIT_FRAC)
    RI = R[R.hidden.isin(isolated)] if len(R) else R
    # ---- cluster-truth (<= 3 km) purity, reported alongside
    CT = cluster_truth(families, events, hidden)
    rows = []
    for h in isolated:
        best = None
        for c, v in clusters.items():
            t = v["times"]
            if not len(t):
                continue
            d = float(S.epi_km(v["lat"], v["lon"], *families[h]))
            mr, _ = GR.one2one(events[h], t, W.SCORE_TOL_S)
            mp, _ = GR.one2one(t, events[h], W.SCORE_TOL_S)
            mc, _ = GR.one2one(t, CT[h], W.SCORE_TOL_S)
            r = dict(hidden=h, cluster=c, dist_km=d, n_cluster=len(t),
                     recall=float(mr.mean()) if len(events[h]) else 0.0,
                     purity=float(mp.mean()), purity_cluster3km=float(mc.mean()))
            if best is None or r["recall"] > best["recall"]:
                best = r
        rows.append(best or dict(hidden=h, cluster=None, recall=0.0, purity=np.nan,
                                 purity_cluster3km=np.nan, dist_km=np.nan, n_cluster=0))
    CTR = pd.DataFrame(rows)
    # ---- gate-PASS clusters matching no Bostock family
    novel = np.nan
    if a.gates:
        G = W.read_csv_safe(a.gates, ["pid", "verdict"])
        if len(G):
            pas = G[G.verdict == "PASS"].pid.astype(str).tolist()
            n = 0
            for c in pas:
                v = clusters.get(c)
                if v is None:
                    continue
                dmin = min(float(S.epi_km(v["lat"], v["lon"], *families[f])) for f in families)
                mm = max(GR.one2one(v["times"], events.get(f, np.zeros(0)),
                                    W.SCORE_TOL_S)[0].mean() for f in families) \
                    if len(v["times"]) else 0.0
                n += int(dmin > W.SCORE_R_KM or mm < 0.10)
            novel = int(n)
    nulls = []
    N = W.read_csv_safe(a.nulls, ["realisation", "n_gate_pass"])
    if len(N):
        nulls = [int(x) for x in N.n_gate_pass]
    out = dict(
        year=int(year), windows=wtxt, bed_mode=W.BED_MODE,
        reference_note=("relative clauses use the frozen template detector's 2012 per-family "
                        "recall/purity (truth_template_reference.csv). For year != 2012 no such "
                        "reference exists -- the frozen templates are built on stations that do "
                        "not exist then -- so the 2012 bar is TRANSFERRED unchanged, fixed before "
                        "scoring; the absolute recall/purity are reported alongside in "
                        "*_isolated_detail.csv so the bar can be re-read."
                        if int(year) != 2012 else "reference measured on the scored year"),
        n_clusters=len(clusters),
        n_isolated=len(isolated),
        isolated=isolated,
        n_isolated_rediscovered=int(RI.rediscovered.sum()) if len(RI) else 0,
        isolated_no_reference=int(RI.no_reference.sum()) if len(RI) and "no_reference" in RI else 0,
        n_hidden_rediscovered_all30=int(R.rediscovered.sum()) if len(R) else 0,
        spurious_split_rate=summ.get("spurious_split_rate"),
        n_kept_split=summ.get("n_kept_split"),
        median_centroid_km=float(np.nanmedian(CTR.dist_km)) if len(CTR) else None,
        median_recall=float(np.nanmedian(CTR.recall)) if len(CTR) else None,
        median_purity=float(np.nanmedian(CTR.purity)) if len(CTR) else None,
        median_purity_cluster3km=float(np.nanmedian(CTR.purity_cluster3km)) if len(CTR) else None,
        gate_pass_no_bostock=novel,
        null_gated_per_realisation=nulls,
        null_mean=float(np.mean(nulls)) if nulls else None,
        ceiling_note="frozen template detector scores 9/30 absolute, 5/10 isolated",
        accept=dict(rediscovered_ge=W.ACCEPT_REDISCOVERED, split_le=W.ACCEPT_SPLIT_RATE,
                    null_gated_le=W.ACCEPT_NULL_GATED),
    )
    nr = out["n_isolated_rediscovered"]
    ss = out["spurious_split_rate"]
    out["verdict"] = (
        "PASS" if (nr >= W.ACCEPT_REDISCOVERED
                   and (ss is None or not np.isfinite(ss) or ss <= W.ACCEPT_SPLIT_RATE)
                   and (not nulls or max(nulls) <= W.ACCEPT_NULL_GATED))
        else ("FA_RULE_BROKEN" if nulls and max(nulls) >= 3
              else ("ENGINE_INADEQUATE" if nr < 2 else "BELOW_BAR")))
    if len(R):
        R.to_csv(a.out + "_rediscovery.csv", index=False)
    if len(SP):
        SP.to_csv(a.out + "_splits.csv", index=False)
    CTR.to_csv(a.out + "_isolated_detail.csv", index=False)
    json.dump(out, open(a.out + ".json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
