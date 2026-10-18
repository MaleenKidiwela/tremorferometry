#!/usr/bin/env python
"""Stage 6 -- LOCATION (the frozen `w_locate_gates` estimator) and GATES 0-2 on the WASHINGTON bed.

The LOCATION half is `w_locate_gates.run(..., gates=False)` unchanged, running under the WA bed
patch: per-site stacks, S on the 3-C envelope at peak/noise >= 3, `g_locate.locate_interface` on
the 2-D Slab2 grid centred 48.0/-123.2 with half-width 62 km, plus the four controls (200 shuffled
picks, 50-realisation synthetic floor, leave-one-site-out, free-depth sanity) and Merlin's
"located" clause with the WA bar of >= 5 picked sites within 60 km.

The GATES half cannot reuse `w_locate_gates`'s own call, because that builds the stock northern
`g_gates.GateContext` (e5r station roster, the 130 Bostock families as the accepted list, the north
cross-family cc table, the PGC/B011 HF veto).  It is rebuilt here for Washington:

  accepted   the WA candidate list ITSELF (leave-one-out by `drop_bkey`), because no Washington
             family has been accepted yet -- Gate 0's distinctness question in this run is "is this
             candidate distinct from the OTHER candidates of this run", which is exactly the
             adopted registering-site rule applied to the only accepted set that exists.  A
             SENTINEL row at (0, 0) with bkey -1 is appended so that a run with a single candidate
             still evaluates (nothing to be confused with -> distinct by distance), and it can
             never be nearest to a real proposal.
  TPL        the iterated 12-s templates from this run's own bank, re-cut to start at S - 1 s so
             that `g_common.cc_template`'s `pre_tpl=1.0` convention is matched exactly (the bank
             stores them from S - 2 s).
  ccp95      the cross-template cc p95 that `w_dedupe` computed on this bank
             (`<dedupe>_crossfam_p95.csv`) -- the WA analogue of `g_crossfam_cc.csv`.
  HF veto    the WA 100-Hz borehole and the cut derived on WA quiet days by `x_hf_cut.py`.

Gate 3 is NOT run (Merlin F6: untestable inside one ETS window).  Gate 1's positive-control table
is supplied only when `--controls` is given -- it is produced by the Step-5b injections, so the
first pass runs without it and a non-PASS Gate 1 is UNTESTABLE by construction (Merlin F4).

    python scripts/x_locate_gates.py --bank results/step5/bank2.npz \
        --indep results/step5/p2_indep.parquet --dedupe results/step5/dedupe.csv \
        --survival results/step5/p2_survival.csv --out results/step5/lg
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
import x_common as X                                                        # noqa: E402


def build_tpl(tids, stas, Wt, keep):
    """{(bkey, sta): {12.0: (3, n) L2-normalised template starting at S - 1 s}}."""
    import w_common as W
    drop = int(round((W.TPL1_PRE_S - 1.0) * W.FS))       # 2 s -> 1 s
    out = {}
    for bkey, tid in enumerate(tids):
        if tid not in keep:
            continue
        bi = tids.index(tid)
        for si, s in enumerate(stas):
            w = Wt[bi, si]
            if not np.isfinite(w).all():
                continue
            seg = np.asarray(w[:, drop:], float) if drop > 0 else np.asarray(w, float)
            nn = np.linalg.norm(seg, axis=1)
            if np.any(nn <= 0):
                continue
            out[(float(bkey), s)] = {12.0: seg / nn[:, None]}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bank", required=True)
    ap.add_argument("--indep", required=True)
    ap.add_argument("--dedupe", default="")
    ap.add_argument("--survival", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--shift", default="")
    ap.add_argument("--controls", default="")
    ap.add_argument("--no-gates", action="store_true")
    a = ap.parse_args()
    W = X.activate()
    import w_locate_gates as LG
    shift = json.loads(a.shift) if a.shift else None
    os.makedirs(os.path.dirname(a.out), exist_ok=True)

    # ---- 1. location (frozen estimator + four controls), gates deferred
    LG.run(a.bank, a.indep, a.dedupe, a.survival, a.out, gates=False, shift=shift)
    C = pd.read_parquet(a.out + "_clusters.parquet")
    L = pd.read_csv(a.out + "_locations.csv")
    if a.no_gates or not len(C):
        json.dump(dict(n_templates=int(len(L)), n_located=int(L.located.sum()) if "located" in L
                       else 0, n_gate_pass=0, gates_run=False),
                  open(a.out + ".json", "w"), indent=1)
        pd.DataFrame(columns=["pid", "verdict"]).to_csv(a.out + "_gates.csv", index=False)
        X.log("[gates] no located proposals -- nothing to gate")
        return

    # ---- 2. the WA gate context
    import g_gates as GG
    import v_common as V
    tids, stas, Wt, off, meta = V.load_bank(a.bank)
    cids = list(dict.fromkeys(C.cid.astype(str)))
    acc = pd.DataFrame([dict(bkey=i, lat=float(C[C.cid.astype(str) == c].lat.iloc[0]),
                             lon=float(C[C.cid.astype(str) == c].lon.iloc[0]), dep=np.nan)
                        for i, c in enumerate(cids)]
                       + [dict(bkey=-1, lat=0.0, lon=0.0, dep=np.nan)])
    tpl = build_tpl(list(tids), list(stas), Wt, set(cids))
    p95p = (a.dedupe.replace(".csv", "_crossfam_p95.csv") if a.dedupe else "")
    ccp95 = {}
    if p95p and os.path.exists(p95p):
        P95 = pd.read_csv(p95p, index_col=0)
        col = P95.columns[0]
        ccp95 = {str(k): float(v) for k, v in P95[col].items() if np.isfinite(v)}
    ctl = pd.read_csv(a.controls) if (a.controls and os.path.exists(a.controls)) else None
    ctx = X.gate_context(acc, tpl, ccp95, controls=ctl)
    X.log(f"[gates] WA context: {len(acc) - 1} accepted candidates, {len(tpl)} templates, "
          f"{len(ccp95)} cross-family p95 entries, HF veto {ctx.cfg.g0_hf_station} "
          f"{ctx.cfg.g0_hf_veto:.3f}, controls={'yes' if ctl is not None else 'no'}")

    days = np.array([np.datetime64(d, "D") for d in W.availability().date])
    ver = []
    for i, cid in enumerate(cids):
        g = C[C.cid.astype(str) == cid]
        p = GG.Proposal(pid=str(cid), times=g.t.to_numpy(float), lat=float(g.lat.iloc[0]),
                        lon=float(g.lon.iloc[0]), dep=np.nan,
                        offset_search_s=W.GATES_OFFSET_SEARCH_S, meta=dict(drop_bkey=i))
        try:
            res = GG.evaluate(ctx, p, gates=W.GATES, evaluated_days=days)
            row = GG.flatten(res)
            row["g1_sites_pass_detail"] = ";".join(
                f"{r['sta']}:{r['verdict']}:{r.get('excess_db', float('nan')):.2f}"
                for r in res["gate1"]["per_station"]) if "gate1" in res else ""
            ver.append(row)
            json.dump(res, open(f"{a.out}_gate_full_{cid}.json", "w"), indent=1, default=str)
        except Exception as e:                                              # noqa: BLE001
            ver.append(dict(pid=str(cid), verdict="ERROR", error=str(e)[:300]))
        X.log(f"[gates] {cid} -> {ver[-1].get('verdict')}")
    Gv = pd.DataFrame(ver)
    Gv.to_csv(a.out + "_gates.csv", index=False)
    npass = int((Gv.verdict == "PASS").sum()) if len(Gv) else 0
    json.dump(dict(n_templates=int(len(L)), n_located=int(L.located.sum()) if "located" in L else 0,
                   n_gate_pass=npass, gates_run=True,
                   n_gate1_2sites=int((Gv.get("g1_sites_pass", pd.Series(dtype=float))
                                       .fillna(0) >= 2).sum()) if len(Gv) else 0),
              open(a.out + ".json", "w"), indent=1)
    X.log(f"[gates] {npass} PASS of {len(Gv)}")


if __name__ == "__main__":
    main()
