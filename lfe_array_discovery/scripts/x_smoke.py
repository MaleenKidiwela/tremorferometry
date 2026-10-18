#!/usr/bin/env python
"""SMOKE ONLY -- exercise the location + Gate 0-2 wiring on a FAKE proposal.

This is a plumbing test, not science: it plants no signal and it proves nothing about the data.
It builds a proposal of `--n` random times on two tremor days at a random active cell and runs the
Washington gate context end to end, so that the grid build, the station table, the corrections, the
DayStore, the HF veto, the Gate-1 nulls at the withheld stations and the Gate-2 location + controls
are all shown to run before the queue is launched.  A FAIL or UNTESTABLE verdict on random times is
the EXPECTED outcome.

    python scripts/x_smoke.py fake-proposal --n 40 --out results/step5_smoke/fake
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x_common as X                                                        # noqa: E402


def fake_proposal(n, out, seed=20260920):
    W = X.activate()
    import g_gates as GG
    cells = pd.read_csv(f"{X.OUT}/active_cells.csv")
    rng = np.random.default_rng(seed)
    c = cells.iloc[int(rng.integers(0, len(cells)))]
    H = W.candidate_hours()
    days = list(H.groupby("date").size().sort_values(ascending=False).head(2).index)
    hh = H[H.date.isin(days)]
    t = np.sort(np.array([pd.Timestamp(r.date).timestamp() + int(r.hour) * 3600
                          + float(rng.uniform(60, 3400))
                          for r in hh.sample(n=min(n, len(hh)), replace=len(hh) < n,
                                             random_state=seed).itertuples()]))
    acc = pd.DataFrame([dict(bkey=0, lat=float(c.lat), lon=float(c.lon), dep=np.nan),
                        dict(bkey=-1, lat=0.0, lon=0.0, dep=np.nan)])
    ctx = X.gate_context(acc, {}, {})
    p = GG.Proposal(pid="FAKE", times=t, lat=float(c.lat), lon=float(c.lon), dep=np.nan,
                    offset_search_s=W.GATES_OFFSET_SEARCH_S, meta=dict(drop_bkey=0))
    days_arr = np.array([np.datetime64(d, "D") for d in W.availability().date])
    t0 = time.time()
    res = GG.evaluate(ctx, p, gates=W.GATES, evaluated_days=days_arr)
    dt = time.time() - t0
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out + ".json", "w"), indent=1, default=str)
    flat = GG.flatten(res)
    pd.DataFrame([flat]).to_csv(out + ".csv", index=False)
    X.log(f"[smoke] fake proposal N={len(t)} at {c.lat:.3f}/{c.lon:.3f} ({c.dep:.1f} km) on "
          f"{days} -> {res['verdict']} in {dt:.1f} s")
    for k in ("gate0", "gate1", "gate2"):
        g = res.get(k, {})
        X.log(f"[smoke]   {k}: {g.get('verdict')} {str(g.get('reason', ''))[:110]}")
    g1 = res.get("gate1", {})
    for q in g1.get("per_station", []):
        X.log(f"[smoke]     g1 {q['sta']}: {q['verdict']} n_used={q.get('n_used')} "
              f"excess={q.get('excess_db')}")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["fake-proposal"])
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--out", default=f"{X.RES}/step5_smoke/fake")
    a = ap.parse_args()
    fake_proposal(a.n, a.out)


if __name__ == "__main__":
    main()
