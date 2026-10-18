#!/usr/bin/env python
"""P3D DISCOVERY -- the pre-registered KILL TEST on KEPT families.

Merlin's Arm B (and, read the same way, Arm A) must first be shown to work on families it is GIVEN:
do loud events of a KEPT (non-hidden) family link to other events of that SAME family?

The test reads an arm's own edge list and asks, of the windows that coincide (+-tol s) with a kept
family's hp event, how many of their edges land on another hp event of the SAME family, against the
rate at which they land on a DIFFERENT family's hp event or on no catalogued event at all.  A
construction that cannot link a family it is handed cannot be believed when it links nothing.

  python scripts/s_discover_killtest.py --tag Bseed_true --arm B
"""
from __future__ import annotations
import os, sys, json, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

OUT = S.P3D


def label_windows(tw, H, tol=2.0):
    """-> family id per window (-1 if it coincides with no hp event)."""
    lab = np.full(len(tw), -1, np.int64)
    ht = H.t.to_numpy()
    hb = H.bkey.to_numpy()
    o = np.argsort(ht)
    ht, hb = ht[o], hb[o]
    idx = np.searchsorted(ht, tw)
    for k, i in enumerate(idx):
        best, bb = tol + 1, -1
        for j in (i - 1, i):
            if 0 <= j < len(ht):
                d = abs(ht[j] - tw[k])
                if d < best:
                    best, bb = d, hb[j]
        if bb >= 0 and best <= tol:
            lab[k] = bb
    return lab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--arm", default="B")
    ap.add_argument("--tol", type=float, default=2.0)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    meta = json.load(open(f"{OUT}/disc_{a.tag}.json"))
    Wdf = pd.read_parquet(f"{OUT}/disc_{a.tag}_windows.parquet")
    E = pd.read_parquet(f"{OUT}/disc_{a.tag}_edges.parquet")
    hid, iso, _ = S.holdouts()
    H = S.hp_events(0.8, episodes=meta["episodes"])
    Hk = H[~H.bkey.isin(hid)]                                  # KEPT families only
    lab = label_windows(Wdf.t.to_numpy(), Hk, a.tol)
    n_lab = int((lab >= 0).sum())
    print(f"[kill] {len(Wdf)} windows; {n_lab} coincide (+-{a.tol} s) with a KEPT family's hp event "
          f"({n_lab/max(len(Wdf),1):.4%} of windows), over {Hk.bkey.nunique()} families", flush=True)
    li, lj = lab[E.i.values], lab[E.j.values]
    both = (li >= 0) & (lj >= 0)
    same = both & (li == lj)
    diff = both & (li != lj)
    none = (li < 0) & (lj < 0)
    half = ~both & ~none
    # chance expectation for same-family among labelled-labelled edges
    p = Hk.bkey.value_counts(normalize=True).to_numpy()
    p_same_chance = float((p ** 2).sum())
    out = dict(tag=a.tag, arm=a.arm, n_windows=int(len(Wdf)), n_edges=int(len(E)),
               n_windows_labelled=n_lab,
               frac_windows_labelled=n_lab / max(len(Wdf), 1),
               n_edges_both_labelled=int(both.sum()), n_edges_same_family=int(same.sum()),
               n_edges_diff_family=int(diff.sum()),
               n_edges_one_labelled=int(half.sum()), n_edges_none_labelled=int(none.sum()),
               frac_same_of_both=(float(same.sum() / both.sum()) if both.sum() else float("nan")),
               chance_frac_same=p_same_chance,
               enrichment=(float(same.sum() / both.sum() / p_same_chance)
                           if both.sum() and p_same_chance > 0 else float("nan")),
               n_families_with_a_same_family_edge=int(pd.Series(li[same]).nunique()))
    out["verdict"] = ("CONSTRUCTION LINKS KEPT FAMILIES" if out["n_edges_same_family"] >= 10
                      else "CONSTRUCTION FAILS ITS OWN KILL TEST -- it cannot link families it is given")
    json.dump(out, open(f"{OUT}/disc_{a.tag}_killtest.json", "w"), indent=1, default=str)
    for k, v in out.items():
        print(f"  {k:34s} {v}")


if __name__ == "__main__":
    main()
