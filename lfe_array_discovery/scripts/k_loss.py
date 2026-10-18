#!/usr/bin/env python
"""Per-stage loss accounting for the 100 kept families."""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matched", default=f"{K.OUT}/auto_full_matched.csv")
    ap.add_argument("--purity", default=f"{K.OUT}/gb_purity.csv")
    ap.add_argument("--tplinfo", default="")
    ap.add_argument("--score", required=True)
    ap.add_argument("--times", default=f"{K.OUT}/s1all_times.csv")
    ap.add_argument("--out", default=f"{K.OUT}/loss.csv")
    a = ap.parse_args()
    C = pd.read_csv(a.times)
    fams = sorted(C.bkey.unique())
    M = pd.read_csv(a.matched)
    P = pd.read_csv(a.purity)
    S = pd.read_csv(a.score)
    npair = M.groupby("bkey").size()
    grp = P[P.purity >= 0.5].groupby("bkey").size()
    grp8 = P[P.purity >= 0.8].groupby("bkey").size()
    gr = S[S.frac >= 0.30].groupby("bkey").size()
    best = S.sort_values("n_match", ascending=False).drop_duplicates("bkey").set_index("bkey")
    rows = []
    for b in fams:
        rows.append(dict(bkey=int(b), n_ev=int((C.bkey == b).sum()),
                         pair_detections=int(npair.get(b, 0)),
                         groups_pure50=int(grp.get(b, 0)), groups_pure80=int(grp8.get(b, 0)),
                         grown_families=int(gr.get(b, 0)),
                         best_frac=float(best.frac.get(b, np.nan)),
                         best_recall_60=float(best.recall_60.get(b, np.nan)),
                         rediscovered=bool(gr.get(b, 0) > 0)))
    R = pd.DataFrame(rows)
    R.to_csv(a.out, index=False)
    n = len(R)
    K.log(f"kept families {n}")
    K.log(f"  seen by the pair stage        : {(R.pair_detections > 0).sum()}")
    K.log(f"  >=1 group at purity >= 0.5    : {(R.groups_pure50 > 0).sum()}")
    K.log(f"  >=1 group at purity >= 0.8    : {(R.groups_pure80 > 0).sum()}")
    K.log(f"  REDISCOVERED (Bostock >=30%)  : {(R.rediscovered).sum()}")
    lost = R[(R.groups_pure50 > 0) & (~R.rediscovered)]
    K.log(f"  lost between grouping and growth: {len(lost)} -> {sorted(lost.bkey)[:30]}")
    lost2 = R[(R.pair_detections > 0) & (R.groups_pure50 == 0)]
    K.log(f"  lost at grouping (pairs but no pure group): {len(lost2)}")


if __name__ == "__main__":
    main()
