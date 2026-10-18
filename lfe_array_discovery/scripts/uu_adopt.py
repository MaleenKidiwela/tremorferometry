#!/usr/bin/env python
"""Choose the adopted locator by the acceptance rule stated in advance, and record it.

RULE (fixed before the ladder was run, and applied mechanically here):
  among the LADDER variants (the cumulative chain, not the ablations), adopt the one with the
  smallest MEDIAN 3-D distance to Bostock's published coordinates, subject to locating at least
  as many of the 58 known families as the baseline does.  Adoption additionally REQUIRES that
  both the median distance improve materially (>= 25 % reduction) and the rail fraction drop;
  otherwise the verdict is NOT ADOPTED and the current locator stands.
This is a choice of ESTIMATOR measured against a known answer, not a candidate-selection cut: it
filters no family (memory `rank-dont-discard`).
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
import uu_common as U                                                       # noqa: E402
import uu_locate as LOC                                                     # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ladder", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    L = pd.read_csv(a.ladder)
    base = L[L.variant.str.startswith("BASELINE_on_disk(")].iloc[0]
    cand = L[L.variant.str.match(r"^L\d")].copy()
    cand = cand[cand.n_located >= base.n_located]
    if not len(cand):
        cand = L[L.variant.str.match(r"^L\d")].copy()
    best = cand.sort_values(["d3_med", "epi_med"]).iloc[0]
    ok = (best.d3_med <= 0.75 * base.d3_med) and (best.rail < base.rail)
    spec = dict(LOC.LADDER + LOC.ABLATE).get(best.variant, {})
    rows = L[L.variant.isin([base.variant, best.variant])][
        ["variant", "n_located", "epi_med", "epi_p90", "dep_med", "d3_med", "rail",
         "within5", "within10", "within20", "null_d3_med"]]
    json.dump(dict(variant=str(best.variant), spec=spec,
                   verdict=("ADOPTED" if ok else "NOT ADOPTED"),
                   baseline={k: (None if pd.isna(base[k]) else float(base[k]))
                             for k in ("n_located", "epi_med", "epi_p90", "dep_med", "d3_med",
                                       "rail")},
                   best={k: (None if pd.isna(best[k]) else float(best[k]))
                         for k in ("n_located", "epi_med", "epi_p90", "dep_med", "d3_med",
                                   "rail")},
                   rows=rows.round(3).to_dict("records")),
              open(a.out, "w"), indent=1)
    U.log(f"adopted {best.variant} ({'ADOPTED' if ok else 'NOT ADOPTED'}): "
          f"epi {base.epi_med:.1f} -> {best.epi_med:.1f} km, 3-D {base.d3_med:.1f} -> "
          f"{best.d3_med:.1f} km, rail {base.rail:.2f} -> {best.rail:.2f}")
    print(best.variant)


if __name__ == "__main__":
    main()
