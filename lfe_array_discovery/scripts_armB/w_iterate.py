#!/usr/bin/env python
"""Stage 4 -- ITERATED STACKS (2-3 passes; 12-s templates from S - 2 s after the first re-pick).

Thin wrapper over `v_iterate.iterate` with `w_common.patch_v` installed, so the iteration uses the
data-driven sites, MIN_SITES 4 and the valid-only station_z.  Nothing in lfe_picker is edited.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bank", required=True)
    ap.add_argument("--indep", required=True)
    ap.add_argument("--survivors", required=True)
    ap.add_argument("--out-bank", required=True)
    ap.add_argument("--out-diag", required=True)
    ap.add_argument("--shift", default="")
    ap.add_argument("--min-indep", type=int, default=W.SURV_MIN_INDEP)
    a = ap.parse_args()
    if os.path.exists(a.out_bank):
        W.log(f"[iterate] {a.out_bank} exists -- skip")
        return
    import json
    shift = json.loads(a.shift) if a.shift else None
    V, E = W.patch_v(shift=shift)
    import v_iterate as I                                                   # noqa: E402
    SV = W.read_csv_safe(a.survivors, W.SURV_COLS)
    if not len(SV) or not SV.survives.astype(bool).any():
        det = W.sites()[3]
        import numpy as np
        V.save_bank(a.out_bank, [], det, np.zeros((0, len(det), 3, 1), np.float32),
                    np.zeros((0, len(det))), dict(kind="iterated_empty"))
        pd.DataFrame().to_csv(a.out_diag, index=False)
        W.log("[iterate] no survivors -- empty bank")
        return
    I.CUT_PRE_S = W.HANDOFF_CUT_PRE_S
    I.CUT_LEN_S = W.HANDOFF_CUT_LEN_S
    I.MAX_STACK = W.MAX_STACK
    I.iterate(a.bank, a.indep, a.survivors, a.out_bank, a.out_diag, use_all=False,
              min_indep=a.min_indep)


if __name__ == "__main__":
    main()
