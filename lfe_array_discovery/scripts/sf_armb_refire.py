#!/usr/bin/env python
"""Re-fire ARM B templates, unchanged, on another window -- through Arm B's OWN detector.

Arm A's re-fire (`n_refire.py`) uses `k_grow.matched_filter`: the network MEAN of the per-station
cc at 8 x MAD.  Arm B's families were found with Chestler's detector instead -- a 2nd-root network
stack thresholded at 45 x MAD on 2-8 Hz member stations -- so a confirmation of an Arm-B family
has to use that same statistic or it is not the same test.  This module calls `o_chestler.mf_nroot`
verbatim on the sealed window; nothing is re-stacked, re-cut or re-thresholded.

  O_BED=2011 python scripts/sf_armb_refire.py --tpl .../tpl_final.npz --grps ... \
      --dates <comma list> --out DIR --tag real_tremor
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import o_common as O                                                        # noqa: E402
import o_chestler as OC                                                     # noqa: E402
import k_grow as KG                                                         # noqa: E402
import sf_common as SF                                                      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--grps", default="")
    ap.add_argument("--dates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    O.activate()
    grps = None
    if a.grps:
        G = pd.read_csv(a.grps)
        grps = set(int(x) for x in G.grp)
    TPL = SF.load_tpl(a.tpl, grps)
    stas = OC.chestler_sites()
    dates = [d for d in a.dates.split(",") if d]
    SF.log(f"[{a.tag}] {len(TPL)} templates, {len(stas)} member stations in the {O.BED} bed, "
           f"{len(dates)} dates {dates[0]} .. {dates[-1]}")
    rows = []
    for g, T in TPL.items():
        s_src = sorted(T)
        s_new = [s for s in s_src if s in set(stas)]
        rows.append(dict(grp=g, n_sta_tpl_src=len(s_src), n_sta_tpl_new=len(s_new),
                         frac_sta_new=len(s_new) / max(len(s_src), 1)))
    pd.DataFrame(rows).to_csv(f"{a.out}/tpl_stations_{a.tag}.csv", index=False)
    dc = OC.BandCache(stas)
    D = OC.mf_nroot(dc, stas, TPL, dates, tpl_len=KG.TPL_LEN)
    D.to_parquet(f"{a.out}/det_{a.tag}.parquet", index=False)
    n = D.groupby("grp").size() if len(D) else pd.Series(dtype=int)
    summ = dict(tag=a.tag, n_templates=len(TPL), n_stations=len(stas), n_dates=len(dates),
                n_det=int(len(D)),
                n_templates_with_det=int(D.grp.nunique()) if len(D) else 0,
                det_per_template_median=float(n.median()) if len(n) else 0.0,
                det_per_template_day=float(len(D) / max(len(TPL) * len(dates), 1)))
    json.dump(summ, open(f"{a.out}/det_{a.tag}_summary.json", "w"), indent=1)
    print(json.dumps(summ, indent=1), flush=True)


if __name__ == "__main__":
    main()
