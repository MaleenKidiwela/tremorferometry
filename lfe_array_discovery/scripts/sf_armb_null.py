#!/usr/bin/env python
"""ARM B's missing control: the whole Chestler chain run on a SITE-SHIFT NULL.

Agent O ran Chestler & Creager's own chain (`o_chestler.py`: 18 x MAD beam peaks -> cross-array
boxcar stack -> 15-s 2-8 Hz seeds -> 2nd-root network stack at 45 x MAD -> 150-window templates,
3 iterations -> R_com > 0.5 dedupe) and got 739 families from 197,566 detections -- but never ran
a null, so that count has no reference.  This module runs the IDENTICAL chain, stage for stage,
on data whose inter-site coherence has been destroyed.

The null is agent O's own (`o_common.activate(null=True)`): each site reads its data from another
day of the same partition (tremor days with tremor days, quiet days with quiet days), restamped
with the target day's t0, so the intra-array beam, the noise, the glitches and the geometry are
untouched and only the inter-site coherence is gone.  One thing has to be added here: Arm B also
correlates MEMBER stations (3 per array plus the singles) through `o_beam.load_member`, which the
site-level patch does not cover, so every member is shifted by ITS OWN SITE's offset -- the array
stays internally coherent, exactly as in the real chain.

  O_BED=2010 O_NULL_SHIFT=1 python scripts/sf_armb_null.py --stage seeds|grow|dedupe|fa --out DIR
"""
from __future__ import annotations

import argparse
import os
import sys
import types

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import o_common as O                                                        # noqa: E402
import o_beam as OB                                                         # noqa: E402
import o_chestler as OC                                                     # noqa: E402


def install_null():
    O.activate(null=True)
    off = O._S.get("OFF")
    if not off:
        raise SystemExit("site-shift offsets were not installed")
    site_of = {}
    for r in O.sites().itertuples():
        for m in str(r.members).split(","):
            site_of[m] = r.site
    td, qd = O.tremor_dates(), O.quiet_days()
    part = {d: ("t", i) for i, d in enumerate(td)}
    part.update({d: ("q", i) for i, d in enumerate(qd)})
    lists = {"t": td, "q": qd}
    base = OB.load_member

    def shifted(sta, date):
        s = site_of.get(sta)
        if date not in part or s is None or s not in off:
            return base(sta, date)
        kind, i = part[date]
        L = lists[kind]
        if len(L) < 2:
            return base(sta, date)
        v = base(sta, L[(i + off[s]) % len(L)])
        if v is None:
            return None
        return float(pd.Timestamp(date).timestamp()), v[1]

    OB.load_member = shifted
    O.log(f"[null] member-level site-shift installed for {len(site_of)} members "
          f"({len(off)} site offsets)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["seeds", "grow", "dedupe", "fa"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-seeds", type=int, default=1200)
    ap.add_argument("--iter", type=int, default=OC.ITERS)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if os.environ.get("O_NULL_SHIFT", "") != "1":
        raise SystemExit("refusing to run: O_NULL_SHIFT=1 must be set for the null chain")
    install_null()
    ns = types.SimpleNamespace(out=a.out, max_seeds=a.max_seeds, iter=a.iter, stage=a.stage)
    {"seeds": OC.stage_seeds, "grow": OC.stage_grow, "dedupe": OC.stage_dedupe,
     "fa": OC.stage_fa}[a.stage](ns)


if __name__ == "__main__":
    main()
