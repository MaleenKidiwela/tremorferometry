#!/usr/bin/env python
"""RE-FIRE the UNCHANGED 2010 candidate templates on the sealed 2011 window.  Prefix `m_`.

`k_grow.matched_filter` is called verbatim: the network mean of the per-station normalised 3-C
cc, threshold = median + 8 x MAD of THAT DAY's own function, 4 s declustering, 12 s templates.
The template waveforms are read straight out of the 2010 run's `grow/tpl_it0.npz` -- no
re-stacking, no re-cutting, no re-thresholding, no re-tuning of any kind.  A template station
that does not exist in 2011 is simply dropped (k_grow already skips stations absent from the day),
and the per-template count of surviving stations is recorded.

  python scripts/m_refire.py --tpl results/l_wa2010/grow/tpl_it0.npz \
      --grps <csv of candidate grps>|<comma list> --dates <comma list> --out <dir> --tag real_tremor
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
import m_common as M                                                        # noqa: E402
import k_grow as KG                                                         # noqa: E402


def load_templates(path, grps=None):
    """{g: {sta: (tpl (3,n) float32, off_s, peak_mad)}} exactly as k_grow saved it."""
    z = np.load(path, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    TPL = {}
    for k, v in meta.items():
        g, s = k.split("_", 1)
        g = int(g)
        if grps is not None and g not in grps:
            continue
        TPL.setdefault(g, {})[s] = (np.asarray(z[k], np.float32), float(v[0]), float(v[1]))
    return TPL


def read_grps(spec):
    if not spec:
        return None
    if os.path.exists(spec):
        V = pd.read_csv(spec)
        if "candidate" in V.columns:
            V = V[V.candidate.astype(bool)]
        return sorted(int(x) for x in V.grp)
    return sorted(int(x) for x in spec.split(","))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--grps", default="")
    ap.add_argument("--dates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    K = M.activate()
    stas = M.stations()
    grps = read_grps(a.grps)
    TPL = load_templates(a.tpl, set(grps) if grps else None)
    if grps and sorted(TPL) != sorted(grps):
        miss = sorted(set(grps) - set(TPL))
        M.log(f"[{a.tag}] WARNING {len(miss)} requested groups have no template: {miss[:20]}")
    dates = [d for d in a.dates.split(",") if d]

    # per-template station bookkeeping: how much of the 2010 template exists in 2011
    keep = set(stas)
    rows = []
    for g, T in TPL.items():
        s2010 = sorted(T)
        s2011 = sorted(s for s in s2010 if s in keep)
        rows.append(dict(grp=g, n_sta_tpl_2010=len(s2010), n_sta_tpl_2011=len(s2011),
                         frac_sta_2011=len(s2011) / max(len(s2010), 1),
                         lost_2011=",".join(s for s in s2010 if s not in keep)))
    ST = pd.DataFrame(rows).sort_values("grp")
    ST.to_csv(f"{a.out}/tpl_stations_{a.tag}.csv", index=False)
    M.log(f"[{a.tag}] {len(TPL)} templates, {len(stas)} stations in the 2011 bed; "
          f"template stations surviving to 2011: median "
          f"{ST.n_sta_tpl_2011.median():.0f}/{ST.n_sta_tpl_2010.median():.0f} "
          f"(min frac {ST.frac_sta_2011.min():.2f})")
    M.log(f"[{a.tag}] {len(dates)} dates: {dates[0]} .. {dates[-1]}")

    dc = K.DayCache(stas)
    D = KG.matched_filter(dc, stas, TPL, dates)
    D.to_parquet(f"{a.out}/det_{a.tag}.parquet", index=False)
    n = D.groupby("grp").size() if len(D) else pd.Series(dtype=int)
    summ = dict(tag=a.tag, n_templates=len(TPL), n_stations=len(stas), n_dates=len(dates),
                n_det=int(len(D)),
                n_templates_with_det=int(D.grp.nunique()) if len(D) else 0,
                det_per_template_median=float(n.median()) if len(n) else 0.0,
                det_per_template_day=float(len(D) / max(len(TPL) * len(dates), 1)))
    json.dump(summ, open(f"{a.out}/det_{a.tag}_summary.json", "w"), indent=1)
    M.log(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
