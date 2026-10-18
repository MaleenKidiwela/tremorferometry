#!/usr/bin/env python
"""STAGE 2 -- CHOOSE THE SITE TREATMENT BY MEASUREMENT, not by assumption.  Prefix `y_`.

Agent O's protocol (`results/WA_2010_ARRAY_RESULTS.md` sec.0d) applied to the 2008 and 2007
Puget / S-Olympic beds: the IDENTICAL literal autocorrelation (`k_auto`: 6 s windows, 1-8 Hz,
per-component cc summed, median + 8 x MAD of the block's own distribution, |lag| >= 6 s, 4 s
decluster) is run three ways --

  member   every admitted station is one term of the network mean (what agent N and agent P did,
           and the configuration that FAILED on the Olympic bed);
  beam     one zero-slowness delay-and-sum beam per array (what Ghosh / Chestler / Sweet did);
  siteavg  every member keeps its own correlation but the network function is the mean over
           SITES of the per-site mean over channels -- one vote per site, member information kept

-- on the SAME ~12 tremor hours and ~12 quiet hours, restricted to the member set present in BOTH
populations, and the winner is the one with the highest tremor/quiet pair-count ratio.  Nothing
else is compared and no threshold is moved.

  stage `keep`  writes the common member set and the two hour lists
  stage `pick`  reads the three _stats.csv and writes mode_smoke.json / mode_chosen.txt
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
import y_common as Y                                                        # noqa: E402


def stage_keep(a):
    T = Y.candidate_hours()
    Q = Y.quiet_hours()
    nt = min(a.hours, len(T))
    nq = min(a.hours, len(Q))
    T = T.iloc[sorted(pd.Series(range(len(T))).sample(nt, random_state=0))]
    Q = Q.iloc[sorted(pd.Series(range(len(Q))).sample(nq, random_state=0))]
    idx = Y.member_cache_index()

    def have(dates):
        s = None
        for d in dates:
            ts = pd.Timestamp(d)
            k = {m for (m, y, doy) in idx if (y, doy) == (ts.year, ts.dayofyear)}
            s = k if s is None else (s | k)
        return s or set()

    keep = have(sorted(set(T.date))) & have(sorted(set(Q.date)))
    os.makedirs(a.out, exist_ok=True)
    open(f"{a.out}/mode_keep.txt", "w").write(",".join(sorted(keep)))
    open(f"{a.out}/mode_hours_tremor.txt", "w").write(
        ";".join(f"{r.date}:{int(r.hour)}" for r in T.itertuples()))
    open(f"{a.out}/mode_hours_quiet.txt", "w").write(
        ";".join(f"{r.date}:{int(r.hour)}" for r in Q.itertuples()))
    print(f"bed {Y.BED}: {len(keep)} member stations common to both populations; "
          f"{len(T)} tremor hours, {len(Q)} quiet hours", flush=True)


def stage_pick(a):
    rows = []
    for m in ("member", "beam", "siteavg"):
        r = {}
        for w in ("tremor", "quiet"):
            d = f"{a.out}/mode/{m}_{w}"
            p = f"{d}/_stats.csv"
            if not os.path.exists(p):
                continue
            S = pd.read_csv(p)
            ok = S[S.nmad_in_sd > 3]
            r[f"n_blocks_{w}"] = int(len(ok))
            r[f"det_per_hour_{w}"] = float(ok.n_det.sum() / max(len(ok), 1))
            r[f"traces_{w}"] = float(ok.n_sta.median()) if len(ok) else np.nan
            r[f"thr_sd_{w}"] = float(ok.nmad_in_sd.median()) if len(ok) else np.nan
        if r.get("det_per_hour_quiet", 0) > 0:
            r["tremor_quiet_ratio"] = r["det_per_hour_tremor"] / r["det_per_hour_quiet"]
        else:
            r["tremor_quiet_ratio"] = (float("inf") if r.get("det_per_hour_tremor", 0) > 0
                                       else 0.0)
        rows.append(dict(mode=m, **r))
    T = pd.DataFrame(rows)
    best = max(rows, key=lambda x: (x["tremor_quiet_ratio"]
                                    if np.isfinite(x["tremor_quiet_ratio"]) else -1))
    why = (f"highest tremor/quiet pair-count ratio at 8 x MAD on {a.hours} + {a.hours} hours "
           f"restricted to the member set present in both populations, bed {Y.BED}; the 2010 "
           f"N-Olympic values were member 1.385 / beam 0.991 / siteavg 15.81 "
           f"(WA_2010_ARRAY_RESULTS sec.0d), and the member-wise whole-bed 2008 value on record "
           f"is 0.89 (WA_2008_RESULTS sec.3)")
    json.dump(dict(bed=Y.BED, table=rows, chosen=best["mode"], why=why),
              open(f"{a.out}/mode_smoke.json", "w"), indent=1, default=float)
    open(f"{a.out}/mode_chosen.txt", "w").write(best["mode"])
    print(T.to_string(), flush=True)
    print("CHOSEN:", best["mode"], flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["keep", "pick"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--hours", type=int, default=12)
    a = ap.parse_args()
    {"keep": stage_keep, "pick": stage_pick}[a.stage](a)


if __name__ == "__main__":
    main()
