#!/usr/bin/env python
"""Re-screen the Washington 2008 (agent N) and 2010 (agent L) templates with the RECOMMENDED
screen measured on the 2005 POLARIS bed, beside the screen as currently coded.  Prefix `q_sv_`.

RECOMMENDED (results/SCREEN_VALIDATION_NORTH.md sec.5):
  tier A (gates)  located on >=5 stations; |dz| <= DZ_MAX of Slab2; >= MIN_DAYS distinct tremor
                  days; quiet-day FA <= FA_MAX per template-day with the tremor-day rate >= 3x it
  tier B          + the shuffled-station control
  tier C          + tremor co-timing above the same-day base rate
The 8-km Slab2 width, the interior clause, the 0.1/day false-alarm cap and the second-half
stability clause were measured to reject more than 15 % of the 58 rediscovered Bostock families
and are dropped or reformulated; the control and the co-timing test are demoted to tiers.

  python scripts/q_sv_wa.py --dz-max 19 --min-days 3 --fa-max 0.25
"""
from __future__ import annotations
import argparse, json, os
import numpy as np, pandas as pd

OUT = "results/q_screenval"


def tiers(V, a, has_fa):
    V = V.copy()
    V["t_c1"] = V.c_c1.astype(bool)
    V["t_c4"] = np.array([bool(np.isfinite(x) and abs(x) <= a.dz_max) for x in V.dz]) & V.t_c1
    V["t_d"] = V.n_days >= a.min_days
    if has_fa:
        V["t_a"] = ((V.fa_quiet_per_tpl_day <= a.fa_max)
                    & (V.rate_tremor_per_tpl_day >= 3 * V.fa_quiet_per_tpl_day) & (V.n_det > 0))
    else:
        V["t_a"] = True
    V["tierA"] = V.t_c1 & V.t_c4 & V.t_d & V.t_a
    V["tierB"] = V.tierA & V.c_c2.astype(bool)
    V["tierC"] = V.tierB & V.c_b.astype(bool)
    return V


def wa2008(a):
    res = {}
    for tag, root in (("real", "results/n_puget2008/2008"),
                      ("null", "results/n_puget2008/2008_null")):
        V = pd.read_csv(f"{root}/screen_final.csv")
        V["current"] = (V.c_a & V.c_b & V.c_c1 & V.c_c2 & V.c_c3 & V.c_c4 & V.c_d & V.c_e)
        V = tiers(V, a, has_fa=True)
        V.to_csv(f"{OUT}/wa2008_{tag}_rescreen.csv", index=False)
        res[tag] = dict(n=int(len(V)),
                        current={c: int(V[c].sum()) for c in
                                 ["c_a", "c_b", "c_c1", "c_c2", "c_c3", "c_c4", "c_d", "c_e",
                                  "current"]},
                        recommended={c: int(V[c].sum()) for c in
                                     ["t_a", "t_c1", "t_c4", "t_d", "tierA", "tierB", "tierC"]})
    return res


def wa2010(a):
    res = {}
    for tag, root in (("real", "results/l_wa2010"), ("null", "results/l_wa2010_null")):
        V = pd.read_csv(f"{root}/verify_it0.csv")
        V["c_c1"] = V.located.astype(bool)
        V["c_c2"] = V.ctrl_pass.fillna(False).astype(bool)
        V["c_c4"] = [bool(np.isfinite(x) and abs(x) <= 8.0) for x in V.dz]
        V["c_d"] = V.n_days >= 3
        V["c_b"] = V.co_timing > V.co_base
        V["current"] = V.c_c1 & V.c_c2 & V.c_c4 & V.c_d & V.c_b
        V = tiers(V, a, has_fa=False)                 # no quiet-day re-fire exists for this bed
        V.to_csv(f"{OUT}/wa2010_{tag}_rescreen.csv", index=False)
        res[tag] = dict(n=int(len(V)),
                        current={c: int(V[c].sum()) for c in
                                 ["c_c1", "c_c2", "c_c4", "c_d", "c_b", "current"]},
                        recommended={c: int(V[c].sum()) for c in
                                     ["t_c1", "t_c4", "t_d", "tierA", "tierB", "tierC"]})
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dz-max", type=float, default=19.0)
    ap.add_argument("--min-days", type=int, default=3)
    ap.add_argument("--fa-max", type=float, default=0.25)
    ap.add_argument("--out", default=f"{OUT}/wa_rescreen.json")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    R = dict(dz_max=a.dz_max, min_days=a.min_days, fa_max=a.fa_max,
             wa2008=wa2008(a), wa2010=wa2010(a))
    json.dump(R, open(a.out, "w"), indent=1)
    print(json.dumps(R, indent=1))
