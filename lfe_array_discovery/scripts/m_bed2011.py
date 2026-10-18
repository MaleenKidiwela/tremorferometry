#!/usr/bin/env python
"""The 2011 bed: station set, tremor hours/days, quiet days.  Prefix `m_`.

Writes (under results/m_2011/):
  bed.json            the whole bed summary (stations kept/dropped and WHY, hours, days)
  tremor_dates.txt    comma-separated tremor days (the re-fire days)
  quiet_dates.txt     comma-separated strictly-quiet days (the false-alarm reference)
  nontremor_dates.txt comma-separated non-tremor days (the wider, weaker FA reference)
  quiet_days.csv      every bed day with its distance to the nearest tremor day
  stations.csv        the per-station admission/health/availability ledger
  candidate_hours.csv the 2011 PNSN candidate hours

Self-check: the identical code on the 2010 window must reproduce the run's 454 hours / 29 days.
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m_common as M                                                        # noqa: E402


def main():
    out = M.OUT
    os.makedirs(out, exist_ok=True)

    H10 = M.candidate_hours("2010-07-20", "2010-09-10")
    M.log(f"[self-check] 2010 window with this code: {len(H10)} candidate hours on "
          f"{H10.date.nunique()} days (the 2010 run recorded 454 / 29)")

    b = M.bed_summary()
    b["selfcheck_2010_hours"] = int(len(H10))
    b["selfcheck_2010_days"] = int(H10.date.nunique())
    H = M.candidate_hours()
    H.to_csv(f"{out}/candidate_hours.csv", index=False)
    Q = M.quiet_frame()
    Q.to_csv(f"{out}/quiet_days.csv", index=False)
    M.station_frame().to_csv(f"{out}/stations.csv", index=False)
    open(f"{out}/tremor_dates.txt", "w").write(",".join(M.tremor_dates()))
    open(f"{out}/quiet_dates.txt", "w").write(",".join(M.quiet_dates()))
    open(f"{out}/nontremor_dates.txt", "w").write(",".join(M.nontremor_dates()))

    # cached station-days actually on disk for the 2011 window
    K = M.activate()
    idx = K.cache_index() if not K._IDX else K._IDX
    days = set(M.bed_days())
    n = 0
    per_day = {}
    for (sta, yr, doy) in idx:
        if yr != M.YEAR:
            continue
        d = (pd.Timestamp(year=yr, month=1, day=1) + pd.Timedelta(days=doy - 1)).strftime(
            "%Y-%m-%d")
        if d in days and sta in set(M.stations()):
            n += 1
            per_day[d] = per_day.get(d, 0) + 1
    b["n_cached_station_days"] = n
    b["cached_stations_per_day_median"] = (int(pd.Series(per_day).median())
                                           if per_day else 0)
    b["cached_stations_on_tremor_days_median"] = (
        int(pd.Series({d: per_day.get(d, 0) for d in M.tremor_dates()}).median())
        if per_day else 0)
    json.dump(b, open(f"{out}/bed.json", "w"), indent=1)
    M.log(json.dumps({k: v for k, v in b.items()
                      if k not in ("stations", "tremor_dates", "quiet_dates")}, indent=1))


if __name__ == "__main__":
    main()
