#!/usr/bin/env python
"""S13 pilot cross-check (PLAN §3, column): our SVI LFEs (P>=0.5, counted once) within +-60 s of the Hutchison & Ghosh
2018 Table 1 2011/2014 template VLFE times vs the day's mean rate. Writes tables/hg2018_crosscheck.csv."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
T = [('1', '2011-08-17 04:30:28', 47.825, -123.075), ('2', '2011-08-21 04:09:06', 48.3, -123.6),
     ('3', '2011-08-21 08:43:29', 48.4, -123.45), ('4', '2011-08-24 10:47:38', 48.35, -123.625),
     ('5', '2011-08-26 13:44:34', 48.35, -123.175), ('6', '2014-12-02 04:58:36', 48.825, -123.725),
     ('7', '2014-12-04 21:02:27', 48.975, -123.725), ('8', '2014-12-05 11:56:58', 48.825, -123.775),
     ('9', '2014-12-11 10:29:35', 48.7, -123.75), ('10', '2014-12-11 12:21:20', 48.725, -123.125),
     ('11', '2014-12-11 19:09:29', 48.775, -123.7), ('12', '2014-12-12 00:13:17', 49.1, -124.05)]
R = C.svi_rows(); R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.5)]
rows = []
for k, t, la, lo in T:
    t0 = pd.Timestamp(t).value / 1e9
    day = R[R.date == t[:10]]
    n60 = ((day.t > t0 - 60) & (day.t < t0 + 60)).sum()
    rows.append(dict(template=k, time_utc=t, lat=la, lon=lo, n_lfe_day=len(day), n_within_60s=int(n60),
                     expected_60s=len(day) * 120 / 86400, n_within_600s=int(((day.t > t0 - 600) & (day.t < t0 + 600)).sum()),
                     expected_600s=len(day) * 1200 / 86400))
D = pd.DataFrame(rows); D.to_csv(f'{C.TB}/hg2018_crosscheck.csv', index=False); print(D.round(2).to_string())
