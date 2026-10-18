#!/usr/bin/env python
"""Digitise H&G 2019 Figure 6 (daily VLFEs, all five 2011 templates summed, red; BGN noise line, cyan) from the 300-dpi
render of p. 561. Calibration (300-dpi pixels): x 873 = 2011-07-01, 2138 = 2012-10-01 (2.762 px/day); y 927 = 0,
508 = 30 VLFE/day (right axis). The red line is a daily polyline; the centre of the red extent in each pixel column is the
line value there, so the column mean over a period = the time-averaged daily count (approximate, ~0.07/day per px).
Legend box (x 1200-1270, y < 600) excluded. Output tables/paper_fig6_daily_digitised.csv and period means."""
import numpy as np, pandas as pd
from PIL import Image

HG = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro'
im = np.array(Image.open(f'{HG}/data/paper_fig6_p561_300dpi.png').convert('RGB')).astype(int)
r, g, b = im[..., 0], im[..., 1], im[..., 2]
red = (r > 200) & (g < 60) & (b < 60)
red[:600, 1195:1275] = False
cy = (r < 80) & (g > 200) & (b > 200)
X0, X1, Y0, Y30 = 873.0, 2138.0, 927.0, 508.0
ppd = (X1 - X0) / (pd.Timestamp('2012-10-01') - pd.Timestamp('2011-07-01')).days
rows = []
for x in range(875, 2306):
    ys = np.where(red[480:935, x])[0] + 480
    if len(ys) == 0:
        continue
    v = (Y0 - 0.5 * (ys.min() + ys.max())) / (Y0 - Y30) * 30
    rows.append((pd.Timestamp('2011-07-01') + pd.Timedelta(days=(x - X0) / ppd), v))
D = pd.DataFrame(rows, columns=['date', 'vlfe_per_day'])
D.to_csv(f'{HG}/tables/paper_fig6_daily_digitised.csv', index=False)
bgn = (Y0 - np.median(np.where(cy[480:935, 900:2300])[0] + 480)) / (Y0 - Y30) * 30
P = {'pre (07-01..07-23-2011)': ('2011-07-01', '2011-07-23'), 'ETS2011 (07-23..09-06)': ('2011-07-23', '2011-09-06'),
     'interETS (09-06-2011..08-30-2012)': ('2011-09-06', '2012-08-30'),
     'interETS_b (10-01-2011..08-30-2012)': ('2011-10-01', '2012-08-30'),
     'ETS2012 (08-30..10-11)': ('2012-08-30', '2012-10-11'), 'post (10-11..12-01-2012)': ('2012-10-11', '2012-12-01'),
     'all (07-01-2011..12-01-2012)': ('2011-07-01', '2012-12-01'),
     '2011 analysis window (07-01..10-01-2011)': ('2011-07-01', '2011-10-01')}
out = [(k, round(D[(D.date >= a) & (D.date < b)].vlfe_per_day.mean(), 2),
        (pd.Timestamp(b) - pd.Timestamp(a)).days) for k, (a, b) in P.items()]
O = pd.DataFrame(out, columns=['period', 'mean_vlfe_per_day_all5', 'days'])
O['implied_count'] = (O.mean_vlfe_per_day_all5 * O.days).round(0)
O.loc[len(O)] = ['BGN cyan line (noise, all 5 summed)', round(bgn, 2), np.nan, np.nan]
O.to_csv(f'{HG}/tables/paper_fig6_period_means.csv', index=False)
print(O.to_string(index=False))
