#!/usr/bin/env python
"""Digitise H&G 2019 Figure 7 (cumulative VLFE function per 2011 template, Jul 2011 - Dec 2012) from the 300-dpi
render of p. 563. Pixel calibration (found by axis/grid-line detection, see NOTES.md 2026-10-15):
  box x: 416 px = 2011-07-01, grid lines 726/1036/1342/1649/1958 = 2011-10-01/2012-01-01/04-01/07-01/10-01
  panels (top=1, bottom=0 on the right axis): rows 391-807, 861-1263, 1317-1748, 1803-2256, 2303-2689
Output: tables/paper_fig7_cumulative.csv (template, date, cumfrac) and tables/paper_fig7_period_fractions.csv.
APPROXIMATE (digitised): ~3.4 px per day, ~0.003 in cumulative fraction per px."""
import numpy as np, pandas as pd
from PIL import Image

HG = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro'
im = np.array(Image.open(f'{HG}/data/paper_fig7_p563_300dpi.png').convert('RGB')).astype(int)
r, g, b = im[..., 0], im[..., 1], im[..., 2]
org = (r > 200) & (g > 60) & (g < 170) & (b < 80)
gx = np.array([416, 726, 1036, 1342, 1649, 1958])
gd = pd.to_datetime(['2011-07-01', '2011-10-01', '2012-01-01', '2012-04-01', '2012-07-01', '2012-10-01'])
gdn = (gd - gd[0]).days.values.astype(float)
a, c = np.polyfit(gx, gdn, 1)
PAN = [(391, 807), (861, 1263), (1317, 1748), (1803, 2256), (2303, 2689)]
rows = []
for k, (y0, y1) in enumerate(PAN, 1):
    prev = y1 - 1.0                      # curve starts at 0 (bottom) on 2011-07-01
    for x in range(418, 2150):           # stop short of the right-axis tick marks
        ys = np.where(org[y0 + 2:y1 - 1, x])[0] + y0 + 2
        ys = ys[np.abs(ys - prev) <= 12]  # track the curve; ignores the legend line and stray orange
        if len(ys) == 0:
            continue
        yv = float(np.median(ys)); prev = yv
        rows.append((k, gd[0] + pd.Timedelta(days=a * x + c), (y1 - yv) / (y1 - y0)))
D = pd.DataFrame(rows, columns=['template', 'date', 'cumfrac'])
# enforce monotone (running max) to remove stray pixels
D['cumfrac_mono'] = D.groupby('template').cumfrac.cummax()
D.to_csv(f'{HG}/tables/paper_fig7_cumulative.csv', index=False)
P = {'pre_2011ETS (07-01..07-23)': ('2011-07-01', '2011-07-23'), 'ETS2011 (07-23..09-06)': ('2011-07-23', '2011-09-06'),
     'interETS_a (09-06-2011..08-30-2012)': ('2011-09-06', '2012-08-30'),
     'interETS_b (10-01-2011..08-30-2012)': ('2011-10-01', '2012-08-30'),
     'ETS2012 (08-30..10-11)': ('2012-08-30', '2012-10-11'), 'post (10-11..12-01)': ('2012-10-11', '2012-12-01')}
out = []
for k, G in D.groupby('template'):
    G = G.set_index('date').cumfrac_mono
    f = lambda t: float(np.interp(pd.Timestamp(t).value, G.index.values.astype('int64'), G.values))
    for name, (t0, t1) in P.items():
        nd = (pd.Timestamp(t1) - pd.Timestamp(t0)).days
        fr = f(t1) - f(t0)
        out.append((k, name, nd, round(fr, 3), round(fr / nd * 520, 2)))
O = pd.DataFrame(out, columns=['template', 'period', 'days', 'frac_of_total', 'rate_rel_to_mean'])
O.to_csv(f'{HG}/tables/paper_fig7_period_fractions.csv', index=False)
print(O.pivot(index='period', columns='template', values='rate_rel_to_mean'))
