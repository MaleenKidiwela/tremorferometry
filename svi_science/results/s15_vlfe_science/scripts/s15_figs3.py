"""S15 Q4 raw figure: GF-injection recovery per node (2012, 2024; Mw 3.8) on a map with the VLFE concurrent-LFE anchors."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import s15_common as C
from s15_figs import BLUE, ORANGE, GRAY, INK, MUTED
T = pd.read_csv(f'{C.OUT_T}/q4_recovery.csv'); T = T[(T.period == 'all') & (T.Mw == 3.8) & ~T.source.str.endswith('_emp')]
Q1 = pd.read_csv(f'{C.OUT_T}/q1_events.csv'); w = Q1.n_lfe_10m >= 1
F = pd.read_csv(f'{C.RES}/s0/families.csv')
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for a, y in zip(ax, (2012, 2024)):
    a.scatter(F.lon, F.lat, s=6, color='#d6d5cf', label='SVI LFE families')
    a.scatter(C.LON0 + Q1.anchor_x[w] / (C.KM * np.cos(np.radians(C.LAT0))), C.LAT0 + Q1.anchor_y[w] / C.KM, s=5, color=GRAY, alpha=0.6,
              label='VLFE concurrent-LFE anchors')
    d = T[T.year == y]
    sc = a.scatter(d.lon, d.lat, c=d.recovery, cmap='Blues', vmin=0, vmax=1, s=180, edgecolor=INK, linewidth=0.8, zorder=3)
    for r in d.itertuples():
        a.annotate(f'{r.source} {r.recovery:.2f}', (r.lon, r.lat), xytext=(7, 4), textcoords='offset points', fontsize=8, color=INK)
    a.set_title(f'{y}: recovery of GF thrust sources, Mw 3.8', loc='left'); a.set_xlabel('longitude')
    a.set_aspect(1 / np.cos(np.radians(48.3))); a.set_xticks([-124.25, -123.75, -123.25, -122.75])
ax[0].set_ylabel('latitude'); ax[0].legend(frameon=False, fontsize=7.5, loc='lower left')
fig.colorbar(sc, ax=ax, shrink=0.8, label='recovered fraction')
fig.savefig(f'{C.OUT_F}/q4_detection_map.png', dpi=150, bbox_inches='tight'); plt.close(fig)
