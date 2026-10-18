#!/usr/bin/env python
"""S13 write-up figure (before gap fill): per-episode 20-50 s moment per LFE vs time, SVI and WA-NW, points coloured
by z vs scrambled control, size ~ sqrt(n LFEs), no error bars. Units with a flat channel (control SD < 1e10 N m;
2019-04 in both regions) are left out and marked with a grey tick on the time axis.
Output: figures/s13_scale_both_m0_per_lfe.png"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

fig, ax = plt.subplots(2, 1, figsize=(10, 6.8), sharex=True)
for a, reg, name in ((ax[0], 'svi', 'S Vancouver Island'), (ax[1], 'wa', 'Washington NW')):
    e = pd.read_csv(f'{C.TB}/scale_{reg}_episodes.csv')
    e = e[e.built].copy()
    e['mid'] = pd.to_datetime(e.start) + (pd.to_datetime(e.end) - pd.to_datetime(e.start)) / 2
    bad = e.thrust_ctrl_sd < 1e10
    g = e[~bad]
    a.axhline(0, color='#9a9a9a', lw=0.6)
    sc = a.scatter(g.mid, g.thrust_m0 / 1e12, c=g.thrust_z.clip(-3, 8), cmap='viridis', vmin=-3, vmax=8,
                   s=10 + 70 * np.sqrt(g.n / 4000).clip(0, 1.5), edgecolor='k', linewidth=0.3, zorder=3)
    for _, r in e[bad].iterrows():
        a.axvline(r.mid, color='#bbbbbb', lw=0.8, ls=':')
    if reg == 'svi':
        s = g[g.start >= '2008']; w = 1 / s.thrust_ctrl_sd ** 2
        mu = (s.thrust_m0 * w).sum() / w.sum() / 1e12
        a.axhline(mu, color='#555555', lw=0.8, ls='--')
        a.text(pd.Timestamp('2003-03-01'), mu + 0.4, f'2008–2026 weighted mean {mu:.1f}e12', fontsize=7, color='#555555')
    a.set_ylabel('VLF moment per LFE (1e12 N m)')
    a.set_ylim(-10, 14)
    a.set_title(f'{name}: {int((g.thrust_z >= 3).sum())} of {len(g)} episodes at z ≥ 3, '
                f'{int((g.thrust_z <= -3).sum())} at z ≤ −3', fontsize=9, loc='left')
fig.colorbar(sc, ax=ax, label='z vs scrambled control (clipped −3…8)')
fig.suptitle('S13 scale-up (before gap fill): 20–50 s moment time-locked to our LFEs, per major episode\n'
             'point size ∝ √(LFEs stacked); dotted line = episode left out (flat channel)', fontsize=9)
fig.savefig(f'{C.FIG}/s13_scale_both_m0_per_lfe.png', dpi=140)
print('ok')
