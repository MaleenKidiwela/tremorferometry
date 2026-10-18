#!/usr/bin/env python
"""S13 Amendment 5: WA pre-registered column grid c6 (lag -14..+2 s, own scrambled-control z) per episode and zone, plus
the zone-1 / pooled grand-stack lag scans (W1 failure branch). Reads tables/v2_wa_units_all.csv, v2_wa_episodes.csv,
v2_wa_{zone1,pooled}_grand_lagscan.csv -> figures/s13_v2_wa_c6.png. No error bars; colour = z."""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
U = pd.read_csv(f'{C.TB}/v2_wa_units_all.csv'); E = pd.read_csv(f'{C.TB}/v2_wa_episodes.csv').set_index('episode')
mid = pd.to_datetime(E.start) + (pd.to_datetime(E.end) - pd.to_datetime(E.start)) / 2
fig, ax = plt.subplots(3, 1, figsize=(10, 9))
P = U[U.unit == 'pooled'].set_index('episode'); m = mid.reindex(P.index)
ax[0].plot(m, P.thrust_c6_m0, color='#9a9a9a', lw=0.8)
sc = ax[0].scatter(m, P.thrust_c6_m0, c=P.thrust_c6_z.clip(-2, 8), cmap='viridis', vmin=-2, vmax=8, s=12 + 60 * np.sqrt(P.n / P.n.max()), zorder=3)
ax[0].axhline(0, color='#9a9a9a', lw=0.6); ax[0].set_ylabel('per LFE, c6 grid (N m)')
ax[0].set_title('S13 WA-NW, column grid c6 (lag −14…+2 s, i.e. origin ≈ time − 5.5 s ± 8 s): pooled per episode (size ∝ √n)', fontsize=9)
for i, z in enumerate(('zone1', 'zone2', 'zone3')):
    Z = U[U.unit == z].set_index('episode')
    ax[1].scatter(mid.reindex(Z.index) + pd.Timedelta(days=20 * (i - 1)), Z.thrust_c6_m0, c=Z.thrust_c6_z.clip(-2, 8), cmap='viridis',
                  vmin=-2, vmax=8, s=18, marker='osD'[i], label=z, zorder=3)
ax[1].axhline(0, color='#9a9a9a', lw=0.6); ax[1].set_ylabel('per LFE by zone, c6 (N m)'); ax[1].legend(fontsize=7, title='o SE … D NW', title_fontsize=7)
for f, c in (('zone1', 'C0'), ('pooled', 'C1')):
    G = pd.read_csv(f'{C.TB}/v2_wa_{f}_grand_lagscan.csv')
    ax[2].plot(G.lag_s, G.m0_per_lfe, color=c, label=f'{f} grand stack (all episodes, joint fit)')
ax[2].axhline(0, color='#9a9a9a', lw=0.6); ax[2].axvspan(-8, 8, color='#eeeeee', zorder=0, label='headline ±8 s grid')
ax[2].set_xlabel('lag relative to origin = time + 0.5 s (s)'); ax[2].set_ylabel('M0 per LFE at that lag (N m)'); ax[2].legend(fontsize=7)
fig.colorbar(sc, ax=ax[:2], label='z vs scrambled control (c6 grid)')
fig.savefig(f'{C.FIG}/s13_v2_wa_c6.png', dpi=130); plt.close(fig)
