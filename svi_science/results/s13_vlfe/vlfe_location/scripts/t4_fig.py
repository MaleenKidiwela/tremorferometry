#!/usr/bin/env python
"""Test 4 raw figure: distance to the LFE centroid per GF variant (ECDF) and the relocated points of the 2011-13 events."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = pd.read_csv(f'{OUT}/tables/t4_reloc.csv')
g = R[R.year.between(2011, 2013) & (R.period == 'ETS') & np.isfinite(R.c10_lat)]
V = [('prod', 'production (catalogue)', '#0b0b0b'), ('ws', '(0) whole-space, my VR', '#2a78d6'), ('ws_st', '(c0) whole-space + stack terms', '#1baf7a'),
     ('ok', '(b) Okada half-space static', '#eb6834'), ('ok_st', '(c) Okada + stack terms', '#e87ba4'), ('wsst', 'whole-space static ×2 (diag.)', '#8a8984')]
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for k, lab, c in V:
    d = np.sort(g[f'{k}_dist'].dropna().values)
    ax[0].step(d, np.arange(1, len(d) + 1) / len(d), where='post', color=c, lw=2 if k in ('prod', 'ws_st') else 1.3,
               label=f'{lab}: median {np.median(d):.1f} km')
ax[0].set_xlabel('distance to LFE centroid ±10 min (km)'); ax[0].set_ylabel('fraction of events'); ax[0].set_xlim(0, 80)
ax[0].legend(fontsize=7.5, frameon=False, loc='lower right'); ax[0].set_title(f'2011–2013 ETS tier VLFEs with LFEs (n={len(g)})', fontsize=9)
a = ax[1]
for k, lab, c in V:
    a.scatter(g[f'{k}_dE'] if k != 'prod' else g.prod_dE, g[f'{k}_dN'] if k != 'prod' else g.prod_dN, s=10, color=c, alpha=0.6, label=lab)
    mx = (g[f'{k}_dE'] if k != 'prod' else g.prod_dE).mean(); my = (g[f'{k}_dN'] if k != 'prod' else g.prod_dN).mean()
    a.plot(mx, my, 'o', color=c, ms=11, mec='white', mew=1.5)
a.axhline(0, color='#8a8984', lw=0.8); a.axvline(0, color='#8a8984', lw=0.8); a.set_aspect('equal'); a.set_xlim(-70, 70); a.set_ylim(-70, 70)
a.set_xlabel('east offset (km)'); a.set_ylabel('north offset (km)'); a.set_title('offsets per variant (big dots = mean vectors)', fontsize=9)
a.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f'{OUT}/figures/t4_gf_variants.png', dpi=140)
