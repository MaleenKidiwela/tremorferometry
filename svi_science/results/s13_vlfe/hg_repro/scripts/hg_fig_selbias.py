#!/usr/bin/env python
"""Raw figure for the selection-bias check (PLAN Amendment 1 A6): Fig.-2-style panels for (a) T2 detections, (b) the
detections of noise template j=0 under the T2 configuration, (c) random times, all at TA.B05D BHN, threshold thr_rate079.
Output figures/fig2_selection_bias.png."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import hg_common as H, hg_analyze as A
TH = pd.read_csv(f'{H.TB}/thresholds.csv'); CAT = pd.read_csv(f'{H.TB}/catalogue_prim_thr_rate079.csv')
a18, b18 = H.per_idx('all18')
thr = float(TH[(TH.template == 2) & (TH.variant == 'prim')].thr_rate079.iloc[0])
z = np.load(f'{H.MF}/noise_000.npz', allow_pickle=True)
pi, pv, _ = A.noise_peaks(z, 'T2_prim')
nsel = pi[(pv > thr) & (np.abs(pi - int(z['i0'])) > H.L) & (pi >= a18) & (pi < b18)]
real = CAT[(CAT.template == 2) & (CAT.i >= a18) & (CAT.i < b18)].i.values
rnd = np.sort(np.random.default_rng(5).integers(a18, b18, len(real)))
Y = np.load(f'{H.MF}/cont_TA.B05D.npy', mmap_mode='r')[0]
fig, ax = plt.subplots(2, 3, figsize=(13, 5.5), gridspec_kw={'height_ratios': [1, 2]})
for c, (lab, idx, ref) in enumerate([('T2 detections (real template)', real, int(H.tidx(H.templates().time[1]))),
                                     ('noise template 0 detections, T2 config', nsel, int(z['i0'])),
                                     ('random times', rnd, None)]):
    W = [np.asarray(Y[i - 205:i + 295], float) for i in idx]
    W = np.array([w / np.max(np.abs(w)) for w in W if w.size == 500 and np.all(np.isfinite(w))])
    st = W.mean(0)
    ax[0, c].plot(st / np.max(np.abs(st)), 'b', lw=1, label='stack')
    if ref is not None:
        w = np.asarray(Y[ref - 205:ref + 295], float); ax[0, c].plot(w / np.nanmax(np.abs(w)), 'k', lw=0.6, alpha=0.6, label='selecting window')
    ax[0, c].set_title(f'{lab} (n={len(W)})', fontsize=9); ax[0, c].legend(fontsize=7)
    ax[1, c].imshow(W, aspect='auto', cmap='gray', vmin=-1, vmax=1, extent=[0, 500, len(W), 0]); ax[1, c].set_xlabel('Time (s)')
fig.suptitle('Selection bias: windows picked by correlation with ANY 90-s window stack into that window (TA.B05D BHN)')
fig.tight_layout(); fig.savefig(f'{H.FIG}/fig2_selection_bias.png', dpi=120)
