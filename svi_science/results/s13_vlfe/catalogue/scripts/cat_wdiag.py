#!/usr/bin/env python
"""Amendment W follow-up diagnostics (descriptive): per-station share of the weighted MF sum at post-W2b control tier
events vs at random times; template weight shares per station; QC zeroing checks. Run with the env of a run
(PROD_YEAR, PROD_QC=1, [PROD_RELAY=1], PILOT_SFX=_v2). -> prints; tables/wdiag_<name><year>.csv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
import cat_pilot as P
from cat_inj import sliding_cc

D = pd.read_parquet(f'{P.PD}/rows_scored{P.SFX}.parquet')
x = D[D.tier & D.source.isin(['A1', 'A2']) & (D.period == 'between')].sort_values('epoch')
B = np.load(f'{P.PD}/bank.npz')
offs, signs = P.draws()
get = P.getter()
rng = np.random.default_rng(7)


def shares(tpl, src, draw, i):
    z = np.load(f'{P.MFD}/{tpl}.npz'); chs = list(z['chans']); wk = z['weights']; off = P.PRE if tpl[0] in 'CR' or tpl[1:].isdigit() and int(tpl[1:]) >= 6 else P.PRE
    p = i - off
    per = {}
    for j, ch in enumerate(chs):
        s = ch.rsplit('.', 1)[0]; si = P.STAS.index(s); k = P.CH.index(ch)
        q = p - (int(offs[draw, si]) * 86400 if src == 'A1' else 0)
        if q < 0:
            q += P.M
        xx = get(ch, q, q + P.L)
        c = sliding_cc(xx, B[tpl][k])[0] if len(xx) == P.L else 0.0
        sg = signs[draw, k] if src == 'A2' else 1.0
        per[s] = per.get(s, 0.0) + wk[j] * sg * c
    return per


rows = []
for r in x.itertuples():
    per = shares(r.template, r.source, int(r.draw), int(r.i))
    tot = sum(per.values())
    top = max(per, key=per.get)
    rows.append(dict(kind='control', day=pd.Timestamp(r.epoch, unit='s').strftime('%Y-%m-%d %H:%M'), template=r.template, src=f'{r.source}_{r.draw}',
                     zmad=r.zmad, top_sta=top, top_share=per[top] / tot if tot else np.nan, **{f's_{k}': v for k, v in per.items()}))
T = pd.DataFrame(rows)
# template weight share per station (C0, C3, T2)
ws = []
for t in ['C0', 'C3', 'T2', 'T4']:
    z = np.load(f'{P.MFD}/{t}.npz'); chs = list(z['chans']); wk = z['weights']
    w = pd.Series(wk, index=[c.rsplit('.', 1)[0] for c in chs]).groupby(level=0).sum()
    ws.append((w / w.sum()).rename(t))
W = pd.concat(ws, axis=1).round(3)
# QC zeroing checks
qc = []
for ch in [f'UW.LRIV.{c}' for c in 'NEZ']:
    if ch in P.CH:
        m = P.qc_mask(ch, P.N)
        if m is not None:
            cc, live = P.ccser(ch, B['C0'][P.CH.index(ch)])
            qc.append(dict(ch=ch, excluded_frac=float(m.mean()), cc_nonzero_in_excluded=int(np.count_nonzero(cc[m[:P.M]])),
                           live_in_excluded=int(live[m[:P.M]].sum())))
nz = []
for nm in [n for n in B.files if n.startswith('N') and not n.startswith('w__')]:
    i0 = int(pd.read_csv(f'{K.TB}/{P.NAME}{P.PILOT}_bank.csv').set_index('name').loc[nm, 'i0'])
    for ch in [f'UW.LRIV.{c}' for c in 'NEZ']:
        if ch in P.CH:
            m = P.qc_mask(ch, P.N)
            if m is not None and m[i0]:
                nz.append((nm, ch, bool(np.isfinite(B[nm][P.CH.index(ch)]).all())))
T.to_csv(f'{K.TB}/wdiag_{P.NAME}{P.PILOT}.csv', index=False)
pd.set_option('display.width', 250)
print(P.NAME, P.PILOT, 'stations', [s.split('.')[1] for s in P.STAS])
print(T[['day', 'template', 'src', 'zmad', 'top_sta', 'top_share']].round(2).to_string(index=False))
print('top station counts', T.top_sta.value_counts().to_dict(), 'median top share', round(float(T.top_share.median()), 2))
print('template weight share per station'); print(W.to_string())
print('QC LRIV checks', qc); print('noise templates with LRIV in excluded month: (name, ch, still finite?)', nz[:6], len(nz))
