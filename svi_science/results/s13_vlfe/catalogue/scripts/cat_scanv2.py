#!/usr/bin/env python
"""Scan extension v2 (PLAN_CAT Amendment R): scan station list = v1 (7 + LZB) + every live broadband regional station of the
year with a bp year file. usage (SCAN_YEAR=Y):
  cat_scanv2.py setup Y     -> data/scan_Y/bp_Yv2/ (symlinks), data/scan_Y/templates_grid_v2.npz (v1 grid + extra stations)
  cat_scanv2.py sources     -> adds the extra stations' C0/C1 displacement stacks to data/inj/sources.npz (copy kept)"""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

EXTRA = ['UW.GNW', 'UW.DOSE', 'UW.SP2', 'TA.D03D', 'TA.B05D', 'UW.WISH', 'US.NLWA', 'TA.A04D', 'UW.OTR', 'UW.HDW'] + \
    [x for x in os.environ.get('Y12_EXTRA', '').split(',') if x]          # Amendment Y12: added near stations (env-controlled)


def _node(args):
    import scan_common as S
    import s13_green as G
    i, src, coords, stas = args
    out = {}
    for m, mech in S.MECHS.items():
        M = G.mt_ned(*mech, 1.0)
        out[m] = np.stack([S.template(M, src, coords[s]) for s in stas])
    out['thrust20'] = np.stack([S.template(G.mt_ned(*S.MECHS['thrust'], 1.0), src, coords[s], tau=20.0) for s in stas])
    out['elem'] = np.stack([np.stack([S.template(E, src, coords[s]) for s in stas]) for E in S.ELEM])
    return i, out


def setup(year):
    import scan_common as S
    import scan_refs as SR
    SD = S.SD
    v1 = np.load(f'{SD}/templates_grid.npz')
    base = list(v1['stations'])
    out = f'{SD}/bp_{year}v2'
    os.makedirs(out, exist_ok=True)
    stas = []
    for s in base:
        for ext in ('.npy', '_daystats.npz'):
            src = f'{SD}/bp_{year}g/{s}{ext}'
            if os.path.isfile(src) and not os.path.exists(f'{out}/{s}{ext}'):
                os.symlink(src, f'{out}/{s}{ext}')
        stas.append(s)
    extra = []
    for s in EXTRA:
        RB = os.environ.get('SCAN_REFBP', f'{SD}/bp2_ref')
        src = f'{RB}/{s}.npy'
        if not os.path.isfile(src):
            continue
        X = np.load(src, mmap_mode='r')
        if np.isfinite(X[2, ::3600]).mean() < 180 / 365:            # >= ~180 live days
            continue
        for ext in ('.npy', '_daystats.npz'):
            if not os.path.exists(f'{out}/{s}{ext}'):
                os.symlink(f'{RB}/{s}{ext}', f'{out}/{s}{ext}')
        extra.append(s)
    cc = SR.all_coords()
    coords = {s: (cc.loc[s, 'lat'], cc.loc[s, 'lon']) for s in extra}
    from multiprocessing import Pool
    res = {}
    nodes = list(zip(v1['lat'], v1['lon'], v1['dep']))
    with Pool(int(os.environ.get('CAT_NW', '8'))) as p:
        for i, o in p.imap_unordered(_node, [(i, n, coords, extra) for i, n in enumerate(nodes)]):
            res[i] = o
    A = {}
    for k in ('thrust', 'strikeslip', 'steepdip', 'thrust20'):
        new = np.stack([res[i][k] for i in range(len(nodes))])
        A[k] = np.concatenate([v1[k], new], axis=1)
    newe = np.stack([res[i]['elem'] for i in range(len(nodes))])                 # (nodes, 5, extra, 3, 181)
    A['elem'] = np.concatenate([v1['elem'], newe], axis=2)
    np.savez_compressed(f'{SD}/templates_grid_v2.npz', stations=np.array(base + extra), tau=v1['tau'], lat=v1['lat'],
                        lon=v1['lon'], dep=v1['dep'], **A)
    K.log(f'scan v2 {year}: stations {base + extra}', 'scanv2')


def sources():
    import cat_inj as I
    import shutil
    if not os.path.isfile(f'{I.INJ}/sources_v1.npz'):
        shutil.copy(f'{I.INJ}/sources.npz', f'{I.INJ}/sources_v1.npz')
    z = dict(np.load(f'{I.INJ}/sources_v1.npz'))
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    EP0 = I.epoch0_H()
    for c in (0, 1):
        g = V[V.cluster == c]
        orig = (g.i + g.shift_s).astype(int).values
        for s in EXTRA:
            ws = []
            for o in orig:
                ep = EP0 + o; y = pd.Timestamp(ep, unit='s').year
                f = f'{K.S13}/data/scan/bp2_ref/{s}.npy' if y == 2011 else f'{K.S13}/data/scan_{y}/bp2_ref/{s}.npy'
                if not os.path.isfile(f):
                    continue
                Xb = np.load(f, mmap_mode='r')
                a0 = int(round(ep - pd.Timestamp(f'{y}-01-01').value / 1e9))
                w = np.array(Xb[:, a0 + I.W0:a0 + I.W1], np.float64)
                if np.isfinite(w).all():
                    ws.append(w - w.mean(1, keepdims=True))
            if len(ws) < 2:
                continue
            sig = np.zeros((3, I.W1 - I.W0))
            for k in range(3):
                U = [w[k] / np.linalg.norm(w[k]) for w in ws if np.linalg.norm(w[k]) > 0]
                st = np.mean(U, 0); st /= np.linalg.norm(st)
                sig[k] = np.mean([w[k] @ st for w in ws]) * st
            z[f'C{c}|sc|{s}'] = sig
            K.log(f'sources C{c}: {s} from {len(ws)} members', 'scanv2')
    np.savez(f'{I.INJ}/sources.npz', **z)


if __name__ == '__main__':
    if sys.argv[1] == 'setup':
        setup(int(sys.argv[2]))
    else:
        sources()
