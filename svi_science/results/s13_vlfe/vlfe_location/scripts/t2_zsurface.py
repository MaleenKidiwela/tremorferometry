#!/usr/bin/env python
"""Test 2 (PLAN.md): production scan z over all grid nodes at each event's best time (cat_inj.local_z, run real_v2),
Δz <= 1 / <= 2 region size; also dumps each event's production data windows (for Test 4).
usage: t2_zsurface.py YEAR   (sets the production scan environment itself; read-only on ../catalogue and ../data)
-> data/zsurf_YEAR.npz, data/win_YEAR.npz, tables/t2_zsurface_YEAR.csv"""
import os, sys
Y = int(sys.argv[1])
S13 = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe'
SD = f'{S13}/data/scan' if Y == 2011 else f'{S13}/data/scan_{Y}'
os.environ.update(OMP_NUM_THREADS='1', SCAN_YEAR=str(Y), SCAN_GRID=f'{SD}/templates_grid_v2.npz', SCAN_GAIN_FALLBACK='1',
                  SCAN_RUN='real_v2', SCAN_BP=f'{SD}/bp_{Y}v2')
os.environ.pop('SCAN_REFBP', None)
if Y in (2011, 2012):
    os.environ['SCAN_REFBP'] = f'{SD}/bp2_ref_v2'
sys.path.insert(0, f'{S13}/catalogue/scripts')
import numpy as np, pandas as pd
import cat_common as K
import cat_inj as I

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.makedirs(f'{OUT}/data', exist_ok=True)


def main():
    E = pd.read_csv(f'{OUT}/tables/events_final.csv')
    E = E[E.year == Y]
    I.bat_init()
    X, SC, S = I._B['X'], I._B['SC'], I._B['S']
    W = SC.W
    zmax, amax = X.EV['zmax'], X.EV['argmax']
    y0 = X.EV['y0']
    N = W['nodes'].iloc[:-1]
    nx, ny = S.kmxy(N.lat.values, N.lon.values)
    rows, Z, WIN = [], [], []
    for r in E.itertuples():
        tsec = int(round(r.epoch - y0))
        lo, hi = max(tsec - 30, 0), min(tsec + 31, len(zmax))
        zw = np.asarray(zmax[lo:hi], np.float64)
        k = int(np.nanargmax(zw)); ts = lo + k; node = int(amax[ts]); zb = float(zw[k])
        z = I.local_z(ts, half=0)[0]                                   # (nodes-1,)
        kb = int(np.nanargmax(z))
        d = np.hypot(nx - nx[kb], ny - ny[kb])
        o = dict(eid=r.eid, year=Y, tsec_best=ts, node_stored=node, z_stored=zb, node_local=kb, z_local=float(z[kb]),
                 z_at_stored_node=float(z[node]) if node < len(z) else np.nan,
                 lat_local=N.lat.values[kb], lon_local=N.lon.values[kb])
        for dz in (1, 2):
            m = z >= z[kb] - dz
            o[f'n_dz{dz}'] = int(m.sum()); o[f'area_dz{dz}_km2'] = 100.0 * m.sum()
            o[f'req_dz{dz}_km'] = float(np.sqrt(100.0 * m.sum() / np.pi)); o[f'maxd_dz{dz}_km'] = float(d[m].max())
            # extent along east and north
            o[f'ext_E_dz{dz}_km'] = float(nx[m].max() - nx[m].min()); o[f'ext_N_dz{dz}_km'] = float(ny[m].max() - ny[m].min())
        if np.isfinite(r.c10_lat):
            cx, cy = S.kmxy(r.c10_lat, r.c10_lon)
            kc = int(np.argmin(np.hypot(nx - cx, ny - cy)))
            o['node_lfe'] = kc; o['z_at_lfe_node'] = float(z[kc]); o['dz_lfe_node'] = float(z[kb] - z[kc])
            o['lfe_node_in_dz1'] = bool(z[kc] >= z[kb] - 1); o['lfe_node_in_dz2'] = bool(z[kc] >= z[kb] - 2)
        rows.append(o); Z.append(z)
        # data windows for Test 4: production cleaned windows, FIT0-15 .. FIT1+15, all stations
        d_, u = divmod(ts, 86400)
        srcday = (d_ + W['offs']) % S.NDAY
        xs, sig, won, g, sv = [], [], [], [], []
        for i, s in enumerate(W['stas']):
            tsrc = int(srcday[i]) * 86400 + u
            xs.append(SC.cleaned_window(i, s, tsrc, ts, S.FIT0 - 15, S.FIT1 + 1 + 15))
            sig.append(W['DS'][s]['sig'][srcday[i]] * 1e9)
            won.append(W['won'][d_, i]); g.append(W['g'][d_, i]); sv.append(W['signs'][0][i])
        WIN.append(dict(x=np.array(xs), sig=np.array(sig), won=np.array(won), g=np.array(g), sv=np.array(sv)))
    T = pd.DataFrame(rows)
    T.to_csv(f'{OUT}/tables/t2_zsurface_{Y}.csv', index=False)
    np.savez_compressed(f'{OUT}/data/zsurf_{Y}.npz', eid=T.eid.values, z=np.array(Z), lat=N.lat.values, lon=N.lon.values,
                        dep=N.dep.values)
    np.savez_compressed(f'{OUT}/data/win_{Y}.npz', eid=T.eid.values, tsec=T.tsec_best.values, stas=np.array(W['stas']),
                        sll=W['sll'], x=np.array([w['x'] for w in WIN]), sig=np.array([w['sig'] for w in WIN]),
                        won=np.array([w['won'] for w in WIN]), g=np.array([w['g'] for w in WIN]),
                        sv=np.array([w['sv'] for w in WIN]), T_thrust=W['T']['thrust'][:-1])
    print(T[['eid', 'node_stored', 'node_local', 'z_stored', 'z_local', 'n_dz1', 'req_dz1_km', 'maxd_dz1_km']].round(2).to_string())


if __name__ == '__main__':
    main()
