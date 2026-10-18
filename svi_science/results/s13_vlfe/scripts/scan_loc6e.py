#!/usr/bin/env python
"""VLFE scan 6e (PLAN_SCAN.md §6e): ETS-stack location test. Grid search of the VR-maximising node (thrust, +-8 s lag scan,
S13 stored weights) for the 2012-08 (refit_svi_2012-08) and 2018-06 (scalef_svi_2018-06) pooled stacks.
Writes tables/scan_loc6e.csv (per node VR), tables/scan_loc6e_summary.csv, figures/scan_loc6e.png."""
import os, sys, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G
import scan_common as S
from s13_stack import TREL, FITM, LAGS, fit

CASES = [('2012-08', 'refit_svi_2012-08'), ('2018-06', 'scalef_svi_2018-06')]


def node_vr(args):
    i, src, stk, wts, coords = args
    M = G.mt_ned(*S.MECHS['thrust'], 1.0)
    T = {s: S.template(M, src, coords[s], out_t=TREL) for s in stk}
    m0, L, vr = fit(stk, T, wts, LAGS)
    return i, m0, L, vr


def main():
    N = pd.read_csv(f'{C.TB}/scan_grid_nodes.csv').iloc[:-1]
    cc = C.station_coords()
    rows, summ = [], []
    from multiprocessing import Pool
    for ep, tag in CASES:
        st = pickle.load(open(f'{C.DATA}/stacks/{tag}.pkl', 'rb'))['units']['pooled']
        U = pd.read_csv(f'{C.TB}/{tag}_units.csv'); u = U[U.unit == 'pooled'].iloc[0]
        stk = {s: st['stk'][s] for s in st['stk'] if st['wts'][s].sum() > 0}
        wts = {s: st['wts'][s] for s in stk}
        coords = {s: (cc.loc[s, 'lat'], cc.loc[s, 'lon']) for s in stk}
        with Pool(7) as p:
            res = p.map(node_vr, [(i, (r.lat, r.lon, r.dep), stk, wts, coords) for i, r in N.iterrows()])
        R = pd.DataFrame(res, columns=['node', 'm0', 'lag', 'vr']).assign(episode=ep)
        R = R.join(N[['lat', 'lon', 'dep']], on='node')
        rows.append(R)
        b = R.loc[R.vr.idxmax()]
        dt = S.dist_km(b.lat, b.lon, u.src_lat, u.src_lon); dm = S.dist_km(b.lat, b.lon, 48.37, -123.61)
        near = R[R.vr >= R.vr.max() - 0.02]
        summ.append(dict(episode=ep, best_lat=b.lat, best_lon=b.lon, best_dep=b.dep, best_vr=b.vr, best_m0=b.m0,
                         best_lag=b.lag, centroid_lat=u.src_lat, centroid_lon=u.src_lon, dist_to_centroid_km=dt,
                         dist_to_merlin_pt_km=dm, n_nodes_within_002_vr=len(near),
                         spread_km_within_002=float(S.dist_km(b.lat, b.lon, near.lat.values, near.lon.values).max()),
                         pass_20km=bool(dt <= 20), fail_30km=bool(dt > 30), stations=','.join(stk)))
    R = pd.concat(rows); R.to_csv(f'{C.TB}/scan_loc6e.csv', index=False)
    Sm = pd.DataFrame(summ); Sm.to_csv(f'{C.TB}/scan_loc6e_summary.csv', index=False)
    print(Sm.round(3).to_string())
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for a_, (ep, tag), sm in zip(ax, CASES, summ):
        r = R[R.episode == ep]
        sc = a_.scatter(r.lon, r.lat, c=r.vr, s=40, cmap='viridis', vmin=max(r.vr.min(), 0))
        a_.plot(sm['centroid_lon'], sm['centroid_lat'], 'r*', ms=14, label='LFE count-weighted centroid')
        a_.plot(sm['best_lon'], sm['best_lat'], 'wx', ms=10, mew=2, label='best node')
        for s in sm['stations'].split(','):
            a_.plot(cc.loc[s, 'lon'], cc.loc[s, 'lat'], 'k^', ms=6); a_.annotate(s.split('.')[1], (cc.loc[s, 'lon'], cc.loc[s, 'lat']), fontsize=7)
        a_.set_title(f'{ep} pooled stack: VR of thrust fit per node\nbest {sm["dist_to_centroid_km"]:.0f} km from centroid', fontsize=9)
        a_.set_xlim(-125.2, -122.3); a_.set_ylim(47.6, 49.5); fig.colorbar(sc, ax=a_, label='VR'); a_.legend(fontsize=7, loc='lower left')
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan_loc6e.png', dpi=130)


if __name__ == '__main__':
    main()
