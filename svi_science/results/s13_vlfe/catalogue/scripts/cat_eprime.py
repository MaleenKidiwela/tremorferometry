#!/usr/bin/env python
"""Merlin step 3 / PLAN_CAT §6: validate e' on the R3 sets (between ETS: 130 candidates, 650 scrambled, 1,300
selection-matched noise; ETS: 115, 575, 1,150) with ratio_far as test (a).
-> tables/eprime_validation.csv, tables/eprime_rows.parquet, figures/eprime_*.png"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

R3 = f'{K.HGR}/data/r3'


def load(sfx):
    R = pd.read_parquet(f'{R3}/rows{sfx}.parquet')
    R = R[R.kind != 'pwbase']
    E = pd.concat([pd.read_parquet(f'{R3}/events{sfx}_{y}.parquet') for y in (2011, 2012)], ignore_index=True)
    D = R.merge(E, on='row_id', how='left')
    N = K.node_table()
    ok = D.node_best.notna()
    nb = D.node_best.fillna(0).astype(int)
    D['dist_family_km'] = np.where(ok, N.dist_family_km.values[nb], np.nan)
    D['in_footprint_poly'] = np.where(ok, N.fp_interior.values[nb], False)
    D['head2_maxmag'] = K.head2_maxmag(D.epoch.values, D.year.values)
    D = K.apply_rules(D)
    D['set'] = 'ETS' if sfx else 'between'
    return D


def main():
    out, rows = [], []
    for sfx in ('', '_ets'):
        D = load(sfx)
        rows.append(D)
        for kind, g in D.groupby('kind'):
            out.append(dict(set=D.set.iloc[0], kind=kind, n=len(g), battery_e=int(g.battery.sum()),
                            battery_eprime=int(g.battery_e2.sum()), pass_e=int(g.pass_e.sum()), pass_eprime=int(g.pass_e2.sum()),
                            in_head2=int(g.in_head2.sum()), head2_ge55=int((g.in_head2 & (g.head2_maxmag >= 5.5)).sum()),
                            head2_far_excited2=int((g.in_head2 & (g.n_far_excited >= 2)).sum()),
                            abcf=int((g.pass_a & g.pass_b & g.pass_c & g.pass_f).sum())))
    T = pd.DataFrame(out)
    A = pd.concat(rows, ignore_index=True)
    A.to_parquet(f'{K.TB}/eprime_rows.parquet')
    ets = T[(T.set == 'ETS') & (T.kind == 'cand')].iloc[0]
    btw = T[(T.set == 'between') & (T.kind == 'cand')].iloc[0]
    scr = T[(T.set == 'between') & (T.kind == 'scr')].iloc[0]
    verdict = dict(ets_e=int(ets.battery_e), ets_eprime=int(ets.battery_eprime), between_eprime=int(btw.battery_eprime),
                   scr_eprime=int(scr.battery_eprime), scr_n=int(scr.n))
    verdict['accept'] = bool(verdict['ets_e'] <= verdict['ets_eprime'] <= 14 and verdict['scr_eprime'] <= 2
                             and verdict['between_eprime'] == 0)
    T.to_csv(f'{K.TB}/eprime_validation.csv', index=False)
    pd.Series(verdict).to_json(f'{K.TB}/eprime_verdict.json')
    print(T.to_string(index=False)); print(verdict)
    # events gained / lost
    g = A[(A.kind == 'cand') & (A.battery != A.battery_e2)]
    cols = ['set', 'cand_id', 'epoch', 'battery', 'battery_e2', 'in_head2', 'head2_maxmag', 'n_far_excited', 'ratio_far',
            'z_best', 'eq_names']
    g = g[cols].copy()
    g['time_pt'] = pd.to_datetime(g.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    g.to_csv(f'{K.TB}/eprime_changed_candidates.csv', index=False)
    print(g.to_string(index=False))
    # raw figure: ratio_far vs z_best for every row, coloured by kind; e' pass/fail markers
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, (st, d) in zip(axs, A.groupby('set')):
        for kind, c, s in (('noise', '0.75', 6), ('scr', 'tab:purple', 8), ('cand', 'tab:red', 22)):
            x = d[d.kind == kind]
            ax.scatter(x.z_best, x.ratio_far, s=s, c=c, alpha=0.6, lw=0, label=f'{kind} ({len(x)})')
            y = x[x.battery_e2]
            ax.scatter(y.z_best, y.ratio_far, s=s * 3, facecolors='none', edgecolors='k', lw=0.8)
        ax.axhline(0.551, color='k', lw=0.6); ax.axvline(5, color='k', lw=0.6)
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('scan z_best'); ax.set_title(f'{st}: black ring = battery with e′')
        ax.legend(fontsize=8)
    axs[0].set_ylabel('ratio_far (far / near amplitude)')
    fig.tight_layout(); fig.savefig(f'{K.FIG}/eprime_ratio_vs_z.png', dpi=130)


if __name__ == '__main__':
    main()
