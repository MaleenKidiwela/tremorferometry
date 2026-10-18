#!/usr/bin/env python
"""S13 pilot figures + pass-rule verdict (PLAN.md §3).
Writes figures/s13_pilot_<ep>_pooled.png, figures/s13_pilot_m0.png, tables/pilot_verdict.json."""
import os, sys, json, pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G

EPS = ['2011-08', '2012-08']
FAMS = ['fam56', 'fam10', 'fam43', 'fam0', 'fam3']
INK, ACC, GREY, SHIFT = '#1f1f1f', '#c0392b', '#9a9a9a', '#2c6fbb'


def stack_fig(ep, unit='pooled'):
    S = pickle.load(open(f'{C.DATA}/stacks/pilot_{ep}.pkl', 'rb'))
    U = pd.read_csv(f'{C.TB}/pilot_{ep}_units.csv').set_index('unit')
    st = S['units'][unit]; t = S['trel']; L = st['lag']; m0 = U.loc[unit, 'thrust_m0']
    co = C.station_coords(); src = (U.loc[unit, 'src_lat'], U.loc[unit, 'src_lon'])
    stas = sorted(st['stk'], key=lambda s: np.hypot(*G.ned_offset(src[0], src[1], co.loc[s, 'lat'], co.loc[s, 'lon'])))
    fig, ax = plt.subplots(len(stas), 3, figsize=(11, 1.55 * len(stas) + 1), sharex=True)
    for i, s in enumerate(stas):
        dkm = np.hypot(*G.ned_offset(src[0], src[1], co.loc[s, 'lat'], co.loc[s, 'lon'])) / 1e3
        for c, cn in enumerate('NEZ'):
            a = ax[i, c]
            a.plot(t, st['c_lo'][s][c] * 1e9, color=GREY, lw=0.6)
            a.plot(t, st['c_hi'][s][c] * 1e9, color=GREY, lw=0.6)
            a.plot(t, st['shp'][s][c] * 1e9, color=SHIFT, lw=0.7, ls='--')
            a.plot(t, st['shm'][s][c] * 1e9, color=SHIFT, lw=0.7, ls=':')
            a.plot(t, np.roll(st['tpl'][s][c], L) * m0 * 1e9, color=ACC, lw=1.2)
            a.plot(t, st['stk'][s][c] * 1e9, color=INK, lw=1.4)
            a.axvline(0, color=GREY, lw=0.5)
            a.set_xlim(-150, 250)
            if c == 0:
                a.set_ylabel(f'{s.split(".")[1]}\n{dkm:.0f} km\nnm', fontsize=8)
            if i == 0:
                a.set_title(cn, fontsize=10)
            a.tick_params(labelsize=7)
    for c in range(3):
        ax[-1, c].set_xlabel('s after LFE origin')
    fig.suptitle(f'S13 {ep} {unit}: 20–50 s displacement stack, n = {U.loc[unit, "n"]} LFEs (P ≥ 0.75). '
                 f'\nBlack = stack; red = thrust template × M̂0 = {m0:.2g} N m per LFE (z = {U.loc[unit, "thrust_z"]:.1f});\n'
                 f'grey = scrambled control 2.5/97.5 %; blue = ±1 day', fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f'{C.FIG}/s13_pilot_{ep}_{unit}.png', dpi=130); plt.close(fig)


def m0_fig():
    D = pd.concat([pd.read_csv(f'{C.TB}/pilot_{ep}_units.csv') for ep in EPS])
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    units = ['pooled', 'pooled_odd', 'pooled_even'] + FAMS
    cm = plt.get_cmap('viridis')
    for k, u in enumerate(units):
        d = D[D.unit == u].sort_values('episode')
        x = np.arange(len(d)) + (k - len(units) / 2) * 0.04
        ax.plot(x, d.thrust_m0, color=GREY, lw=0.8)
        sc = ax.scatter(x, d.thrust_m0, c=d.thrust_z.clip(0, 8), cmap=cm, vmin=0, vmax=8, s=40, zorder=3)
        ax.annotate(u, (x[-1] + 0.03, d.thrust_m0.iloc[-1]), fontsize=7, va='center')
        ax.plot(x, d.thrust_ctrl_p99, color=GREY, lw=0.6, ls=':')
    ax.axhline(0, color=GREY, lw=0.6)
    ax.set_xticks(range(len(EPS))); ax.set_xticklabels(EPS)
    ax.set_ylabel('time-locked 20–50 s moment per LFE (N m)')
    ax.set_title('S13 pilot: thrust-template moment per LFE (colour = z vs scrambled control; dotted = control 99 %)',
                 fontsize=8)
    fig.colorbar(sc, ax=ax, label='z')
    fig.tight_layout(); fig.savefig(f'{C.FIG}/s13_pilot_m0.png', dpi=130); plt.close(fig)


def verdict():
    V = {}
    D = {ep: pd.read_csv(f'{C.TB}/pilot_{ep}_units.csv').set_index('unit') for ep in EPS}
    p1 = {}; p2 = {}; p3 = {}
    for ep, U in D.items():
        r = U.loc['pooled']
        p1[ep] = bool(r.thrust_m0 > r.thrust_ctrl_p99 and r.thrust_m0 > max(r.thrust_m0_shift_p1, r.thrust_m0_shift_m1))
        p2[ep] = bool(r.thrust_m0 > 0 and r.thrust_vr >= 0.2 and r.thrust_vr > r.strikeslip_vr and r.thrust_vr > r.steepdip_vr
                      and r.thrust_sign_agree_big >= 0.7)
        p3[ep] = bool(all(U.loc[h].thrust_m0 > 0 and U.loc[h].thrust_z >= 2 for h in ('pooled_odd', 'pooled_even')))
    V['P1_detection'] = p1; V['P2_thrust'] = p2; V['P3_halves'] = p3
    fam_ok = [f for f in FAMS if any(D[ep].loc[f].thrust_z >= 3 for ep in EPS if f in D[ep].index)]
    V['P4_families_z3'] = fam_ok; V['P4'] = len(fam_ok) >= 2
    good = [ep for ep in EPS if p1[ep] and p2[ep]]
    V['pilot_works'] = bool(any(p1[ep] and p2[ep] and p3[ep] for ep in EPS))
    V['episodes_P1P2'] = good
    for ep, U in D.items():
        r = U.loc['pooled']
        V[f'{ep}_pooled'] = dict(n=int(r.n), m0=float(r.thrust_m0), Mw=float((np.log10(max(r.thrust_m0, 1)) - 9.1) / 1.5),
                                 z=float(r.thrust_z), p=float(r.thrust_p), vr=float(r.thrust_vr),
                                 min_detectable_m0=float(3 * r.thrust_ctrl_sd))
    json.dump(V, open(f'{C.TB}/pilot_verdict.json', 'w'), indent=1)
    print(json.dumps(V, indent=1))


if __name__ == '__main__':
    for ep in EPS:
        stack_fig(ep, 'pooled')
        for f in FAMS:
            stack_fig(ep, f)
    m0_fig()
    verdict()
