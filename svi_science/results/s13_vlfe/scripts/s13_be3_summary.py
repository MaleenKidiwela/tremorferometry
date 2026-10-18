#!/usr/bin/env python
"""S13 Amendment 7 summary: between-ETS stacks at candidate times. Reads tables/be3_{hi,c030,c010,pnsn}_units.csv and
data/stacks/be3_*.pkl. Writes tables/be3_summary.csv, tables/be3_verdict.json, figures/s13_be3_classes.png (rho per unit,
colour = z, no error bars), figures/s13_be3_stacks_<unit>.png (raw stacks of the main units)."""
import os, sys, json, pickle
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G

ETS = 5.39e12                    # SVI 2008-26 A2 weighted mean per LFE (Amendment 5 results)
LFE_M0 = 4.5e11                  # Bostock mag_med 1.7 (stated, not measured)
TAGS = ['be3_hi', 'be3_c030', 'be3_c010', 'be3_pnsn']
MAIN = ['C075', 'C050', 'C030', 'C010', 'PNSN_noLFE', 'PNSN_noCand', 'BE_broad']
INK, ACC, GREY, SHIFT = '#1f1f1f', '#c0392b', '#9a9a9a', '#2c6fbb'


def word(r, s):
    if not (np.isfinite(r) and np.isfinite(s)):
        return 'n/a'
    if abs(r - 1) < 2 * s and s <= 0.35:
        return 'ETS-like'
    if r + 2 * s < 1:
        return 'below ETS'
    if r - 2 * s > 1:
        return 'above ETS'
    return 'undetermined'


def main():
    U = pd.concat([pd.read_csv(f'{C.TB}/{t}_units.csv').assign(tag=t) for t in TAGS if os.path.isfile(f'{C.TB}/{t}_units.csv')],
                  ignore_index=True)
    Pm = {}
    for t in TAGS:
        f = f'{C.DATA}/rows_be3/{t}.parquet'
        if os.path.isfile(f):
            R = pd.read_parquet(f); Pm.update(R.groupby('unit').P.mean().to_dict())
    out = []
    for r in U.itertuples():
        sd = r.thrust_ctrl_sd
        zp1 = (r.thrust_m0_shift_p1 - r.thrust_ctrl_mean) / sd; zm1 = (r.thrust_m0_shift_m1 - r.thrust_ctrl_mean) / sd
        rho, srho = r.thrust_m0 / ETS, sd / ETS
        out.append(dict(unit=r.unit, tag=r.tag, n=r.n, mean_P=Pm.get(r.unit, np.nan), date0=r.date0, date1=r.date1, n_sta=r.n_sta,
                        stations=r.stations, m0_per_cand=r.thrust_m0, z=r.thrust_z, p=r.thrust_p, ctrl_sd=sd, lag=r.thrust_lag,
                        vr_thrust=r.thrust_vr, vr_free=getattr(r, 'vr_free', np.nan), vr_ss=getattr(r, 'strikeslip_vr', np.nan),
                        vr_steep=getattr(r, 'steepdip_vr', np.nan), m0_day_p1=r.thrust_m0_shift_p1, m0_day_m1=r.thrust_m0_shift_m1,
                        z_day_p1=zp1, z_day_m1=zm1, wide_lag=getattr(r, 'thrust_wide_lag', np.nan),
                        wide_z=getattr(r, 'thrust_wide_z', np.nan), m0_Pw=r.thrust_m0_Pw,
                        present=bool(r.thrust_z >= 3 and abs(zp1) < 3 and abs(zm1) < 3),
                        rho=rho, sigma_rho=srho, R_vlf_lfe=12.0 * rho, R_2sig_lo=12.0 * (rho - 2 * srho), R_2sig_hi=12.0 * (rho + 2 * srho),
                        wording=word(rho, srho), max_share=r.max_share, max_share_ch=r.max_share_ch,
                        loso_min=r.loso_m0_min, loso_min_without=r.loso_min_without, n_nb10=r.n_nb10, nb_kern=r.nb_kern,
                        guard_dropped=r.guard_dropped))
    D = pd.DataFrame(out)
    D.to_csv(f'{C.TB}/be3_summary.csv', index=False)
    with pd.option_context('display.width', 250):
        print(D[['unit', 'n', 'mean_P', 'm0_per_cand', 'z', 'z_day_p1', 'z_day_m1', 'present', 'rho', 'sigma_rho', 'R_vlf_lfe',
                 'wording', 'lag', 'vr_thrust', 'n_sta']].round(3).to_string())
    V = {r.unit: dict(n=int(r.n), m0=float(r.m0_per_cand), z=float(r.z), present=bool(r.present), rho=float(r.rho),
                      sigma_rho=float(r.sigma_rho), R=float(r.R_vlf_lfe), wording=r.wording, mean_P=float(r.mean_P))
         for r in D.itertuples()}
    json.dump(dict(ETS_ref=ETS, LFE_M0=LFE_M0, units=V), open(f'{C.TB}/be3_verdict.json', 'w'), indent=1)
    # class figure
    order = [u for u in MAIN if u in set(D.unit)]
    extra = [u for u in D.unit if u not in order]
    fig, ax = plt.subplots(figsize=(12, 4.8))
    allu = order + sorted(extra)
    x = np.arange(len(allu))
    d = D.set_index('unit').loc[allu]
    sc = ax.scatter(x, d.rho, c=d.z.clip(-3, 8), cmap='viridis', vmin=-3, vmax=8, s=[90 if u in order else 30 for u in allu],
                    edgecolor='k', lw=0.4, zorder=3)
    ax.axhline(1, color=ACC, lw=0.8, ls='--', label='ETS value (5.39e12 N m per LFE; VLF:LFE ≈ 12)')
    ax.axhline(0, color=GREY, lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels([f'{u}\nn={int(n)}' for u, n in zip(allu, d.n)], rotation=60, fontsize=7)
    ax.set_ylabel('VLF moment per candidate / ETS value (ρ)')
    ax.set_title('Between-ETS 20–50 s stacks at candidate times (SVI, 2003–2026). Colour = z vs scrambled control; no error bars.\n'
                 'Big dots = main units; small = era (e1 2003–08, e2 2009–13, e3 2014–26) and tremor-day splits', fontsize=8)
    ax2 = ax.secondary_yaxis('right', functions=(lambda r: 12 * r, lambda R: R / 12)); ax2.set_ylabel('VLF:LFE moment ratio (×4.5e11)')
    fig.colorbar(sc, ax=ax, label='z', pad=0.08); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/s13_be3_classes.png', dpi=130); plt.close(fig)
    # raw stacks of the main units
    co = C.station_coords()
    for t in TAGS:
        f = f'{C.DATA}/stacks/{t}.pkl'
        if not os.path.isfile(f):
            continue
        S = pickle.load(open(f, 'rb'))
        Ut = U[U.tag == t].set_index('unit')
        for unit in [u for u in MAIN if u in S['units']]:
            st = S['units'][unit]; tt = S['trel']; L = st['lag']; m0 = Ut.loc[unit, 'thrust_m0']
            src = (Ut.loc[unit, 'src_lat'], Ut.loc[unit, 'src_lon'])
            stas = sorted(st['stk'], key=lambda s: np.hypot(*G.ned_offset(src[0], src[1], co.loc[s, 'lat'], co.loc[s, 'lon'])))
            fig, ax = plt.subplots(len(stas), 3, figsize=(11, 1.55 * len(stas) + 1.2), sharex=True, squeeze=False)
            for i, s in enumerate(stas):
                dkm = np.hypot(*G.ned_offset(src[0], src[1], co.loc[s, 'lat'], co.loc[s, 'lon'])) / 1e3
                for c, cn in enumerate('NEZ'):
                    a = ax[i, c]
                    a.plot(tt, st['c_lo'][s][c] * 1e9, color=GREY, lw=0.6); a.plot(tt, st['c_hi'][s][c] * 1e9, color=GREY, lw=0.6)
                    a.plot(tt, st['shp'][s][c] * 1e9, color=SHIFT, lw=0.7, ls='--'); a.plot(tt, st['shm'][s][c] * 1e9, color=SHIFT, lw=0.7, ls=':')
                    a.plot(tt, np.roll(st['tpl'][s][c], L) * m0 * 1e9, color=ACC, lw=1.2)
                    a.plot(tt, st['stk'][s][c] * 1e9, color=INK, lw=1.4)
                    a.axvline(0, color=GREY, lw=0.5); a.set_xlim(-150, 250)
                    if c == 0:
                        a.set_ylabel(f'{s.split(".")[1]}\n{dkm:.0f} km\nnm', fontsize=8)
                    if i == 0:
                        a.set_title(cn, fontsize=10)
                    a.tick_params(labelsize=7)
            for c in range(3):
                ax[-1, c].set_xlabel('s after candidate origin')
            fig.suptitle(f'S13 between-ETS {unit}: 20–50 s stack, n = {int(Ut.loc[unit, "n"])}. Black = stack; red = thrust × '
                         f'M̂0 = {m0:.2g} N m per candidate (z = {Ut.loc[unit, "thrust_z"]:.1f}); grey = scrambled control '
                         f'2.5/97.5 %; blue = ±1 day', fontsize=8)
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fig.savefig(f'{C.FIG}/s13_be3_stacks_{unit}.png', dpi=110); plt.close(fig)


if __name__ == '__main__':
    main()
