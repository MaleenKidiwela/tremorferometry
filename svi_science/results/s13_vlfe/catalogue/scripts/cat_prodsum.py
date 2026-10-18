#!/usr/bin/env python
"""Production 2010-2018 summary (Amendment S): per-year and per-episode tables, ledger, raw figure.
-> tables/prod_years.csv, prod_episodes_all.csv, prod_ledger_all.csv, prod_between_all.csv, figures/prod_episodes.png"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

YEARS = range(2010, 2019)
FINAL = len(sys.argv) > 1 and sys.argv[1] == 'final'          # Amendment U: 2016-2018 from the relay re-run (prodr)
NM = lambda y: 'prodr' if (FINAL and y >= 2016) else 'prod'
OUT = '_final' if FINAL else ''


def main():
    import s13_common as C
    import cat_inj as I
    S_, SM = I.load_sources()
    F = K.families(); F = F[F.region == 'S_Vancouver_Is'].copy()
    F['family_id'] = pd.to_numeric(F.family_id, errors='coerce'); F = F[F.family_id.notna()]
    near = set(F.family_id[(K.dkm(F.lat.values, F.lon.values, SM.loc['C0', 'lat'], SM.loc['C0', 'lon']) <= 10)
                           | (K.dkm(F.lat.values, F.lon.values, SM.loc['C1', 'lat'], SM.loc['C1', 'lon']) <= 10)].astype(int))
    R = C.svi_rows(); R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.75)]
    Rp = R[R.fam.isin(near)]
    years, eps, led, btw = [], [], [], []
    for y in YEARS:
        s = json.load(open(f'{K.TB}/{NM(y)}{y}_summary_v2.json'))
        cp = pd.read_csv(f'{K.TB}/{NM(y)}{y}_R_completeness_v2.csv')
        o = dict(year=y, n_stations=len(s['stations']), stations=' '.join(x.split('.')[1] for x in s['stations']),
                 tier_events_episodes=s['tier_events_episodes'], between_raw=s['between_raw'], between_tierA=s['between_tierA'],
                 ctrl_per_yr=s['controls']['ctrl_rate_between_per_yr'], noise_pair_per_yr=s['controls']['noise_pair_rate_between_per_yr'],
                 Mw_c_ep=s['completeness']['Mw_c'].get('ETS'), Mw_c_btw=s['completeness']['Mw_c'].get('between'),
                 T1_n=s['T1']['n_ETS'], T2b=f"{s['T2']['T2b_obs']} vs {s['T2']['T2b_null_mean']:.1f}", T3=f"{s['T3']['n_pass']}/{s['T3']['n']}",
                 flag_relay=s['flag_relay'], ledger=s['ledger_found'])
        for r in cp[~cp.period.str.contains('\\|')].itertuples():
            k = 'ep' if r.period == 'ETS' else 'btw'
            o[f'{k}_{r.Mw}'] = f'{r.recovered:.2f} ({r.ci_lo:.2f}-{r.ci_hi:.2f}) MF {r.mf_only:.2f} B|MF {r.battery_given_mf:.2f}'
        o['bank'] = '2011+relay2018' if NM(y) == 'prodr' else '2011'
        years.append(o)
        E = pd.read_csv(f'{K.TB}/{NM(y)}{y}_episodes_v2.csv')
        for r in E.itertuples():
            a, b = pd.Timestamp(r.start).value / 1e9, pd.Timestamp(r.end).value / 1e9
            n_patch = int(((Rp.t >= a) & (Rp.t < b)).sum()); n_svi = int(((R.t >= a) & (R.t < b)).sum())
            eps.append(dict(year=y, episode=r.episode, major=r.major, start=r.start, end=r.end, tier=r.tier_events, with_LFEs=r.with_LFEs,
                            sum_M0=r.sum_M0_Nm, Mw_sum=r.Mw_of_sum, ctrl_mean=r.ctrl_mean,
                            compl_34=E.loc[r.Index, 'compl_3.4'], compl_36=E.loc[r.Index, 'compl_3.6'], compl_38=E.loc[r.Index, 'compl_3.8'],
                            n_inj_34=E.loc[r.Index, 'n_inj_3.4'], svi_lfe075=n_svi, patch_lfe075=n_patch))
        f = f'{K.TB}/{NM(y)}{y}_ledger_v2.csv'
        if os.path.getsize(f) > 10:
            led.append(pd.read_csv(f))
        f = f'{K.TB}/{NM(y)}{y}_R_between_v2.csv'
        if os.path.isfile(f) and os.path.getsize(f) > 10:
            b_ = pd.read_csv(f); b_['year'] = y; btw.append(b_)
    Y = pd.DataFrame(years); Y.to_csv(f'{K.TB}/prod_years{OUT}.csv', index=False)
    E = pd.DataFrame(eps)
    E['corrected'] = np.where(E.compl_34 >= 0.3, E.tier / E.compl_34, np.nan)
    E.to_csv(f'{K.TB}/prod_episodes_all{OUT}.csv', index=False)
    L = pd.concat(led, ignore_index=True); L.to_csv(f'{K.TB}/prod_ledger_all{OUT}.csv', index=False)
    B = pd.concat(btw, ignore_index=True) if btw else pd.DataFrame(); B.to_csv(f'{K.TB}/prod_between_all{OUT}.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
    print(Y.drop(columns=['stations']).to_string(index=False))
    print(E.round(3).to_string(index=False))
    print(L.groupby('known').found_tier.agg(['sum', 'count']))
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t = pd.to_datetime(E.start) + (pd.to_datetime(E.end) - pd.to_datetime(E.start)) / 2
    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    sc = axs[0].scatter(t, E.tier, c=E.compl_36, cmap='viridis', vmin=0, vmax=1, s=np.where(E.major, 70, 30), edgecolors='k', lw=0.4)
    axs[0].set_ylabel('tier VLFE events per episode'); axs[0].set_title('Production 2010–2018: SVI episodes (large = major), colour = completeness at Mw 3.6')
    m = E.sum_M0 > 0
    axs[1].scatter(t[m], E.sum_M0[m], c=E.compl_36[m], cmap='viridis', vmin=0, vmax=1, s=np.where(E.major[m], 70, 30), edgecolors='k', lw=0.4)
    axs[1].set_yscale('log'); axs[1].set_ylabel('summed scan M0 (N m)')
    fig.colorbar(sc, ax=axs, label='completeness at Mw 3.6')
    fig.savefig(f'{K.FIG}/prod_episodes{OUT}.png', dpi=130)


if __name__ == '__main__':
    main()
