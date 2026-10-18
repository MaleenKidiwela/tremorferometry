"""S5 step 3 POST-HOC (added after seeing step-3 numbers; see PLAN_S5.md appendix A):
era and depth comparisons STANDARDISED for within-episode timing. Raw era/depth contrasts mix the tidal response with
how many events fall early (first day, ~no response) vs late (days 2-6, strong response) in each era/depth.
Strata = days since the family's own onset (10th pct) in major episodes, 7 bins as PLAN_S5.md section 4.
Standardised excess C = sum_s w_s Cex_s with common weights w_s (pooled event share of both sides); surrogate SD with the
same weights; two-sided p from the surrogate distribution of the difference; MDD = 2.8 SD.
"""
import json, pickle, numpy as np, pandas as pd
from pathlib import Path
import s5_lib as L
from s5_step3_ours import COMPS, codes_dt, R, OUT

DTN = ['<1d', '1-2d', '2-3d', '3-4d', '4-6d', '6-10d', '>=10d']


def strat(o, s, codes_by_stratum):
    """codes_by_stratum: list (per stratum) of code index. returns per-stratum n, Cex, surrogate deviations."""
    n = np.array([o['n'][c] for c in codes_by_stratum])
    ce = np.array([o['cos_m2'][c] / max(o['n'][c], 1) for c in codes_by_stratum])
    cs = np.array([s['cos_m2'][:, c] / max(o['n'][c], 1) for c in codes_by_stratum])   # (nstrata, nsurr)
    return n, ce - cs.mean(1), cs - cs.mean(1, keepdims=True)


def compare(A, B, min_n=20):
    nA, eA, dA = A; nB, eB, dB = B
    ok = (nA >= min_n) & (nB >= min_n)
    w = (nA + nB) * ok; w = w / w.sum()
    sA = (w * eA).sum(); sB = (w * eB).sum()
    zA = (w[:, None] * dA).sum(0); zB = (w[:, None] * dB).sum(0)
    d = sA - sB; dz = zA - zB
    return dict(std_Cex_A=sA, std_Cex_B=sB, a_A_pct=200 * sA, a_B_pct=200 * sB, diff=d, sd=dz.std(),
                p_two=(1 + np.sum(np.abs(dz) >= abs(d))) / (1 + len(dz)), MDD_C=2.8 * dz.std(), MDD_a_pct=200 * 2.8 * dz.std(),
                strata_used=[DTN[i] for i in np.where(ok)[0]], weights=w[ok].round(3).tolist())


def main():
    F = pd.read_csv(R / 's0/families.csv')
    q1, q2 = F.downdip_km.quantile([1 / 3, 2 / 3])
    F['d3'] = np.where(F.downdip_km >= q2, 2, np.where(F.downdip_km <= q1, 0, 1))
    A = pd.read_parquet(R / 's0/analysis_table.parquet', columns=['series', 'fam', 'P_same', 'keep_xfam8', 'time_utc', 'episode', 'major'])
    A['episode'] = A.episode.astype(str).replace('nan', '')
    ev = A[(A.series == 'data') & (A.P_same >= 0.75) & A.keep_xfam8].rename(columns={'fam': 'fid', 'time_utc': 't'}).reset_index(drop=True)
    S = L.Stress(pd.read_parquet(R / 's5_setup/tidal_stress_svi_hourly.parquet', columns=['fid', 'time'] + COMPS), COMPS)
    es = L.EventSet(ev, S, COMPS); ev = es.ev
    mj = ev.major.values & (ev.episode.values != '')
    o = ev[mj].groupby(['fid', 'episode']).t.quantile(0.10).rename('t0').reset_index()
    ev = ev.merge(o, on=['fid', 'episode'], how='left')
    dt = codes_dt(np.where(mj, (ev.t - ev.t0).dt.total_seconds().values / 86400, np.nan))
    half = np.where(ev.t.dt.year.values <= 2019, 0, 1)
    d3 = ev.fid.map(F.set_index('fid').d3).values.astype(int)
    labels = {'dt_half': (np.where(dt >= 0, dt * 2 + half, -1), 14),
              'dt_depth': (np.where(dt >= 0, dt * 3 + d3, -1), 21),
              'dt_depth_half': (np.where(dt >= 0, dt * 6 + d3 * 2 + half, -1), 42),
              'dt': (dt, 7)}
    obs, surr, _, _ = L.run(es, labels, nsurr=1000, seed=5)
    res = {}
    o_, s_ = obs['dt_half'], surr['dt_half']
    res['ours_2020_26_vs_2014_19'] = compare(strat(o_, s_, [i * 2 + 1 for i in range(7)]), strat(o_, s_, [i * 2 for i in range(7)]))
    o_, s_ = obs['dt_depth'], surr['dt_depth']
    res['ours_deep_vs_shallow'] = compare(strat(o_, s_, [i * 3 + 2 for i in range(7)]), strat(o_, s_, [i * 3 for i in range(7)]))
    o_, s_ = obs['dt_depth_half'], surr['dt_depth_half']
    for d, nm in enumerate(('shallow', 'mid', 'deep')):
        res[f'ours_{nm}_2020_26_vs_2014_19'] = compare(strat(o_, s_, [i * 6 + d * 2 + 1 for i in range(7)]), strat(o_, s_, [i * 6 + d * 2 for i in range(7)]))
    # Bostock (onset-based dtbin, dt_depth) from step-2 surrogates, same standardisation
    bo, bs, _, _ = pickle.load(open(OUT / 'step2_surrogates.pkl', 'rb'))
    res['bostock_deep_vs_shallow'] = compare(strat(bo['dt_depth'], bs['dt_depth'], [i * 3 + 2 for i in range(7)]),
                                             strat(bo['dt_depth'], bs['dt_depth'], [i * 3 for i in range(7)]))
    # Bostock vs ours (amplitude NOT interpretable as change: detectors differ) — reported for the record + shape
    res['ours_vs_bostock_amplitude_not_interpretable'] = compare(strat(obs['dt'], surr['dt'], list(range(7))),
                                                                 strat(bo['dtbin'], bs['dtbin'], list(range(7))))
    nb, eb, _ = strat(bo['dtbin'], bs['dtbin'], list(range(7))); no, eo, dzo = strat(obs['dt'], surr['dt'], list(range(7)))
    res['shape'] = dict(bins=DTN, bostock_Cex=eb.tolist(), ours_Cex=eo.tolist(), ours_sd=dzo.std(1).tolist(),
                        bostock_norm=(eb / eb[2:5].mean()).tolist(), ours_norm=(eo / eo[2:5].mean()).tolist(),
                        corr_shape=float(np.corrcoef(eb, eo)[0, 1]))
    json.dump(res, open(OUT / 'step3_posthoc_standardised.json', 'w'), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
