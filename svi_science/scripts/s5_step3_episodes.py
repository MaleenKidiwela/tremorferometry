"""S5 POST-HOC (PLAN_S5.md appendix B): is the tidal response different from episode to episode by more than tidal
sampling noise? If yes, era contrasts must be tested against between-EPISODE scatter (episode-label permutation), not
only against the surrogate (tidal-sampling) noise.
Per major episode (ours: major, >= 100 once-counted P>=0.75 events; Bostock: his 10 windows): excess C_M2 for
(i) all events, (ii) 'active' events = 1 <= days since family's own onset < 10 (the strong-response days).
Heterogeneity: chi2 = sum ((Cex_e - mean_w)/sd_e)^2 vs df (sd_e from the surrogates).
Era test: weighted mean Cex 2020-26 minus 2014-19 with p from permuting episode labels between halves (all splits).
"""
import json, itertools, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
import s5_lib as L
import s5_step2_bostock as B
from s5_step3_ours import COMPS, R, OUT


def per_episode(es, ep_codes, dt_days, names, nsurr=1000):
    act = np.where((dt_days >= 1) & (dt_days < 10), ep_codes, -1)
    ne = len(names)
    obs, surr, _, _ = L.run(es, {'ep': (ep_codes, ne), 'ep_act': (act, ne)}, nsurr=nsurr, seed=5)
    rows = []
    for ln in ('ep', 'ep_act'):
        o, s = obs[ln], surr[ln]
        for i, nm in enumerate(names):
            n = o['n'][i]
            if n == 0:
                continue
            cs = s['cos_m2'][:, i] / n
            rows.append(dict(set=ln, episode=nm, n=int(n), Cex=o['cos_m2'][i] / n - cs.mean(), sd=cs.std()))
    return pd.DataFrame(rows)


def hetero(d):
    w = 1 / d.sd ** 2; m = (w * d.Cex).sum() / w.sum()
    chi2 = (((d.Cex - m) / d.sd) ** 2).sum(); df = len(d) - 1
    tau2 = max(0.0, (chi2 - df) / (w.sum() - (w ** 2).sum() / w.sum()))   # DerSimonian-Laird between-episode variance
    return dict(n_episodes=len(d), mean=m, chi2=chi2, df=df, p=float(stats.chi2.sf(chi2, df)), between_sd=float(np.sqrt(tau2)),
                median_sampling_sd=float(d.sd.median()))


def era_perm(d, is_late):
    w = 1 / (d.sd.values ** 2 + d.attrs.get('tau2', 0))
    x = d.Cex.values; g = np.asarray(is_late)
    def stat(gg):
        return (w[gg] * x[gg]).sum() / w[gg].sum() - (w[~gg] * x[~gg]).sum() / w[~gg].sum()
    obs = stat(g)
    k = g.sum(); null = []
    for c in itertools.combinations(range(len(x)), k):
        gg = np.zeros(len(x), bool); gg[list(c)] = True; null.append(stat(gg))
    null = np.array(null)
    return dict(diff=obs, p_two=float(np.mean(np.abs(null) >= abs(obs) - 1e-12)), n_splits=len(null),
                n_late=int(k), n_early=int(len(x) - k), MDD_C=float(2.8 * null.std()), MDD_a_pct=float(200 * 2.8 * null.std()))


def main():
    res = {}
    # ---- ours
    A = pd.read_parquet(R / 's0/analysis_table.parquet', columns=['series', 'fam', 'P_same', 'keep_xfam8', 'time_utc', 'episode', 'major'])
    A['episode'] = A.episode.astype(str).replace('nan', '')
    ev = A[(A.series == 'data') & (A.P_same >= 0.75) & A.keep_xfam8 & A.major].rename(columns={'fam': 'fid', 'time_utc': 't'})
    cnt = ev.episode.value_counts(); keep = cnt[cnt >= 100].index
    ev = ev[ev.episode.isin(keep)].reset_index(drop=True)
    S = L.Stress(pd.read_parquet(R / 's5_setup/tidal_stress_svi_hourly.parquet', columns=['fid', 'time'] + COMPS), COMPS)
    es = L.EventSet(ev, S, COMPS); ev = es.ev
    o = ev.groupby(['fid', 'episode']).t.quantile(0.10).rename('t0').reset_index()
    ev = ev.merge(o, on=['fid', 'episode'], how='left')
    names = sorted(ev.episode.unique()); code = ev.episode.map({n: i for i, n in enumerate(names)}).values
    dt = (ev.t - ev.t0).dt.total_seconds().values / 86400
    D = per_episode(es, code, dt, names); D['catalogue'] = 'ours'
    # ---- Bostock
    bev, W = B.load_bostock(); bev = B.onsets(bev)
    SB = L.Stress(pd.read_parquet(OUT / 'tidal_stress_bostock_hourly.parquet', columns=['fid', 'time'] + B.COMPS), B.COMPS)
    bes = L.EventSet(bev, SB, B.COMPS); bev = bes.ev
    bdt = (bev.t - bev.t0_p10).dt.total_seconds().values / 86400
    DB = per_episode(bes, bev.ep.values.astype(int), bdt, W.episode.tolist()); DB['catalogue'] = 'bostock'
    T = pd.concat([D, DB], ignore_index=True)
    T.to_csv(OUT / 'step3_per_episode.csv', index=False)
    for cat in ('ours', 'bostock'):
        for st in ('ep', 'ep_act'):
            d = T[(T.catalogue == cat) & (T.set == st)].reset_index(drop=True)
            h = hetero(d); res[f'{cat}_{st}_heterogeneity'] = h
            if cat == 'ours':
                d.attrs['tau2'] = h['between_sd'] ** 2
                res[f'ours_{st}_era_episode_perm'] = era_perm(d, d.episode.str[:4].astype(int) >= 2020)
            else:
                d.attrs['tau2'] = h['between_sd'] ** 2
                res[f'bostock_{st}_2009_13_vs_2003_08_episode_perm'] = era_perm(d, d.episode.str[:4].astype(int) >= 2009)
    json.dump(res, open(OUT / 'step3_episode_heterogeneity.json', 'w'), indent=1, default=float)
    print(T.round(3).to_string()); print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
