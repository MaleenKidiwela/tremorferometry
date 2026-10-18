"""S5 self-tests of s5_lib on the setup's 2014-2026 hourly stress (no LFE data touched):
 1. phase convention: |hc| f cos(theta) reproduces the M2 part of tau (corr ~1)
 2. shifted-time phase formula == direct computation at the shifted time
 3. synthetic injection: events drawn with rate 1 + 0.3 cos(theta_M2) on random days -> recovered a ~ 30 %, theta0 ~ 0;
    and with a -90 deg lag -> theta0 ~ -90
"""
import json, numpy as np, pandas as pd
from pathlib import Path
import s5_lib as L
import pyTMD.constituents as PC

OUT = Path('/home/jovyan/tremorferometry/svi_science/results/s5')
S0 = Path('/home/jovyan/tremorferometry/svi_science/results/s5_setup/tidal_stress_svi_hourly.parquet')
COMPS = ['tau_Pa', 'sigma_n_Pa']
df = pd.read_parquet(S0, columns=['fid', 'time'] + COMPS, filters=[('fid', 'in', [10, 40, 77])])
S = L.Stress(df, COMPS)
res = {}
# 1
i = 0; t = df.time[df.fid == S.fids[i]].values
pu, pf, G = PC.arguments(L.mjd(t), ['m2'], corrections='FES')
hc = S.hc['tau_Pa']['m2'][i]
m2_series = (pf[:, 0] * (hc.real * np.cos(np.radians(G[:, 0]) + pu[:, 0]) - hc.imag * np.sin(np.radians(G[:, 0]) + pu[:, 0])))
ev = pd.DataFrame({'fid': S.fids[i], 't': t[::7]})
es = L.EventSet(ev, S, COMPS)
q = es.quantities(np.zeros(len(ev), np.int64))
res['corr_cos_theta_vs_M2_series'] = float(np.corrcoef(q['cos_m2'], m2_series[::7])[0, 1])
res['corr_cos_theta_vs_total_tau'] = float(np.corrcoef(q['cos_m2'], S.A['tau_Pa'][i, ::7])[0, 1])
# 2
rng = np.random.default_rng(1)
d = rng.integers(-3000, 3000, len(ev))
q2 = es.quantities(d)
th_direct = L.astro_th((es.tns + d * L.NS_D).astype('datetime64[ns]'))
ph_direct = th_direct[:, 0] + es.arg['m2']
res['max_abs_diff_shifted_phase_deg'] = float(np.degrees(np.abs(np.angle(np.exp(1j * (np.arctan2(q2['sin_m2'], q2['cos_m2']) - ph_direct))))).max())
ph_o1 = th_direct[:, 1] + es.arg['o1']
res['max_abs_diff_shifted_phase_o1_deg'] = float(np.degrees(np.abs(np.angle(np.exp(1j * (np.arctan2(q2['sin_o1'], q2['cos_o1']) - ph_o1))))).max())
# 3 synthetic injection with bursts (whole hours of events) on random active days
for lag in (0, -90):
    rows = []
    for f in S.fids:
        days = rng.choice(np.arange(np.datetime64('2015-01-01'), np.datetime64('2025-01-01')), 40, replace=False)
        cand = (days[:, None].astype('datetime64[ns]') + (np.arange(0, 24 * 60, 10) * 60 * 10**9).astype('timedelta64[ns]')).ravel()
        e0 = L.EventSet(pd.DataFrame({'fid': f, 't': cand}), S, COMPS)
        qq = e0.quantities(np.zeros(len(cand), np.int64))
        ph = np.arctan2(qq['sin_m2'], qq['cos_m2'])
        p = 0.05 * (1 + 0.3 * np.cos(ph - np.radians(lag)))
        keep = rng.random(len(cand)) < p
        rows.append(pd.DataFrame({'fid': f, 't': cand[keep]}))
    E = pd.concat(rows, ignore_index=True)
    es3 = L.EventSet(E, S, COMPS)
    obs, surr, _, _ = L.run(es3, {'all': (np.zeros(len(E), int), 1)}, nsurr=200, seed=3)
    T = L.summarize(obs, surr, 'all', ['all'], const=['m2', 'o1'])
    res[f'injection_lag{lag}'] = T[['n', 'a_m2_pct', 'theta0_m2_deg', 'p_C_m2', 'a_o1_pct', 'p_C_o1']].iloc[0].to_dict()
json.dump(res, open(OUT / 'selftest.json', 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
