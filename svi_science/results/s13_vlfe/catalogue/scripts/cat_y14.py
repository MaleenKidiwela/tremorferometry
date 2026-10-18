#!/usr/bin/env python
"""Amendment Y14 summary: per year, f (pre-Y14) vs f' -> tables/y14_summary.csv"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
RUNS = ['prod2010', 'prod2011', 'prod2012', 'prodq2013', 'prod2014', 'prod2015', 'prodq2016', 'prodq2017', 'prodq2018', 'prodrq2018',
        'prodf2019', 'prodf2020', 'prodfy2021', 'prodfy2022', 'prodfy2023']


def main():
    L1 = pd.read_csv(f'{K.TB}/y9_ladder_table.csv'); L0 = pd.read_csv(f'{K.TB}/preY14/y9_ladder_table.csv')
    rows = []
    for RY in RUNS:
        y = int(RY[-4:])
        o = dict(run=RY, year=y)
        for lab, d in (('f', f'{K.TB}/preY14'), ('fp', K.TB)):
            v = json.load(open(f'{d}/{RY}_R_verdict_v2.json')); o[f'ETS_tier_{lab}'] = v['T1'].get('n_ETS'); o[f'GO_{lab}'] = v['GO']
            o[f'T1_{lab}'] = v['T1']['pass']; o[f'T3_{lab}'] = f"{v['T3']['n_pass']}/{v['T3']['n']}"
            E = pd.read_csv(f'{d}/{RY}_events_v2.csv')
            o[f'relay_vetted_{lab}'] = int(((E.period == 'ETS') & (E.lfe_status == 'with LFEs') & (E.vr_thrust >= 0.12)).sum())
        J = pd.read_parquet(f'{K.DATA}/{RY}/inj_list.parquet').merge(pd.read_parquet(f'{K.DATA}/{RY}/inj_v2.parquet'), on='inj_id')
        J = K.rescore_inj(J); J = J[(J.Mw == 3.8) & J.source.isin(['C0', 'C1'])]
        mf2 = (J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5)
        o['C0C1_38_f'] = float((mf2 & J.battery_e2_preY14.fillna(False).astype(bool)).mean()); o['C0C1_38_fp'] = float((mf2 & J.battery_e2).mean())
        for lab, Lt in (('f', L0), ('fp', L1)):
            x = Lt[(Lt.year == str(y)) & (Lt.period == 'between') & (Lt.level_ge == 5)]
            o[f'between_L5_{lab}'] = int(x.real.iloc[0]) if len(x) and RY not in ('prodrq2018',) else np.nan
        rows.append(o)
    T = pd.DataFrame(rows); T.to_csv(f'{K.TB}/y14_summary.csv', index=False)
    pd.set_option('display.width', 250)
    print(T[['run', 'ETS_tier_f', 'ETS_tier_fp', 'between_L5_f', 'between_L5_fp', 'C0C1_38_f', 'C0C1_38_fp', 'T1_fp', 'T3_fp', 'GO_f', 'GO_fp',
             'relay_vetted_f', 'relay_vetted_fp']].round(2).to_string(index=False))
    for p in ('ETS', 'between'):
        for lab, Lt in (('f', L0), ('fp', L1)):
            x = Lt[(Lt.year.isin([str(v) for v in range(2010, 2019)])) & (Lt.period == p) & (Lt.level_ge == 5)]
            print(p, lab, '2010-18 level5 real', int(x.real.sum()), 'noise_x11', round(x.noise_x11.sum(), 2), 'A1', x.A1_per_draw.sum(), 'A2', x.A2_per_draw.sum())


if __name__ == '__main__':
    main()
