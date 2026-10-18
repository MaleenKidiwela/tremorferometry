#!/usr/bin/env python
"""Amendment Y11 (Merlin review 6): check f diagnosis and f' re-scoring.
  inj38 (prodf2021 env)  -> data/prodf2021/inj_y11_m38.parquet (Mw 3.8 injections with loso_dropped / max_share_ch)
  step1                  -> tables/y11_step1_loso_dropped.csv
  rescore                -> tables/y11_*.csv/json (known VLFEs, level-5 real vs controls under f and f', completeness)
  promoted               -> figures/y11_prom_<k>.png, tables/y11_promoted_between.csv
  fdsn / gobb / report   -> tables/y11_fdsn_stations.csv, y11_gobb_psd.csv, y11_episodes_2021_2026.csv"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K


def inj38():
    import cat_pilot as P
    from multiprocessing import Pool
    rows = pd.read_parquet(f'{P.PD}/inj_list.parquet')
    rows = rows[rows.Mw == 3.8].to_dict('records')
    nch = P.NW * 6; out = []
    with Pool(P.NW, initializer=P.inj_init) as p:
        for r in p.imap_unordered(P.inj_chunk, [rows[k::nch] for k in range(nch)]):
            out += r
    pd.DataFrame(out).to_parquet(f'{P.PD}/inj_y11_m38.parquet')
    K.log(f'Y11 inj38 {P.NAME}{P.PILOT}: {len(out)}', 'y11')


def fprime(zl, zb):
    zl = np.asarray(zl, float); zb = np.asarray(zb, float)
    return ((zl >= 3) & (zl >= 0.5 * zb)) | (zl >= 5)


def fstrict(zl, zb):
    zl = np.asarray(zl, float); zb = np.asarray(zb, float)
    return (zl >= 3) & (zl >= 0.5 * zb)


RUNS = {2010: 'prod', 2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq',
        2018: 'prodq', 2019: 'prodf', 2020: 'prodf', 2021: 'prodf'}


def step1():
    J = pd.read_parquet(f'{K.DATA}/prodf2021/inj_list.parquet').merge(pd.read_parquet(f'{K.DATA}/prodf2021/inj_y11_m38.parquet'), on='inj_id')
    J['month'] = pd.to_datetime(J.epoch, unit='s').dt.month
    J['grp'] = np.where(J.source.isin(['C0', 'C1']), J.source, 'relay')
    J['mf2'] = (J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5)
    rows = []
    for (g, m), x in J.groupby(['grp', 'month']):
        vc = x.loso_dropped.value_counts(normalize=True)
        mc = x.max_share_ch.astype(str).str.rsplit('.', n=1).str[0].value_counts(normalize=True)
        rows.append(dict(source=g, month=m, n=len(x), loso_dropped_top=vc.index[0] if len(vc) else '', frac_top=float(vc.iloc[0]) if len(vc) else np.nan,
                         frac_LRIV=float((x.loso_dropped == 'UW.LRIV').mean()), frac_PGC=float((x.loso_dropped == 'CN.PGC').mean()),
                         max_share_sta_top=mc.index[0] if len(mc) else '', max_share_med=float(x.max_share.median()),
                         zloso_over_zbest_med=float((x.z_loso / x.z_best).median()), z_loso_med=float(x.z_loso.median()), z_best_med=float(x.z_best.median()),
                         f=float(fstrict(x.z_loso, x.z_best).mean()), fprime=float(fprime(x.z_loso, x.z_best).mean())))
    T = pd.DataFrame(rows); T.to_csv(f'{K.TB}/y11_step1_loso_dropped.csv', index=False)
    jd = J[(J.month >= 7)]
    S = dict(C0_JulDec_frac_drop_LRIV=float((jd[jd.source == 'C0'].loso_dropped == 'UW.LRIV').mean()),
             C1_JulDec_frac_drop_LRIV=float((jd[jd.source == 'C1'].loso_dropped == 'UW.LRIV').mean()),
             relay_JulDec_frac_drop_PGC=float((jd[jd.grp == 'relay'].loso_dropped == 'CN.PGC').mean()),
             C0_JulDec_loso_dropped=jd[jd.source == 'C0'].loso_dropped.value_counts().to_dict(),
             relay_JulDec_loso_dropped=jd[jd.grp == 'relay'].loso_dropped.value_counts().to_dict())
    S['accept'] = bool(S['C0_JulDec_frac_drop_LRIV'] >= 0.8 and S['relay_JulDec_frac_drop_PGC'] >= 0.8)
    json.dump(S, open(f'{K.TB}/y11_step1.json', 'w'), indent=1, default=str)
    pd.set_option('display.width', 250); print(T.round(2).to_string(index=False)); print(json.dumps(S, indent=1, default=str))


def levels(A, fcol):
    """A: y9_all rows (agreement rows have chk_* and z_loso_b/z_best_b). returns level with f replaced by fcol (bool array)."""
    base = sum(A[f'chk_{k}'].fillna(False).astype(bool).astype(int) for k in 'abce')
    return np.where(A.agree_bin == 1, base + (fcol & A.z_loso_b.notna()).astype(int), -1)


def rescore():
    import cat_ychecks as YC
    out = {}
    # (i) known VLFEs
    import cat_y9 as Y9
    known = Y9.known_lists()
    ALL = []
    for y in range(2010, 2019):
        A = pd.read_parquet(f'{K.DATA}/ranked/y9_all_{y}.parquet'); A['year'] = y
        A['f_old'] = fstrict(A.z_loso_b, A.z_best_b); A['f_new'] = fprime(A.z_loso_b, A.z_best_b)
        A['level_f'] = levels(A, A.f_old); A['level_fp'] = levels(A, A.f_new)
        ALL.append(A)
    ALL = pd.concat(ALL, ignore_index=True)
    R = ALL[ALL.ctrl == 'real']; kn = R[Y9.near(R.epoch.values, known)]
    out['known'] = dict(n=len(kn), level5_f=int((kn.level_f == 5).sum()), level5_fprime=int((kn.level_fp == 5).sum()),
                        f_pass=int(kn.f_old.sum()), fprime_pass=int(kn.f_new.sum()))
    # (ii) pooled level-5 between, real vs controls
    nn = ALL[ALL.ctrl == 'noise_template'].groupby('year').template.nunique()
    rows = []
    for lab, col in (('f', 'level_f'), ('fprime', 'level_fp')):
        b = ALL[(ALL.period == 'between') & (ALL[col] == 5)]
        r = dict(rule=lab, real=int((b.ctrl == 'real').sum()))
        r['noise_x11'] = float(sum(11 * ((b.ctrl == 'noise_template') & (b.year == y)).sum() / nn[y] for y in nn.index))
        r['noise_raw'] = int((b.ctrl == 'noise_template').sum())
        r['A1_per_draw'] = float(((b.ctrl == 'scrambled') & (b.source == 'A1')).sum() / 2)
        r['A2_per_draw'] = float(((b.ctrl == 'scrambled') & (b.source == 'A2')).sum() / 2)
        rows.append(r)
    # matched banks (accepted, all years available), scaled to the 2010-2018 between live time
    S = pd.read_csv(f'{K.TB}/y9nb_summary.csv'); S['accepted'] = S.agree_ratio_to_real.between(1 / 1.5, 1.5)
    nb = []
    for y in sorted(S.year.unique()):
        run = {2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq'}[y]
        d = f'{K.DATA}/y9nb/{run}{y}'
        A = pd.read_parquet(f'{d}/rows.parquet'); b = K.apply_rules(pd.read_parquet(f'{d}/bat.parquet'))
        A = A.merge(b[['i', 'pass_a', 'pass_b', 'pass_c', 'pass_e2', 'a_testable', 'travelling', 'z_best', 'z_loso']].rename(columns={'z_best': 'zb'}), on='i', how='left')
        base = (A.a_testable.fillna(False).astype(bool) & A.pass_a.fillna(False).astype(bool)).astype(int) + (A.travelling.notna() & A.pass_b.fillna(False).astype(bool)).astype(int) \
            + (A.zb.notna() & A.pass_c.fillna(False).astype(bool)).astype(int) + A.pass_e2.fillna(False).astype(bool).astype(int)
        A['l_f'] = base + (A.z_loso.notna() & fstrict(A.z_loso, A.zb)).astype(int)
        A['l_fp'] = base + (A.z_loso.notna() & fprime(A.z_loso, A.zb)).astype(int)
        A['year'] = y; nb.append(A)
    NB = pd.concat(nb, ignore_index=True)
    acc = S[S.accepted][['year', 'bank', 'between_yr']]
    NBa = NB.merge(acc, on=['year', 'bank'])
    T = float(acc.between_yr.sum()); live = 5.535934291581109
    for r, col in zip(rows, ('l_f', 'l_fp')):
        n = int(((NBa.period == 'between') & (NBa[col] == 5)).sum())
        r['matched_banks_level5'] = n; r['matched_bank_years'] = T; r['matched_scaled_2010_18'] = n / T * live
    P = pd.DataFrame(rows); P.to_csv(f'{K.TB}/y11_level5_f_vs_fprime.csv', index=False)
    rise = {c: float(P.loc[1, c] - P.loc[0, c]) for c in ('noise_x11', 'A1_per_draw', 'A2_per_draw', 'matched_scaled_2010_18')}
    out['level5_between'] = P.round(3).to_dict('records'); out['control_rise'] = rise
    # (iii) C0/C1 completeness at Mw 3.8 by year and month, 2013-2021
    cr = []
    for y in range(2013, 2022):
        run = RUNS[y]
        J = pd.read_parquet(f'{K.DATA}/{run}{y}/inj_list.parquet').merge(pd.read_parquet(f'{K.DATA}/{run}{y}/inj_v2.parquet'), on='inj_id')
        J = J[(J.Mw == 3.8) & J.source.isin(['C0', 'C1'])]
        J['month'] = pd.to_datetime(J.epoch, unit='s').dt.month
        J['mf2'] = (J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5)
        base = J.mf2 & J.pass_a.fillna(False).astype(bool) & J.pass_b.fillna(False).astype(bool) & J.pass_c.fillna(False).astype(bool) & J.pass_e2.fillna(False).astype(bool)
        J['rec_f'] = base & fstrict(J.z_loso, J.z_best); J['rec_fp'] = base & fprime(J.z_loso, J.z_best)
        J['f_given_mf2'] = fstrict(J.z_loso, J.z_best)
        for (src, m), x in J.groupby(['source', 'month']):
            cr.append(dict(year=y, run=run, source=src, month=m, n=len(x), compl_f=float(x.rec_f.mean()), compl_fprime=float(x.rec_fp.mean()),
                           f_pass_among_mf2=float(x[x.mf2].f_given_mf2.mean()) if x.mf2.any() else np.nan))
    C = pd.DataFrame(cr); C['rise'] = C.compl_fprime - C.compl_f; C.to_csv(f'{K.TB}/y11_completeness_by_month.csv', index=False)
    c0 = C[C.source == 'C0']; low = c0[c0.f_pass_among_mf2 < 0.5]
    w = low.n
    out['C0_low_f_months'] = dict(n_months=len(low), mean_rise=float((low.rise * w).sum() / w.sum()) if len(low) else None,
                                  min_rise=float(low.rise.min()) if len(low) else None, frac_months_rise_ge_0p15=float((low.rise >= 0.15).mean()) if len(low) else None)
    Y = C.groupby(['year', 'source']).apply(lambda d: pd.Series(dict(compl_f=(d.compl_f * d.n).sum() / d.n.sum(), compl_fprime=(d.compl_fprime * d.n).sum() / d.n.sum()))).reset_index()
    Y.to_csv(f'{K.TB}/y11_completeness_by_year.csv', index=False)
    adopt_ctrl = all(v <= 1 for v in rise.values())
    adopt_c0 = out['C0_low_f_months']['mean_rise'] is not None and out['C0_low_f_months']['mean_rise'] >= 0.15
    out['adopt'] = dict(controls_rise_le_1=adopt_ctrl, c0_rise_ge_0p15=adopt_c0, adopt=bool(adopt_ctrl and adopt_c0))
    # (iv) promoted real between rows
    pr = R[(R.period == 'between') & (R.level_f == 4) & (R.level_fp == 5)]
    pr.to_csv(f'{K.TB}/y11_promoted_between.csv', index=False)
    out['promoted_between_real'] = len(pr)
    json.dump(out, open(f'{K.TB}/y11_rescore.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250)
    print(json.dumps(out, indent=1, default=float)); print(Y.round(2).pivot(index='year', columns='source').to_string())
    print(c0[c0.year >= 2017].round(2).to_string(index=False))


if __name__ == '__main__':
    {'inj38': inj38, 'step1': step1, 'rescore': rescore}[sys.argv[1]]()
