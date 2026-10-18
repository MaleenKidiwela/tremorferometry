#!/usr/bin/env python
"""Amendment Y12 (stations first).
  qc Y STA,STA..   -> tables/y12_qc_<Y>.csv (daily 20-50 s Z level vs CN.PGC; largest-teleseism amplitude vs PGC; pass flag)
  sources          -> data/inj/sources_y12.npz (sources.npz + Green's-function scan sources for the Y12 stations)
  mkrun Y          -> data/prodfy<Y>/ seeded from prodf<Y> (bank, mf, rows, inj_list reused) + STEP markers
  rescore Y        -> tables/y12_fvsfprime_<Y>.csv/json (C0/C1 completeness by month, f and f', pre-Y12 vs Y12)"""
import os, sys, json, glob, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

NEW = ['CN.SYMB', 'UW.BUCKS', 'UW.OSQM', 'UW.RNWY', 'UW.SLDQ', 'UW.ODUC']
DISP = f'{K.S13}/data/disp1sps_v3'


def band(x):
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, [0.02, 0.05], btype='band', fs=1, output='sos')
    x = np.nan_to_num(x - np.nanmean(x))
    return sosfiltfilt(sos, x)


def dayz(s, y, doy):
    f = f'{DISP}/{s}/{y}/{doy:03d}.npz'
    if not os.path.isfile(f):
        return None
    d = np.load(f)['d'][2].astype(float)
    return d if np.isfinite(d).sum() > 0.8 * len(d) else None


def qc(y, stas):
    rows = []
    for s in stas:
        lv, nlive, nflat = [], 0, 0
        for doy in range(1, 367):
            a = dayz(s, y, doy)
            if a is None:
                continue
            nlive += 1
            if np.nanstd(a) == 0:
                nflat += 1; continue
            p = dayz('CN.PGC', y, doy)
            if p is None or np.nanstd(p) == 0:
                continue
            lv.append(np.median(np.abs(band(a))) / np.median(np.abs(band(p))))
        # largest teleseism of the year (global list, > 3000 km)
        EQ = K.eq_table(y); EQ = EQ[(EQ.dist_km > 3000) & (pd.to_datetime(EQ.t, unit='s').dt.year == y)].sort_values('mag', ascending=False)
        tele = np.nan; tdesc = ''
        for e in EQ.head(5).itertuples():
            d = pd.Timestamp(e.t, unit='s'); doy = d.dayofyear
            a = dayz(s, y, doy); p = dayz('CN.PGC', y, doy)
            if a is None or p is None:
                continue
            i0 = int(e.t - pd.Timestamp(f'{d.date()}').value / 1e9); i1 = min(86400, i0 + int(e.dist_km / 3.0) + 1800)
            if i0 < 0 or i1 - i0 < 600:
                continue
            tele = float(np.max(np.abs(band(a)[i0:i1])) / np.max(np.abs(band(p)[i0:i1]))); tdesc = f'{e.time_pt} M{e.mag:.1f} {e.dist_km:.0f} km'
            break
        med = float(np.median(lv)) if lv else np.nan
        ok = bool(lv and 1 / 3 <= med <= 3 and nflat <= 0.1 * max(nlive, 1))   # Y12b: teleseism ratio is a column, not a gate
        rows.append(dict(year=y, sta=s, live_days=nlive, flat_days=nflat, n_ratio_days=len(lv), level_ratio_med=med,
                         level_ratio_p10=float(np.percentile(lv, 10)) if lv else np.nan, level_ratio_p90=float(np.percentile(lv, 90)) if lv else np.nan,
                         tele_ratio=tele, tele=tdesc, qc_pass=ok))
    T = pd.DataFrame(rows); T.to_csv(f'{K.TB}/y12_qc_{y}.csv', index=False)
    print(T.round(2).to_string(index=False))
    print('PASS', ','.join(T[T.qc_pass].sta))


def sources():
    import cat_inj as I
    z = dict(np.load(f'{I.INJ}/sources.npz')); M = pd.read_csv(f'{K.TB}/inj_sources.csv').set_index('source')
    for Y in range(2021, 2027):
        f = f'{K.S13}/data/scan_{Y}/templates_grid_v2.npz'
        if not os.path.isfile(f):
            continue
        G = np.load(f); st = list(G['stations'])
        for s in NEW:
            if s not in st:
                continue
            for src in M.index:
                k = f'{src}|sc|{s}'
                if k in z:
                    continue
                m0 = 10 ** (1.5 * M.loc[src, 'Mw_ref'] + 9.1)
                z[k] = m0 * G['thrust'][int(M.loc[src, 'node']), st.index(s)]          # (3, 181) m on W0..W1
                print('GF source', k, '%.2e' % np.abs(z[k][2]).max())
    np.savez(f'{I.INJ}/sources_y12.npz', **z)
    for fac, tag in ((0.5, 'x0p5'), (2.0, 'x2')):                 # sensitivity: only the Green's-function (new-station) keys scaled
        zz = {k: (v * fac if k.split('|')[1] == 'sc' and k.split('|')[2] in NEW else v) for k, v in z.items()}
        np.savez(f'{I.INJ}/sources_y12_{tag}.npz', **zz)


def injvar(tag):
    """C0/C1 injections (all Mw) with INJ_SOURCES set to a sensitivity variant -> {PD}/inj_v2_<tag>.parquet"""
    import cat_pilot as P
    from multiprocessing import Pool
    rows = pd.read_parquet(f'{P.PD}/inj_list.parquet'); rows = rows[rows.source.isin(['C0', 'C1'])].to_dict('records')
    nch = P.NW * 6; out = []
    with Pool(P.NW, initializer=P.inj_init) as p:
        for r in p.imap_unordered(P.inj_chunk, [rows[k::nch] for k in range(nch)]):
            out += r
    pd.DataFrame(out).to_parquet(f'{P.PD}/inj_v2_{tag}.parquet')


def mkrun(y):
    src, dst = f'{K.DATA}/prodf{y}', f'{K.DATA}/prodfy{y}'
    os.makedirs(dst, exist_ok=True)
    for f in ('bank.npz', 'rows.parquet', 'inj_list.parquet'):
        if not os.path.exists(f'{dst}/{f}') and os.path.exists(f'{src}/{f}'):
            shutil.copy(f'{src}/{f}', f'{dst}/{f}')
    if not os.path.exists(f'{dst}/mf'):
        os.symlink(f'{src}/mf', f'{dst}/mf')
    for t in ('bank', 'noise_threshold'):
        a = f'{K.TB}/prodf{y}_{t}.csv'
        if os.path.isfile(a):
            shutil.copy(a, f'{K.TB}/prodfy{y}_{t}.csv')
    for m in ('stations', 'cont', 'diag', 'qcx', 'bank', 'mf', 'detect', 'inj_list'):
        if os.path.exists(f'{src}/STEP_{m}_DONE') and (m != 'inj_list' or os.path.exists(f'{dst}/inj_list.parquet')):
            open(f'{dst}/STEP_{m}_DONE', 'w').close()
    print('seeded', dst)


def rescore(y):
    from cat_y11 import fstrict, fprime
    rows = []
    for run, fn in ((f'prodf{y}', 'inj_v2'), (f'prodfy{y}', 'inj_v2_x0p5'), (f'prodfy{y}', 'inj_v2'), (f'prodfy{y}', 'inj_v2_x2')):
        if not os.path.isfile(f'{K.DATA}/{run}/{fn}.parquet'):
            continue
        J = pd.read_parquet(f'{K.DATA}/{run}/inj_list.parquet').merge(pd.read_parquet(f'{K.DATA}/{run}/{fn}.parquet'), on='inj_id')
        run = run if run.startswith('prodf2') else f'{run}{fn[6:] or "_x1"}'
        J = J[J.source.isin(['C0', 'C1'])]
        J['month'] = pd.to_datetime(J.epoch, unit='s').dt.month
        J['mf2'] = (J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5)
        base = J.mf2 & J.pass_a.fillna(False).astype(bool) & J.pass_b.fillna(False).astype(bool) & J.pass_c.fillna(False).astype(bool) & J.pass_e2.fillna(False).astype(bool)
        J['rec_f'] = base & fstrict(J.z_loso, J.z_best); J['rec_fp'] = base & fprime(J.z_loso, J.z_best)
        for (s, mw, m), x in J.groupby(['source', 'Mw', 'month']):
            rows.append(dict(run=run, source=s, Mw=mw, month=m, n=len(x), compl_f=float(x.rec_f.mean()), compl_fprime=float(x.rec_fp.mean()),
                             f_pass_among_mf2=float(fstrict(x[x.mf2].z_loso, x[x.mf2].z_best).mean()) if x.mf2.any() else np.nan,
                             zloso_over_zbest_med=float((x.z_loso / x.z_best).median()),
                             loso_dropped_top=(x.loso_dropped.value_counts().index[0] if 'loso_dropped' in x and x.loso_dropped.notna().any() else '')))
    T = pd.DataFrame(rows); T.to_csv(f'{K.TB}/y12_fvsfprime_{y}.csv', index=False)
    t = T[(T.Mw == 3.8) & (T.source == 'C0') & (T.month >= 7)]
    S = {r: dict(C0_JulDec_f_pass_among_mf2=float(g.f_pass_among_mf2.mean()), C0_JulDec_compl_f=float((g.compl_f * g.n).sum() / g.n.sum()),
                 C0_JulDec_compl_fprime=float((g.compl_fprime * g.n).sum() / g.n.sum())) for r, g in t.groupby('run')}
    yy = f'prodfy{y}_x1'
    S['f_still_leaks'] = bool(S.get(yy, {}).get('C0_JulDec_f_pass_among_mf2', 0) < 0.5)
    o = T[(T.Mw == 3.8) & (T.source == 'C0') & (T.month == 10)].set_index('run')
    S['C0_Oct_Mw38_compl_f'] = o.compl_f.round(3).to_dict(); S['C0_Oct_Mw38_compl_fprime'] = o.compl_fprime.round(3).to_dict()
    json.dump(S, open(f'{K.TB}/y12_fvsfprime_{y}.json', 'w'), indent=1)
    pd.set_option('display.width', 250)
    print(T[T.Mw == 3.8].pivot_table(index=['source', 'month'], columns='run', values=['compl_f', 'compl_fprime', 'f_pass_among_mf2']).round(2).to_string())
    print(json.dumps(S, indent=1))


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'qc':
        qc(int(sys.argv[2]), sys.argv[3].split(','))
    elif m == 'mkrun':
        mkrun(int(sys.argv[2]))
    elif m == 'injvar':
        injvar(sys.argv[2])
    elif m == 'rescore':
        rescore(int(sys.argv[2]))
    else:
        sources()
