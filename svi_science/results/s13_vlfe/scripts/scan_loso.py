#!/usr/bin/env python
"""POST-HOC diagnostic (not a pre-registered rule; written after the pilot verdict): the high-z false alarms that pass
the headline filter are single-station events (max_share median 0.99 at z >= 15) at nodes under that station. Here: for
every headline-filtered positive cluster with z >= 7 (headline mask), recompute the statistic at its best node WITHOUT
the station that contributes most to sum w T' d ("z_loso"), normalised by the same daily MAD scale. Runs: real (incl. the
50 polarity draws), ctrlA_s*, inj2. Writes data/scan/runs/RUN/loso.parquet and tables/scan_loso_separation.csv."""
import os, sys
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scan_common as S
import scan_columns as SC
from scan_validate import headline_ok, between_mask, ets_vec, live_years, R_

USE = ['n_sta', 'n_Z_agree', 'share_excess', 'vr_ratio', 'vr_gt_ss', 'vr_gt_steep', 'isolation', 'not_hf_loud']
W = SC.W


def zloso(r):
    stas, offs = W['stas'], W['offs']
    v, sgn, tsec, node = int(r.variant), int(r.sign), int(r.tsec), int(r.node)
    d, u = divmod(tsec, 86400)
    sv = W['signs'][v] * sgn
    srcday = (d + offs) % S.NDAY
    g = W['g'][d]
    T = W['T']['thrust'][node] * g[:, :, None]
    num = np.zeros(len(stas)); den = np.zeros(len(stas))
    for i, s in enumerate(stas):
        ts = int(srcday[i]) * 86400 + u
        x = SC.window(W['BP'][s], ts + S.FIT0, ts + S.FIT1 + 1) * 1e9
        sig = W['DS'][s]['sig'][srcday[i]] * 1e9
        for c in range(3):
            if W['won'][d, i, c] and np.isfinite(x[c]).all() and sig[c] > 0:
                w = 1 / sig[c] ** 2
                num[i] += w * (T[i, c] * x[c] * sv[i, c]).sum(); den[i] += w * (T[i, c] ** 2).sum()
    mad = W['mad'][v, node, d]
    k = int(np.argmax(num))
    keep = np.arange(len(stas)) != k
    zf = num.sum() / np.sqrt(den.sum()) / mad if den.sum() > 0 else np.nan
    zl = num[keep].sum() / np.sqrt(den[keep].sum()) / mad if den[keep].sum() > 0 else np.nan
    return dict(z_full_chk=zf, z_loso=zl, loso_dropped=stas[k], dom_den_share=den[k] / den.sum() if den.sum() > 0 else np.nan)


def run_all(rows):
    return [dict(idx=i, **zloso(r)) for i, r in rows]


def one(run, bp):
    D = pd.read_parquet(f'{R_}/{run}/cand_cols.parquet')
    mh = S.mask_vector('head')
    D = D[(~mh[D.tsec.values]) & (D.sign > 0) & (D.z >= 7)]
    D = D[headline_ok(D, USE)]
    from multiprocessing import Pool
    rows = list(D.iterrows()); chunks = [rows[k::32] for k in range(32)]
    res = []
    with Pool(8, initializer=SC.init, initargs=(run, bp)) as p:
        for x in p.imap_unordered(run_all, chunks):
            res += x
    D = D.join(pd.DataFrame(res).set_index('idx'))
    D.to_parquet(f'{R_}/{run}/loso.parquet')
    return D


if __name__ == '__main__':
    out = {}
    out['real'] = one('real', f'{S.SD}/bp_2011g')
    for s in (3, 5, 7, 11, 13):
        out[f'a{s}'] = one(f'ctrlA_s{s}', f'{S.SD}/bp_2011g')
    out['inj'] = one('inj2', f'{S.SD}/bp_2011inj2')
    mh = S.mask_vector('head'); btw = between_mask(); ets = ets_vec()
    nav = np.load(f'{R_}/real/navail.npy'); ly = live_years(nav, mh); lyb = live_years(nav, mh, btw)
    r = out['real']; r0 = r[r.variant == 0]; rb = r[r.variant > 0]; nd = 50
    J = pd.read_parquet(f'{S.SD}/inject_2011.parquet')
    inj = out['inj']
    rows = []
    for thr in (4, 5, 6, 7, 8, 10, 12, 15):
        sel = lambda D: D[(D.z_loso >= thr)]                                  # noqa: E731
        a = np.mean([len(sel(out[f'a{s}'])) / live_years(np.load(f'{R_}/ctrlA_s{s}/navail.npy'), mh) for s in (3, 5, 7, 11, 13)])
        b = sel(rb); b = b[btw[b.tsec.values]]
        p = sel(r0)
        ii = sel(inj); t = ii.tsec.values
        rec = {}
        for mw in (3.4, 3.6, 3.8, 4.0):
            jj = J[(J.gf == 'standard') & (J.mw == mw)]
            rec[mw] = float(np.mean([np.any(np.abs(t - x) <= 60) for x in jj.tsec]))
        rows.append(dict(thr_z_loso=thr, real_pos_per_yr=len(p) / ly, real_between_per_yr=btw[p.tsec.values].sum() / lyb,
                         real_ets=int(ets[p.tsec.values].sum()), ctrlA_per_yr=a, ctrlB_between_per_yr=len(b) / nd / lyb,
                         **{f'inj_rec_mw{k}': v for k, v in rec.items()}))
    T = pd.DataFrame(rows); T.to_csv(f'{S.C.TB}/scan_loso_separation.csv', index=False)
    print(T.round(2).to_string())
