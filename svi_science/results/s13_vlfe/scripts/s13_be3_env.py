#!/usr/bin/env python
"""S13 Amendment 7: origins of the between-ETS PNSN tremor windows from the network 2-8 Hz envelope.
1) 1-s RMS envelopes of the 2-8 Hz vertical (scan_hfenv.process, unchanged; data/scan/hfenv/NET.STA/YYYY/DOY.npy) for
   CN.PGC, CN.VGZ, CN.SNB, CN.NLLB, UW.LRIV on the PNSN days (tables/be3_env_days.csv).
2) per PNSN epoch t (tables/be3_pnsn_epochs.csv): each station envelope / its day median, shifted by -r/3.95 km/s (r = hypocentral
   distance from the PNSN epicentre at Slab2 depth), averaged over the >= 2 finite stations; origin = argmax in [t-150, t+300] s.
Writes tables/be3_pnsn_origins.csv. usage: s13_be3_env.py [--workers 3]"""
import os, sys, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G
import scan_hfenv as H

STAS = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'UW.LRIV']
W0, W1 = -150, 300


def _job(a):
    try:
        return a, H.process(*a)
    except Exception as e:  # noqa
        return a, f'ERR {type(e).__name__}: {str(e)[:80]}'


def env(sta, day):
    d = pd.Timestamp(day)
    f = f'{H.OUTD}/{sta}/{d.year}/{d.dayofyear:03d}.npy'
    if not os.path.isfile(f):
        return None
    e = np.load(f).astype(np.float64)
    m = np.nanmedian(e) if np.isfinite(e).any() else np.nan
    return e / m if np.isfinite(m) and m > 0 else None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--test-date', default='', help='code test: one date, output to --test-out only')
    ap.add_argument('--test-out', default='')
    a = ap.parse_args()
    days = pd.read_csv(f'{C.TB}/be3_env_days.csv').date.tolist()
    if a.test_date:
        days = [a.test_date]
    jobs = [(s, d) for d in days for s in STAS]
    from multiprocessing import Pool
    with Pool(a.workers) as p:
        for (s, d), r in p.imap_unordered(_job, jobs):
            C.log(f'be3 env {s} {d} {r}', f'{C.LOGD}/scan_hfenv.log')
    T = pd.read_csv(f'{C.TB}/be3_pnsn_epochs.csv')
    T = T[T.noLFE.astype(bool)].copy()
    if a.test_date:
        T = T[T.date == a.test_date]
    cc = C.station_coords()
    rows = []
    cache = {}
    for r in T.itertuples():
        d0 = pd.Timestamp(r.date)
        t0day = d0.value / 1e9
        acc, n = None, 0
        for s in STAS:
            if s not in cc.index:
                continue
            nn, ee = G.ned_offset(r.lat, r.lon, cc.loc[s, 'lat'], cc.loc[s, 'lon'])
            tt = np.sqrt(nn ** 2 + ee ** 2 + (r.dep * 1e3) ** 2) / 1e3 / C.VS_KMS
            u0 = int(round(r.t - t0day + W0 + tt)); u1 = int(round(r.t - t0day + W1 + tt)) + 1
            seg = np.full(u1 - u0, np.nan)
            for dd in range(u0 // 86400, (u1 - 1) // 86400 + 1):        # windows can cross midnight
                k = (s, (d0 + pd.Timedelta(days=dd)).strftime('%Y-%m-%d'))
                if k not in cache:
                    cache[k] = env(*k)
                e = cache[k]
                if e is None:
                    continue
                a0, a1 = max(u0, dd * 86400), min(u1, (dd + 1) * 86400)
                seg[a0 - u0:a1 - u0] = e[a0 - dd * 86400:a1 - dd * 86400]
            if np.isfinite(seg).mean() < 0.9:
                continue
            acc = np.where(np.isfinite(seg), seg, 0.0) if acc is None else acc + np.where(np.isfinite(seg), seg, 0.0)
            n += 1
        if n < 2:
            rows.append(dict(id=r.id, n_env=n)); continue
        net = acc / n
        k = int(np.argmax(net))
        rows.append(dict(id=r.id, n_env=n, to=r.t + W0 + k, env_peak=float(net[k]), env_med=float(np.median(net))))
        if len(cache) > 400:
            cache.clear()
    O = T.merge(pd.DataFrame(rows), on='id', how='left')
    if a.test_date:
        O.to_csv(a.test_out, index=False); print(O[['time', 'n_env', 'to', 'env_peak', 'env_med']].head(12).to_string()); return
    O.to_csv(f'{C.TB}/be3_pnsn_origins_all.csv', index=False)
    O = O[O.n_env >= 2]
    O.to_csv(f'{C.TB}/be3_pnsn_origins.csv', index=False)
    print('PNSN noLFE epochs', len(T), 'with >= 2 envelopes', len(O), '; origin - epoch time (s) quantiles',
          np.round(np.quantile(O.to - O.t, [.05, .25, .5, .75, .95]), 1).tolist())


if __name__ == '__main__':
    main()
