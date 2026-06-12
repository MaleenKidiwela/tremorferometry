#!/usr/bin/env python
"""Family Trust Test — Tier 1 (decisive waveform battery) for ONE station with waveforms on disk.
v2 after adversarial review (notes/2026-06-11): glitch-window QC + robust stack, nwin-matched null,
nboot>=200, per-window hour retained -> TRUE day/night T1b, tri-state verdicts.

Per family:
  T1a stack-vs-random : robust (trimmed) stack of its sampled detections vs an nwin-MATCHED bootstrap
                        null of random-time stacks; score = coda(+2..+4 s, outside the 2-s selection
                        window) amplitude in null-sigma units.
  T1b day/night       : local-day vs local-night sub-stacks; coda cross-correlation (source-independence).
Calibration: T1c reversed-template fakes through the same scoring -> per-station FAIL distribution.
Verdict: TRUSTED (> fake 95th, glitch-clean, daynight ok) / FAIL (< fake median) / UNDETERMINED (between).
"""
import os, sys, glob, argparse
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import numpy as np, pandas as pd

FS = 40.0
PRE, POST = 3.0, 11.0
NW = int((PRE + POST) * FS)
S_SL = slice(int(2.0*FS), int(4.0*FS))    # t = -1..+1 (selection window)
C_SL = slice(int(5.0*FS), int(7.0*FS))    # t = +2..+4 (coda, outside selection)
WMAX = 15.0                                # reject windows with |w|>WMAX after S-window normalization
ZRUN = int(0.5*FS)                         # reject windows containing >=0.5 s of exact zeros (gap fill)


def _harvest(args):
    path, jobs = args
    import obspy
    from scipy.signal import resample_poly, butter, sosfiltfilt
    from math import gcd
    out = []
    try:
        st = obspy.read(path); st.merge(fill_value=0); tr = st[0]
        sr = int(round(tr.stats.sampling_rate))
        x = tr.data.astype(float)
        if sr != 40:
            g = gcd(sr, 40); x = resample_poly(x, 40//g, sr//g)
        from scipy.signal import butter as _b
        sos = butter(4, [2/20., 8/20.], btype='band', output='sos')
        xf = sosfiltfilt(sos, x - x.mean())
        t0 = tr.stats.starttime
        for key, ts in jobs:
            tt = obspy.UTCDateTime(pd.Timestamp(ts))
            i0 = int((tt - t0) * 40) - int(PRE*FS)
            if i0 < 0 or i0 + NW > len(x): continue
            raw = x[i0:i0+NW]
            # gap QC on the UNfiltered segment: runs of exact zeros
            z = (raw == 0)
            if z.any():
                runs = np.diff(np.flatnonzero(np.diff(np.r_[0, z.view(np.int8), 0])) .reshape(-1, 2), axis=1)
                if len(runs) and runs.max() >= ZRUN: continue
            w = xf[i0:i0+NW]
            s = w[S_SL].std()
            if s <= 0: continue
            w = w / s
            if np.abs(w).max() > WMAX: continue      # glitch amplitude QC
            out.append((key, pd.Timestamp(ts).hour, w.astype(np.float32)))
    except Exception:
        pass
    return out


def collect(net, sta, jobs, workers):
    byfile = {}
    for key, ts in jobs:
        ts = pd.Timestamp(ts)
        f = f"data/waveforms/{net}.{sta}/{ts.year}/{ts.dayofyear:03d}.mseed"
        byfile.setdefault(f, []).append((key, ts))
    byfile = {f: v for f, v in byfile.items() if os.path.exists(f)}
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    acc = {}
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        for res in ex.map(_harvest, list(byfile.items()), chunksize=4):
            for key, hr, w in res:
                acc.setdefault(key, []).append((hr, w))
    return acc


def tstack(ws):
    """Robust stack: 2% trimmed mean along axis 0."""
    from scipy.stats import trim_mean
    return trim_mean(np.asarray(ws), 0.02, axis=0)


def coda_amp(stack): return float(np.abs(stack[C_SL]).max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--net', required=True); ap.add_argument('--sta', required=True)
    ap.add_argument('--mf', default=None); ap.add_argument('--rev-mf', default=None)
    ap.add_argument('--nsamp', type=int, default=250); ap.add_argument('--nrand', type=int, default=2500)
    ap.add_argument('--nboot', type=int, default=200); ap.add_argument('--workers', type=int, default=16)
    a = ap.parse_args()
    s = a.sta.lower()
    rng = np.random.RandomState(42)

    print(f'[{a.sta}] loading detections...', flush=True)
    m = pd.read_csv(a.mf or f'data/mf_{s}_all.csv', usecols=['template','time'])
    m['t'] = pd.to_datetime(m.time)
    jobs = []
    for fam, g in m.groupby('template'):
        for ts in g.sample(min(a.nsamp, len(g)), random_state=1).t: jobs.append((('F', fam), ts))
    rfams = 0
    if a.rev_mf:
        rv = pd.concat([pd.read_csv(f, usecols=['template','time']) for f in sorted(glob.glob(a.rev_mf))],
                       ignore_index=True)
        rv['t'] = pd.to_datetime(rv.time); rfams = rv.template.nunique()
        for fam, g in rv.groupby('template'):
            for ts in g.sample(min(a.nsamp, len(g)), random_state=2).t: jobs.append((('R', fam), ts))
    t0, t1 = m.t.min(), m.t.max()
    for ts in pd.to_datetime(rng.uniform(t0.value, t1.value, a.nrand).astype('int64')):
        jobs.append((('N', 'rand'), ts))
    print(f'[{a.sta}] harvesting {len(jobs)} windows...', flush=True)
    acc = collect(a.net, a.sta, jobs, a.workers)
    rand = np.asarray([w for _, w in acc.get(('N','rand'), [])])
    print(f'[{a.sta}] random windows kept: {len(rand)} (QC-passed)', flush=True)
    if len(rand) < 500: sys.exit('not enough random windows')

    def null_for(nwin):
        boots = np.empty(a.nboot)
        for k in range(a.nboot):
            idx = rng.randint(0, len(rand), nwin)
            boots[k] = coda_amp(tstack(rand[idx]))
        return boots.mean(), boots.std()

    null_cache = {}
    rows = []
    for kind, fam in [k for k in acc if k[0] in ('F','R')]:
        pairs = acc[(kind, fam)]
        if len(pairs) < 30: continue
        hrs = np.array([h for h, _ in pairs]); A = np.asarray([w for _, w in pairs])
        nb = min(len(A) // 10 * 10, 250)            # bucket nwin for null cache
        if nb not in null_cache: null_cache[nb] = null_for(nb)
        nmu, nsd = null_cache[nb]
        st = tstack(A)
        # TRUE day/night T1b (local day 15-23 UTC / night 5-13 UTC)
        dmask = (hrs >= 15) & (hrs <= 23); nmask = (hrs >= 5) & (hrs <= 13)
        dn_cc = np.nan
        if dmask.sum() >= 15 and nmask.sum() >= 15:
            sd_, sn_ = tstack(A[dmask]), tstack(A[nmask])
            dn_cc = float(np.corrcoef(sd_[C_SL], sn_[C_SL])[0,1])
        rows.append(dict(kind=kind, fam=fam, nwin=len(A),
                         coda=coda_amp(st), coda_sigma=(coda_amp(st)-nmu)/nsd,
                         s_amp=float(np.abs(st[S_SL]).max()), daynight_cc=dn_cc))
    R = pd.DataFrame(rows)
    rv_ = R[R.kind=='R']
    if len(rv_):
        f50, f95 = rv_.coda_sigma.median(), rv_.coda_sigma.quantile(0.95)
    else:
        f50, f95 = 1.0, 6.0
    def verdict(r):
        if r.kind != 'F': return ''
        if r.coda_sigma > f95 and (np.isnan(r.daynight_cc) or r.daynight_cc > 0.3): return 'TRUSTED'
        if r.coda_sigma < f50: return 'FAIL'
        return 'UNDETERMINED'
    R['verdict'] = R.apply(verdict, axis=1)
    out = f'data/family_trust_tier1_{a.sta}.csv'
    R.to_csv(out, index=False)
    rr = R[R.kind=='F']
    print(f'\n[{a.sta}] REAL fams (n={len(rr)}): coda_sigma median {rr.coda_sigma.median():.1f} '
          f'max {rr.coda_sigma.max():.1f} | daynight_cc median {rr.daynight_cc.median():.2f}')
    if len(rv_):
        print(f'[{a.sta}] FAKES (n={len(rv_)}): median {f50:.1f}, 95th {f95:.1f}, max {rv_.coda_sigma.max():.1f}')
    print(f'[{a.sta}] verdicts: {rr.verdict.value_counts().to_dict()}')
    print(f'-> {out}')
    print('TIER1_DONE')

if __name__ == '__main__':
    main()
