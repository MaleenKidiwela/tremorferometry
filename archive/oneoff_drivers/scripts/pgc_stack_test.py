#!/usr/bin/env python
"""Merlin step 1 (decisive): does the LFE CODA survive at the surface (PGC)? Detection at PGC is already
proven (RESULTS.md R5 AUC 0.978); dv/v needs the coda. Stack PGC BHZ at B011's certified-family detection
times (B011 & PGC are ~290 m apart -> negligible lag) and measure the causality ratio RMS(coda 2-4s)/
RMS(mirror -2..0s) at the surface. Acceptance: >=30% of families reach ratio>1.5 (borehole-like) -> proceed;
<10-15% -> broadband coda dv/v not viable, PGC is detection-only. Uses on-disk data only."""
import os, numpy as np, pandas as pd
from obspy import read, UTCDateTime
from scipy.signal import butter, sosfiltfilt
from scipy.stats import trim_mean

FS = 40.0; PRE, POST = 5.0, 15.0; NW = int((PRE + POST) * FS)
SOS = butter(4, [2/20., 8/20.], btype='band', output='sos')
T_AXIS = np.arange(NW)/FS - PRE


def cut_day(item):
    """Read one PGC day-file ONCE, cut+normalize windows for all its (fam,time) jobs."""
    f, jobs = item
    out = []
    try:
        st = read(f).select(component='Z'); st.merge(fill_value=0); tr = st[0]
        if int(round(tr.stats.sampling_rate)) != 40:
            return out
        x = tr.data.astype(float); xf = sosfiltfilt(SOS, x - x.mean()); t0 = tr.stats.starttime
        for fam, t in jobs:
            i0 = int((UTCDateTime(t) - t0) * 40) - int(PRE*FS)
            if i0 < 0 or i0 + NW > len(x):
                continue
            w = xf[i0:i0+NW]; s = w.std()
            if s <= 0:
                continue
            w = w/s
            if np.abs(w).max() > 15:
                continue
            out.append((fam, w.astype(np.float32)))
    except Exception:
        pass
    return out


def main():
    cert = set(pd.read_csv('data/b011_fwd_vs_rev_coda.csv').query('ratio>1.5').fam.astype(str))
    print(f"[pgc-stack] {len(cert)} certified B011 families; collecting high-cc(>=0.85) times 2010-2017...", flush=True)
    times = {}
    for y in range(2010, 2018):
        f = f'data/mf_b011p90f40_{y}.csv'
        if not os.path.exists(f):
            continue
        for ch in pd.read_csv(f, usecols=['template', 'time', 'cc'], chunksize=4_000_000):
            ch['template'] = ch.template.astype(str)
            ch = ch[(ch.cc >= 0.85) & (ch.template.isin(cert))]
            for fam, g in ch.groupby('template'):
                times.setdefault(fam, []).extend(g.time.tolist())
        print(f"  {y} done ({len(times)} families)", flush=True)

    jobs_byfile = {}
    for fam, ts in times.items():
        ts = pd.to_datetime(ts)
        if len(ts) > 800:
            ts = pd.Series(ts).sample(800, random_state=1)
        for t in ts:
            t = pd.Timestamp(t)
            jobs_byfile.setdefault(f'data/waveforms/CN.PGC/{t.year}/{t.dayofyear:03d}.mseed', []).append((fam, t))
    jobs_byfile = {f: v for f, v in jobs_byfile.items() if os.path.exists(f)}
    print(f"[pgc-stack] harvesting {sum(len(v) for v in jobs_byfile.values())} windows from "
          f"{len(jobs_byfile)} PGC day-files (parallel)...", flush=True)

    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    acc = {}
    with ProcessPoolExecutor(max_workers=12, mp_context=mp.get_context('spawn')) as ex:
        for res in ex.map(cut_day, list(jobs_byfile.items()), chunksize=8):
            for fam, w in res:
                acc.setdefault(fam, []).append(w)
    print(f"[pgc-stack] harvested; stacking {len(acc)} families...", flush=True)

    rows = []
    for fam, wins in acc.items():
        if len(wins) < 30:
            continue
        stack = trim_mean(np.asarray(wins), 0.02, axis=0)
        # FIXED-ORIGIN windows (finalize_causality convention: detection origin at t=0, S ~ +1s):
        # coda 2-4s (post-S), mirror -2..0s (pre-arrival noise, well before the S).
        coda = (T_AXIS >= 2) & (T_AXIS <= 4); mir = (T_AXIS >= -2) & (T_AXIS < 0)
        tpk = float(T_AXIS[np.argmax(np.abs(stack) * ((T_AXIS >= -1) & (T_AXIS <= 2)))])  # arrival-lag diagnostic
        rms = lambda m: float(np.sqrt(np.mean(stack[m]**2)))
        rm = rms(mir)
        rows.append((fam, len(wins), round(tpk, 2), rms(coda)/rm if rm > 0 else np.nan))

    R = pd.DataFrame(rows, columns=['fam', 'n', 'tpk', 'ratio']).dropna()
    R.to_csv('data/pgc_stack_test.csv', index=False)
    n15 = int((R.ratio > 1.5).sum())
    print(f"\n[RESULT] PGC stack test: {len(R)} families stacked; {n15} reach causality ratio>1.5 "
          f"({100*n15/max(1,len(R)):.0f}%)")
    print(f"  ratio: p50={R.ratio.median():.2f} p75={R.ratio.quantile(.75):.2f} "
          f"p90={R.ratio.quantile(.9):.2f} max={R.ratio.max():.2f}")
    print(f"  arrival lag tpk: median {R.tpk.median():.2f}s (expect ~0 for co-located 290 m)")
    print(f"  ACCEPTANCE >=30% viable; <10-15% -> broadband coda dv/v NOT viable")


if __name__ == '__main__':
    main()
