#!/usr/bin/env python
"""Segment-2 physics gate (Merlin step 3): does the LFE coda survive at PGC in the HHZ (100 Hz) era
2018-2026, decimated 100->40 Hz? Same as T1/T3 but on the new instrument era. Stack PGC HHZ(->40) at B011
certified-family cc>=0.85 times 2018-2026; per family compute overall causality ratio (T1) and per 30-cal-day
block ratio + coda-cc vs the family all-time coda (T3). Accept (Merlin): T1 >=30% families ratio>1.5 AND
T3 >=20 families with >=50% of blocks passing (ratio>1.5 & coda_cc>=0.6). Fail -> HHZ-era surface coda
unusable -> ship PGC as BHZ-only 2010-2017 segment."""
import os, numpy as np, pandas as pd
from obspy import read, UTCDateTime
from scipy.signal import butter, sosfiltfilt, resample_poly
from scipy.stats import trim_mean

FS = 40.0; PRE, POST = 5.0, 15.0; NW = int((PRE + POST) * FS)
SOS = butter(4, [2/20., 8/20.], btype='band', output='sos')
T = np.arange(NW)/FS - PRE
CODA = (T >= 2) & (T <= 4); MIR = (T >= -2) & (T < 0)


def cut_day(item):
    f, jobs = item
    out = []
    try:
        st = read(f).select(component='Z'); st.merge(fill_value=0); tr = st[0]
        sr = int(round(tr.stats.sampling_rate))
        x = tr.data.astype(float)
        if sr == 100:
            x = resample_poly(x, 2, 5)            # 100 -> 40 Hz polyphase (anti-aliased) — the fleet standard
        elif sr != 40:
            return out
        xf = sosfiltfilt(SOS, x - x.mean()); t0 = tr.stats.starttime
        for fam, blk, ep in jobs:
            i0 = int((UTCDateTime(ep) - t0) * 40) - int(PRE*FS)
            if i0 < 0 or i0 + NW > len(x):
                continue
            w = xf[i0:i0+NW]; s = w.std()
            if s <= 0 or np.abs(w/s).max() > 15:
                continue
            out.append((fam, blk, (w/s).astype(np.float32)))
    except Exception:
        pass
    return out


def main():
    cert = set(pd.read_csv('data/b011_fwd_vs_rev_coda.csv').query('ratio>1.5').fam.astype(str))
    print(f"[hhz] {len(cert)} B011 certified families; collecting cc>=0.85 times 2018-2026...", flush=True)
    times = {}
    for y in range(2018, 2027):
        f = f'data/mf_b011p90f40_{y}.csv'
        if not os.path.exists(f):
            continue
        for ch in pd.read_csv(f, usecols=['template', 'time', 'cc'], chunksize=4_000_000):
            ch['template'] = ch.template.astype(str)
            ch = ch[(ch.template.isin(cert)) & (ch.cc >= 0.85)]
            for fam, g in ch.groupby('template'):
                times.setdefault(fam, []).extend(pd.to_datetime(g.time).astype('int64').values / 1e9)
        print(f"  {y}", flush=True)

    rng = np.random.RandomState(3)
    jobs_byfile = {}
    for fam, eps in times.items():
        eps = np.array(eps); blk = (eps // (30*86400)).astype(np.int64)
        for b, e in pd.DataFrame({'ep': eps, 'blk': blk}).groupby('blk'):
            ev = e.ep.values
            if len(ev) > 300:
                ev = rng.choice(ev, 300, replace=False)
            for ep in ev:
                t = pd.Timestamp(ep, unit='s')
                jobs_byfile.setdefault(f'data/waveforms/CN.PGC/{t.year}/{t.dayofyear:03d}.mseed', []).append((fam, int(b), float(ep)))
    jobs_byfile = {f: v for f, v in jobs_byfile.items() if os.path.exists(f)}
    print(f"[hhz] harvesting {sum(len(v) for v in jobs_byfile.values())} windows from {len(jobs_byfile)} day-files (decimating 100->40)...", flush=True)

    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    acc = {}
    with ProcessPoolExecutor(max_workers=12, mp_context=mp.get_context('spawn')) as ex:
        for res in ex.map(cut_day, list(jobs_byfile.items()), chunksize=8):
            for fam, blk, w in res:
                acc.setdefault((fam, blk), []).append(w)

    fam_all = {}
    for (fam, blk), ws in acc.items():
        fam_all.setdefault(fam, []).extend(ws)
    ref = {fam: trim_mean(np.asarray(ws), 0.02, axis=0) for fam, ws in fam_all.items() if len(ws) >= 30}

    # T1: overall per-family causality ratio
    t1 = []
    for fam, stk in ref.items():
        rm = np.sqrt(np.mean(stk[MIR]**2))
        t1.append((fam, np.sqrt(np.mean(stk[CODA]**2))/rm if rm > 0 else np.nan))
    T1 = pd.DataFrame(t1, columns=['fam', 'ratio']).dropna()

    # T3: per 30-day block
    rows = []
    for (fam, blk), ws in acc.items():
        if len(ws) < 20 or fam not in ref:
            continue
        stk = trim_mean(np.asarray(ws), 0.02, axis=0)
        rm = np.sqrt(np.mean(stk[MIR]**2))
        ratio = np.sqrt(np.mean(stk[CODA]**2))/rm if rm > 0 else np.nan
        cc = np.corrcoef(stk[CODA], ref[fam][CODA])[0, 1]
        rows.append((fam, blk, round(float(ratio), 2), round(float(cc), 2)))
    R = pd.DataFrame(rows, columns=['fam', 'blk', 'ratio', 'coda_cc']).dropna()
    R.to_csv('data/pgc_hhz_stack_test.csv', index=False)
    R['ok'] = (R.ratio > 1.5) & (R.coda_cc >= 0.6)
    nfam = int((R.groupby('fam').ok.mean() >= 0.5).sum())
    n1 = int((T1.ratio > 1.5).sum())

    print(f"\n[RESULT] PGC HHZ-era (2018-2026) segment-2 gate:")
    print(f"  T1 (overall): {n1}/{len(T1)} families causal ratio>1.5 ({100*n1/max(1,len(T1)):.0f}%)  [BHZ era: 59%]")
    print(f"  T3 (30-day): {nfam} families >=50% blocks pass (ratio>1.5 & coda_cc>=0.6)  [BHZ era: 43]")
    print(f"  median block ratio {R.ratio.median():.2f}, coda_cc {R.coda_cc.median():.2f}")
    print(f"  ACCEPT: T1>=30% AND T3>=20 -> HHZ era viable (2-segment PGC); else BHZ-only")


if __name__ == '__main__':
    main()
