#!/usr/bin/env python
"""Merlin prescription step 1 (DECISIVE for PGC dv/v): the 8-year stack test (T1) passed, but dv/v runs on
30-CAL-DAY rolling stacks (FINAL_PIPELINE step 6/8). Does the coda survive at 30-day resolution at the
surface? For each of the 70 T1-surviving families, bin its B011 cc>=0.85 detection times into calendar
30-day blocks, stack PGC BHZ per block, and per block compute (a) causality ratio RMS(coda 2-4s)/RMS(mirror
-2..0s) and (b) correlation of the block coda (2-4s) vs the family's all-time coda (the dv/v cc_max analog).
Acceptance (Merlin): >=20 families for which >=50% of ETS-SEASON blocks have ratio>1.5 AND coda-cc>=0.6.
Failure: <10 families qualify, or coda-cc collapses outside ETS peaks -> PGC dv/v episodic-only or not viable.
"""
import os, numpy as np, pandas as pd
from obspy import read, UTCDateTime
from scipy.signal import butter, sosfiltfilt
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
        if int(round(tr.stats.sampling_rate)) != 40:
            return out
        x = tr.data.astype(float); xf = sosfiltfilt(SOS, x - x.mean()); t0 = tr.stats.starttime
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
    surv = set(pd.read_csv('data/pgc_stack_test.csv').query('ratio>1.5').fam.astype(str))   # 70 T1 families
    print(f"[30d] {len(surv)} T1-surviving families; collecting B011 cc>=0.85 times...", flush=True)
    times = {}
    for y in range(2010, 2018):
        f = f'data/mf_b011p90f40_{y}.csv'
        if not os.path.exists(f):
            continue
        for ch in pd.read_csv(f, usecols=['template', 'time', 'cc'], chunksize=4_000_000):
            ch['template'] = ch.template.astype(str)
            ch = ch[(ch.template.isin(surv)) & (ch.cc >= 0.85)]
            for fam, g in ch.groupby('template'):
                times.setdefault(fam, []).extend(pd.to_datetime(g.time).astype('int64').values / 1e9)
        print(f"  {y}", flush=True)

    # assign each detection a calendar 30-day block index; build harvest jobs (cap per family-block)
    rng = np.random.RandomState(3)
    jobs_byfile = {}
    for fam, eps in times.items():
        eps = np.array(eps); blk = (eps // (30*86400)).astype(np.int64)
        df = pd.DataFrame({'ep': eps, 'blk': blk})
        for b, g in df.groupby('blk'):
            e = g.ep.values
            if len(e) > 300:
                e = rng.choice(e, 300, replace=False)
            for ep in e:
                t = pd.Timestamp(ep, unit='s')
                jobs_byfile.setdefault(f'data/waveforms/CN.PGC/{t.year}/{t.dayofyear:03d}.mseed', []).append((fam, int(b), float(ep)))
    jobs_byfile = {f: v for f, v in jobs_byfile.items() if os.path.exists(f)}
    print(f"[30d] harvesting {sum(len(v) for v in jobs_byfile.values())} windows from {len(jobs_byfile)} day-files...", flush=True)

    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    acc = {}   # (fam,blk) -> [windows]
    with ProcessPoolExecutor(max_workers=12, mp_context=mp.get_context('spawn')) as ex:
        for res in ex.map(cut_day, list(jobs_byfile.items()), chunksize=8):
            for fam, blk, w in res:
                acc.setdefault((fam, blk), []).append(w)

    # per family: all-time coda reference, then per-block ratio + coda-cc
    fam_ref = {}
    for (fam, blk), ws in acc.items():
        fam_ref.setdefault(fam, []).extend(ws)
    ref_coda = {fam: trim_mean(np.asarray(ws), 0.02, axis=0)[CODA] for fam, ws in fam_ref.items() if len(ws) >= 30}

    rows = []
    for (fam, blk), ws in acc.items():
        if len(ws) < 20 or fam not in ref_coda:
            continue
        stk = trim_mean(np.asarray(ws), 0.02, axis=0)
        rm = np.sqrt(np.mean(stk[MIR]**2))
        ratio = np.sqrt(np.mean(stk[CODA]**2)) / rm if rm > 0 else np.nan
        cc = np.corrcoef(stk[CODA], ref_coda[fam])[0, 1]
        rows.append((fam, blk, len(ws), round(float(ratio), 2), round(float(cc), 2)))
    R = pd.DataFrame(rows, columns=['fam', 'blk', 'n', 'ratio', 'coda_cc']).dropna()
    R.to_csv('data/pgc_30day_stack_test.csv', index=False)

    # a block "qualifies" if ratio>1.5 AND coda_cc>=0.6; family qualifies if >=50% of its blocks (with >=20 det) qualify
    R['ok'] = (R.ratio > 1.5) & (R.coda_cc >= 0.6)
    fam_ok = R.groupby('fam').ok.mean()
    nfam = int((fam_ok >= 0.5).sum())
    print(f"\n[RESULT] PGC 30-day-stack test: {R.fam.nunique()} families with usable blocks; "
          f"{nfam} families have >=50% of 30-day blocks passing (ratio>1.5 & coda_cc>=0.6)")
    print(f"  blocks: {len(R)} total, {int(R.ok.sum())} pass ({100*R.ok.mean():.0f}%)")
    print(f"  median block: ratio {R.ratio.median():.2f}, coda_cc {R.coda_cc.median():.2f}, n {int(R.n.median())}")
    print(f"  ACCEPTANCE >=20 families -> PGC dv/v viable year-round; 10-20 -> episodic; <10 -> not viable")


if __name__ == '__main__':
    main()
