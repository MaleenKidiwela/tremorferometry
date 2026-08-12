#!/usr/bin/env python
"""Broadband fingerprint training data via CO-LOCATED labels (Merlin plan, 2026-07-14).
PGC (BHZ 40 Hz) is ~290 m from borehole B011. B011's causality-certified LFE families give precise,
abundant detection times that are TRUE labels for PGC's surface record (no circularity: labels come from
the borehole waveforms, features from PGC). Extract Stage-1 features at:
  - COLOC : B011 certified-family high-cc detections -> PGC windows (positives)
  - TNOISE: random times on tremor-ACTIVE days, >=60 s from any certified detection -> PGC windows
            (the hard negative that decides whether the picker works DURING tremor)
Window geometry matches feat_ref_pgc (PRE=10,POST=30, S~anchor, t_pre=10) so all sets combine.
Output: lfe_features/data/feat_coloc_pgc.parquet (labels COLOC / TNOISE).
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from feature_defs import compute_all, FEATURE_ORDER

PRE, POST = 10.0, 30.0
NET, STA, FS = 'CN', 'PGC', 40.0


def cut_day(item):
    from obspy import read, UTCDateTime
    f, jobs = item          # jobs: list of (label, fam, epoch)
    out = []
    try:
        st = read(f).select(component='Z'); st.merge(fill_value=0); tr = st[0]
        if int(round(tr.stats.sampling_rate)) != 40:
            return out
        for label, fam, ep in jobs:
            a = UTCDateTime(ep)
            z = tr.slice(a - PRE, a + POST).data.astype(float)
            if len(z) < int((PRE + POST) * FS) - 2:
                continue
            try:
                feats = compute_all(z, FS, PRE)
            except Exception:
                continue
            feats['label'] = label; feats['fam'] = fam
            feats['year'] = a.year; feats['time'] = ep
            out.append(feats)
    except Exception:
        pass
    return out


def main():
    cert = set(pd.read_csv('data/b011_fwd_vs_rev_coda.csv').query('ratio>1.5').fam.astype(str))
    print(f"[coloc] {len(cert)} certified B011 families; collecting positives + building TNOISE...", flush=True)

    rng = np.random.RandomState(7)
    # POSITIVES = per-family HIGH-CONFIDENCE detections (cc>=0.92 = cleaner individual-window LFE labels
    # than the 0.85 densify threshold, which is noise-dominated at the window level).
    CC_POS = 0.92
    pos = {}             # fam -> np.array of epochs (reservoir for COLOC positives)
    pos_ep_all = []      # ALL cc>=CC_POS epochs (for fine ±5 s TNOISE exclusion)
    for y in range(2010, 2018):
        f = f'data/mf_b011p90f40_{y}.csv'
        if not os.path.exists(f):
            continue
        for ch in pd.read_csv(f, usecols=['template', 'time', 'cc'], chunksize=4_000_000):
            ch['template'] = ch.template.astype(str)
            ch = ch[(ch.template.isin(cert)) & (ch.cc >= CC_POS)]
            if not len(ch):
                continue
            ep = pd.to_datetime(ch.time).astype('int64').values / 1e9
            pos_ep_all.append(ep)
            for fam, g in ch.groupby('template'):
                e = pd.to_datetime(g.time).astype('int64').values / 1e9
                cur = pos.get(fam)
                allv = e if cur is None else np.concatenate([cur, e])
                pos[fam] = allv if len(allv) <= 4000 else rng.choice(allv, 4000, replace=False)
        print(f"  {y} done ({len(pos)} fams)", flush=True)
    pos_ep_all = np.sort(np.concatenate(pos_ep_all))

    jobs_byfile = {}
    for fam, eps in pos.items():
        if len(eps) > 400:
            eps = rng.choice(eps, 400, replace=False)
        for ep in eps:
            t = pd.Timestamp(ep, unit='s')
            jobs_byfile.setdefault(f'data/waveforms/{NET}.{STA}/{t.year}/{t.dayofyear:03d}.mseed', []).append(('COLOC', fam, float(ep)))
    npos = sum(len(v) for v in jobs_byfile.values())

    # TNOISE = genuine IN-TREMOR noise: random times inside PNSN tremor windows near B011, >=60 s from a
    # real LFE positive. This is the HARD negative (the continuous rumble during ETS, not a discrete LFE).
    trem = pd.read_csv('catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time', 'lat', 'lon'])
    trem = trem[(trem.lat.between(48.15, 49.15)) & (trem.lon.between(-124.1, -122.8))]
    tep = pd.to_datetime(trem.time).astype('int64').values / 1e9
    tep = tep[(tep >= pd.Timestamp('2010-01-01').timestamp()) & (tep < pd.Timestamp('2018-01-01').timestamp())]
    rng.shuffle(tep)
    ntn = 0
    for te in tep:
        if ntn >= npos:
            break
        cand = te + rng.uniform(-120, 120)      # random moment within the tremor window
        i = np.searchsorted(pos_ep_all, cand)   # exclude only if an LFE is at the window CENTER (±5 s)
        near = min(abs(pos_ep_all[min(i, len(pos_ep_all)-1)] - cand),
                   abs(pos_ep_all[max(0, i-1)] - cand))
        if near < 5.0:
            continue
        t = pd.Timestamp(cand, unit='s')
        jobs_byfile.setdefault(f'data/waveforms/{NET}.{STA}/{t.year}/{t.dayofyear:03d}.mseed', []).append(('TNOISE', 'trem', float(cand)))
        ntn += 1
    print(f"[coloc] placed {ntn} TNOISE (in-tremor, >5 s from any LFE center)", flush=True)
    jobs_byfile = {f: v for f, v in jobs_byfile.items() if os.path.exists(f)}
    print(f"[coloc] harvesting {npos} COLOC + {ntn} TNOISE windows from {len(jobs_byfile)} PGC day-files...", flush=True)

    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    rows = []
    with ProcessPoolExecutor(max_workers=12, mp_context=mp.get_context('spawn')) as ex:
        for res in ex.map(cut_day, list(jobs_byfile.items()), chunksize=8):
            rows.extend(res)
    df = pd.DataFrame(rows)
    df.to_parquet('lfe_features/data/feat_coloc_pgc.parquet', index=False)
    print(f"[coloc] wrote feat_coloc_pgc.parquet: {df.shape}, labels {df.label.value_counts().to_dict()}")


if __name__ == '__main__':
    main()
