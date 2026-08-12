#!/usr/bin/env python
"""Family Trust Test — Tier 0 + Tier 2b, margin-wide, from EXISTING files (no waveforms needed).
Per (station, family): detection-list flags (rate, cap saturation, day/night, cc-shape, local-tremor corr),
template shape (centroid/kurtosis + station-health gate), and product-level coda coherence (2-4s cc).
Output: data/family_trust_provisional.csv + per-station summary. Flags + provisional class only —
Tier 1 (stack-vs-random, day/night split) upgrades these to verdicts where waveforms exist.
"""
import glob, os, re, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.stats import kurtosis as kurt_fn

STA = {'B927':(49.2188,-124.8113),'NLLB':(49.2271,-123.9882),'B928':(48.834,-125.134),'PGC':(48.6498,-123.4521),
 'B011':(48.65,-123.448),'B004':(48.202,-124.427),'B013':(47.813,-122.9108),'HDW':(47.6490,-123.0530),
 'GNW':(47.5641,-122.8250),'B014':(47.5133,-123.8125),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B020':(46.3827,-123.8445),'B201':(46.3033,-122.2648),'B204':(46.136,-122.169),'B023':(46.1112,-123.0787),
 'B022':(45.9546,-123.931),'B026':(45.3094,-123.8231),'COLT':(45.17044,-122.438152),'COR':(44.5855,-123.3046),
 'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),
 'B036':(42.5058,-123.3817),'B040':(41.8308,-122.4205),'B039':(41.4667,-122.4847),'B935':(40.4787,-123.5732),
 'B017':(46.9960,-123.5575),'B001':(48.0431,-123.1314),'B005':(48.0596,-123.5034),'B003':(48.0623,-124.1416),
 'B045':(40.4360,-124.0008),'B932':(40.2825,-124.2245),'B049':(40.2403,-123.8225)}

# local tremor daily series per station (60 km)
cat = pd.read_csv('catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time','lat','lon'])
cat['t'] = pd.to_datetime(cat.time, errors='coerce'); cat = cat.dropna(subset=['t'])

def tremor_daily(la, lo):
    dkm = np.hypot((cat.lat-la)*111, (cat.lon-lo)*111*np.cos(np.radians(la)))
    sub = cat[dkm < 60]
    return sub.groupby(sub.t.dt.floor('D')).size()

def template_feats(sta):
    s = sta.lower()
    f = f'data/{s}_pnsn_families_100km.npz'
    if not os.path.exists(f): return {}
    z = np.load(f, allow_pickle=True)
    out = {}
    for fam in z.files:
        w = np.asarray(z[fam], float); w = w - w.mean()
        if w.std() == 0: continue
        W = np.abs(np.fft.rfft(w)); fr = np.fft.rfftfreq(len(w), 1/40.)
        out[fam] = ((W*fr).sum()/W.sum(), float(kurt_fn(w, fisher=True)))
    return out

rows = []
for sta, (la, lo) in STA.items():
    mf = f'data/mf_{sta.lower()}_all.csv'
    if not os.path.exists(mf):
        print(f'[{sta}] no mf_all — skip', flush=True); continue
    print(f'[{sta}] loading detections...', flush=True)
    m = pd.read_csv(mf, usecols=['template','time','cc'])
    m['t'] = pd.to_datetime(m.time); m['hr'] = m.t.dt.hour; m['day'] = m.t.dt.floor('D')
    tr = tremor_daily(la, lo)
    tf = template_feats(sta)
    # product-level coda cc from the 2-4s cal product
    try:
        dv = pd.read_csv(f'data/daily_dvv_{sta}_2to4_cal.csv', usecols=['patch','cc_max'])
        codacc = dv.groupby('patch').cc_max.median()
    except Exception:
        codacc = pd.Series(dtype=float)
    # station-health gate for templates
    cens = pd.Series({k: v[0] for k, v in tf.items()})
    cen_cut = 4.3 if (len(cens) and cens.median() > 4.3) else max(4.3, cens.quantile(0.85) if len(cens) else 4.3)
    for fam, g in m.groupby('template'):
        dd = g.groupby('day').size()
        span = pd.date_range(dd.index.min(), dd.index.max(), freq='D')
        a = dd.reindex(span, fill_value=0); b = tr.reindex(span, fill_value=0)
        rcorr = float(np.corrcoef(a, b)[0,1]) if (a.std() > 0 and b.std() > 0) else np.nan
        dayn = int(g.hr.between(12,21).sum()); nightn = int(g.hr.between(0,9).sum())
        cen, kt = tf.get(fam, (np.nan, np.nan))
        rows.append(dict(station=sta, fam=fam, n=len(g), rate=len(g)/len(dd),
            cap_frac=float((dd >= 100).mean()), dn=dayn/max(nightn,1),
            cc_med=float(g.cc.median()), cc_hi_frac=float((g.cc > 0.92).mean()),
            tremor_r=rcorr, tpl_cen=cen, tpl_kurt=kt,
            tpl_spiky=bool(cen > cen_cut and kt > 4) if np.isfinite(cen) else None,
            coda_cc=float(codacc.get(fam, np.nan))))
    print(f'[{sta}] {m.template.nunique()} families done', flush=True)

A = pd.DataFrame(rows)
# provisional class (FLAGS, not verdicts; Tier-1 upgrades)
sta_coda_med = A.groupby('station').coda_cc.transform('median')
A['flag_codacollapse'] = A.coda_cc < (sta_coda_med - 0.15)
A['flag_spiky'] = A.tpl_spiky == True
A['flag_diurnal'] = A.dn > 1.5
A['flag_capsat'] = A.cap_frac > 0.8
A['flag_cchug'] = (A.cc_med < 0.87) & (A.cc_hi_frac < 0.02)
A['n_flags'] = A[['flag_codacollapse','flag_spiky','flag_diurnal','flag_capsat','flag_cchug']].sum(axis=1)
A['provisional'] = np.select(
    [A.flag_codacollapse | A.flag_spiky, A.n_flags >= 3, A.n_flags <= 1],
    ['SUSPECT-strong', 'SUSPECT', 'CLEANISH'], default='MIXED')
A.to_csv('data/family_trust_provisional.csv', index=False)
print('\n=== PER-STATION SUMMARY ===')
s = A.groupby('station').agg(fams=('fam','size'), cleanish=('provisional', lambda x:(x=='CLEANISH').sum()),
    suspect=('provisional', lambda x: x.str.startswith('SUSPECT').sum()),
    dn_med=('dn','median'), cap_med=('cap_frac','median'), tremor_med=('tremor_r','median'),
    codacc_med=('coda_cc','median'))
print(s.round(2).to_string())
print(f'\nTOTAL: {len(A)} families | {(A.provisional=="CLEANISH").sum()} CLEANISH | '
      f'{A.provisional.str.startswith("SUSPECT").sum()} SUSPECT | -> data/family_trust_provisional.csv')
print('TIER0_DONE')
