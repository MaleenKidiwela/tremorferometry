#!/usr/bin/env python
"""rf06_align_check (2026-10-14, quick check on existing outputs, not a queue stage): for the episodes that have a local
tremor time in tables/slipfront_<ets>.csv (2009-06+, tremor within 30 km), recompute Delta CC_SV (15-day post - pre)
with the same RFs aligned on (a) the GPS midpoint, (b) the median local tremor time, (c) the first local tremor time
(2nd percentile). Scrambled control: CC values shuffled among all RFs of the run (4000 draws).
Output: tables/align_check_15d.csv
"""
import pandas as pd, numpy as np
W="/home/jovyan/tremorferometry/svi_science/results/s12_lvl/rf_gosselin"
rng=np.random.default_rng(1)
ROWS=[]
def delta(x,r,Wd=15):
    pre=(r>=-Wd)&(r<0); post=(r>=0)&(r<Wd)
    return x[post].mean()-x[pre].mean(), pre.sum(), post.sum()
for tag,sta,ets,sf in [('full','CN.PGC','SVI_ALBH','SVI_ALBH'),('full','CN.VGZ','SVI_ALBH','SVI_ALBH'),('full','CN.SNB','SVI_ALBH','SVI_ALBH'),
                       ('ext_svi','CN.NLLB','SVI_ALBH','SVI_ALBH'),('ext_svi','CN.YOUB','SVI_ALBH','SVI_ALBH'),('ext_svi','CN.GOBB','SVI_ALBH','SVI_ALBH'),('ext_svi','CN.MGB','SVI_ALBH','SVI_ALBH'),
                       ('ext_wa','UW.GNW','WA_SEAT','WA_SEAT'),('ext_wa','UW.SP2','WA_SEAT','WA_SEAT')]:
    M=pd.read_csv(f"{W}/tables/{tag}/{sta}_events.csv",parse_dates=['t']); M=M[np.isfinite(M.cc_sv_clean)]
    E=pd.read_csv(f"{W}/tables/ets_mid_{ets}.csv",parse_dates=['mid'])
    F=pd.read_csv(f"{W}/tables/slipfront_{sf}.csv",parse_dates=['mid']); F=F[(F.station==sta)&F.median_d.notna()]
    mids=F.mid.values.astype('datetime64[ns]')
    tt=M.t.values.astype('datetime64[ns]')
    # all midpoints of the run for nearest assignment, then keep RFs whose nearest midpoint has a front value
    allm=np.sort(E.mid.values.astype('datetime64[ns]'))
    rel=(tt[:,None]-allm[None,:])/np.timedelta64(1,'D'); j=np.argmin(np.abs(rel),1); r=rel[np.arange(len(tt)),j]
    near=allm[j]
    off_med=dict(zip(F.mid.values.astype('datetime64[ns]'),F.median_d)); off_first=dict(zip(F.mid.values.astype('datetime64[ns]'),F.first_d))
    keep=np.array([m in off_med for m in near])
    x=M.cc_sv_clean.values
    xs=x[keep]; rm=r[keep]; rmed=rm-np.array([off_med[m] for m in near[keep]]); rfst=rm-np.array([off_first[m] for m in near[keep]])
    out=[]; rec=dict(station=sta, tag=tag, n_episodes=len(F))
    for lab,rr in (('midpoint',rm),('median tremor',rmed),('first tremor',rfst)):
        d,npre,npost=delta(xs,rr)
        B=np.array([delta(rng.permutation(x)[keep],rr)[0] for _ in range(4000)])
        out.append(f"{lab}: {d:+.3f} (pct {np.mean(B<=d):.3f}, n {npre}/{npost})")
        k=lab.split()[0]; rec.update({f"d_{k}":round(d,4), f"pct_{k}":round(float(np.mean(B<=d)),4), f"npre_{k}":int(npre), f"npost_{k}":int(npost)})
    ROWS.append(rec)
    print(f"{sta:8s} episodes {len(F)}  " + " | ".join(out))

pd.DataFrame(ROWS).to_csv(f"{W}/tables/align_check_15d.csv", index=False)
