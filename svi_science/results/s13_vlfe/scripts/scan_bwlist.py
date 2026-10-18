#!/usr/bin/env python
"""VLFE scan pilot 2 (PLAN_SCAN.md S7.4): reconstruct the known body-wave list among the 47 pilot-1 survivors
(real, variant 0, headline mask, S5 filter, z >= 7, z_loso >= 7). A survivor is a known body-wave cluster if its time is in
[t_P - 60 s, t_S + 600 s] of a ComCat M >= 5.0 event at D >= 1,500 km (bw5 windows of tables/scan2_eq_events.csv).
Writes tables/scan2_survivors47.csv (all 47 with the matching event, if any) and tables/scan2_bw_known_list.csv."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S

D = pd.read_parquet(f'{S.SD}/runs/real/loso.parquet')
D = D[(D.variant == 0) & (D.z_loso >= 7)].copy()
E = pd.read_csv(f'{C.TB}/scan2_eq_events.csv')
E = E[E.bw5_t0.notna()].drop_duplicates(subset=['t', 'mag']).copy()
y0 = S.Y0.value / 1e9
rows = []
for i, r in D.iterrows():
    t = y0 + r.tsec
    m = E[(E.bw5_t0 <= t) & (E.bw5_t1 >= t)]
    rec = dict(tsec=int(r.tsec), time_utc=str(r.time_utc), time_pt=r.time_pt, node=int(r.node), lat=r.lat, lon=r.lon,
               z=r.z, z_loso=r.z_loso, mw=r.mw, loso_dropped=r.loso_dropped, max_share_ch=r.max_share_ch,
               in_mask_m45=bool(r.in_mask_m45), n_bw_events=len(m))
    if len(m):
        b = m.sort_values('mag', ascending=False).iloc[0]
        rec.update(bw_event_time=b.time, bw_event_mag=b.mag, bw_event_magtype=b.magType, bw_event_dist_km=b.dist_km,
                   bw_event_depth=b.depth, bw_event_lat=b.lat, bw_event_lon=b.lon, t_after_P_s=t - (b.t + b.tP_s),
                   t_after_S_s=t - (b.t + b.tS_s), bw_event_id=b.id)
    rows.append(rec)
A = pd.DataFrame(rows).sort_values('tsec')
A['known_body_wave'] = A.n_bw_events > 0
A.to_csv(f'{C.TB}/scan2_survivors47.csv', index=False)
L = A[A.known_body_wave]
L.to_csv(f'{C.TB}/scan2_bw_known_list.csv', index=False)
print(len(A), 'survivors;', len(L), 'known body-wave clusters (Merlin by eye: 11)')
print(L[['time_utc', 'z', 'z_loso', 'lat', 'lon', 'bw_event_time', 'bw_event_mag', 'bw_event_dist_km', 't_after_P_s', 't_after_S_s']].round(1).to_string())
