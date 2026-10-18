#!/usr/bin/env python
"""VLFE scan pilot 2 (PLAN_SCAN.md S7.6, S7.4): mask tiers head2 (Ms = mb - 0.3), bw5 (body-wave windows of M >= 5.0 at
D >= 1,500 km: [t_P - 60 s, t_S + 600 s], iasp91 first P-type / S-type arrivals at (48.37, -123.61)), reg3 (ComCat
M >= 3.0 within 2,500 km: [origin, origin + D/2.5 km/s + 600 s]). Headline mask v2 = head2 | bw5 | reg | reg3 (reg3 only
if it covers <= 15 % of 2011). Writes tables/scan2_eq_events.csv, tables/scan2_mask_fraction.csv,
data/scan/masks2_2011.npz (bool per second of 2011: head, m45, t58, reg, head2, bw5, reg3, v1, v2)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S

LAT0, LON0 = 48.37, -123.61
Y0 = pd.Timestamp(f'{S.YEAR}-01-01').value / 1e9                    # SCAN_YEAR (Amendment S8; 2011 unchanged)
Y1 = pd.Timestamp(f'{S.YEAR + 1}-01-01').value / 1e9
OUT = f'{S.SD}/masks2_{S.YEAR}.npz'


def gc_km(la, lo):
    la0, lo0 = np.radians(LAT0), np.radians(LON0)
    la, lo = np.radians(la), np.radians(lo)
    return np.degrees(np.arccos(np.clip(np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0), -1, 1))) * 111.19


def tsec(s):
    return pd.to_datetime(s, utc=True, format='mixed').dt.tz_localize(None).astype('datetime64[ns]').astype('int64') / 1e9


def ps_times(dep, deg, model):
    dep = float(np.clip(dep if np.isfinite(dep) else 10.0, 0.0, 700.0))
    tp = ts = np.nan
    try:
        a = model.get_travel_times(source_depth_in_km=dep, distance_in_degree=deg, phase_list=['ttp'])
        if a:
            tp = min(x.time for x in a)
        b = model.get_travel_times(source_depth_in_km=dep, distance_in_degree=deg, phase_list=['tts'])
        if b:
            ts = min(x.time for x in b)
    except Exception:
        pass
    return tp, ts


def interval_mask(a, b, y0=Y0, y1=Y1):
    m = np.zeros(int(y1 - y0), bool)
    ok = np.isfinite(a) & np.isfinite(b) & (b > y0) & (a < y1)
    for s, e in zip(a[ok], b[ok]):
        m[int(max(s - y0, 0)):int(min(e - y0, len(m)))] = True
    return m


def build():
    from obspy.taup import TauPyModel
    model = TauPyModel('iasp91')
    E = pd.read_csv(f'{C.TB}/scan_eq_events.csv')
    E = E[(E.t > Y0 - 2 * 86400) & (E.t < Y1)].copy()
    mb = E.magType.astype(str).str.lower().str.startswith('mb')
    E['ms_est2'] = np.where(mb, E.mag - 0.3, E.mag)
    dd = np.clip(E.dist_km / 111.19, 20, 160)
    E['a20_pred2_nm'] = 20 * 10 ** (E.ms_est2 - 1.66 * np.log10(dd) - 3.3) * 1e3
    tele = E.dist_km > 1500
    extra = np.select([E.mag < 6, E.mag < 7, E.mag < 8], [1200, 3600, 3 * 3600], 86400)
    full_end = E.t + E.dist_km / 2.5 + extra
    hd = tele & ((E.a20_pred2_nm >= 20) | (E.mag >= 6))
    E['head2_t0'] = np.where(hd, np.where(E.mag < 6, E.t + E.dist_km / 4.5 - 300, E.t), np.nan)
    E['head2_t1'] = np.where(hd, np.where(E.mag < 6, E.t + E.dist_km / 2.8 + 900, full_end), np.nan)
    bw = tele & (E.mag >= 5.0)
    tp = np.full(len(E), np.nan); ts = tp.copy()
    for j, (i, r) in enumerate(E.iterrows()):
        if bw.loc[i]:
            tp[j], ts[j] = ps_times(r.depth, r.dist_km / 111.19, model)
    E['tP_s'] = tp; E['tS_s'] = ts
    E['bw5_t0'] = np.where(bw & np.isfinite(tp), E.t + tp - 60, np.nan)
    E['bw5_t1'] = np.where(bw & np.isfinite(ts), E.t + ts + 600, np.where(bw & np.isfinite(tp), E.t + tp + 900, np.nan))
    E['src_tab'] = 'scan_eq_events'
    R = pd.read_csv(f'{C.DATA}/eq/comcat_reg2500_M30_{S.YEAR}.csv', usecols=['time', 'latitude', 'longitude', 'depth', 'mag', 'magType', 'id'])
    R = R.rename(columns=dict(latitude='lat', longitude='lon')).assign(src='comcat_reg2500_M30', src_tab=f'comcat_reg2500_M30_{S.YEAR}')
    R['t'] = tsec(R.time); R['dist_km'] = gc_km(R.lat.values, R.lon.values)
    R['reg3_t0'] = R.t; R['reg3_t1'] = R.t + R.dist_km / 2.5 + 600
    A = pd.concat([E, R], ignore_index=True).sort_values('t').reset_index(drop=True)
    A['time_pt'] = pd.to_datetime(A.t, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    A.to_csv(S.tab('scan2_eq_events.csv'), index=False)
    M = {}
    for tier in ('head', 'm45', 't58', 'reg', 'head2', 'bw5', 'reg3'):
        M[tier] = interval_mask(A[f'{tier}_t0'].values.astype(float), A[f'{tier}_t1'].values.astype(float))
    M['v1'] = M['head'] | M['reg']
    reg3_ok = M['reg3'].mean() <= 0.15
    M['v2'] = M['head2'] | M['bw5'] | M['reg'] | (M['reg3'] if reg3_ok else False)
    np.savez_compressed(OUT, **M, reg3_in_v2=reg3_ok)
    F = pd.DataFrame([dict(tier=k, **{f'frac_{S.YEAR}': float(v.mean())}) for k, v in M.items()])
    F.loc[len(F)] = dict(tier='reg3_in_v2', **{f'frac_{S.YEAR}': float(reg3_ok)})
    F.to_csv(S.tab('scan2_mask_fraction.csv'), index=False)
    print(F.round(4).to_string())
    return A, M


_CACHE = {}


def mask(tier):
    if tier not in _CACHE:
        if not os.path.isfile(OUT):
            build()
        Z = np.load(OUT)
        for k in Z.files:
            _CACHE[k] = Z[k]
    return _CACHE[tier]


if __name__ == '__main__':
    build()
