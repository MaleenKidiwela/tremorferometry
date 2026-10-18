#!/usr/bin/env python
"""VLFE scan (PLAN_SCAN.md §1): truth tables and earthquake masks. No detector numbers are used or produced here.
Writes
  tables/hg2018_table1.csv           Hutchison & Ghosh (2019, JGR 124, 554-565; accepted 2018) Table 1, transcribed
  tables/scan_footprint_families.csv SVI families with >= 30 catalogued P >= 0.75 counted-once LFEs 2003-2026
  tables/scan_footprint_polygon.csv  10-km-buffered convex hull of those families (lat, lon vertices)
  tables/station_epochs.csv          BH/HH channel epochs from data/resp/*.xml
  tables/scan_eq_events.csv          every masking event with its mask windows (tiers: head, merlin45, tele58, regional)
  tables/scan_mask_fraction.csv      masked fraction per year and tier
"""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

LAT0, LON0 = 48.37, -123.61          # SVI count-weighted family centroid (Merlin step 6e)


def hg_table():
    src = ('Hutchison, A. A. & Ghosh, A. (2019). Repeating VLFEs during ETS events in Cascadia track slow slip and continue '
           'throughout inter-ETS period. JGR Solid Earth 124, 554-565, doi:10.1029/2018JB016138 (accepted 27 Nov 2018), '
           'Table 1, p. 557; transcribed 2026-10-14 PT from notes/papers/Hutchison_Ghosh_2018_JGR_VLFE.pdf. Lon given as deg W in '
           'the paper -> negative here. Date in the paper is dd/mm/yyyy.')
    rows = [('17/08/2011', '04:30:28', 123.075, 47.825, 53, 3.4, 1), ('21/08/2011', '04:09:06', 123.6, 48.3, 54, 3.7, 2),
            ('21/08/2011', '08:43:29', 123.45, 48.4, 52, 3.3, 3), ('24/08/2011', '10:47:38', 123.625, 48.35, 30, 3.5, 4),
            ('26/08/2011', '13:44:34', 123.175, 48.35, 55, 3.4, 5), ('02/12/2014', '04:58:36', 123.725, 48.825, 31, 3.6, 6),
            ('04/12/2014', '21:02:27', 123.725, 48.975, 47, 3.7, 7), ('05/12/2014', '11:56:58', 123.775, 48.825, 34, 3.6, 8),
            ('11/12/2014', '10:29:35', 123.75, 48.7, 36, 3.9, 9), ('11/12/2014', '12:21:20', 123.125, 48.725, 33, 3.8, 10),
            ('11/12/2014', '19:09:29', 123.7, 48.775, 36, 4.1, 11), ('12/12/2014', '00:13:17', 124.05, 49.1, 60, 4.1, 12),
            ('14/12/2014', '05:23:51', 123.875, 48.65, 30, 3.7, 13)]
    D = pd.DataFrame(rows, columns=['date_ddmmyyyy', 'time_utc', 'lon_W', 'lat', 'depth_km', 'Mw', 'template'])
    D['time'] = pd.to_datetime(D.date_ddmmyyyy + ' ' + D.time_utc, format='%d/%m/%Y %H:%M:%S')
    D['lon'] = -D.lon_W
    D['time_PT'] = D.time.dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    D['M0_Nm'] = 10 ** (1.5 * D.Mw + 9.1)
    with open(f'{C.TB}/hg2018_table1.csv', 'w') as fh:
        fh.write(f'# source: {src}\n')
        fh.write('# text statements (no numeric between-ETS rate is given): 2011 ETS window used by H&G 2011-07-23..09-06 '
                 '(PNSN); 1,394 new matched-filter detections in the 2011 ETS from 5 templates; inter-ETS 2011-10-01..2012-08-30 '
                 '"a fair amount of VLFE activity" for all five 2011 templates, sustained above their noise-template background '
                 '(Fig. 6-7, counts only in figures); 3 of 5 templates rise again in the 2012 ETS; threshold = mean 3rd-highest '
                 'summed CC over 100 noise templates (~6 MAD, ~0.6 false detections per template per day quoted); VLFE '
                 'magnitudes in Cascadia Mw 3.1-4.3; signal duration ~90 s; matched-filter detections are SMALLER than the '
                 'templates (below the CMT grid-search detection level, p. 9).\n')
        D[['template', 'time', 'time_PT', 'lat', 'lon', 'depth_km', 'Mw', 'M0_Nm']].to_csv(fh, index=False)
    return D


def footprint():
    from shapely.geometry import MultiPoint
    R = C.svi_rows()
    R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.75)]
    n = R.groupby('fam').size().rename('n_P075')
    F = C.svi_families().join(n, how='left').fillna({'n_P075': 0})
    F['in_footprint_set'] = F.n_P075 >= 30
    F.to_csv(f'{C.TB}/scan_footprint_families.csv')
    S = F[F.in_footprint_set]
    kx = 111.19 * np.cos(np.radians(LAT0)); ky = 111.19
    pts = MultiPoint([((lo - LON0) * kx, (la - LAT0) * ky) for la, lo in zip(S.lat, S.lon)])
    poly = pts.convex_hull.buffer(10.0, quad_segs=16)
    xy = np.array(poly.exterior.coords)
    P = pd.DataFrame(dict(lat=LAT0 + xy[:, 1] / ky, lon=LON0 + xy[:, 0] / kx))
    P.to_csv(f'{C.TB}/scan_footprint_polygon.csv', index=False)
    return F, P, poly.area


def epochs():
    from obspy import read_inventory
    rows = []
    for f in sorted(glob.glob(f'{C.RESP}/*.xml')):
        inv = read_inventory(f)
        for net in inv:
            for sta in net:
                for ch in sta:
                    if ch.code[:2] not in ('BH', 'HH'):
                        continue
                    sens = ch.sensor.description if ch.sensor is not None else ''
                    rows.append(dict(netsta=f'{net.code}.{sta.code}', loc=ch.location_code, cha=ch.code,
                                     start=str(ch.start_date)[:19], end=str(ch.end_date)[:19] if ch.end_date else '',
                                     sr=ch.sample_rate, az=ch.azimuth, dip=ch.dip, sensor=sens, lat=ch.latitude,
                                     lon=ch.longitude, file=os.path.basename(f)))
    E = pd.DataFrame(rows).drop_duplicates(['netsta', 'loc', 'cha', 'start'])
    E.to_csv(f'{C.TB}/station_epochs.csv', index=False)
    return E


def dist_km(la, lo):
    la0, lo0 = np.radians(LAT0), np.radians(LON0)
    la, lo = np.radians(la), np.radians(lo)
    return np.degrees(np.arccos(np.clip(np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0), -1, 1))) * 111.19


def tsec(s):
    return pd.to_datetime(s, utc=True, format='mixed').dt.tz_localize(None).astype('datetime64[ns]').astype('int64') / 1e9


def eq_events():
    G = pd.concat([pd.read_csv(f, usecols=['time', 'latitude', 'longitude', 'depth', 'mag', 'magType', 'id'])
                   for f in sorted(glob.glob(f'{C.DATA}/eq/comcat_global_M45_*.csv'))], ignore_index=True)
    G = G.rename(columns=dict(latitude='lat', longitude='lon')).assign(src='comcat_global_M45')
    T = pd.read_csv(f'{C.ROOT}/data/rf_teleseis/events_M58_1993_2026.csv').rename(
        columns=dict(latitude='lat', longitude='lon', evid='id')).assign(src='rf_teleseis_M58')
    T = T[T.time >= '2002-12-01']
    R3 = pd.read_csv(f'{C.DATA}/eq/comcat_reg300_M25_2003_2026.csv', usecols=['time', 'latitude', 'longitude', 'depth', 'mag', 'magType', 'id'])
    R3 = R3.rename(columns=dict(latitude='lat', longitude='lon')).assign(src='comcat_reg300_M25')
    R15 = pd.read_csv(f'{C.DATA}/eq/comcat_reg1500_M40_2003_2026.csv', usecols=['time', 'latitude', 'longitude', 'depth', 'mag', 'magType', 'id'])
    R15 = R15.rename(columns=dict(latitude='lat', longitude='lon')).assign(src='comcat_reg1500_M40')
    K = pd.read_csv(f'{C.ROOT}/lfe_extension/results/trackC/eq_catalog.csv').rename(columns=dict(ot='time'))
    K = K.assign(magType='', id='', src='trackC_' + K.src.astype(str))
    E = pd.concat([G, T, R3, R15, K[['time', 'lat', 'lon', 'depth', 'mag', 'magType', 'id', 'src']]], ignore_index=True)
    E['t'] = tsec(E.time)
    E['dist_km'] = dist_km(E.lat.values, E.lon.values)
    mb = E.magType.astype(str).str.lower().str.startswith('mb')
    E['ms_est'] = np.where(mb, (E.mag - 2.5) / 0.63, E.mag)            # Gutenberg-Richter 1956 mb = 0.63 Ms + 2.5
    dd = np.clip(E.dist_km / 111.19, 20, 160)
    E['a20_pred_nm'] = 20 * 10 ** (E.ms_est - 1.66 * np.log10(dd) - 3.3) * 1e3   # Ms (Prague) at T = 20 s
    tele = E.dist_km > 1500
    extra = np.select([E.mag < 6, E.mag < 7, E.mag < 8], [1200, 3600, 3 * 3600], 86400)
    full_end = E.t + E.dist_km / 2.5 + extra
    # tier 'merlin45': Merlin literal, every teleseism M >= 4.5: origin -> origin + D/2.5 + extra(M)
    E['m45_t0'] = np.where(tele, E.t, np.nan); E['m45_t1'] = np.where(tele, full_end, np.nan)
    # tier 'tele58': Merlin literal on the M >= 5.8 list only
    t58 = tele & (E.mag >= 5.8)
    E['t58_t0'] = np.where(t58, E.t, np.nan); E['t58_t1'] = np.where(t58, full_end, np.nan)
    # tier 'head' (PLAN_SCAN amendment S1): teleseisms with predicted 20-s amplitude >= 20 nm or M >= 6.
    #   M < 6: surface-wave group window origin + D/4.5 - 300 s .. origin + D/2.8 + 900 s; M >= 6: Merlin's full window.
    hd = tele & ((E.a20_pred_nm >= 20) | (E.mag >= 6))
    E['head_t0'] = np.where(hd, np.where(E.mag < 6, E.t + E.dist_km / 4.5 - 300, E.t), np.nan)
    E['head_t1'] = np.where(hd, np.where(E.mag < 6, E.t + E.dist_km / 2.8 + 900, full_end), np.nan)
    # regional (all tiers): M >= 2.5 within 300 km +-5 min; M >= 4 within 1500 km +-20 min
    r1 = (E.dist_km <= 300) & (E.mag >= 2.5); r2 = (E.dist_km <= 1500) & (E.mag >= 4.0)
    E['reg_t0'] = np.where(r2, E.t - 1200, np.where(r1, E.t - 300, np.nan))
    E['reg_t1'] = np.where(r2, E.t + 1200, np.where(r1, E.t + 300, np.nan))
    keep = E[['m45_t0', 'head_t0', 'reg_t0', 't58_t0']].notna().any(axis=1)
    E = E[keep].sort_values('t').reset_index(drop=True)
    E.to_csv(f'{C.TB}/scan_eq_events.csv', index=False)
    rows = []
    for y in range(2003, 2027):
        y0 = pd.Timestamp(f'{y}-01-01').value / 1e9; y1 = pd.Timestamp(f'{y + 1}-01-01').value / 1e9
        if y == 2026:
            y1 = pd.Timestamp('2026-09-15').value / 1e9
        row = dict(year=y)
        for tier in ('head', 'm45', 't58', 'reg'):
            m = interval_mask(E[f'{tier}_t0'].values, E[f'{tier}_t1'].values, y0, y1)
            row[f'frac_{tier}'] = m.mean()
        mh = interval_mask(E.head_t0.values, E.head_t1.values, y0, y1) | interval_mask(E.reg_t0.values, E.reg_t1.values, y0, y1)
        row['frac_headline_total'] = mh.mean()
        rows.append(row)
    M = pd.DataFrame(rows); M.to_csv(f'{C.TB}/scan_mask_fraction.csv', index=False)
    return E, M


def interval_mask(a, b, y0, y1):
    m = np.zeros(int(y1 - y0), bool)
    ok = np.isfinite(a) & np.isfinite(b) & (b > y0) & (a < y1)
    for s, e in zip(a[ok], b[ok]):
        m[int(max(s - y0, 0)):int(min(e - y0, len(m)))] = True
    return m


if __name__ == '__main__':
    D = hg_table(); print(D)
    F, P, area = footprint(); print('footprint families', int(F.in_footprint_set.sum()), 'of', len(F), 'area km2', round(area))
    E = epochs(); print(E[E.netsta.isin(['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV', 'CN.LZB'])]
                        [['netsta', 'cha', 'start', 'end', 'sensor']].to_string())
    EQ, M = eq_events(); print(M.round(3).to_string())
