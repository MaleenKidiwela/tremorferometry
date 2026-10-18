#!/usr/bin/env python
"""Catalogue step 5 (metadata only; PLAN_CAT.md §8): candidate far LHZ references 2003-2026 at 250-800 km from the
SVI scan centre (48.37 N, -123.61 E; scan_common LAT0/LON0). FDSN station metadata (level=channel, text) from
EarthScope/IRIS (IU, US, UW, TA, BK, CC, CN, II), NCEDC (BK), and NRCan (CN); EarthScope availability extents (time
spans, metadata only) for each candidate LHZ channel. NO waveforms are requested.
Plain urllib GET, no identifying headers beyond the default Python user agent; nothing personal is sent.
Outputs: tables/farref_channels.csv (one row per channel epoch), tables/farref_availability.csv (spans),
tables/farref_per_year.csv (per station-year: metadata-open days, available days, azimuth, distance),
tables/farref_count_per_year.csv, tables/farref_download_estimate.csv."""
import os, io, sys, time, urllib.request, urllib.parse
import numpy as np, pandas as pd

CAT = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/catalogue'
TB = f'{CAT}/tables'
LAT0, LON0 = 48.37, -123.61
RMIN, RMAX = 250.0, 800.0
NETS = 'IU,US,UW,TA,BK,CC,CN,II'
T0, T1 = '2003-01-01', '2026-10-15'
SERV = {'IRIS': 'https://service.iris.edu/fdsnws/station/1/query',
        'NCEDC': 'https://service.ncedc.org/fdsnws/station/1/query',
        'NRCAN': 'https://earthquakescanada.nrcan.gc.ca/fdsnws/station/1/query'}
AVAIL = 'https://service.iris.edu/mustang/measurements/1/query'   # fdsnws/availability is retired (2026-10-15 test)
MB_PER_SY = {'IU.COR': 80.5, 'US.NEW': 47.6, 'US.HAWA': 50.3, 'UW.OMAK': 46.8, 'TA.I02D': 57.1}   # on-disk 2011 LHZ


def log(m):
    line = f'{pd.Timestamp.now(tz="America/Los_Angeles"):%Y-%m-%d %H:%M:%S} PT farref {m}'
    print(line, flush=True)
    open(f'{CAT}/logs/farref.log', 'a').write(line + '\n')


def get(url, tries=3):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode('utf-8', 'replace'), r.status
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                return '', e.code
            err = f'HTTP {e.code}'
        except Exception as e:  # noqa
            err = f'{type(e).__name__}: {e}'
        time.sleep(3 * (k + 1))
    return None, err


def dist_az(la, lo):
    la0, lo0, la1, lo1 = map(np.radians, (LAT0, LON0, np.asarray(la, float), np.asarray(lo, float)))
    d = np.arccos(np.clip(np.sin(la0) * np.sin(la1) + np.cos(la0) * np.cos(la1) * np.cos(lo1 - lo0), -1, 1))
    az = np.degrees(np.arctan2(np.sin(lo1 - lo0) * np.cos(la1),
                               np.cos(la0) * np.sin(la1) - np.sin(la0) * np.cos(la1) * np.cos(lo1 - lo0))) % 360
    return np.degrees(d) * 111.19, az


def stations():
    rows = []
    for name, base in SERV.items():
        nets = NETS if name == 'IRIS' else ('BK' if name == 'NCEDC' else 'CN')
        q = dict(net=nets, cha='LHZ', starttime=T0 + 'T00:00:00', endtime=T1 + 'T00:00:00', latitude=LAT0, longitude=LON0,
                 minradius=round(RMIN / 111.19 - 0.05, 3), maxradius=round(RMAX / 111.19 + 0.05, 3),
                 level='channel', format='text')
        txt, st = get(base + '?' + urllib.parse.urlencode(q))
        if not txt:
            log(f'{name}: no text ({st})')
            continue
        D = pd.read_csv(io.StringIO(txt), sep='|')
        D.columns = [c.strip().lstrip('#').strip() for c in D.columns]
        D['source'] = name
        rows.append(D)
        log(f'{name}: {len(D)} channel epochs')
    D = pd.concat(rows, ignore_index=True)
    D['netsta'] = D.Network + '.' + D.Station
    D['Location'] = D.Location.fillna('').astype(str).str.replace('nan', '')
    D['dist_km'], D['az_deg'] = dist_az(D.Latitude, D.Longitude)
    D = D[(D.dist_km >= RMIN) & (D.dist_km <= RMAX)].copy()
    D['StartTime'] = pd.to_datetime(D.StartTime, utc=True, errors='coerce').dt.tz_localize(None)
    D['EndTime'] = pd.to_datetime(D.EndTime.replace('', np.nan), utc=True, errors='coerce').dt.tz_localize(None).fillna(pd.Timestamp('2599-01-01'))
    D = D.drop_duplicates(['netsta', 'Location', 'Channel', 'StartTime'])
    D.to_csv(f'{TB}/farref_channels.csv', index=False)
    return D


def _avail_one(ns):
    net, sta = ns.split('.')
    q = dict(metric='percent_availability', net=net, sta=sta, cha='LHZ',
             timewindow=f'{T0}T00:00:00,{T1}T00:00:00', format='text')
    txt, st = get(AVAIL + '?' + urllib.parse.urlencode(q))
    if txt is None:
        return [dict(netsta=ns, avail_status=f'error {st}')]
    if not txt.strip() or txt.count('\n') < 2:
        return [dict(netsta=ns, avail_status='no MUSTANG rows')]
    A = pd.read_csv(io.StringIO(txt), skiprows=1)
    A['day'] = pd.to_datetime(A.start.str[:10], format='%Y/%m/%d')
    A = A.groupby('day').value.max().reset_index()           # any location code
    return [dict(netsta=ns, avail_status='ok', day=r.day, pct=float(r.value)) for r in A.itertuples()]


def availability(D):
    """MUSTANG percent_availability (daily, QC metadata; no waveforms) for every candidate LHZ station."""
    from concurrent.futures import ThreadPoolExecutor
    out = []
    with ThreadPoolExecutor(4) as ex:
        for k, r in enumerate(ex.map(_avail_one, sorted(D.netsta.unique()))):
            out += r
            if k % 20 == 0:
                log(f'availability: {k + 1} stations done')
    A = pd.DataFrame(out)
    A.to_parquet(f'{TB}/farref_availability_daily.parquet')
    log(f'availability: {len(A)} rows over {A.netsta.nunique()} stations; status {A.drop_duplicates("netsta").avail_status.value_counts().to_dict()}')
    return A


def per_year(D, A):
    years = range(2003, 2027)
    rows = []
    A = A[A.avail_status == 'ok'].copy()
    for ns, g in D.groupby('netsta'):
        a = A[A.netsta == ns]
        for y in years:
            y0, y1 = pd.Timestamp(f'{y}-01-01'), min(pd.Timestamp(f'{y + 1}-01-01'), pd.Timestamp(T1))
            if y0 >= y1:
                continue
            iv = sorted((max(r.StartTime, y0), min(r.EndTime, y1)) for r in g.itertuples() if r.EndTime > y0 and r.StartTime < y1)
            md, cur = 0.0, None                                   # union of channel epochs (all location codes)
            for a0, a1 in iv:
                if cur is None or a0 > cur[1]:
                    if cur:
                        md += (cur[1] - cur[0]).total_seconds() / 86400
                    cur = [a0, a1]
                else:
                    cur[1] = max(cur[1], a1)
            if cur:
                md += (cur[1] - cur[0]).total_seconds() / 86400
            if md <= 0:
                continue
            ad = float(((a.day >= y0) & (a.day < y1) & (a.pct >= 50)).sum()) if len(a) and 'day' in a else 0.0
            r0 = g.iloc[0]
            rows.append(dict(netsta=ns, year=y, dist_km=round(r0.dist_km, 1), az_deg=round(r0.az_deg, 1),
                             lat=r0.Latitude, lon=r0.Longitude, locs=','.join(sorted(set(g.Location.fillna('').astype(str)))),
                             metadata_days=round(md, 1), avail_days_mustang=round(ad, 1),
                             sensor=str(g.iloc[-1].get('SensorDescription', ''))[:40]))
    P = pd.DataFrame(rows)
    P.to_csv(f'{TB}/farref_per_year.csv', index=False)
    C = []
    for y in years:
        p = P[P.year == y]
        ok = p[p.avail_days_mustang >= 180]
        okm = p[p.metadata_days >= 180]
        C.append(dict(year=y, n_metadata_ge180d=len(okm), n_avail_ge180d=len(ok),
                      stations_avail=' '.join(f'{r.netsta}({r.dist_km:.0f}km,{r.az_deg:.0f}°)' for r in ok.sort_values('az_deg').itertuples()),
                      az_gap_max_deg=float(np.max(np.diff(np.r_[np.sort(ok.az_deg.values), np.sort(ok.az_deg.values)[:1] + 360]))) if len(ok) >= 2 else 360.0))
    C = pd.DataFrame(C)
    C.to_csv(f'{TB}/farref_count_per_year.csv', index=False)
    mb = float(np.median(list(MB_PER_SY.values())))
    E = dict(mb_per_station_year_median=mb, mb_per_station_year_range=f'{min(MB_PER_SY.values())}-{max(MB_PER_SY.values())}',
             basis='on-disk 2011 LHZ STEIM2 day files of the 5 R1 far references (data/raw_ref2011)')
    pd.DataFrame([E]).to_csv(f'{TB}/farref_download_estimate.csv', index=False)
    log('per-year tables written')
    return P, C


def cn_broadband():
    """secondary class (Merlin step 5 fallback): CN broadband BHZ/HHZ at 250-800 km (NRCan FDSN; CN has no LHZ)."""
    q = dict(net='CN', cha='BHZ,HHZ', starttime=T0 + 'T00:00:00', latitude=LAT0, longitude=LON0,
             minradius=round(RMIN / 111.19 - 0.05, 3), maxradius=round(RMAX / 111.19 + 0.05, 3), level='channel', format='text')
    txt, st = get(SERV['NRCAN'] + '?' + urllib.parse.urlencode(q))
    if not txt:
        log(f'NRCAN broadband: none ({st})'); return None
    D = pd.read_csv(io.StringIO(txt), sep='|')
    D.columns = [c.strip().lstrip('#').strip() for c in D.columns]
    D = D[D.SensorDescription.str.contains('CMG3|Trillium|STS|CMG-3|T120|T240', case=False, na=False)].copy()
    D['netsta'] = D.Network + '.' + D.Station
    D['dist_km'], D['az_deg'] = dist_az(D.Latitude, D.Longitude)
    D = D[(D.dist_km >= RMIN) & (D.dist_km <= RMAX)]
    D['StartTime'] = pd.to_datetime(D.StartTime, errors='coerce')
    D['EndTime'] = pd.to_datetime(D.EndTime, errors='coerce').fillna(pd.Timestamp('2599-01-01'))
    D.to_csv(f'{TB}/farref_cn_broadband_channels.csv', index=False)
    rows = []
    for y in range(2003, 2027):
        y0, y1 = pd.Timestamp(f'{y}-01-01'), pd.Timestamp(f'{y + 1}-01-01')
        on = D[(D.StartTime < y1 - pd.Timedelta(days=180)) & (D.EndTime > y0 + pd.Timedelta(days=180))]
        u = on.drop_duplicates('netsta')
        rows.append(dict(year=y, n_cn_broadband=len(u), stations=' '.join(f'{r.netsta}({r.dist_km:.0f}km,{r.az_deg:.0f}°,{r.Channel})'
                                                                          for r in u.sort_values('az_deg').itertuples())))
    pd.DataFrame(rows).to_csv(f'{TB}/farref_cn_broadband_per_year.csv', index=False)
    log(f'NRCAN broadband: {D.netsta.nunique()} stations')


def proposal(P):
    """download estimate for a fixed long-running 5-station LHZ set (Merlin step 5) + per-year totals."""
    S5 = ['IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']
    mb = pd.read_csv(f'{TB}/farref_download_estimate.csv').iloc[0]
    q = P[P.netsta.isin(S5)].copy()
    q['station_years'] = np.maximum(q.metadata_days, q.avail_days_mustang) / 365.25
    on_disk = {('IU.COR', 2011), ('US.NEW', 2011), ('US.HAWA', 2011), ('IU.COR', 2012), ('US.NEW', 2012), ('US.HAWA', 2012)}
    q['on_disk'] = [(r.netsta, r.year) in on_disk for r in q.itertuples()]
    sy = q[~q.on_disk].station_years.sum()
    out = dict(stations=' '.join(S5), station_years_needed=round(sy, 1),
               gb_at_median=round(sy * mb.mb_per_station_year_median / 1000, 1),
               gb_range=f'{sy * 46.8 / 1000:.1f}-{sy * 80.5 / 1000:.1f}')
    pd.DataFrame([out]).to_csv(f'{TB}/farref_download_proposal.csv', index=False)
    q.to_csv(f'{TB}/farref_proposal_station_years.csv', index=False)
    log(f'proposal: {out}')


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'tables':
    D = pd.read_csv(f'{TB}/farref_channels.csv', parse_dates=['StartTime', 'EndTime'])
    A = pd.read_parquet(f'{TB}/farref_availability_daily.parquet')
    P, C = per_year(D, A)
    proposal(P)
elif __name__ == '__main__':
    cn_broadband()
    D = stations()
    A = availability(D)
    P, C = per_year(D, A)
    proposal(P)
    open(f'{CAT}/FARREF_META_DONE', 'w').write(pd.Timestamp.now(tz='America/Los_Angeles').strftime('%Y-%m-%d %H:%M PT\n'))
