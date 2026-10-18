#!/usr/bin/env python
"""VLFE scan pilot 2 (PLAN_SCAN.md Amendment S7.1): reference-station table, coordinates and response lookup.

Reference stations (role 'veto' if >= 100 km from (48.37 N, -123.61), else 'column'); the pending far stations (download
asked, not approved) are listed with role 'veto' and status 'pending' so that they drop in once their data exist.
Coordinates: lfe_extension_wa/results/wa_temp_networks/stations.csv where listed, else the StationXML in
lfe_array_discovery/results/stationxml/{NET}.xml, else data/resp_scan/{NET}.xml.
Writes tables/scan2_ref_stations.csv when run as a script."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

ROOT = C.ROOT
XML_ALT = f'{ROOT}/lfe_array_discovery/results/stationxml'
XML_SCAN = f'{C.DATA}/resp_scan'
TEMP_CSV = f'{ROOT}/lfe_extension_wa/results/wa_temp_networks/stations.csv'
REFS = [  # (netsta, set)
    ('UW.GNW', 'full2011'), ('UW.DOSE', 'full2011'),
    ('CC.OBSR', 'ets_days'), ('UW.LEBA', 'ets_days'), ('UW.RATT', 'ets_days'), ('UW.ELW', 'ets_days'),
    ('TA.D04D', 'ets_days'), ('UW.WISH', 'ets_days'), ('UW.SP2', 'ets_days'), ('TA.B05D', 'ets_days'),
    ('US.NLWA', 'ets_days'), ('TA.D03D', 'ets_days'),
    ('IU.COR', 'far_lhz'), ('US.NEW', 'far_lhz'), ('US.HAWA', 'far_lhz'), ('UW.OMAK', 'far_lhz'),
    ('TA.I02D', 'far_lhz'),
    ('US.MSO', 'far_lhz'), ('US.WVOR', 'far_lhz')]           # catalogue Amendment Q0 (2026-10-15): 2003-2026 download
# far set fetched after user approval (2026-10-15 ~0:30 am PT); CN.PNT and CN.BBB have no LH channel (NRCan/EarthScope) ->
# substitutes UW.OMAK (for PNT, same azimuth/distance) and TA.I02D (for BBB; no LH station north of the network in 2011)
ETS_DAYS = ('2011-07-30', '2011-08-27')          # DOY 211-239
LAT0, LON0 = 48.37, -123.61
_INV = {}


def _xml(path):
    from obspy import read_inventory
    if path not in _INV:
        _INV[path] = read_inventory(path) if os.path.isfile(path) else None
    return _INV[path]


def xml_candidates(net, sta=None):
    first = [f'{XML_SCAN}/{net}.{sta}.xml'] if sta and os.path.isfile(f'{XML_SCAN}/{net}.{sta}.xml') else []
    if sta and os.path.isfile(f'{XML_SCAN}/{net}.{sta}.2012.xml'):          # 2012 far-reference StationXML (Amendment S8.2)
        first.append(f'{XML_SCAN}/{net}.{sta}.2012.xml')
    if sta and os.path.isfile(f'{XML_SCAN}/{net}.{sta}.all.xml'):           # 2003-2026 StationXML (catalogue Amendment P)
        first.append(f'{XML_SCAN}/{net}.{sta}.all.xml')
    return first + [f'{C.RESP}/{net}_bb.xml', f'{XML_ALT}/{net}.xml', f'{XML_SCAN}/{net}.xml']


def inventory_for(netsta, t=None):
    """First StationXML (in the S7.1 order) holding the station (with an epoch at UTCDateTime t if given), restricted to
    that station."""
    net, sta = netsta.split('.')
    for f in xml_candidates(net, sta):
        inv = _xml(f)
        if inv is None:
            continue
        s = inv.select(station=sta, time=t) if t is not None else inv.select(station=sta)
        if len(s) and any(len(st.channels) for n in s for st in n):
            return s, f
    return None, None


def coords(netsta):
    net, sta = netsta.split('.')
    T = pd.read_csv(TEMP_CSV)
    r = T[(T.net == net) & (T.sta.astype(str) == sta)]
    if len(r):
        return float(r.lat.iloc[0]), float(r.lon.iloc[0]), 'stations.csv'
    cands = xml_candidates(net, sta)
    for f in cands[:1] + [c for c in cands[1:] if 'stationxml' in c] + [c for c in cands[1:] if 'stationxml' not in c]:
        inv = _xml(f)
        if inv is None:
            continue
        s = inv.select(station=sta)
        if len(s):
            st = s[0][0]
            return float(st.latitude), float(st.longitude), os.path.relpath(f, ROOT)
    return np.nan, np.nan, 'missing'


def dist0(lat, lon):
    kx = 111.19 * np.cos(np.radians(LAT0))
    return float(np.hypot((lon - LON0) * kx, (lat - LAT0) * 111.19))


def table():
    rows = []
    for k, (ns, st) in enumerate(REFS):
        la, lo, src = coords(ns)
        d = dist0(la, lo) if np.isfinite(la) else np.nan
        bp = os.path.isfile(f'{C.DATA}/scan/bp2_ref/{ns}.npy')
        rows.append(dict(netsta=ns, set=st, k_ctrl=7 + k, lat=la, lon=lo, coord_src=src, dist_km=d,
                         role=('column' if (np.isfinite(d) and d < 100) else 'veto'),
                         status=('fetched after user approval 2026-10-15' if st == 'far_lhz' else 'on disk'),
                         bp_exists=bp))
    return pd.DataFrame(rows)


def load_table():
    f = f'{C.TB}/scan2_ref_stations.csv'
    T = pd.read_csv(f) if os.path.isfile(f) else table()
    sd = f'{C.DATA}/scan' if os.environ.get('SCAN_YEAR', '2011') == '2011' else f"{C.DATA}/scan_{os.environ['SCAN_YEAR']}"
    T['bp_exists'] = [os.path.isfile(f'{sd}/bp2_ref/{s}.npy') for s in T.netsta]          # SCAN_YEAR (Amendment S8)
    return T


def all_coords():
    """Pilot + reference coordinates as a DataFrame indexed by NET.STA (lat, lon)."""
    cc = C.station_coords()[['lat', 'lon']]
    T = load_table().set_index('netsta')[['lat', 'lon']]
    T = T[T.lat.notna()]
    return pd.concat([cc[~cc.index.isin(T.index)], T])


if __name__ == '__main__':
    T = table()
    T.to_csv(f'{C.TB}/scan2_ref_stations.csv', index=False)
    print(T.to_string())
