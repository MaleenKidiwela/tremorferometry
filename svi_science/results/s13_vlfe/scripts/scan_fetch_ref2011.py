#!/usr/bin/env python
"""VLFE scan pilot 2 (S7.1, user approval 2026-10-15 ~0:30 am PT): 1-sps LHZ for calendar 2011 of the far references from
EarthScope FDSN, monthly requests, <= 2 in parallel, polite retries; split into day files
data/raw_ref2011/NET.STA/2011/DOY.mseed; StationXML (response level) -> data/resp_scan/NET.STA.xml.
CN.PNT and CN.BBB have no LH channel at NRCan or EarthScope (BH 40 Hz only; a year of BHZ is ~2 GB per station, far over
the approved ~150 MB) -> substitutes with LHZ at EarthScope: UW.OMAK for CN.PNT (same azimuth ~E-NE, 315 km) and TA.I02D for
CN.BBB (no LH station exists north of the network in 2011; I02D = coastal Oregon, ~470 km S, adds the coastal azimuth).
LHN/LHE not fetched (3x the volume; not 'trivially small'). Log: logs/ref_download.log."""
import os, sys, time, io, datetime as dt
from concurrent.futures import ThreadPoolExecutor
import urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

DS = 'https://service.earthscope.org/fdsnws/dataselect/1/query'
ST = 'https://service.earthscope.org/fdsnws/station/1/query'
STAS = [('IU', 'COR', '00'), ('US', 'NEW', '--,00'), ('US', 'HAWA', '--,00'), ('UW', 'OMAK', '--'), ('TA', 'I02D', '--')]
OUT = f'{C.DATA}/raw_ref2011'
RESP = f'{C.DATA}/resp_scan'
LOG = f'{C.LOGD}/ref_download.log'


def log(m):
    s = f'{dt.datetime.now(dt.timezone.utc).astimezone(dt.timezone(dt.timedelta(hours=-7))):%Y-%m-%d %H:%M:%S} PT {m}'
    print(s, flush=True)
    open(LOG, 'a').write(s + '\n')


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code == 204:
                return 204, b''
            log(f'HTTP {e.code} {url} (try {k + 1})')
        except Exception as e:  # noqa
            log(f'ERR {type(e).__name__} {url} (try {k + 1})')
        time.sleep(10 * (k + 1))
    return None, b''


def month(job):
    from obspy import read
    net, sta, loc, m = job
    t0 = dt.datetime(2011, m, 1); t1 = dt.datetime(2011 + (m == 12), m % 12 + 1, 1)
    url = f'{DS}?net={net}&sta={sta}&loc={loc}&cha=LHZ&starttime={t0:%Y-%m-%dT%H:%M:%S}&endtime={t1:%Y-%m-%dT%H:%M:%S}'
    code, b = get(url)
    if not b:
        log(f'{net}.{sta} 2011-{m:02d}: no data (code {code})'); return 0
    st = read(io.BytesIO(b))
    d = f'{OUT}/{net}.{sta}/2011'; os.makedirs(d, exist_ok=True)
    nday = 0
    day = t0
    from obspy import UTCDateTime
    while day < t1:
        a = UTCDateTime(day); s = st.slice(a, a + 86400 - 1e-3)
        if len(s):
            s.write(f'{d}/{day.timetuple().tm_yday:03d}.mseed', format='MSEED'); nday += 1
        day += dt.timedelta(days=1)
    gaps = len(st.get_gaps())
    log(f'{net}.{sta} 2011-{m:02d}: {len(b) / 1e6:.2f} MB, {len(st)} traces, {gaps} gaps, {nday} day files')
    return len(b)


def main():
    os.makedirs(RESP, exist_ok=True)
    log('start: approved by the user 2026-10-15 ~0:30 am PT; EarthScope FDSN, LHZ 2011; CN.PNT/CN.BBB -> UW.OMAK/TA.I02D (no LH)')
    for net, sta, loc in STAS:
        f = f'{RESP}/{net}.{sta}.xml'
        if not os.path.isfile(f):
            code, b = get(f'{ST}?net={net}&sta={sta}&cha=LH?&starttime=2010-12-01&endtime=2012-02-01&level=response')
            open(f, 'wb').write(b); log(f'StationXML {net}.{sta}: {len(b) / 1e3:.0f} kB (code {code})')
    jobs = [(n, s, l, m) for n, s, l in STAS for m in range(1, 13)]
    with ThreadPoolExecutor(2) as ex:
        tot = sum(ex.map(month, jobs))
    log(f'done: {tot / 1e6:.1f} MB total')
    open(f'{C.OUT}/REF2011_DOWNLOAD_DONE', 'w').write('ok\n')


if __name__ == '__main__':
    main()
