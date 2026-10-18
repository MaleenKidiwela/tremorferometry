#!/usr/bin/env python
"""Fetch the missing 1 Jun 2011 - 1 Dec 2013 days for the H&G (2019) Fig.-1 stations from EarthScope FDSN (PLAN.md step 1).
- TA.A04D, TA.B05D, TA.D03D: LH? (1 sps), stored as is.
- UW.SP2, UW.LRIV, UW.DOSE, UW.GNW: BH? (40 sps) requested with 600 s padding, decimated in memory (8 x 5, obspy
  anti-alias) to 1 sps and trimmed to the day; only the 1-sps file is stored (channel codes kept, BH response applies:
  the prep v3 recipe also decimates 40 -> 1 before remove_response).
- CN.NLLB, CN.VGZ: not at EarthScope (CN at NRCan); on-disk coverage is 576-580 / 580 days, so nothing fetched.
A day is fetched only if no raw day file exists in s13_common.TRACE_ROOTS. Output data/raw/NET.STA/YYYY/DOY.mseed.
StationXML (response) -> data/resp/NET.STA.xml. Bytes transferred are logged; stops if the running total > 5 GB.
usage: hg_fetch.py [--dry]"""
import os, sys, io, time, datetime as dt, urllib.request, urllib.error, threading
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
sys.path.insert(0, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts')
import s13_common as C

HG = f'{C.OUT}/hg_repro'
RAW, RESP, LOG = f'{HG}/data/raw', f'{HG}/data/resp', f'{HG}/logs/fetch.log'
DS = 'https://service.earthscope.org/fdsnws/dataselect/1/query'
STQ = 'https://service.earthscope.org/fdsnws/station/1/query'
STAS = {'TA.A04D': 'LH?', 'TA.B05D': 'LH?', 'TA.D03D': 'LH?', 'UW.SP2': 'BH?', 'UW.LRIV': 'BH?', 'UW.DOSE': 'BH?',
        'UW.GNW': 'BH?'}
T0, T1 = '2011-06-01', '2013-12-01'
CAP = 5e9
_lock = threading.Lock()
_bytes = [0]


def log(m):
    s = f'{dt.datetime.now(dt.timezone.utc).astimezone(dt.timezone(dt.timedelta(hours=-7))):%Y-%m-%d %H:%M:%S} PT {m}'
    print(s, flush=True)
    with _lock:
        open(LOG, 'a').write(s + '\n')


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=600) as r:
                b = r.read()
                with _lock:
                    _bytes[0] += len(b)
                return r.status, b
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                return e.code, b''
            log(f'HTTP {e.code} {url} (try {k + 1})')
        except Exception as e:  # noqa
            log(f'ERR {type(e).__name__} {url} (try {k + 1})')
        time.sleep(10 * (k + 1))
    return None, b''


def day_job(job):
    from obspy import read, UTCDateTime
    netsta, cha, day = job
    if _bytes[0] > CAP:
        return netsta, day, 'cap'
    net, sta = netsta.split('.')
    out = f'{RAW}/{netsta}/{day.year}/{day.dayofyear:03d}.mseed'
    if os.path.isfile(out):
        return netsta, day, 'exists'
    a = UTCDateTime(str(day.date())); b = a + 86400
    pad = 600 if cha.startswith('BH') else 0
    url = (f'{DS}?net={net}&sta={sta}&loc=*&cha={cha}&starttime={(a - pad).isoformat()[:19]}'
           f'&endtime={(b + pad).isoformat()[:19]}')
    code, raw = get(url)
    if not raw:
        return netsta, day, f'nodata {code}'
    st = read(io.BytesIO(raw))
    if cha.startswith('BH'):
        st.merge(method=1, fill_value=None)
        st = st.split()
        keep = []
        for tr in st:
            if tr.stats.npts < 40 * 1200:
                continue
            sr = tr.stats.sampling_rate
            if abs(sr - 40) > 0.01:
                continue
            tr.stats.sampling_rate = 40.0
            tr.data = tr.data.astype('float64')
            tr.detrend('demean')
            tr.decimate(8, no_filter=False, strict_length=False)
            tr.decimate(5, no_filter=False, strict_length=False)
            tr.data = tr.data.astype('float32')
            keep.append(tr)
        from obspy import Stream
        st = Stream(keep).trim(a, b - 1e-3, nearest_sample=False)
        enc = 'FLOAT32'
    else:
        st = st.trim(a, b - 1e-3, nearest_sample=False)
        enc = None
    if len(st) == 0:
        return netsta, day, 'empty'
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + '.part'
    st.write(tmp, format='MSEED', **({'encoding': enc} if enc else {}))
    os.replace(tmp, out)
    return netsta, day, f'ok {len(st)}tr {len(raw) / 1e6:.1f}MB'


def main(dry=False):
    os.makedirs(RESP, exist_ok=True)
    jobs = []
    for netsta, cha in STAS.items():
        net, sta = netsta.split('.')
        x = f'{RESP}/{netsta}.xml'
        if not dry and not os.path.isfile(x):
            code, b = get(f'{STQ}?net={net}&sta={sta}&cha=BH?,LH?&starttime=2011-01-01&endtime=2014-01-01&level=response')
            if b:
                open(x, 'wb').write(b)
                log(f'resp {netsta} {len(b)} B')
        for d in pd.date_range(T0, T1, freq='D', inclusive='left'):
            if not C.day_files(netsta, d.strftime('%Y-%m-%d')):
                jobs.append((netsta, cha, d))
    cnt = pd.Series([j[0] for j in jobs]).value_counts()
    log(f'missing days to fetch: {len(jobs)} ' + ' '.join(f'{k}:{v}' for k, v in cnt.items()))
    if dry:
        return
    with ThreadPoolExecutor(3) as ex:
        for netsta, day, r in ex.map(day_job, jobs):
            log(f'{netsta} {day:%Y-%m-%d} {r} | total {_bytes[0] / 1e9:.2f} GB')
    log(f'FETCH END total transferred {_bytes[0] / 1e9:.2f} GB')


if __name__ == '__main__':
    main('--dry' in sys.argv)
