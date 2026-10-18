#!/usr/bin/env python
"""Amendment Y12 fetch: CN.SYMB 2021-01-01 .. 2026-10-16 and UW.BUCKS 2023-09-17 .. 2026-10-16, HH? (100 sps) decimated in
memory to 1 sps (user approval relayed 2026-10-16 ~6 pm PT). Adapted from fetch_2011carriers.py (same rules below). EarthScope FDSN dataselect, <= 2 parallel, polite retries, default
Python user agent, no personal data.
- A day is fetched only if no file exists for it in s13_common.TRACE_ROOTS, hg_repro/data/raw, data/disp1sps_v3 or
  hg_repro/data/disp1sps.
- TA.A04D, TA.D03D: LH? (1 sps) stored as delivered.
- UW.SP2, UW.WISH, US.NLWA: BH? requested with 600 s padding, decimated in memory to 1 sps (obspy anti-alias, 8 x 5
  for 40 sps, 10 x 5 for 50 sps), rounded to integer counts, trimmed to the day, stored STEIM2 (the prep recipe also
  decimates to 1 sps before response removal). Only the 1-sps file is stored.
Output: data/raw/NET.STA/YYYY/DOY.mseed (a TRACE_ROOT). Log catalogue/logs/fetch_2011carriers.log. Marker
catalogue/FETCH2011_DONE. Stops if stored bytes > 2 GB."""
import os, sys, io, time, datetime as dt, urllib.request, urllib.error, threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
sys.path.insert(0, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts')
import s13_common as C

CAT = f'{C.OUT}/catalogue'
RAW = f'{C.DATA}/raw'
LOG = f'{CAT}/logs/fetch_y12.log'
DS = 'https://service.earthscope.org/fdsnws/dataselect/1/query'
STAS = {'CN.SYMB': ('HH?', '2021-01-01'), 'UW.BUCKS': ('HH?', '2023-09-17')}
END = '2026-10-16'
CAP = 3e9
_lock = threading.Lock()
_b = {'xfer': 0, 'stored': 0}
STOP = threading.Event()


def log(m):
    s = f'{dt.datetime.now(dt.timezone.utc).astimezone(dt.timezone(dt.timedelta(hours=-7))):%Y-%m-%d %H:%M:%S} PT {m}'
    print(s, flush=True)
    with _lock:
        open(LOG, 'a').write(s + '\n')


def have(netsta, d):
    if C.day_files(netsta, d.strftime('%Y-%m-%d')):
        return True
    for p in (f'{C.OUT}/hg_repro/data/raw/{netsta}/{d.year}/{d.dayofyear:03d}.mseed',
              f'{C.DATA}/disp1sps_v3/{netsta}/{d.year}/{d.dayofyear:03d}.npz',
              f'{C.OUT}/hg_repro/data/disp1sps/{netsta}/{d.year}/{d.dayofyear:03d}.npz'):
        if os.path.isfile(p):
            return True
    return False


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=600) as r:
                b = r.read()
                with _lock:
                    _b['xfer'] += len(b)
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
    from obspy import read, UTCDateTime, Stream
    netsta, cha, day = job
    if STOP.is_set():
        return netsta, day, 'stopped'
    net, sta = netsta.split('.')
    out = f'{RAW}/{netsta}/{day.year}/{day.dayofyear:03d}.mseed'
    if os.path.isfile(out) or os.path.isfile(out + '.nd'):
        return netsta, day, 'exists'
    a = UTCDateTime(str(day.date())); b = a + 86400
    pad = 600 if cha[:2] in ('BH', 'HH') else 0
    url = f'{DS}?net={net}&sta={sta}&loc=*&cha={cha}&starttime={(a - pad).isoformat()[:19]}&endtime={(b + pad).isoformat()[:19]}'
    code, raw = get(url)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not raw:
        if code in (204, 404):
            open(out + '.nd', 'w').close()
        return netsta, day, f'nodata {code}'
    st = read(io.BytesIO(raw))
    if cha[:2] in ('BH', 'HH'):
        for tr in st:                       # snap near-identical rates (100.0000076 etc.)
            sr = tr.stats.sampling_rate
            if abs(sr - round(sr)) < 1e-3 * sr:
                tr.stats.sampling_rate = float(round(sr))
        # 2026-10-15: some days mix 40 and 50 Hz pieces (sensor/digitizer change) -> merge per rate, never across rates
        parts = []
        for sr in sorted({tr.stats.sampling_rate for tr in st}):
            sub = Stream([tr for tr in st if tr.stats.sampling_rate == sr])
            sub.merge(method=1, fill_value=None); parts += sub.split().traces
        st = Stream(parts)
        keep = []
        for tr in st:
            sr = tr.stats.sampling_rate
            f = {40: (8, 5), 50: (10, 5), 20: (4, 5), 100: (10, 10), 200: (10, 20)}.get(int(round(sr)))
            if f is None or tr.stats.npts < sr * 1200:
                continue
            tr.data = tr.data.astype('float64'); tr.detrend('demean')
            tr.decimate(f[0], no_filter=False, strict_length=False); tr.decimate(f[1], no_filter=False, strict_length=False)
            tr.data = np.round(tr.data).astype('int32')
            keep.append(tr)
        st = Stream(keep).trim(a, b - 1e-3, nearest_sample=False)
    else:
        st = st.trim(a, b - 1e-3, nearest_sample=False)
    if len(st) == 0:
        return netsta, day, 'empty'
    tmp = out + '.part'
    st.write(tmp, format='MSEED', encoding='STEIM2')
    os.replace(tmp, out)
    n = os.path.getsize(out)
    with _lock:
        _b['stored'] += n
        if _b['stored'] > CAP:
            STOP.set()
    return netsta, day, f'ok {len(st)}tr {len(raw) / 1e6:.1f}MB->{n / 1e6:.2f}MB'


def main():
    jobs = []
    for netsta, (cha, start) in STAS.items():
        for d in pd.date_range(start, END, freq='D'):
            if not have(netsta, d):
                jobs.append((netsta, cha, d))
    jobs.sort(key=lambda j: (j[2], j[0]))                     # oldest days first, so 2021 is complete first
    cnt = pd.Series([j[0] for j in jobs]).value_counts()
    log('start (user approval 2026-10-16 ~6 pm PT): missing days ' + ' '.join(f'{k}:{v}' for k, v in cnt.items()) + f' total {len(jobs)}')
    with ThreadPoolExecutor(2) as ex:
        for netsta, day, r in ex.map(day_job, jobs):
            log(f'{netsta} {day:%Y-%m-%d} {r} | stored {_b["stored"] / 1e9:.3f} GB, transferred {_b["xfer"] / 1e9:.2f} GB')
    log(f'END stored {_b["stored"] / 1e9:.3f} GB, transferred {_b["xfer"] / 1e9:.2f} GB' + (' — STOPPED at 3 GB cap' if STOP.is_set() else ''))
    if not STOP.is_set():
        open(f'{CAT}/FETCH_Y12_DONE', 'w').write(f'{dt.datetime.now(dt.timezone.utc).astimezone(dt.timezone(dt.timedelta(hours=-7))):%Y-%m-%d %H:%M} PT; '
                                                  f'stored {_b["stored"] / 1e9:.3f} GB\n')


if __name__ == '__main__':
    main()
