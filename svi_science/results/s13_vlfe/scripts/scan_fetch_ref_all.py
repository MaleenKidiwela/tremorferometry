#!/usr/bin/env python
"""Far-reference download 2003-2026 (user approval relayed by the coordinator 2026-10-15 ~1:20 pm PT): 1-sps LHZ
(vertical only) of a fixed long-running far-reference set chosen from catalogue/tables/farref_per_year.csv
(250-800 km from 48.37 N, 123.61 W; >= 3 per year in every year 2003-2026):
  IU.COR (421 km, az 177), US.NEW (480 km, 89), US.HAWA (378 km, 124; to 2024), US.MSO (744 km, 100), US.WVOR (765 km, 148).
No year has < 2 LHZ far stations (catalogue RESULTS §3), so no other data class is fetched. All references lie at azimuth
89-177 (no LHZ station exists north or west of the network).
EarthScope FDSN dataselect, monthly requests, <= 2 in parallel, polite retries (10/20/30/40 s), default Python user
agent, no personal data. Split into day files data/raw_ref_all/NET.STA/YYYY/DOY.mseed (same layout as raw_ref2011, so the
scan_prep_v3 / scan_prep_ref path can read it via TRACE_ROOTS). Month marker .mMM.done per year dir (resumable).
Station-year days already on disk in raw_ref2011/raw_ref2012 (COR, NEW, HAWA 2011-12) are not fetched again.
StationXML (response level, all epochs) -> data/resp_scan/NET.STA.all.xml. Budget: stop if the total exceeds 10 GB.
Log: logs/ref_all_download.log. Marker: REF_ALL_DONE (s13_vlfe/)."""
import os, sys, time, io, datetime as dt, threading
from concurrent.futures import ThreadPoolExecutor
import urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

DS = 'https://service.earthscope.org/fdsnws/dataselect/1/query'
ST = 'https://service.earthscope.org/fdsnws/station/1/query'
STAS = [('IU', 'COR', '00'), ('US', 'NEW', '--,00'), ('US', 'HAWA', '--,00'), ('US', 'MSO', '--,00'), ('US', 'WVOR', '--,00')]
OUT = f'{C.DATA}/raw_ref_all'
RESP = f'{C.DATA}/resp_scan'
LOG = f'{C.LOGD}/ref_all_download.log'
END = dt.datetime(2026, 10, 15)
ON_DISK = {('IU.COR', 2011), ('US.NEW', 2011), ('US.HAWA', 2011), ('IU.COR', 2012), ('US.NEW', 2012), ('US.HAWA', 2012)}
BUDGET = 10e9
TOT = [0]
LOCK = threading.Lock()
STOP = threading.Event()


def log(m):
    s = f'{dt.datetime.now(dt.timezone.utc).astimezone(dt.timezone(dt.timedelta(hours=-7))):%Y-%m-%d %H:%M:%S} PT {m}'
    print(s, flush=True)
    with LOCK:
        open(LOG, 'a').write(s + '\n')


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=600) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                return e.code, b''
            log(f'HTTP {e.code} {url} (try {k + 1})')
        except Exception as e:  # noqa
            log(f'ERR {type(e).__name__} {url} (try {k + 1})')
        time.sleep(10 * (k + 1))
    return None, b''


def month(job):
    from obspy import read, UTCDateTime
    net, sta, loc, y, m = job
    if STOP.is_set():
        return 0
    d = f'{OUT}/{net}.{sta}/{y}'
    mk = f'{d}/.m{m:02d}.done'
    if os.path.isfile(mk):
        return 0
    t0 = dt.datetime(y, m, 1); t1 = min(dt.datetime(y + (m == 12), m % 12 + 1, 1), END)
    url = f'{DS}?net={net}&sta={sta}&loc={loc}&cha=LHZ&starttime={t0:%Y-%m-%dT%H:%M:%S}&endtime={t1:%Y-%m-%dT%H:%M:%S}&nodata=404'
    code, b = get(url)
    os.makedirs(d, exist_ok=True)
    if code is None:
        log(f'{net}.{sta} {y}-{m:02d}: FAILED after retries (no marker; rerun resumes)'); return 0
    if not b:
        log(f'{net}.{sta} {y}-{m:02d}: no data (code {code})')
        open(mk, 'w').write('nodata\n'); return 0
    st = read(io.BytesIO(b))
    nday, nb, day = 0, 0, t0
    while day < t1:
        a = UTCDateTime(day); s = st.slice(a, a + 86400 - 1e-3)
        if len(s):
            f = f'{d}/{day.timetuple().tm_yday:03d}.mseed'
            s.write(f + '.part', format='MSEED', encoding='STEIM2' if all(tr.data.dtype.kind == 'i' for tr in s) else None)
            os.replace(f + '.part', f); nday += 1; nb += os.path.getsize(f)
        day += dt.timedelta(days=1)
    open(mk, 'w').write(f'{nday} days {nb} bytes\n')
    with LOCK:
        TOT[0] += nb
        tot = TOT[0]
    log(f'{net}.{sta} {y}-{m:02d}: {len(b) / 1e6:.2f} MB requested, {len(st)} traces, {nday} day files, '
        f'{nb / 1e6:.2f} MB stored; total {tot / 1e9:.2f} GB')
    if tot > BUDGET:
        STOP.set(); log(f'BUDGET: total {tot / 1e9:.2f} GB > 10 GB -> stopping (report to the user)')
    return nb


def main():
    os.makedirs(RESP, exist_ok=True); os.makedirs(OUT, exist_ok=True)
    # bytes already stored (resume)
    for root, _, files in os.walk(OUT):
        TOT[0] += sum(os.path.getsize(os.path.join(root, f)) for f in files if f.endswith('.mseed'))
    log(f'start (user approval relayed 2026-10-15 ~1:20 pm PT): {[f"{n}.{s}" for n, s, _ in STAS]} LHZ 2003-01-01 - '
        f'{END:%Y-%m-%d}; already stored {TOT[0] / 1e9:.2f} GB')
    for net, sta, loc in STAS:
        f = f'{RESP}/{net}.{sta}.all.xml'
        if not os.path.isfile(f):
            code, b = get(f'{ST}?net={net}&sta={sta}&cha=LH?&starttime=2003-01-01T00:00:00&level=response')
            open(f, 'wb').write(b); log(f'StationXML {net}.{sta}: {len(b) / 1e3:.0f} kB (code {code})')
    jobs = [(n, s, l, y, m) for y in range(2003, 2027) for n, s, l in STAS for m in range(1, 13)
            if (f'{n}.{s}', y) not in ON_DISK and dt.datetime(y, m, 1) < END]
    log(f'{len(jobs)} station-months to request')
    with ThreadPoolExecutor(2) as ex:
        list(ex.map(month, jobs))
    nmiss = sum(not os.path.isfile(f'{OUT}/{n}.{s}/{y}/.m{m:02d}.done') for n, s, l, y, m in jobs)
    log(f'finished: total stored {TOT[0] / 1e9:.2f} GB; station-months without marker {nmiss}'
        f'{"; STOPPED at budget" if STOP.is_set() else ""}')
    if nmiss == 0 and not STOP.is_set():
        open(f'{C.OUT}/REF_ALL_DONE', 'w').write(f'{dt.datetime.now(dt.timezone.utc).astimezone(dt.timezone(dt.timedelta(hours=-7))):%Y-%m-%d %H:%M} PT; {TOT[0] / 1e9:.2f} GB\n')


if __name__ == '__main__':
    main()
