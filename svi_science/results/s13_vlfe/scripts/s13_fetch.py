#!/usr/bin/env python
"""S13 fetch: day-file downloads from EarthScope FDSN dataselect, cached and resumable.

usage: s13_fetch.py LIST.csv [--parallel 3] [--tag NAME]
LIST.csv columns: netsta, date (YYYY-MM-DD), cha (e.g. 'BH?,HH?' or 'LH?').
Writes data/raw/NET.STA/YYYY/DOY.mseed (or DOY.mseed.nd when the service has no data). Existing files are skipped.
Writes data/raw/FETCH_<tag>_DONE when every row has a file or an .nd marker. Max 3 parallel requests."""
import os, sys, time, argparse, subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

RAW = f'{C.DATA}/raw'
URL = 'https://service.earthscope.org/fdsnws/dataselect/1/query'


def one(row):
    net, sta = row.netsta.split('.')
    d = pd.Timestamp(row.date)
    out = f'{RAW}/{row.netsta}/{d.year}/{d.dayofyear:03d}.mseed'
    if os.path.isfile(out) or os.path.isfile(out + '.nd'):
        return 'skip'
    os.makedirs(os.path.dirname(out), exist_ok=True)
    t0 = d.strftime('%Y-%m-%dT00:00:00'); t1 = (d + pd.Timedelta(days=1)).strftime('%Y-%m-%dT00:00:00')
    q = f'{URL}?net={net}&sta={sta}&loc=*&cha={row.cha}&start={t0}&end={t1}&nodata=404'
    tmp = out + '.part'
    for attempt in range(4):
        r = subprocess.run(['curl', '-sL', '--max-time', '600', '-o', tmp, '-w', '%{http_code}', q],
                           capture_output=True, text=True)
        code = r.stdout.strip()
        if code == '200' and os.path.getsize(tmp) > 0:
            os.replace(tmp, out); return f'ok {os.path.getsize(out)/1e6:.1f} MB'
        if code in ('204', '404'):
            if os.path.exists(tmp):
                os.remove(tmp)
            open(out + '.nd', 'w').close(); return 'nodata'
        time.sleep(10 * (attempt + 1))
    return f'fail {code}'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('list'); ap.add_argument('--parallel', type=int, default=3); ap.add_argument('--tag', default='list')
    a = ap.parse_args()
    L = pd.read_csv(a.list, dtype=str)
    lf = f'{C.LOGD}/fetch_{a.tag}.log'
    n_fail = 0
    with ThreadPoolExecutor(min(a.parallel, 3)) as ex:
        for row, res in zip(L.itertuples(), ex.map(one, L.itertuples())):
            C.log(f'{row.netsta} {row.date} {res}', lf)
            n_fail += res.startswith('fail')
    if n_fail == 0:
        open(f'{RAW}/FETCH_{a.tag}_DONE', 'w').write(f'{len(L)} rows\n')
    C.log(f'finished {len(L)} rows, {n_fail} failed', lf)
