#!/usr/bin/env python
"""S13 step 1: inventory of continuous waveform data already on disk.

Walks every trace root (ROOT/NET.STA/YEAR/DOY.mseed), counts day files per
station-year (skips *.nd no-data markers), sums bytes, and reads the header of
the first and a middle file per station-year to record channel codes, sampling
rate, sample dtype/encoding and covered seconds (full day or cut window).

Output: tables/inventory_station_year.csv  (one row per root x station x year)
        tables/inventory_station.csv       (merged over roots, best channel)
"""
import os, sys, glob
import pandas as pd
import obspy

REPO = '/home/jovyan/tremorferometry'
ROOTS = ['lfe_extension/traces', 'lfe_extension/traces_backward',
         'lfe_extension_wa/traces', 'lfe_extension_wa/traces_dl',
         'lfe_extension_wa/traces_dl2', 'lfe_catalog_repair/traces',
         'lfe_array_discovery/traces', 'lfe_array_discovery/traces_eq',
         'lfe_array_discovery/traces_sealed', 'data/waveforms']
OUT = os.path.join(REPO, 'svi_science/results/s13_vlfe/tables')


def head(path):
    try:
        st = obspy.read(path, headonly=True)
    except Exception as e:  # noqa
        return dict(err=str(e)[:60])
    chans = sorted({tr.stats.channel for tr in st})
    sr = sorted({float(tr.stats.sampling_rate) for tr in st})
    secs = max((tr.stats.endtime - tr.stats.starttime) for tr in st) if len(st) else 0
    enc = sorted({str(tr.stats.mseed.encoding) for tr in st if 'mseed' in tr.stats})
    locs = sorted({tr.stats.location for tr in st})
    return dict(chans=','.join(chans), sr=','.join(f'{s:g}' for s in sr),
                secs=round(secs), enc=','.join(enc), locs=','.join(locs),
                ntr=len(st))


rows = []
for root in ROOTS:
    rp = os.path.join(REPO, root)
    if not os.path.isdir(rp):
        continue
    for sta in sorted(os.listdir(rp)):
        sp = os.path.join(rp, sta)
        if not os.path.isdir(sp):
            continue
        for yr in sorted(os.listdir(sp)):
            yp = os.path.join(sp, yr)
            if not os.path.isdir(yp):
                continue
            files, nbytes = [], 0
            with os.scandir(yp) as it:
                for e in it:
                    if e.is_file() and not e.name.endswith('.nd') and e.name.endswith(('.mseed', '.ms', '.MSEED')):
                        sz = e.stat().st_size
                        if sz > 0:
                            files.append(e.name)
                            nbytes += sz
            if not files:
                continue
            files.sort()
            h1 = head(os.path.join(yp, files[len(files) // 2]))
            rows.append(dict(root=root, sta=sta, year=yr, n_files=len(files),
                             gb=nbytes / 1e9, mb_per_file=nbytes / 1e6 / len(files),
                             first=files[0], last=files[-1], **h1))
    print(root, len(rows), flush=True)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, 'inventory_station_year.csv'), index=False)
print('wrote', len(df))
open(os.path.join(REPO, 'svi_science/results/s13_vlfe/INVENTORY_DONE'), 'w').write('ok\n')
