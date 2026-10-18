#!/usr/bin/env python
"""E5 step 5 STEP A verification -- per station: days with data, no-data days, gap fraction, sample rate, channel codes and the
Z/N/E resolution by e4_chanorder (PB EH1/EH2 -> N/E).  Reads the mseed headers only.  -> results/e5r_download_table.csv"""
import glob, os, sys
import numpy as np, pandas as pd
from obspy import read
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e4_chanorder as co
ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"
S = pd.read_csv(f"{ROOT}/results/e5r_stations.csv")
rows = []
for r in S.itertuples():
    d = f"{ROOT}/traces/{r.netsta}/2012"
    files = sorted(glob.glob(f"{d}/*.mseed")); nd = sorted(glob.glob(f"{d}/*.mseed.nd"))
    chans, fs, cov, nfrag, mb = set(), set(), [], [], 0.0
    for f in files:
        st = read(f, headonly=True); mb += os.path.getsize(f) / 1e6
        chans |= {tr.stats.channel for tr in st}; fs |= {float(tr.stats.sampling_rate) for tr in st}
        per = {}
        for tr in st:
            per[tr.stats.channel] = per.get(tr.stats.channel, 0) + tr.stats.npts / tr.stats.sampling_rate
        cov.append(min(per.values()) / 86400.0); nfrag.append(len(st) / max(len(per), 1))
    try:
        order = co.resolve_order(chans, station=r.netsta)
    except Exception as e:
        order = f"ERR {e}"
    rows.append(dict(station=r.netsta, chan_band=r.chan, fs_hz="/".join(f"{x:g}" for x in sorted(fs)), channels=",".join(sorted(chans)),
                     order_ZNE=order if isinstance(order, str) else "/".join(order), days_data=len(files), days_nodata=len(nd),
                     days_expected=52, coverage_median=float(np.median(cov)) if cov else np.nan,
                     coverage_min=float(np.min(cov)) if cov else np.nan, days_cov_lt_0p9=int(sum(c < 0.9 for c in cov)),
                     frag_per_chan_median=float(np.median(nfrag)) if nfrag else np.nan, size_MB=round(mb, 1),
                     nodata_doys=",".join(os.path.basename(x)[:3] for x in nd)))
T = pd.DataFrame(rows)
T.to_csv(f"{ROOT}/results/e5r_download_table.csv", index=False)
pd.set_option("display.width", 250)
print(T.to_string(index=False))
print(f"TOTAL {T.size_MB.sum()/1e3:.2f} GB, {T.days_data.sum()} station-days with data, {T.days_nodata.sum()} no-data")
