#!/usr/bin/env python
"""Round-4 artifact: per-station USED-BAND physical instrument-era boundary table. For each station's used band
(pick_band logic: high-gain seismometer 2nd-char H, most span, prefer BH/HH native-40), fetch FDSN level=response
and declare a boundary ONLY where the response POLES/ZEROS physically differ between consecutive epochs (sensor/
datalogger change). Gain/sensitivity-only revisions are NON-events (pipeline is gain-immune). Merge boundaries
<30 days. Output res_catalog/era_table.csv. Everything downstream (refined scan, era-split inversion) uses this."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from obspy.clients.fdsn import Client
from collections import defaultdict
OUT = "fault_tomography/inversion/res_catalog"
pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5]
tags = sorted(pm.tag.unique())
fo = pd.read_csv("data/broadband_fleet_order.csv"); net_of = dict(zip(fo.sta, fo.net))
cand = pd.read_csv("data/candidate_stations_post2020.csv")
for _, r in cand.iterrows(): net_of.setdefault(r.sta, r.get("net", "PB"))
def sta_net(tag):
    low = tag.lower()
    for pre, s in [("pgc", "PGC"), ("shb", "SHB"), ("clrs", "CLRS")]:
        if low.startswith(pre): return s, "CN"
    sta = tag.replace("p90f40", "").upper(); return sta, net_of.get(sta, "PB" if sta[0] == "B" and sta[1:2].isdigit() else "UW")
def prov(net): return "NCEDC" if net in ("BK", "NC") else "IRIS"
def paz(c):
    try:
        p = c.response.get_paz(); return (tuple(np.round(np.sort_complex(np.asarray(p.poles)), 2)),
                                          tuple(np.round(np.sort_complex(np.asarray(p.zeros)), 2)))
    except Exception:
        return None
rows = []
for tag in tags:
    sta, net = sta_net(tag)
    try:
        inv = Client(prov(net), timeout=45).get_stations(network=net, station=sta, level="response")
    except Exception:
        rows.append(dict(tag=tag, sta=sta, net=net, used_band="FDSN_ERR", n_bnd=-1, boundaries="")); continue
    chs = [c for c in inv[0][0].channels if c.code.endswith("Z") and len(c.code) == 3 and c.code[1] == "H"]
    span = defaultdict(float); n40 = {}
    for c in chs:
        s = c.start_date.year; e = c.end_date.year if c.end_date else 2026
        b = c.code[:2]; span[b] += max(0, min(e, 2026)-max(s, 2010)); n40[b] = abs(c.sample_rate-40) < 1
    if not span:
        rows.append(dict(tag=tag, sta=sta, net=net, used_band="NO_H_BAND", n_bnd=-1, boundaries="")); continue
    used = max(span, key=lambda b: span[b] + (2 if n40.get(b) else 0) + (100 if b[0] in "BH" else 0))
    uch = sorted([c for c in chs if c.code[:2] == used], key=lambda c: c.start_date)
    bset = []; prev_pz = None
    for c in uch:
        pz = paz(c)
        if prev_pz is not None and pz is not None and pz != prev_pz:
            bd = pd.Timestamp(c.start_date.datetime)
            if 2009 < bd.year < 2026: bset.append(bd)
        if pz is not None: prev_pz = pz
    merged = []
    for d in sorted(set(bset)):
        if not merged or (d - merged[-1]).days > 30: merged.append(d)
    rows.append(dict(tag=tag, sta=sta, net=net, used_band=used, n_bnd=len(merged),
                     boundaries=";".join(str(d.date()) for d in merged)))
    if merged: print(f"{tag} ({used}): {len(merged)} physical boundaries {[str(d.date()) for d in merged]}")
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/era_table.csv", index=False)
good = R[R.n_bnd >= 0]
print(f"\nstations resolved: {len(good)}/{len(R)} (FDSN/band errors {int((R.n_bnd<0).sum())})")
print(f"stations with >=1 physical mid-record boundary in USED band: {int((good.n_bnd>0).sum())}")
print(f"total physical boundaries: {int(good.n_bnd.sum())} (vs 360 raw all-band epoch starts)")
yc = pd.to_datetime(pd.Series([b for bs in good.boundaries for b in bs.split(';') if b])).dt.year.value_counts().sort_index()
print("USED-BAND physical-boundary histogram by year: " + " ".join(f"{y}:{yc[y]}" for y in yc.index))
print(f"wrote {OUT}/era_table.csv")
