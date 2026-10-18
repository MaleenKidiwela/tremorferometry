#!/usr/bin/env python
"""Stage 0: families, permanent stations, sites, station-days, array-era truth.

  -> results/setup/families.csv     329 dup_group representatives of the WA regions
  -> results/setup/stations.csv     permanent stations (anchors of NW_WA/SW_WA/OLY2010 + bridge-chain
                                    carriers + permanent clear stations), with 3-km sites
  -> results/setup/station_days.csv one row per (station, day) with the raw file on disk if any
  -> results/setup/truth.parquet    the families' own array-era detections inside the windows
  -> results/setup/scan_days.csv    days each source run scanned (any family detected that day)
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

os.makedirs(W.SETUP, exist_ok=True)
F = pd.read_csv(f"{W.ROOT}/data/catalogue/families_v1_1.csv")
F = F[F.region.isin(W.REGIONS) & (F.representative != False)].copy()   # noqa: E712
F["family_id"] = F.family_id.astype(int)
F["lat_use"] = F.pos_lat
F["lon_use"] = F.pos_lon
F["val_windows"] = F.region.map({"Royer_NW_WA": "V2008,V2010", "Royer_SW_WA": "V2008",
                                 "Olympic_2010": "V2010", "Olympic_2010_armB": "V2010"})
F["fid"] = np.arange(len(F))
F[["fid", "family_uid", "region", "family_id", "det_source", "template_file", "n_detections",
   "era_start", "era_end", "lat_use", "lon_use", "located", "location_quality", "perm_carriers",
   "n_perm_carriers", "perm_weak", "val_windows"]].to_csv(f"{W.SETUP}/families.csv", index=False)
print("families", len(F), F.region.value_counts().to_dict())

# ---------------------------------------------------------------- stations
anch = set()
for r in ["NW_WA", "SW_WA", "OLY2010"]:
    anch |= set(pd.read_csv(f"{W.EXT}/{r}_anchors.csv").site)
extra = set()
for c in ["perm_carriers", "perm_weak", "clear_stations"]:
    src = pd.read_csv(f"{W.ROOT}/data/catalogue/families_v1_1.csv")
    src = src[src.region.isin(W.REGIONS) & (src.representative != False)]   # noqa: E712
    for v in src[c].dropna():
        extra |= {x.split(":")[0] for x in v.split(";")}
C = pd.concat([pd.read_csv(f"{W.EXT}/{r}_channels.csv") for r in ["NW_WA", "SW_WA", "OLY2010"]])
C = C[~C.net.isin(W.TEMP_NETS)]
known = C.drop_duplicates(["net", "sta"])[["net", "sta", "lat", "lon"]]
rows = []
for r in known.itertuples():
    ns = f"{r.net}.{r.sta}"
    if ns in anch or r.sta in extra:
        rows.append(dict(netsta=ns, net=r.net, sta=r.sta, lat=r.lat, lon=r.lon,
                         why="anchor" if ns in anch else "carrier/clear"))
S = pd.DataFrame(rows).drop_duplicates("netsta").sort_values("netsta").reset_index(drop=True)
# 3-km sites, single linkage
site = list(range(len(S)))
for i in range(len(S)):
    for j in range(i + 1, len(S)):
        if W.km(S.lat[i], S.lon[i], S.lat[j], S.lon[j]) < W.SITE_KM:
            a, b = site[i], site[j]
            site = [a if x == b else x for x in site]
S["site_ix"] = site
S["site"] = S.groupby("site_ix").sta.transform(lambda v: "+".join(sorted(v)))
S.drop(columns="site_ix").to_csv(f"{W.SETUP}/stations.csv", index=False)
print("stations", len(S), "sites", S.site.nunique())
print(S[S.site.str.contains(r"\+")][["netsta", "site"]].to_string())

# ---------------------------------------------------------------- station-days
N = pd.concat([pd.read_csv(f"{W.EXT}/{r}_bridge_nodes_strict.csv") for r in ["NW_WA", "SW_WA", "OLY2010"]]
              ).drop_duplicates("node")
rows = []
for w, d, k in W.all_days():
    for ns in S.netsta:
        net, sta = ns.split(".")
        n = N[(N.net == net) & (N.sta == sta) & (N.start <= d) & (N.end > d)]
        rp = W.raw_path(ns, d)
        rows.append(dict(window=w, date=d, kind=k, netsta=ns, meta_live=len(n) > 0,
                         meta_nodes=";".join(n.node), raw=rp or ""))
SD = pd.DataFrame(rows)
SD.to_csv(f"{W.SETUP}/station_days.csv", index=False)
live = SD[SD.meta_live]
print("station-days", len(SD), "metadata-live", len(live), "on disk", (live.raw != "").sum(),
      "to fetch", (live.raw == "").sum())

# ---------------------------------------------------------------- truth
days = sorted({d for _, d, _ in W.all_days()})
T, SC = [], []
for src, g in F.groupby("det_source"):
    D = pd.read_parquet(f"{W.ROOT}/{src}", columns=["grp", "t", "date"])
    D["date"] = D.date.astype(str)
    SC.append(pd.DataFrame(dict(det_source=src, date=sorted(set(D.date)))))
    D = D[D.date.isin(days)]
    m = dict(zip(g.family_id, g.fid))
    D = D[D.grp.isin(m)].copy()
    D["fid"] = D.grp.map(m)
    D["det_source"] = src
    T.append(D[["fid", "t", "date", "det_source"]])
T = pd.concat(T, ignore_index=True)
T["half"] = W.half_of(T.t.to_numpy())
T.to_parquet(f"{W.SETUP}/truth.parquet", index=False)
pd.concat(SC).to_csv(f"{W.SETUP}/scan_days.csv", index=False)
print("truth rows", len(T), T.groupby(T.date.str[:4]).size().to_dict(), "halves", T.half.value_counts().to_dict())
