#!/usr/bin/env python
"""EXTENSION (user-approved, clearly separate from the literal reproduction): the station pool for
growth and location = every public station of the listed networks operating within EXT_RMAX km of
the subarray (the centroid of its 7 sec.3.1 autocorrelation stations).  Radii 75 / 100 / 150 km are
run as separate columns; this table holds the 150 km superset.

Networks: UW (PNSN), PB (PBO), TA (USArray), XU (CAFE), CN (CNSN) -- the coordinator's list -- plus
the other public permanent/regional networks present in the FDSN inventory: US, CC, C8, GS, YW.
EXCLUDED, stated: dense small-aperture temporary arrays XG (Array of Arrays), ZK, Z3 (their
elements are not independent stations and would dominate every ranking) and ocean-bottom 7D / ZU.

  python scripts/ab_ext_pool.py  ->  results/ab_royer/ext/pool_{NW,SW}.csv, station_days_{NW,SW}.csv
"""
from __future__ import annotations

import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

NETS = ["UW", "PB", "TA", "XU", "CN", "US", "CC", "C8", "GS", "YW"]
EXCLUDED = ["XG", "ZK", "Z3", "7D", "ZU"]
EXT_RMAX = 150.0
RADII = [75.0, 100.0, 150.0]


def main():
    os.makedirs(A.EXT_DIR, exist_ok=True)
    F = pd.read_csv(f"{A.STUDY}/results/fdsn_stations_40-50N_2003-2014.txt", sep="|")
    F.columns = [c.strip("# ").strip() for c in F.columns]
    F = F[F.Network.isin(NETS)].copy()
    F["start"] = F.StartTime.astype(str).str[:10]
    F["end"] = F.EndTime.fillna("2599-12-31").astype(str).str[:10].replace("nan", "2599-12-31")
    for sub in ("NW", "SW", "NVI"):
        ss = A.SUBARRAY[sub]["stas"]
        la = sum(A.COORD[s][0] for s in ss) / len(ss)
        lo = sum(A.COORD[s][1] for s in ss) / len(ss)
        G = F.copy()
        G["dist_km"] = A.epi_km(la, lo, G.Latitude, G.Longitude)
        days = A.scan_days(sub)
        d0, d1 = min(days), max(days)
        G = G[(G.dist_km <= EXT_RMAX) & (G.start <= d1) & (G.end >= d0)]
        # one network per station code: the literal pool's network wins, else the nearest
        G = G.sort_values("dist_km")
        rows = []
        for sta, gg in G.groupby("Station"):
            if any(sta in P_ for P_ in A.POOL.values()) or sta in A.SUBARRAY[sub]["stas"]:
                gg = gg[gg.Network == A.NET[sta]] if (gg.Network == A.NET[sta]).any() else gg
            net = gg.iloc[0].Network
            for r in gg[gg.Network == net].itertuples():
                rows.append(dict(sta=sta, net=net, lat=r.Latitude, lon=r.Longitude,
                                 elev=r.Elevation, start=r.start, end=r.end,
                                 dist_km=round(float(r.dist_km), 2),
                                 in_literal_pool=sta in A.POOL[sub],
                                 dup_codes=",".join(sorted(set(gg.Network) - {net}))))
        E = pd.DataFrame(rows)
        E.to_csv(f"{A.EXT_DIR}/pool_{sub}.csv", index=False)
        sd = []
        for d in days:
            ok = E[(E.start <= d) & (E.end >= d)].sta.unique()
            sd += [dict(sta=s, date=d) for s in ok]
        SD = pd.DataFrame(sd)
        SD.to_csv(f"{A.EXT_DIR}/station_days_{sub}.csv", index=False)
        u = E.drop_duplicates("sta")
        A.log(f"[ext pool {sub}] centroid {la:.3f} {lo:.3f}: "
              + ", ".join(f"<= {int(R)} km: {int((u.dist_km <= R).sum())} stations"
                          for R in RADII)
              + f"; {len(SD)} operating station-days over {len(days)} scan days; networks "
              + str(u.net.value_counts().to_dict()))


if __name__ == "__main__":
    main()
