#!/usr/bin/env python
"""Step 2 - the array-era SITE table (1 km single-linkage), the unit of geometry for this study.

results/sites_array_era_final.csv: site_id, members, n_members, n_members_3c, lat, lon,
slab2_depth_km, nets, fs (set), n_components, sensor_class, chain, ehz_only/rational/short_period
flags, and the DELIVERED alive-days of the site in each array-era window (2007-01, 2008-05,
2010-08, 2011-08) taken from results/download_manifest.csv.
Also results/sites_by_window.csv - 3-C sites per network per window, which REPLACES PLAN.md
section 1's station counts.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_common as D                                                   # noqa: E402

WIN = {"w2007_01": ("2007-01-01", "2007-02-28"),
       "w2008_05": ("2008-04-15", "2008-06-10"),
       "w2010_08": ("2010-07-20", "2010-09-10"),
       "w2011_08": ("2011-07-20", "2011-09-10")}
BROADBAND = ("TRILLIUM", "STS-", "STS1", "STS2", "CMG", "GURALP", "KS-54000", "BROADBAND",
             "VELOCITY TRANSDUCER", "T120", "TRILLIUM120", "NANOMETRICS")
ACCEL = ("EPISENSOR", "ACCELEROMETER", "FBA", "KINEMETRICS ES", "TSA-")


def sensor_class(sensor, band, short_period):
    s = (sensor or "").upper()
    if any(k in s for k in ACCEL) and not any(k in s for k in BROADBAND):
        return "accelerometer"
    if short_period:
        return "short_period"
    if any(k in s for k in BROADBAND):
        return "broadband"
    return "unknown"


def main():
    F = D.station_flags()
    F["key"] = F.net + "." + F.sta
    M = pd.read_csv(D.MANIFEST, keep_default_na=False, low_memory=False)
    M = M[(M.status == "data")].drop_duplicates(["net", "sta", "day"])
    M["key"] = M.net + "." + M.sta
    # DATA-level 3-C audit: the StationXML can declare EH1/EH2/EHZ while the archive only ever
    # returns EHZ (ZK: 1604 of 1642 delivered station-days carry ONE channel).  Everything that
    # counts 3-C below counts DELIVERED 3-C, not declared 3-C.
    M["nc"] = pd.to_numeric(M.n_channels, errors="coerce")
    M["d3c"] = M.nc >= 3
    # --- 1 km single-linkage over EVERY array-era station in the box (EHZ-only included, but
    #     every 3-C count below uses only the 3-C members).
    sm = D.cluster_sites(F[["key", "lat", "lon"]], key="key")
    F["site_name"] = [sm.get(k, k.split(".")[-1]) for k in F.key]
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    sl = pd.read_csv(D.SLAB2, usecols=["lat", "lon", "depth", "etype"], low_memory=False)
    sl = sl[~sl.etype.isin(["BA", "TO"])]
    sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna(subset=["lat", "lon", "depth"])
    L = LinearNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
    N = NearestNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
    rows = []
    for site, g in F.groupby("site_name"):
        g3 = g[g.three_c]
        keys3 = set(g3.key)
        rec = dict(site_name=site, members=";".join(sorted(g.key)), n_members=len(g),
                   n_members_3c=len(g3),
                   members_3c=";".join(sorted(g3.key)),
                   lat=float(g.lat.mean()), lon=float(g.lon.mean()),
                   nets=",".join(sorted(set(g.net))),
                   net3c=",".join(sorted(set(g3.net))) if len(g3) else "",
                   fs=",".join(str(int(v)) for v in sorted(set(g3.fs))) if len(g3) else
                      ",".join(str(int(v)) for v in sorted(set(g.fs))),
                   n_components=int(g3.n_comp.max()) if len(g3) else int(g.n_comp.max()),
                   chain=",".join(sorted(set(g3.chain))) if len(g3) else "",
                   ehz_only=bool(len(g3) == 0),
                   rational=bool(len(g3) and g3.rational_rate.any()),
                   short_period=bool(g.short_period.all()),
                   sensor=(g3.sensor.mode().iloc[0] if len(g3) and len(g3.sensor.mode())
                           else (g.sensor.mode().iloc[0] if len(g.sensor.mode()) else "")))
        rec["sensor_class"] = sensor_class(rec["sensor"], "", rec["short_period"])
        for w, (a, b) in WIN.items():
            mm = M[(M.day >= a) & (M.day <= b)]
            rec[f"{w}_days_any"] = int(mm[mm.key.isin(set(g.key))].day.nunique())
            m3 = mm[mm.key.isin(keys3) & mm.d3c]
            rec[f"{w}_days_3c"] = int(m3.day.nunique())
            rec[f"{w}_members_3c"] = int(m3.key.nunique())
            rec[f"{w}_members_declared3c_zonly_data"] = int(
                mm[mm.key.isin(keys3) & ~mm.d3c].key.nunique()
                - mm[mm.key.isin(keys3) & mm.d3c].key.nunique() if False else
                len(set(mm[mm.key.isin(keys3) & ~mm.d3c].key) - set(m3.key)))
        rows.append(rec)
    S = pd.DataFrame(rows)
    S["slab2_depth_km"] = np.where(np.isfinite(L(np.column_stack([S.lon, S.lat]))),
                                   L(np.column_stack([S.lon, S.lat])),
                                   N(np.column_stack([S.lon, S.lat])))
    S = S.sort_values(["lat", "lon"]).reset_index(drop=True)
    S.insert(0, "site_id", [f"S{i:03d}" for i in range(len(S))])
    cols = ["site_id", "site_name", "members", "members_3c", "n_members", "n_members_3c",
            "lat", "lon", "slab2_depth_km", "nets", "net3c", "fs", "n_components",
            "sensor", "sensor_class", "chain", "ehz_only", "rational", "short_period"] + \
           [f"{w}_{k}" for w in WIN
            for k in ("days_any", "days_3c", "members_3c", "members_declared3c_zonly_data")]
    S[cols].to_csv(os.path.join(D.RES, "sites_array_era_final.csv"), index=False)
    # ---- per-window 3-C site counts by network (REPLACES PLAN section 1)
    out = []
    for w in WIN:
        alive = S[(S[f"{w}_days_3c"] > 0)]
        for net, g in alive.groupby("net3c"):
            out.append(dict(window=w, net=net, n_sites_3c=len(g),
                            n_stations_3c=int(g[f"{w}_members_3c"].sum()),
                            median_days=float(g[f"{w}_days_3c"].median())))
    W = pd.DataFrame(out)
    W.to_csv(os.path.join(D.RES, "sites_by_window.csv"), index=False)
    print(S[cols[:13]].head(3).to_string(index=False))
    print("\n=== 3-C SITES per window per network (replaces PLAN section 1 station counts) ===")
    print(W.pivot_table(index="net", columns="window", values="n_sites_3c",
                        aggfunc="sum").fillna(0).astype(int).to_string())
    print("\n=== 3-C STATIONS per window per network ===")
    print(W.pivot_table(index="net", columns="window", values="n_stations_3c",
                        aggfunc="sum").fillna(0).astype(int).to_string())
    xg = S[S.net3c == "XG"]
    print(f"\nXG: {int(xg.n_members_3c.sum())} 3-C stations -> {len(xg)} sites; "
          f"fs {sorted(set(','.join(xg.fs).split(',')))}; chain {sorted(set(xg.chain))}")
    print(xg[["site_id", "site_name", "lat", "lon", "n_members_3c", "fs",
              "w2010_08_days_3c", "w2011_08_days_3c"]].to_string(index=False))


if __name__ == "__main__":
    main()
