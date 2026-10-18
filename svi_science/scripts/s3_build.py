"""S3 build: per-event tables for SVI and NW big ETS episodes + PNSN tremor in the same boxes.

Outputs (svi_science/results/s3/):
  events_svi.parquet, events_nw.parquet  -- one row per candidate in a padded major-episode window
  tremor_svi.parquet, tremor_nw.parquet  -- PNSN epicentres near the families, same frame
  sources_nw.csv, families_svi_frame.csv -- positions in the along-strike frame
  episodes_all.csv                       -- episode windows (padded) and counts
"""
import glob
import json
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s3"
KM = 111.195
SLAB = f"{ROOT}/data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz"
PAD = 3.0
T0_NW, T1_NW = pd.Timestamp("2015-01-01"), pd.Timestamp("2026-09-16")


def slab_grid():
    S = pd.read_csv(SLAB, header=None, names=["lon", "lat", "z"])
    S["lon"] = S.lon - 360.0
    lons = np.sort(S.lon.unique()); lats = np.sort(S.lat.unique())
    Z = S.pivot_table(index="lat", columns="lon", values="z", dropna=False).reindex(index=lats, columns=lons).to_numpy()
    return lons, lats, Z


def _bilin(lons, lats, Z, lo, la):
    i = np.searchsorted(lons, lo) - 1; j = np.searchsorted(lats, la) - 1
    fx = (lo - lons[i]) / (lons[i + 1] - lons[i]); fy = (la - lats[j]) / (lats[j + 1] - lats[j])
    return ((1 - fx) * (1 - fy) * Z[j, i] + fx * (1 - fy) * Z[j, i + 1] + (1 - fx) * fy * Z[j + 1, i] + fx * fy * Z[j + 1, i + 1])


def frame_from(lat, lon):
    """Same rule as s0_build.frame: strike = circular mean Slab2 down-dip az - 90; origin = mean position."""
    lons, lats, Z = slab_grid(); h = 0.05; az = []
    for lo, la in zip(lon, lat):
        dzdx = (_bilin(lons, lats, Z, lo + h, la) - _bilin(lons, lats, Z, lo - h, la)) / (2 * h * KM * np.cos(np.radians(la)))
        dzdy = (_bilin(lons, lats, Z, lo, la + h) - _bilin(lons, lats, Z, lo, la - h)) / (2 * h * KM)
        az.append(np.degrees(np.arctan2(-dzdx, -dzdy)) % 360)
    dd = np.radians(np.array(az))
    mdd = np.degrees(np.arctan2(np.sin(dd).mean(), np.cos(dd).mean())) % 360
    return dict(origin_lat=float(np.mean(lat)), origin_lon=float(np.mean(lon)), downdip_az=float(mdd), strike_az=float((mdd - 90) % 360))


def project(fr, lat, lon):
    x = (np.asarray(lon) - fr["origin_lon"]) * KM * np.cos(np.radians(fr["origin_lat"]))
    y = (np.asarray(lat) - fr["origin_lat"]) * KM
    s = np.radians(fr["strike_az"]); d = np.radians(fr["downdip_az"])
    return x * np.sin(s) + y * np.cos(s), x * np.sin(d) + y * np.cos(d), x, y


def padded_windows(E):
    E = E.sort_values("start").reset_index(drop=True).copy()
    s = pd.to_datetime(E.start); e = pd.to_datetime(E.end) + pd.Timedelta(days=1)   # end date inclusive
    lo, hi = [], []
    for i in range(len(E)):
        gp = (s[i] - e[i - 1]).total_seconds() / 86400 if i > 0 else 99
        gn = (s[i + 1] - e[i]).total_seconds() / 86400 if i < len(E) - 1 else 99
        lo.append(s[i] - pd.Timedelta(days=min(PAD, max(gp, 0) / 2)))
        hi.append(e[i] + pd.Timedelta(days=min(PAD, max(gn, 0) / 2)))
    E["w0"] = lo; E["w1"] = hi
    return E


def assign(t, E):
    ep = np.full(len(t), None, dtype=object)
    tv = t.values
    for _, r in E.iterrows():
        m = (tv >= np.datetime64(r.w0)) & (tv < np.datetime64(r.w1))
        ep[m] = r.episode
    return ep


def near_mask(plat, plon, flat, flon, km=10.0):
    out = np.zeros(len(plat), bool)
    for la, lo in zip(flat, flon):
        d = np.hypot((plon - lo) * KM * np.cos(np.radians(la)), (plat - la) * KM)
        out |= d <= km
    return out


def main():
    info = {}
    # ------------------------------------------------------------------ SVI
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    frs = json.load(open(f"{ROOT}/svi_science/results/s0/frame_groups.json"))
    fr_svi = dict(origin_lat=frs["origin_lat"], origin_lon=frs["origin_lon"], downdip_az=frs["mean_downdip_az"], strike_az=frs["strike_az"])
    F[["fid", "lat", "lon", "dep", "along_strike_km", "downdip_km", "group", "n_carriers"]].to_csv(f"{OUT}/families_svi_frame.csv", index=False)
    Esvi = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/episode_rates/episodes_svi_final.csv")
    Esvi = padded_windows(Esvi)
    T = pd.read_parquet(f"{ROOT}/svi_science/results/s0/analysis_table.parquet",
                        columns=["series", "fam", "k", "P_same", "keep_xfam8", "time_utc", "n_sta_day", "net_mad",
                                 "flag_ep_2017_03", "flag_ep_2025_10"])
    T = T[T.P_same >= 0.5].copy()
    T["ep"] = assign(T.time_utc, Esvi[Esvi.major])
    T = T[T.ep.notna()].copy()
    T = T.merge(F[["fid", "along_strike_km", "downdip_km", "dep", "group", "lat", "lon"]].rename(columns={"fid": "fam"}), on="fam", how="left")
    T = T.rename(columns={"time_utc": "time", "fam": "unit", "dep": "depth_km", "along_strike_km": "as_km", "downdip_km": "dd_km"})
    T["counted_once"] = T.keep_xfam8.astype(bool)
    T["region"] = "SVI"
    T.to_parquet(f"{OUT}/events_svi.parquet")
    info["svi_rows"] = int(len(T))
    # tremor
    P = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon", "depth"])
    P["time"] = pd.to_datetime(P.time).astype("datetime64[ns]")
    m = near_mask(P.lat.values, P.lon.values, F.lat.values, F.lon.values, 10.0)
    Ps = P[m].copy(); Ps["ep"] = assign(Ps.time, Esvi[Esvi.major]); Ps = Ps[Ps.ep.notna()]
    a, d, _, _ = project(fr_svi, Ps.lat, Ps.lon); Ps["as_km"] = a; Ps["dd_km"] = d
    Ps.to_parquet(f"{OUT}/tremor_svi.parquet")
    info["svi_tremor_rows"] = int(len(Ps))

    # ------------------------------------------------------------------ NW sources
    SG = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3.csv")
    AC = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3_array_confirmed.csv").rename(columns={"group": "source", "group_size": "source_size"})
    SG = SG.merge(AC, on="fid")
    best = SG.sort_values(["source", "erh"]).groupby("source").head(1).set_index("source")
    S = pd.DataFrame({"source": best.index, "lat": best.lat_use.values, "lon": best.lon_use.values, "depth_km": best.depth_km.values,
                      "erh": best.erh.values, "quality": best.location_quality.values, "best_fid": best.fid.values,
                      "source_size": best.source_size.values})
    S["members"] = S.source.map(SG.groupby("source").fid.apply(lambda x: ";".join(map(str, sorted(x)))))
    S["member_span_km"] = S.source.map(SG.groupby("source").apply(
        lambda g: float(np.max(np.hypot((g.lon_use.values[:, None] - g.lon_use.values[None]) * KM * np.cos(np.radians(47.9)),
                                        (g.lat_use.values[:, None] - g.lat_use.values[None]) * KM))), include_groups=False))
    S["N_E3_members"] = S.source.map(SG.groupby("source").N_E3.sum())
    good = S.quality.isin(["A", "B"])
    fr_nw = frame_from(S.lat[good].values, S.lon[good].values)
    a, d, x, y = project(fr_nw, S.lat, S.lon); S["as_km"] = a; S["dd_km"] = d
    S.to_csv(f"{OUT}/sources_nw.csv", index=False)
    info["nw_frame"] = fr_nw; info["svi_frame"] = fr_svi
    fid2src = SG.set_index("fid").source

    # ------------------------------------------------------------------ NW events
    C = pd.read_parquet(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet",
                        columns=["det_id", "time", "region", "family", "P_own_pnet", "selected_v2", "n_sta", "net_mad"],
                        filters=[("region", "==", "NW")])
    C = C[(C.time >= T0_NW) & (C.time < T1_NW) & (C.P_own_pnet >= 0.5)].copy()
    fl = []
    for f in sorted(glob.glob(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/dedup/flags/flags_20*.parquet")):
        ym = f.split("_")[-1][:7]
        if ym < "2015-01":
            continue
        q = pd.read_parquet(f, columns=["det_id", "series", "dup_A", "dup_B", "dup_D", "win_km_D"])
        fl.append(q[(q.series == "data") & q.det_id.isin(C.det_id)].drop(columns="series"))
    fl = pd.concat(fl)
    C = C.merge(fl, on="det_id", how="left")
    info["nw_rows_without_flag"] = int(C.dup_D.isna().sum())
    info["nw_selected_without_flag"] = int((C.selected_v2 & C.dup_D.isna()).sum())
    Enw = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/episode_rates/episodes_nw_final.csv")
    Enw = padded_windows(Enw)
    Enw = Enw[(pd.to_datetime(Enw.start) >= T0_NW)]
    C["ep"] = assign(C.time, Enw[Enw.major]); C = C[C.ep.notna()].copy()
    C["source"] = C.family.map(fid2src)
    C = C.merge(S[["source", "as_km", "dd_km", "depth_km", "quality", "lat", "lon"]], on="source", how="left")
    C = C.rename(columns={"family": "unit_fid"})
    C["unit"] = C.source
    C["counted_once"] = ~C.dup_D.fillna(False).astype(bool)
    C["P_same"] = C.P_own_pnet
    C["series"] = "data"
    C["region"] = "NW"
    C.to_parquet(f"{OUT}/events_nw.parquet")
    info["nw_rows"] = int(len(C))
    m = near_mask(P.lat.values, P.lon.values, S.lat[good].values, S.lon[good].values, 10.0)
    Pn = P[m].copy(); Pn["ep"] = assign(Pn.time, Enw[Enw.major]); Pn = Pn[Pn.ep.notna()]
    a, d, _, _ = project(fr_nw, Pn.lat, Pn.lon); Pn["as_km"] = a; Pn["dd_km"] = d
    Pn.to_parquet(f"{OUT}/tremor_nw.parquet")
    info["nw_tremor_rows"] = int(len(Pn))

    # ------------------------------------------------------------------ episode table
    rows = []
    for reg, E, ev in [("SVI", Esvi[Esvi.major], T), ("NW", Enw[Enw.major], C)]:
        for _, r in E.iterrows():
            g = ev[ev.ep == r.episode]
            hd = g[(g.series == "data") & (g.P_same >= 0.75)]
            if reg == "NW":
                hd = hd[hd.selected_v2]
            rows.append(dict(region=reg, episode=r.episode, start=r.start, end=r.end, w0=r.w0, w1=r.w1,
                             n_headline=int(hd.counted_once.sum()), n_raw075=int(len(hd)),
                             n_once050=int(g[(g.series == "data") & g.counted_once].shape[0]),
                             n_c1_075=int(((g.series == "c1") & (g.P_same >= 0.75)).sum())))
    EP = pd.DataFrame(rows); EP["big"] = EP.n_headline >= 300
    EP.to_csv(f"{OUT}/episodes_all.csv", index=False)
    info["big"] = {reg: EP[(EP.region == reg) & EP.big].episode.tolist() for reg in ["SVI", "NW"]}
    json.dump(info, open(f"{OUT}/build_info.json", "w"), indent=1, default=str)
    print(json.dumps(info, indent=1, default=str))
    print(EP[EP.big].to_string())


if __name__ == "__main__":
    main()
