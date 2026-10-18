#!/usr/bin/env python
"""g2_local_build (2026-10-11) -- Study G2 step 1: merged local earthquake catalogue (PLAN.md §1).

ComCat (dl/comcat/*.csv) + NRCan (dl/nrcan_events.csv) -> duplicates removed (|dt| <= 5 s and epicentre <= 25 km;
border authority: US side keeps ComCat, Canadian side keeps NRCan; a fixed/zero depth is replaced by the other
catalogue's free depth), non-earthquake types removed, 15-km volcano circles removed, distance to the families
(SVI 115; WA-NW quality A/B primaries), Slab2 (and McCrory 2012 column) depth at the epicentre, class
crustal / near-interface / intraslab / overriding-deep / offshore / unknown-depth, down-dip / along-strike
coordinates per region, Gardner-Knopoff declustering, stress estimates.  Writes tables/eq_local.parquet.

    python -I svi_science/scripts/g2_local_build.py
"""
import glob
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator

ROOT = "/home/jovyan/tremorferometry"
G2 = f"{ROOT}/svi_science/results/g2_local"
KM = 111.195
VOLC = {"Baker": (48.777, -121.813), "GlacierPeak": (48.112, -121.114), "Rainier": (46.853, -121.760),
        "StHelens": (46.191, -122.194), "Adams": (46.206, -121.490)}


def families():
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    F = F.assign(region="SVI", fid=F.fid, dep=F.dep)[["region", "fid", "lat", "lon", "dep"]]
    W = pd.read_csv(f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv")
    W = W[W.ab & W.copy_of.isna()]
    W = W.assign(region="WA")[["region", "fid", "lat", "lon", "dep"]]
    return pd.concat([F, W], ignore_index=True)


def border_canada(lat, lon):
    """Crude US-Canada border: 49 N east of Boundary Bay, Haro Strait / mid Juan de Fuca Strait west of it."""
    lat = np.asarray(lat); lon = np.asarray(lon)
    jdf = np.interp(lon, [-124.75, -123.25, -123.15, -123.05], [48.50, 48.25, 48.70, 49.0])  # Strait -> Haro Strait
    return np.where(lon >= -123.05, lat >= 49.0, lat >= jdf)


def slab2():
    S = pd.read_csv(f"{ROOT}/data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz", header=None, names=["lon", "lat", "z"])
    S["lon"] -= 360
    lons = np.sort(S.lon.unique()); lats = np.sort(S.lat.unique())
    Z = S.pivot_table(index="lat", columns="lon", values="z", dropna=False).reindex(index=lats, columns=lons).to_numpy()
    return RegularGridInterpolator((lats, lons), -Z, bounds_error=False, fill_value=np.nan)


def mccrory():
    M = pd.read_csv(f"{ROOT}/data/slab_models/mccrory2012/Drafted_slab_XYZ_Grid_GCS_WGS84.txt")
    M = M.dropna()
    return LinearNDInterpolator(M[["Latitude", "Longitude"]].to_numpy(), -M.Depth_km.to_numpy())


def gk_decluster(t_days, mag, lat, lon):
    """Gardner & Knopoff (1974) windows, largest-first; returns mainshock flag and cluster id."""
    n = len(t_days); cl = -np.ones(n, int); order = np.argsort(-mag, kind="stable")
    for i in order:
        if cl[i] >= 0:
            continue
        m = mag[i]
        T = 10 ** (0.032 * m + 2.7389) if m >= 6.5 else 10 ** (0.5409 * m - 0.547)
        L = 10 ** (0.1238 * m + 0.983)
        dt = t_days - t_days[i]
        d = KM * np.hypot(lat - lat[i], (lon - lon[i]) * np.cos(np.radians(lat[i])))
        sel = (dt >= -1e-9) & (dt <= T) & (d <= L) & (cl < 0)        # aftershocks (forward window)
        sel |= (dt < 0) & (dt >= -T) & (d <= L) & (cl < 0)          # foreshocks (same window backward)
        cl[sel] = i; cl[i] = i
    main = cl == np.arange(n)
    return main, cl


def stress(M, r_km):
    """Order-of-magnitude stress proxies at hypocentral distance r (PLAN §1.6). Returns kPa (static, dynamic)."""
    M0 = 10 ** (1.5 * M + 9.1)
    r = np.maximum(r_km, 1.0) * 1e3
    stat = M0 / (4 * np.pi * r ** 3)
    beta, rho, mu, Q, dsig = 3700.0, 2900.0, 3e10, 300.0, 3e6
    fc = 0.37 * beta * (16 * dsig / (7 * M0)) ** (1 / 3)
    wc = 2 * np.pi * fc
    pgv = 0.6 * M0 * wc ** 2 / (4 * np.pi * rho * beta ** 3 * r) * np.exp(-np.pi * fc * r / (Q * beta))
    return stat / 1e3, mu * pgv / beta / 1e3, pgv


def frame(F, lat, lon):
    """Along-strike / down-dip coordinates (km) about the region's family mean; down-dip azimuth from Slab2."""
    sl = slab2(); la0, lo0 = F.lat.mean(), F.lon.mean(); h = 0.1
    dzdx = (sl([la0, lo0 + h])[0] - sl([la0, lo0 - h])[0]) / (2 * h * KM * np.cos(np.radians(la0)))
    dzdy = (sl([la0 + h, lo0])[0] - sl([la0 - h, lo0])[0]) / (2 * h * KM)
    az = np.arctan2(dzdx, dzdy)                    # depth increases toward az (positive-down depths)
    x = (np.asarray(lon) - lo0) * KM * np.cos(np.radians(la0)); y = (np.asarray(lat) - la0) * KM
    dd = x * np.sin(az) + y * np.cos(az); al = x * np.cos(az) - y * np.sin(az)
    return al, dd, np.degrees(az) % 360


def main():
    C = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{G2}/dl/comcat/comcat_*.csv"))], ignore_index=True)
    C = C.drop_duplicates("id")
    C = pd.DataFrame(dict(id=C.id, src="ComCat", net=C.net, type=C.type, t=pd.to_datetime(C.time, utc=True),
                          lat=C.latitude, lon=C.longitude, depth=C.depth, mag=C.mag, magType=C.magType,
                          depth_fixed=C.depthError.isna() | (C.depthError == 0), place=C.place))
    N = pd.read_csv(f"{G2}/dl/nrcan_events.csv").drop_duplicates("id")
    dep = N.depth.astype(float)
    # NRCan's QuakeML depth is in km (checked: 2015-12-30 07:39 UTC M4.7 Sidney 57.5 in both QuakeML and text)
    N = pd.DataFrame(dict(id="nrcan:" + N.id.astype(str), src="NRCan", net="CN", type=N.type,
                          t=pd.to_datetime(N.time, utc=True), lat=N.lat, lon=N.lon, depth=dep,
                          mag=N.mag, magType=N.magType,
                          depth_fixed=N.depthType.fillna("").str.contains("operator"),
                          place=N.place))
    print("raw rows: ComCat", len(C), C.type.value_counts().to_dict())
    print("raw rows: NRCan", len(N), N.type.value_counts().to_dict())
    A = pd.concat([C, N], ignore_index=True)
    A["t"] = A.t.dt.tz_convert(None).astype("datetime64[ns]")
    A = A.sort_values("t").reset_index(drop=True)
    # ---- duplicates (any type), pairing each ComCat row with the nearest NRCan row
    tn = A.t.astype("int64").to_numpy() / 1e9
    isC = (A.src == "ComCat").to_numpy()
    iN = np.where(~isC)[0]; tN = tn[iN]
    pair = -np.ones(len(A), int)
    for i in np.where(isC)[0]:
        a, b = np.searchsorted(tN, tn[i] - 5), np.searchsorted(tN, tn[i] + 5)
        best, bd = -1, 1e9
        for j in iN[a:b]:
            d = KM * np.hypot(A.lat[i] - A.lat[j], (A.lon[i] - A.lon[j]) * np.cos(np.radians(A.lat[i])))
            if d <= 25 and abs(tn[i] - tn[j]) + d / 10 < bd and pair[j] < 0:
                best, bd = j, abs(tn[i] - tn[j]) + d / 10
        if best >= 0:
            pair[i] = best; pair[best] = i
    keep = np.ones(len(A), bool)
    can = border_canada(A.lat, A.lon)
    A["dup_of"] = ""; A["matched"] = pair >= 0
    A["mag_other"] = np.nan; A["type_other"] = ""
    for i in np.where(isC & (pair >= 0))[0]:
        j = pair[i]
        win, lose = (j, i) if can[i] else (i, j)
        keep[lose] = False
        A.loc[win, "dup_of"] = A.id[lose]; A.loc[win, "mag_other"] = A.mag[lose]; A.loc[win, "type_other"] = A.type[lose]
        if A.depth_fixed[win] and not A.depth_fixed[lose]:
            A.loc[win, "depth"] = A.depth[lose]; A.loc[win, "depth_fixed"] = False
    print("matched pairs", int((isC & (pair >= 0)).sum()), "kept", int(keep.sum()), "of", len(A))
    A = A[keep].reset_index(drop=True)
    # ---- type: keep earthquakes; if either catalogue calls it a blast/explosion -> removed
    nonq = (A.type != "earthquake") | A.type_other.isin(["explosion", "quarry blast", "mining explosion",
                                                          "chemical explosion", "experimental explosion",
                                                          "not existing", "other event"])
    print("removed non-earthquake types:", A.loc[nonq, "type"].value_counts().to_dict(),
          "other-catalogue veto", int((nonq & (A.type == "earthquake")).sum()))
    A["removed"] = np.where(nonq, "type", "")
    # ---- volcano circles
    for v, (la, lo) in VOLC.items():
        d = KM * np.hypot(A.lat - la, (A.lon - lo) * np.cos(np.radians(la)))
        A.loc[(d <= 15) & (A.removed == ""), "removed"] = f"volcano_{v}"
    # ---- geometry
    F = families()
    sl = slab2(); mc = mccrory()
    A["slab2_km"] = sl(np.c_[A.lat, A.lon])
    A["mccrory_km"] = mc(np.c_[A.lat, A.lon])
    A["dz_slab2"] = A.depth - A.slab2_km
    A["dz_mccrory"] = A.depth - A.mccrory_km
    for reg in ("SVI", "WA"):
        f = F[F.region == reg]
        dh = KM * np.hypot(A.lat.to_numpy()[:, None] - f.lat.to_numpy()[None],
                           (A.lon.to_numpy()[:, None] - f.lon.to_numpy()[None]) * np.cos(np.radians(f.lat.mean())))
        dv = np.clip(A.depth.to_numpy(), 0, None)[:, None] - f.dep.to_numpy()[None]
        r = np.sqrt(dh ** 2 + dv ** 2)
        A[f"{reg}_dist_km"] = r.min(1); A[f"{reg}_epi_km"] = dh.min(1)
        A[f"{reg}_nearest_fid"] = f.fid.to_numpy()[r.argmin(1)]
        al, dd, az = frame(f, A.lat, A.lon)
        fa, fd, _ = frame(f, f.lat, f.lon)
        A[f"{reg}_along_km"] = al; A[f"{reg}_downdip_km"] = dd
        A[f"{reg}_in_strike"] = (al >= np.quantile(fa, 0.0) - 15) & (al <= np.quantile(fa, 1.0) + 15)
        z5, z95 = np.nanquantile(sl(np.c_[f.lat, f.lon]), [0.05, 0.95])
        A[f"{reg}_zone"] = np.where(A.slab2_km < z5, "updip", np.where(A.slab2_km > z95, "downdip", "lfe_zone"))
        A.loc[A.slab2_km.isna(), f"{reg}_zone"] = "no_slab"
        print(reg, "down-dip az", round(az, 1), "family slab2 depth 5-95 %", round(z5, 1), round(z95, 1))
        s, dy, pgv = stress(A.mag.to_numpy(), A[f"{reg}_dist_km"].to_numpy())
        A[f"{reg}_static_kpa"] = s; A[f"{reg}_dyn_kpa"] = dy
    # ---- class (PLAN §1.4), Slab2 primary; McCrory column
    def cls(dz, depth, slab, fixed):
        c = np.full(len(dz), "unknown_depth", dtype=object)
        ok = ~fixed
        c[ok & np.isnan(slab)] = "offshore_noslab"
        c[ok & (dz > 5)] = "intraslab"
        c[ok & (np.abs(dz) <= 5)] = "near_interface"
        c[ok & (dz < -5) & (depth < 20)] = "crustal"
        c[ok & (dz < -5) & (depth >= 20)] = "overriding_deep"
        return c
    fx = A.depth_fixed.to_numpy().astype(bool)
    A["cls"] = cls(A.dz_slab2.to_numpy(), A.depth.to_numpy(), A.slab2_km.to_numpy(), fx)
    A["cls_mccrory"] = cls(A.dz_mccrory.to_numpy(), A.depth.to_numpy(), A.mccrory_km.to_numpy(), fx)
    A["dist_min_km"] = A[["SVI_dist_km", "WA_dist_km"]].min(1)
    A["t_PT"] = A.t.dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M %Z")
    # ---- declustering on the analysis set (earthquakes, not removed, <= 150 km of either region, M >= 1.0)
    use = (A.removed == "") & (A.dist_min_km <= 150) & A.mag.notna()
    B = A[use].reset_index()
    main, cl = gk_decluster(B.t.astype("int64").to_numpy() / 86400e9, B.mag.to_numpy(), B.lat.to_numpy(), B.lon.to_numpy())
    A["gk_main"] = False; A["gk_cluster"] = -1
    A.loc[B["index"], "gk_main"] = main; A.loc[B["index"], "gk_cluster"] = B["index"].to_numpy()[cl]
    A["in_analysis"] = use
    print("analysis set", int(use.sum()), "GK mainshocks", int(main.sum()))
    print(A[use].cls.value_counts().to_dict())
    print("classes (McCrory)", A[use].cls_mccrory.value_counts().to_dict())
    A.to_parquet(f"{G2}/tables/eq_local.parquet", index=False)


if __name__ == "__main__":
    main()
