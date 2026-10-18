"""Q6 build: PNSN tremor on one margin axis, space-time components (episodes), front runs, starts and stops.

Rules: svi_science/results/q6_front_arrest/PLAN.md (pre-registered 2026-10-07 4:55 pm PT).
Writes results/q6_front_arrest/:
  axis_nodes.csv       node lat, lon, s_km
  epicentres.parquet   time (UTC), time_pt, s_km, d_km, bin, day, comp (-1 = not in an active bin-day)
  episodes_pnsn.csv    one row per component (>= 3 epicentres): n, times, start, run ends, big flags
  bin_episode.parquet  comp, bin, n, t_med (days since 2010-01-01 PT), active (>=3)
  build_info.json
"""
import json
import numpy as np
import pandas as pd
from scipy import ndimage

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q6_front_arrest"
TZ = "America/Los_Angeles"
BIN = 5.0
T0 = pd.Timestamp("2010-01-01")          # PT day origin
STRAIT = (48.30, -123.40)


_TR = None


def xy(lat, lon):
    """Azimuthal equidistant projection (km) centred on the Strait reference point."""
    global _TR
    if _TR is None:
        import pyproj
        _TR = pyproj.Transformer.from_crs("EPSG:4326", f"+proj=aeqd +lat_0={STRAIT[0]} +lon_0={STRAIT[1]} +units=km +ellps=WGS84",
                                          always_xy=True)
    x, y = _TR.transform(np.asarray(lon, float), np.asarray(lat, float))
    return np.asarray(x), np.asarray(y)


def build_axis(P):
    rows = []
    for lo in np.arange(39.5, 51.0, 0.2):
        m = (P.lat.values >= lo) & (P.lat.values < lo + 0.2)
        if m.sum() >= 2000:
            rows.append((np.median(P.lat.values[m]), np.median(P.lon.values[m])))
    N = pd.DataFrame(rows, columns=["lat", "lon"])
    lon_s = N.lon.rolling(3, center=True, min_periods=1).mean()
    N["lon"] = lon_s.values
    x, y = xy(N.lat, N.lon)
    seg = np.hypot(np.diff(x), np.diff(y))
    N["s_raw"] = np.r_[0, np.cumsum(seg)]
    return N


def project(N, lat, lon, chunk=200000):
    nx, ny = xy(N.lat.values, N.lon.values)
    ax, ay = nx[:-1], ny[:-1]; dx, dy = np.diff(nx), np.diff(ny); L2 = dx ** 2 + dy ** 2
    s0 = N.s_raw.values[:-1]
    S = np.empty(len(lat)); D = np.empty(len(lat))
    for i in range(0, len(lat), chunk):
        px, py = xy(lat[i:i + chunk], lon[i:i + chunk])
        u = ((px[:, None] - ax) * dx + (py[:, None] - ay) * dy) / L2
        uc = np.clip(u, 0, 1)
        qx = ax + uc * dx; qy = ay + uc * dy
        dist = np.hypot(px[:, None] - qx, py[:, None] - qy)
        k = dist.argmin(1); r = np.arange(len(k))
        S[i:i + chunk] = s0[k] + uc[r, k] * np.sqrt(L2[k])
        # signed: + = east of the line (down-dip side)
        cross = dx[k] * (py - ay[k]) - dy[k] * (px - ax[k])
        D[i:i + chunk] = -np.sign(cross) * dist[r, k]
    return S, D


def link_cells(act, rd, rb):
    """Connected components of active cells, linking cells with |dday| <= rd and |dbin| <= rb. Labels 1..n, 0 = inactive."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    nd, nb = act.shape
    ids = -np.ones(act.shape, np.int64); ids[act] = np.arange(act.sum())
    I, J = [], []
    for dd in range(0, rd + 1):
        for db in range(-rb, rb + 1):
            if dd == 0 and db <= 0:
                continue
            a = act[:nd - dd, max(0, -db):nb - max(0, db)]
            b = act[dd:, max(0, db):nb + min(0, db)]
            m = a & b
            ia = ids[:nd - dd, max(0, -db):nb - max(0, db)][m]
            ib = ids[dd:, max(0, db):nb + min(0, db)][m]
            I.append(ia); J.append(ib)
    I = np.concatenate(I); J = np.concatenate(J); n = int(act.sum())
    g = coo_matrix((np.ones(len(I)), (I, J)), shape=(n, n))
    nc, lab = connected_components(g, directed=False)
    out = np.zeros(act.shape, np.int64); out[act] = lab + 1
    return out, nc


def two_class_split(v):
    v = np.sort(np.asarray(v, float)); best = (np.inf, None)
    for i in range(1, len(v)):
        a, b = v[:i], v[i:]
        w = ((a - a.mean()) ** 2).sum() + ((b - b.mean()) ** 2).sum()
        if w < best[0]:
            best = (w, (v[i - 1] + v[i]) / 2)
    return best[1]


def front_run(active_bins, start_bin, lo_bin, hi_bin):
    """Walk from start bin NW (+) and SE (-) allowing single-bin gaps. Returns (lo_end, hi_end) bins."""
    A = set(int(b) for b in active_bins)
    if start_bin not in A:
        arr = np.array(sorted(A)); start_bin = int(arr[np.argmin(np.abs(arr - start_bin))])
    ends = {}
    for step in (+1, -1):
        b = start_bin; last = start_bin
        while True:
            if (b + step) in A:
                b += step; last = b
            elif (b + 2 * step) in A:
                b += 2 * step; last = b
            else:
                break
            if b <= lo_bin or b >= hi_bin:
                break
        ends[step] = last
    return ends[-1], ends[+1], start_bin


def main():
    info = {}
    P = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    P["time"] = pd.to_datetime(P.time).astype("datetime64[ns]")
    P["time_pt"] = P.time.dt.tz_localize("UTC").dt.tz_convert(TZ).dt.tz_localize(None)
    N = build_axis(P)
    ss, _ = project(N, np.array([STRAIT[0]]), np.array([STRAIT[1]]))
    N["s_km"] = N.s_raw - ss[0]
    N.to_csv(f"{OUT}/axis_nodes.csv", index=False)
    S, D = project(N, P.lat.values, P.lon.values)
    P["s_km"] = S - ss[0]; P["d_km"] = D
    info["n_raw"] = int(len(P)); info["n_dropped_d60"] = int((np.abs(D) > 60).sum())
    P = P[np.abs(P.d_km) <= 60].reset_index(drop=True)
    smin, smax = float(N.s_km.min()), float(N.s_km.max())
    info["domain_s"] = [smin, smax]
    P = P[(P.s_km >= smin) & (P.s_km <= smax)].reset_index(drop=True)
    edges0 = np.floor(smin / BIN) * BIN
    P["bin"] = ((P.s_km - edges0) // BIN).astype(int)
    P["tday"] = (P.time_pt - T0) / pd.Timedelta(days=1)
    P["day"] = np.floor(P.tday).astype(int)
    nb, nd = P.bin.max() + 1, P.day.max() + 1
    G = np.zeros((nd, nb), np.int32)
    np.add.at(G, (P.day.values, P.bin.values), 1)
    act = G >= 3
    lab, ncomp = link_cells(act, 2, 3)
    info["n_active_bin_days"] = int(act.sum()); info["n_components"] = int(ncomp)
    P["comp"] = lab[P.day.values, P.bin.values] - 1
    P.loc[~act[P.day.values, P.bin.values], "comp"] = -1
    info["frac_epicentres_in_components"] = float((P.comp >= 0).mean())
    P[["time", "time_pt", "s_km", "d_km", "bin", "tday", "day", "comp"]].to_parquet(f"{OUT}/epicentres.parquet")

    C = P[P.comp >= 0].sort_values(["comp", "tday"])
    be = C.groupby(["comp", "bin"]).agg(n=("tday", "size"), t_med=("tday", "median")).reset_index()
    be["active"] = be.n >= 3
    be.to_parquet(f"{OUT}/bin_episode.parquet")
    lo_bin = 0; hi_bin = nb - 1
    rows = []
    actd = {c: g.bin.values for c, g in be[be.active].groupby("comp")}
    for c, g in C.groupby("comp", sort=True):
        n = len(g); t = g.tday.values; s = g.s_km.values
        k10 = max(1, int(np.ceil(0.1 * n)))
        start = float(np.median(s[:k10]))
        start2d = float(np.median(s[t <= t[0] + 2]))
        ab = actd.get(c, np.array([], int)) if n >= 10 else np.array([], int)
        r = dict(comp=int(c), n=n, t_first=float(t[0]), t_last=float(t[-1]), dur=float(t[-1] - t[0]),
                 t_first_PT=str(g.time_pt.iloc[0])[:16], start_s=start, start2d_s=start2d,
                 s_p05=float(np.quantile(s, 0.05)), s_p95=float(np.quantile(s, 0.95)), n_active_bins=len(ab))
        if len(ab):
            sb = int((start - edges0) // BIN)
            lo_e, hi_e, sb2 = front_run(ab, sb, lo_bin, hi_bin)
            r.update(start_bin=sb2, run_lo_bin=lo_e, run_hi_bin=hi_e,
                     stop_se_s=edges0 + (lo_e + 0.5) * BIN, stop_nw_s=edges0 + (hi_e + 0.5) * BIN)
        rows.append(r)
    E = pd.DataFrame(rows)
    big_thr = two_class_split(np.log10(E.n[E.n >= 30]))
    info["big_split_log10"] = float(big_thr); info["big_split_n"] = float(10 ** big_thr)
    E["big"] = (E.n >= 30) & (np.log10(E.n) > big_thr)
    E["big300"] = E.n >= 300; E["big1000"] = E.n >= 1000
    E["start_date_PT"] = [str((T0 + pd.Timedelta(days=x)).date()) for x in E.t_first]
    E.to_csv(f"{OUT}/episodes_pnsn.csv", index=False)
    info["edges0"] = float(edges0); info["n_bins"] = int(nb)
    info["n_comp_ge10"] = int((E.n >= 10).sum()); info["n_comp_ge30"] = int((E.n >= 30).sum())
    info["n_big"] = int(E.big.sum()); info["n_big300"] = int(E.big300.sum()); info["n_big1000"] = int(E.big1000.sum())
    # where named places fall on the axis
    places = {"Victoria": (48.43, -123.37), "Port Angeles": (48.12, -123.43), "Sequim": (48.08, -123.10),
              "Lake Cowichan": (48.82, -124.05), "Hood Canal S (Union)": (47.36, -123.10), "Olympia": (47.04, -122.90),
              "Columbia R. (46.2N)": (46.2, -123.1), "Port Alberni": (49.23, -124.81), "Mendocino (40.4N)": (40.4, -122.6)}
    ps, _ = project(N, np.array([v[0] for v in places.values()]), np.array([v[1] for v in places.values()]))
    info["places_s_km"] = {k: float(v - ss[0]) for k, v in zip(places, ps)}
    json.dump(info, open(f"{OUT}/build_info.json", "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
