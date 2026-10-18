"""s11_dist_common -- shared helpers for the S11 distant-station lever (svi_science/results/s11_depth/distant/PLAN.md).
Region = S11_REGION (SVI or WA_NW, as s11_common). Outputs under s11_depth/distant/<svi|wa_nw>/. Read-only on every
earlier S11 / S10 / WA table. eps = -dv/v."""
from __future__ import annotations
import os, sys, io, urllib.request
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/scripts")
import s10_common as S
import s11_common as C

REG = C.REGION
SUB = {"SVI": "svi", "WA_NW": "wa_nw"}[REG]
S11 = f"{C.ROOT}/svi_science/results/s11_depth"
DROOT = f"{S11}/distant"
DOUT = f"{DROOT}/{SUB}"
DTB, DLOG, DFIG = f"{DOUT}/tables", f"{DOUT}/logs", f"{DROOT}/figures"
NEAR_TB = f"{S11}/tables" if REG == "SVI" else f"{S11}/wa_nw/tables"      # near-run tables (cells, store, kernels)
TEMP_NETS = {"XG", "XU", "ZK", "X4", "ZZ", "YW", "YG", "XD", "7D", "Z2", "X9"}
HMIN, HMAX = 80.0, 150.0
MIN_DATES = 150
NEAR_GRID, INV_GRID, FAR_GRID = (-14.0, 10.0), (-30.0, 12.0), (-28.0, 10.0)
VEL = f"{C.ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/{'south_VI' if REG == 'SVI' else 'royer_NW_WA'}.csv"
SITE_KM = 3.0


def ensure():
    for d in (DTB, DLOG, DFIG):
        os.makedirs(d, exist_ok=True)


def set_grid(g):
    """switch the s10_common window grid (module globals read at call time by s10_windows / s10_measure / s11_measure)"""
    t0, t1 = g
    S.T0W, S.T1W = float(t0), float(t1)
    S.NW = int(round((t1 - t0) * S.FS)) + 1
    S.TAU = t0 + np.arange(S.NW) / S.FS
    S.CUT_PRE = -t0 + 1.0
    S.CUT_POST = 30.0


def km(la1, lo1, la2, lo2):
    return np.hypot((np.asarray(lo1) - lo2) * 111.195 * np.cos(np.radians((np.asarray(la1) + la2) / 2)),
                    (np.asarray(la1) - la2) * 111.195)


# ------------------------------------------------------------------------------ raw roots and reader
def roots():
    import wa_prod_common as WP
    r = list(WP.prod_roots())
    for x in (S.RAW_BW, S.RAW_FW):
        if x not in r:
            r.append(x)
    return r


def raw_path(ns, date, rts):
    t = pd.Timestamp(date)
    for base in rts:
        p = f"{base}/{ns}/{t.year}/{t.dayofyear:03d}.mseed"
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    return None


# ------------------------------------------------------------------------------ geometry
def station_coords():
    W = pd.read_csv(f"{C.ROOT}/lfe_extension_wa/results/wa_prod/setup/stations.csv")[["netsta", "lat", "lon"]]
    B = pd.read_csv(f"{C.ROOT}/svi_science/results/backward_scope/inventory_stations.csv")
    B["netsta"] = B.net + "." + B.sta
    out = pd.concat([W, B[["netsta", "lat", "lon"]]], ignore_index=True).drop_duplicates("netsta")
    return out.set_index("netsta")


def sensor_depths(stas, log=print):
    T = pd.concat([pd.read_csv(f"{S11}/tables/sensor_depths.csv"), pd.read_csv(f"{S11}/wa_nw/tables/sensor_depths.csv")])
    T["netsta"] = T.Network + "." + T.Station
    d = T.drop_duplicates("netsta").set_index("netsta").sensor_depth_m.to_dict()
    p = f"{DTB}/sensor_depths_far.csv"
    old = pd.read_csv(p).set_index("netsta").sensor_depth_m.to_dict() if os.path.exists(p) else {}
    rows = []
    for ns in stas:
        if ns in d and np.isfinite(d[ns]):
            rows.append(dict(netsta=ns, sensor_depth_m=float(d[ns]), src="s11 tables")); continue
        if ns in old:
            rows.append(dict(netsta=ns, sensor_depth_m=float(old[ns]), src="cached FDSN")); continue
        net, sta = ns.split(".")
        url = f"https://service.earthscope.org/fdsnws/station/1/query?net={net}&sta={sta}&level=channel&format=text"
        try:
            txt = urllib.request.urlopen(url, timeout=60).read().decode()
            D = pd.read_csv(io.StringIO(txt), sep="|"); D.columns = [c.strip().lstrip("#").strip() for c in D.columns]
            D = D[D.Channel.str[:2].isin(["HH", "EH", "BH"])]
            rows.append(dict(netsta=ns, sensor_depth_m=float(D.Depth.max()), src=f"EarthScope FDSN {C.pt_now()}"))
        except Exception as e:
            log(f"FDSN {ns}: {e} -> 0 m (flagged)")
            rows.append(dict(netsta=ns, sensor_depth_m=0.0, src=f"FDSN failed -> 0 (flagged): {e}"))
    T = pd.DataFrame(rows); T.to_csv(p, index=False)
    return T.set_index("netsta").sensor_depth_m.to_dict()


def epochs_table(ns):
    """[(epoch, start, end)] : s10_common.EPOCHS if present, else the WA node / epoch-change logic (s11_wa_inventory)"""
    if ns in S.EPOCHS:
        return [(e, a, b) for e, a, b in S.EPOCHS[ns]]
    WA = f"{C.ROOT}/lfe_extension_wa"
    roots_ = [f"{C.ROOT}/data/extension_stations/{r}_bridge_nodes_strict.csv" for r in ("NW_WA", "SW_WA", "OLY2010")]
    N = pd.concat([pd.read_csv(p) for p in roots_] + [pd.read_csv(f"{WA}/results/step4b/setup/extra_nodes.csv")])
    N["netsta"] = N.net + "." + N.sta
    E = pd.read_csv(f"{WA}/results/wa_pre/w2/epoch_changes.csv")
    b = sorted(set(N[N.netsta == ns].start.astype(str)) | set(E[E.netsta == ns].date.astype(str)))
    if not b:
        return [("E1", "1990-01-01", "2100-01-01")]
    out = []
    for i, d in enumerate(b):
        out.append((f"E{i + 1}", d if i else "1990-01-01", b[i + 1] if i + 1 < len(b) else "2100-01-01"))
    return out


def epoch_of(ep, dates):
    d = pd.Series(np.asarray(dates, dtype=str)); out = np.full(len(d), "", dtype=object)
    for e, a, b in ep:
        out[((d >= a) & (d < b)).to_numpy()] = e
    return out


def dead_days():
    bad = pd.read_csv(f"{C.ROOT}/lfe_extension_wa/results/wa_pre/w2/r2_bad_for_dvv.csv")
    bad = bad[bad.reason == "dead"]
    s = set()
    for r in bad.itertuples():
        for dd in pd.date_range(r.start, r.end, freq="D"):
            s.add((r.netsta, dd.strftime("%Y-%m-%d")))
    return s


# ------------------------------------------------------------------------------ region near tables
def near_tables():
    """-> dict(F families [fam, lat, lon, depth_km], X rows, O near offsets [fam, netsta, offset], A active days,
    near_pairs set, near_stations list)"""
    if REG == "SVI":
        cat = pd.read_csv(f"{C.ROOT}/data/catalogue/families_v1_3.csv")
        cat = cat[cat.region == "S_Vancouver_Is"].copy(); cat["fam"] = pd.to_numeric(cat.family_id, errors="coerce")
        f0 = set(pd.read_csv(f"{C.MTB}/families.csv").fam); f1 = set(pd.read_csv(f"{S11}/tables/ext/families.csv").fam)
        F = cat[cat.fam.isin(f0 | f1)][["fam", "lat", "lon", "depth_km"]].copy(); F["fam"] = F.fam.astype(int)
        X = pd.concat([pd.read_parquet(f"{C.MTB}/rows.parquet"), pd.read_parquet(f"{S11}/tables/ext/rows.parquet")], ignore_index=True)
        O = pd.concat([pd.read_csv(f"{C.MTB}/offsets_sp.csv"), pd.read_csv(f"{S11}/tables/ext/offsets_sp.csv")], ignore_index=True)
        A = pd.concat([pd.read_csv(f"{C.MTB}/active_days.csv"), pd.read_csv(f"{S11}/tables/ext/active_days.csv")], ignore_index=True)
        near = list(C.STATIONS)
    else:
        F = pd.read_csv(f"{C.MTB}/families.csv")[["fam", "lat", "lon", "depth_km"]]
        X = pd.read_parquet(f"{C.MTB}/rows.parquet")
        O = pd.read_csv(f"{C.MTB}/offsets_sp.csv")
        A = pd.read_csv(f"{C.MTB}/active_days.csv")
        near = sorted(O.netsta.unique())
    O = O[np.isfinite(O.offset)][["fam", "netsta", "offset"]]
    X = X[X.fam.isin(set(F.fam))]
    near_pairs = set(zip(O.fam.astype(int), O.netsta)) if REG != "SVI" else {(int(f), ns) for f in F.fam for ns in near}
    return dict(F=F.reset_index(drop=True), X=X.reset_index(drop=True), O=O.reset_index(drop=True), A=A, near_pairs=near_pairs,
                near_stations=near)


# ------------------------------------------------------------------------------ 1-D flat-layer travel times
class TT:
    def __init__(self, csv=VEL, zmax=200.0):
        P = pd.read_csv(csv)
        self.top = P.ztop_km.to_numpy(float)
        self.bot = np.where(np.isfinite(P.zbot_km.to_numpy(float)), P.zbot_km.to_numpy(float), zmax)
        self.v = {"S": P.vs.to_numpy(float), "P": P.vp.to_numpy(float)}

    def _above(self, zs, ph):
        h, v = [], []
        for t, b, vv in zip(self.top, self.bot, self.v[ph]):
            if t >= zs:
                break
            h.append(min(b, zs) - t); v.append(vv)
        return np.array(h), np.array(v)

    def direct(self, x, zs, ph):
        h, v = self._above(zs, ph)
        if x < 1e-6:
            return float(np.sum(h / v))
        vm = v.max(); lo, hi = 0.0, 1.0 / vm * (1 - 1e-12)
        for _ in range(80):
            p = 0.5 * (lo + hi); q = np.sqrt(1 - (p * v) ** 2)
            X = np.sum(h * p * v / q)
            lo, hi = (p, hi) if X < x else (lo, p)
        p = 0.5 * (lo + hi); q = np.sqrt(1 - (p * v) ** 2)
        return float(np.sum(h / (v * q)))

    def heads(self, x, zs, ph):
        best = np.inf
        v_all = self.v[ph]
        for k in range(len(self.top)):
            zk = self.top[k]
            if zk <= zs:
                continue
            vk = v_all[k]
            if vk <= v_all[:k].max():
                continue
            p = 1.0 / vk; T = x * p; xc = 0.0
            for t, b, vv in zip(self.top[:k], self.bot[:k], v_all[:k]):
                hh_up = max(0.0, min(b, zs) - t)                       # above the source: one leg
                hh_dn = max(0.0, min(b, zk) - max(t, zs))              # between source and refractor: two legs
                eta = np.sqrt(1 / vv ** 2 - p ** 2); tn = p * vv / np.sqrt(1 - (p * vv) ** 2)
                T += (hh_up + 2 * hh_dn) * eta; xc += (hh_up + 2 * hh_dn) * tn
            if x >= xc:
                best = min(best, T)
        return best

    def first(self, x, zs, ph):
        return min(self.direct(x, zs, ph), self.heads(x, zs, ph))
