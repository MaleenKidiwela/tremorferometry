#!/usr/bin/env python
"""P4 quiet flag, attempt 2, step 1: envelope cross-correlation (Wech & Creager 2008 style) window table per UTC day.
Rules: svi_science/results/backward/p4/envflag/ENVFLAG.md (written before computing).  CPU only.

  python bw_p4env_detect.py --years 2010,2011 [--workers 8]  -> p4/envflag/win/env_<Y>.parquet
  python bw_p4env_detect.py --dates 2010-08-20 --dry          (one day, prints)

Per day one job: envelopes of every roster station with a cache, then the grid search for the main set (CN4) and the
variant (all, one per site), each on the data and on the scrambled control.
"""
import os, sys, argparse, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

ROOT = "/home/jovyan/tremorferometry"
LFE = f"{ROOT}/lfe_extension"
CACHES = [f"{LFE}/cache20hz_backward", f"{LFE}/cache20hz_r1", f"{LFE}/cache20hz",
          f"{ROOT}/lfe_catalog_repair/cache20hz", f"{ROOT}/lfe_picker/cache20hz"]
OUT = f"{ROOT}/svi_science/results/backward/p4/envflag/win"
VEL = f"{ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/south_VI.csv"
FAMS = f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv"
STA_TAB = f"{ROOT}/svi_science/results/backward_scope/inventory_stations.csv"
CN4 = ["PGC", "VGZ", "NLLB", "SNB"]
ROSTER = {"PGC": "CN", "VGZ": "CN", "NLLB": "CN", "SNB": "CN", "LZB": "CN",
          "TWKB": "PO", "TSJB": "PO", "KLNB": "PO", "SILB": "PO", "SSIB": "PO", "TWGB": "PO",
          "SHVB": "C8", "MGCB": "C8", "TWBB": "C8", "SHDB": "C8", "KHVB": "C8",
          "B001": "PB", "B003": "PB", "B004": "PB", "B005": "PB", "B006": "PB", "B007": "PB",
          "B009": "PB", "B010": "PB", "B011": "PB", "B926": "PB", "LRIV": "UW"}
SITES = [["PGC", "B011", "B010", "B009"], ["LRIV", "B005", "B006", "B007"]]
BAND = (2.0, 8.0)
WIN, STEP = 300, 150
NWIN = (86400 - WIN) // STEP + 1                     # 575
DEPTH = 35.0
AREA_KM = 30.0
GAP_PAD = 10


def epi_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def coords():
    S = pd.read_csv(STA_TAB).drop_duplicates("sta").set_index("sta")
    return {s: (float(S.loc[s, "lat"]), float(S.loc[s, "lon"])) for s in ROSTER}


def grid():
    glat = np.round(np.arange(47.50, 49.55 + 1e-9, 0.05), 3)
    glon = np.round(np.arange(-125.10, -122.30 + 1e-9, 0.05), 3)
    G = np.array([(a, b) for a in glat for b in glon])
    F = pd.read_csv(FAMS)
    dmin = epi_km(G[:, 0:1], G[:, 1:2], F.lat.to_numpy()[None], F.lon.to_numpy()[None]).min(1)
    return G, dmin


def tt_table():
    sys.path.insert(0, f"{ROOT}/lfe_picker/scripts")
    import p_vel1d as V
    M = pd.read_csv(VEL)
    V.ZTOP = M.ztop_km.to_numpy(float)
    V.VP = M.vs.to_numpy(float)                       # S speeds in the "P" slot -> S times with vpvs = 1
    xs = np.arange(0, 260.5, 0.5)
    ts = np.array([V.tt(float(x), DEPTH, "P", 1.0) for x in xs])
    return xs, ts


_G = _D = _T = _C = None


def setup():
    global _G, _D, _T, _C
    if _G is None:
        _C = coords(); _G, _D = grid(); xs, ts = tt_table()
        _T = {s: np.interp(epi_km(_G[:, 0], _G[:, 1], *_C[s]), xs, ts) for s in ROSTER}


def cache_path(sta, date):
    d = pd.Timestamp(date)
    for c in CACHES:
        p = f"{c}/{ROSTER[sta]}.{sta}/{d.year}/{d.dayofyear:03d}.npz"
        if os.path.exists(p):
            return p
    return None


def envelope(sta, date):
    """1 sps 2-8 Hz 3-C envelope of one UTC day (86400 values, NaN where no data)."""
    from scipy.signal import butter, sosfiltfilt, hilbert
    from scipy.fft import next_fast_len
    p = cache_path(sta, date)
    if p is None:
        return None
    try:
        z = np.load(p); x = z["x"].astype(np.float64); t0 = float(z["t0"]); fs = float(z["fs"]) if "fs" in z.files else 20.0
    except Exception:
        return None
    day0 = pd.Timestamp(date).value / 1e9
    n = int(round(86400 * fs))
    X = np.zeros((x.shape[0], n)); valid = np.zeros(n, bool)
    i0 = int(round((t0 - day0) * fs)); a, b = max(i0, 0), min(i0 + x.shape[1], n)
    if b <= a:
        return None
    X[:, a:b] = x[:, a - i0:b - i0]
    X[~np.isfinite(X)] = 0.0
    valid[a:b] = np.any(X[:, a:b] != 0, axis=0)
    if valid.sum() < 3600 * fs:
        return None
    sos = butter(4, [BAND[0] / (fs / 2), min(BAND[1], 0.95 * fs / 2) / (fs / 2)], "band", output="sos")
    nf = next_fast_len(n)
    e2 = np.zeros(n)
    for c in range(X.shape[0]):
        y = sosfiltfilt(sos, X[c])
        e2 += np.abs(hilbert(y, N=nf)[:n]) ** 2
    env = np.sqrt(e2)
    k = int(round(fs))
    eb = env[: 86400 * k].reshape(86400, k).mean(1)
    vb = valid[: 86400 * k].reshape(86400, k).all(1)
    bad = ~vb
    if bad.any():                                         # pad gaps by +-GAP_PAD s
        bad = np.convolve(bad.astype(int), np.ones(2 * GAP_PAD + 1, int), mode="same") > 0
    lp = butter(2, 0.1 / 0.5, "low", output="sos")
    f = np.where(bad, np.median(eb[~bad]) if (~bad).any() else 0.0, eb)
    out = sosfiltfilt(lp, f)
    out[bad] = np.nan
    return out.astype(np.float64)


def windows(e):
    idx = np.arange(NWIN)[:, None] * STEP + np.arange(WIN)[None]
    return e[idx]                                          # (NWIN, WIN)


def search(E, stas):
    """E: dict sta -> 86400 envelope. Returns per-window arrays for the station list `stas` (those present)."""
    stas = [s for s in stas if s in E]
    out = dict(n_valid=np.zeros(NWIN, int), coh=np.full(NWIN, np.nan), s3=np.full(NWIN, np.nan),
               node=np.full(NWIN, -1), imp=np.full(NWIN, np.nan))
    if len(stas) < 3:
        return out
    W = np.stack([windows(E[s]) for s in stas])            # (S, NWIN, WIN)
    ok = np.isfinite(W).all(2)                             # (S, NWIN)
    Wn = np.nan_to_num(W - np.nanmean(W, 2, keepdims=True))
    Wn = Wn / np.maximum(np.linalg.norm(Wn, axis=2, keepdims=True), 1e-12)
    F = np.fft.rfft(Wn, n=2 * WIN, axis=2)
    nS, nG = len(stas), len(_G)
    acc = np.zeros((NWIN, nG)); cnt = np.zeros((NWIN, 1))
    pairs = [(i, j) for i in range(nS) for j in range(i + 1, nS)]
    lagi = {}
    for i, j in pairs:
        cc = np.fft.irfft(F[i] * F[j].conj(), n=2 * WIN, axis=1)       # cc[l] = sum x_i[t+l] x_j[t]
        lag = np.round(_T[stas[i]] - _T[stas[j]]).astype(int)           # i later than j by lag -> peak at +lag
        lag = np.clip(lag, -(WIN - 1), WIN - 1)
        v = cc[:, lag % (2 * WIN)]
        okp = (ok[i] & ok[j]).astype(float)[:, None]
        acc += v * okp; cnt += okp
        lagi[(i, j)] = (cc, lag)
    nval = ok.sum(0)
    use = nval >= 3
    coh = acc / np.maximum(cnt, 1)
    node = coh.argmax(1)
    out["n_valid"] = nval
    out["coh"] = np.where(use, coh[np.arange(NWIN), node], np.nan)
    out["node"] = np.where(use, node, -1)
    # per-station mean cc at the located node -> s3
    sm = np.zeros((nS, NWIN)); sc = np.zeros((nS, NWIN))
    for (i, j), (cc, lag) in lagi.items():
        v = cc[np.arange(NWIN), lag[node] % (2 * WIN)]
        okp = (ok[i] & ok[j]).astype(float)
        sm[i] += v * okp; sm[j] += v * okp; sc[i] += okp; sc[j] += okp
    smean = np.where(ok & (sc > 0), sm / np.maximum(sc, 1), -np.inf)
    s3 = -np.sort(-smean, axis=0)[2] if nS >= 3 else np.full(NWIN, -np.inf)
    out["s3"] = np.where(use & np.isfinite(s3), s3, np.nan)
    net = np.nanmean(np.where(np.isfinite(W), W, np.nan), 0)            # (NWIN, WIN)
    with np.errstate(all="ignore"):
        out["imp"] = np.nanmax(net, 1) / np.nanmedian(net, 1)
    return out


def variant_set(present):
    s = [x for x in ROSTER if x in present]
    for site in SITES:
        keep = next((x for x in site if x in present), None)
        s = [x for x in s if x not in site or x == keep]
    return s


def shift_control(E):
    out = {}
    for i, s in enumerate(ROSTER):                         # fixed roster order -> a station's shift never depends on the set
        if s in E:
            sh = (i * (4 * 3600 + 13 * 60) + 17 * 60) % 86400
            out[s] = np.roll(E[s], sh)
    return out


def one_day(date):
    warnings.filterwarnings("ignore")
    setup()
    E = {}
    for s in ROSTER:
        e = envelope(s, date)
        if e is not None and np.isfinite(e).sum() > 3600:
            E[s] = e
    rows = []
    vset = variant_set(set(E))
    for series, EE in (("data", E), ("c1", shift_control(E))):
        for sset, stas in (("CN4", CN4), ("all", vset)):
            r = search(EE, stas)
            nd = r["node"]
            lat = np.where(nd >= 0, _G[np.maximum(nd, 0), 0], np.nan)
            lon = np.where(nd >= 0, _G[np.maximum(nd, 0), 1], np.nan)
            dfam = np.where(nd >= 0, _D[np.maximum(nd, 0)], np.nan)
            rows.append(pd.DataFrame(dict(date=date, series=series, sset=sset, win=np.arange(NWIN, dtype=np.int16),
                                          n_valid=r["n_valid"].astype(np.int8), coh=r["coh"].astype(np.float32),
                                          s3=r["s3"].astype(np.float32), lat=lat.astype(np.float32),
                                          lon=lon.astype(np.float32), dfam_km=dfam.astype(np.float32),
                                          imp=r["imp"].astype(np.float32),
                                          stas=",".join(s for s in stas if s in EE))))
    D = pd.concat(rows, ignore_index=True)
    return D, dict(date=date, stations=",".join(sorted(E)), n_sta=len(E))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", default=""); ap.add_argument("--dates", default="")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.dates:
        for d in a.dates.split(","):
            D, m = one_day(d)
            print(m)
            x = D[(D.series == "data")]
            print(x.groupby("sset")[["coh", "s3", "n_valid"]].describe().T.to_string())
            if not a.dry:
                D.to_parquet(f"{OUT}/env_day_{d}.parquet", index=False)
        return
    for y in [int(v) for v in a.years.split(",") if v]:
        f = f"{OUT}/env_{y}.parquet"
        if os.path.exists(f) and not a.force:
            print(f"[env] {y} exists, skip", flush=True); continue
        days = [str(d.date()) for d in pd.date_range(f"{y}-01-01", f"{y}-12-31")]
        parts, meta = [], []
        with ProcessPoolExecutor(a.workers) as ex:
            for k, (D, m) in enumerate(ex.map(one_day, days, chunksize=2), 1):
                parts.append(D); meta.append(m)
                if k % 50 == 0:
                    print(f"[env] {y}: {k}/{len(days)} days", flush=True)
        A = pd.concat(parts, ignore_index=True)
        A.to_parquet(f + ".tmp", index=False); os.replace(f + ".tmp", f)
        pd.DataFrame(meta).to_csv(f"{OUT}/env_{y}_stations.csv", index=False)
        print(f"[env] {y}: {len(A)} rows written", flush=True)


if __name__ == "__main__":
    main()
