"""S6 step 1: per detection x station 3-C 2-8 Hz mean-square envelopes (0.1 s bins, -15..+30 s about the expected S)
from the SVI cut-window store, in ground velocity (counts / channel-epoch sensitivity). Rules: results/s6/PLAN_S6.md.
Reads the store lazily (index rows filtered by det_id, one h5 slice per window). Nothing is written outside results/s6.

  python s6_extract.py select          # build results/s6/targets.parquet
  python s6_extract.py run --workers 6 # envelopes per month -> results/s6/env/env_<YYYY-MM>.npz
"""
import os, sys, glob, json, argparse, time
import numpy as np, pandas as pd

ROOT = "/home/jovyan/tremorferometry"
CUT = f"{ROOT}/lfe_extension/results/cut/svi_ext"
OUT = f"{ROOT}/svi_science/results/s6"
ENV = f"{OUT}/env"
NB = 450; BIN = 0.1; PRE = 15.0
STAS = ["CN.NLLB", "CN.PFB", "CN.PGC", "CN.SNB", "CN.VGZ", "PB.B001", "PB.B003", "PB.B004", "PB.B005", "PB.B006", "PB.B007",
        "PB.B009", "PB.B010", "PB.B011", "PB.B926", "UW.LRIV"]
SIDX = {s: i for i, s in enumerate(STAS)}


def select():
    t = pd.read_parquet(f"{ROOT}/svi_science/results/s0/analysis_table.parquet",
                        columns=["fam", "series", "P_same", "keep_xfam8", "time_utc", "net_mad"])
    d = t[t.series == "data"].copy()
    hi = d[d.P_same >= 0.5].assign(kind="event")
    lo = d[d.P_same < 0.05].sample(6000, random_state=0).assign(kind="nonevent")
    x = pd.concat([hi, lo]).sort_values("time_utc")
    D = pd.read_parquet(f"{CUT}/detections.parquet", columns=["det_id", "fam", "t", "t_epoch"])
    D["fam"] = D.fam.astype(int); D = D.sort_values("t")
    m = pd.merge_asof(x, D, left_on="time_utc", right_on="t", by="fam", tolerance=pd.Timedelta("0.1s"), direction="nearest")
    print("matched", m.det_id.notna().mean(), len(m))
    m = m.dropna(subset=["det_id"]); m["det_id"] = m.det_id.astype(np.int64)
    m["month"] = m.time_utc.dt.strftime("%Y-%m")
    m[["det_id", "fam", "time_utc", "t_epoch", "P_same", "keep_xfam8", "net_mad", "kind", "month"]].to_parquet(f"{OUT}/targets.parquet")
    print(m.kind.value_counts())


def gain_table():
    E = pd.read_csv(f"{CUT}/channel_epochs.csv", dtype={"loc": str})
    E["loc"] = E["loc"].fillna("")
    return {f"{r.netsta}.{r.loc}.{r.cha}@{r.start}": float(r.gain) for r in E.itertuples()}


def do_month(m):
    import h5py, pyarrow.parquet as pq
    from scipy.signal import butter, sosfiltfilt
    fo = f"{ENV}/env_{m}.npz"
    if os.path.exists(fo): return m, "done"
    T = pd.read_parquet(f"{OUT}/targets.parquet", filters=[("month", "==", m)])
    if not len(T): return m, "empty"
    fi = f"{CUT}/index_{m}.parquet"
    if not os.path.exists(fi): return m, "noindex"
    cols = ["det_id", "netsta", "loc", "cha", "fs", "t0", "offset_s", "n_samples", "n_gap", "dataset", "h5_start", "h5_file",
            "resp_epoch", "offset_source"]
    I = pq.read_table(fi, columns=cols, filters=[("det_id", "in", T.det_id.tolist())]).to_pandas()
    I = I[I.netsta.isin(SIDX)]
    G = gain_table()
    te = dict(zip(T.det_id, T.t_epoch))
    sos = {}
    keys = I.groupby(["det_id", "netsta"]).size().index.tolist()
    kpos = {k: i for i, k in enumerate(keys)}
    env = np.full((len(keys), NB), np.nan, np.float32)
    ncomp = np.zeros(len(keys), np.int8)
    gmiss = np.zeros(len(keys), bool)
    acc = np.zeros((len(keys), NB)); cnt = np.zeros((len(keys), NB))
    offs = {}; osrc = {}; repochs = {}
    I = I.sort_values(["h5_file", "dataset", "h5_start"])
    for h5f, g in I.groupby("h5_file"):
        with h5py.File(f"{CUT}/{h5f}", "r") as h:
            for r in g.itertuples():
                k = kpos[(r.det_id, r.netsta)]
                x = h[r.dataset][int(r.h5_start):int(r.h5_start) + int(r.n_samples)].astype(np.float64)
                if x.size < 50: continue
                fs = float(r.fs)
                if fs not in sos:
                    sos[fs] = butter(4, [2.0 / (fs / 2), 8.0 / (fs / 2)], btype="band", output="sos")
                good = x != 0 if r.n_gap > 0 else np.ones(x.size, bool)
                if good.sum() < 50: continue
                x = np.where(good, x - x[good].mean(), 0.0)
                y = sosfiltfilt(sos[fs], x)
                gain = G.get(r.resp_epoch)
                if gain is None or not np.isfinite(gain):
                    gmiss[k] = True; gain = np.nan
                y = y / gain
                if r.n_gap > 0:   # mask gaps dilated by 0.5 s
                    bad = ~good; w = int(0.5 * fs)
                    bad = np.convolve(bad.astype(float), np.ones(2 * w + 1), "same") > 0
                    y = np.where(bad, np.nan, y)
                ws = te[r.det_id] + float(r.offset_s) - PRE
                tt = float(r.t0) + np.arange(x.size) / fs - ws
                b = np.floor(tt / BIN).astype(int)
                ok = (b >= 0) & (b < NB) & np.isfinite(y)
                acc[k] += np.bincount(b[ok], weights=y[ok] ** 2, minlength=NB)[:NB]
                cnt[k] += np.bincount(b[ok], minlength=NB)[:NB]
                ncomp[k] += 1
                offs[k] = float(r.offset_s); osrc[k] = r.offset_source; repochs[k] = r.resp_epoch
    with np.errstate(invalid="ignore", divide="ignore"):
        # mean square per bin summed over components = (sum y^2 / samples per comp) ; cnt counts samples over all comps
        nc = np.maximum(ncomp, 1)[:, None].astype(float)
        env = (acc / (cnt / nc)).astype(np.float32)      # = sum over comps of per-comp mean square (if all comps present)
        env[cnt == 0] = np.nan
    det = np.array([k[0] for k in keys], np.int64); sta = np.array([SIDX[k[1]] for k in keys], np.int8)
    off = np.array([offs.get(i, np.nan) for i in range(len(keys))])
    src = np.array([osrc.get(i, "") for i in range(len(keys))])
    rep = np.array([repochs.get(i, "") for i in range(len(keys))])
    np.savez_compressed(fo + ".tmp.npz", det_id=det, sta=sta, env=env, ncomp=ncomp, gain_missing=gmiss, offset_s=off,
                        offset_source=src, resp_epoch=rep)
    os.replace(fo + ".tmp.npz", fo)
    return m, f"{len(keys)} det-stations, {len(T)} targets"


WT0, WN = -14.5, 380      # waveform window: -14.5..+4.5 s about expected S at 20 Hz
CORD = {"Z": 0, "N": 1, "1": 1, "E": 2, "2": 2}


def do_month_wave(m):
    """3-C 2-8 Hz velocity waveforms at 20 Hz, -14.5..+4.5 s about the expected S -> wave_<m>.npz"""
    import h5py, pyarrow.parquet as pq
    from scipy.signal import butter, sosfiltfilt
    from scipy.interpolate import CubicSpline
    fo = f"{ENV}/wave_{m}.npz"
    if os.path.exists(fo): return m, "done"
    T = pd.read_parquet(f"{OUT}/targets.parquet", filters=[("month", "==", m)])
    fi = f"{CUT}/index_{m}.parquet"
    if not len(T) or not os.path.exists(fi): return m, "empty"
    cols = ["det_id", "netsta", "cha", "fs", "t0", "offset_s", "n_samples", "n_gap", "dataset", "h5_start", "h5_file", "resp_epoch"]
    I = pq.read_table(fi, columns=cols, filters=[("det_id", "in", T.det_id.tolist())]).to_pandas()
    I = I[I.netsta.isin(SIDX)]
    G = gain_table(); te = dict(zip(T.det_id, T.t_epoch)); sos = {}
    keys = I.groupby(["det_id", "netsta"]).size().index.tolist(); kpos = {k: i for i, k in enumerate(keys)}
    W = np.full((len(keys), 3, WN), np.nan, np.float32)
    I = I.sort_values(["h5_file", "dataset", "h5_start"])
    for h5f, g in I.groupby("h5_file"):
        with h5py.File(f"{CUT}/{h5f}", "r") as h:
            for r in g.itertuples():
                c = CORD.get(r.cha[-1])
                if c is None: continue
                k = kpos[(r.det_id, r.netsta)]
                x = h[r.dataset][int(r.h5_start):int(r.h5_start) + int(r.n_samples)].astype(np.float64)
                if x.size < 50: continue
                fs = float(r.fs)
                if fs not in sos: sos[fs] = butter(4, [2.0 / (fs / 2), 8.0 / (fs / 2)], btype="band", output="sos")
                good = x != 0 if r.n_gap > 0 else np.ones(x.size, bool)
                if good.sum() < 50: continue
                x = np.where(good, x - x[good].mean(), 0.0)
                y = sosfiltfilt(sos[fs], x) / G.get(r.resp_epoch, np.nan)
                tt = float(r.t0) + np.arange(x.size) / fs - (te[r.det_id] + float(r.offset_s))   # time about expected S
                tg = WT0 + np.arange(WN) / 20.0
                if tg[0] < tt[0] or tg[-1] > tt[-1]: continue
                if r.n_gap > 0:
                    i0, i1 = np.searchsorted(tt, tg[0]) - 2, np.searchsorted(tt, tg[-1]) + 2
                    if not good[max(i0, 0):i1].all(): continue
                j0, j1 = max(np.searchsorted(tt, tg[0]) - 10, 0), min(np.searchsorted(tt, tg[-1]) + 10, tt.size)
                W[k, c] = CubicSpline(tt[j0:j1], y[j0:j1])(tg).astype(np.float32)
    det = np.array([k[0] for k in keys], np.int64); sta = np.array([SIDX[k[1]] for k in keys], np.int8)
    np.savez_compressed(fo + ".tmp.npz", det_id=det, sta=sta, W=W)
    os.replace(fo + ".tmp.npz", fo)
    return m, f"{len(keys)} det-stations"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what"); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--months", default="")
    a = ap.parse_args()
    os.makedirs(ENV, exist_ok=True)
    if a.what == "select": select()
    elif a.what in ("run", "wave"):
        T = pd.read_parquet(f"{OUT}/targets.parquet", columns=["month"])
        ms = sorted(T.month.unique()) if not a.months else a.months.split(",")
        import multiprocessing as mp
        t00 = time.time()
        with mp.Pool(a.workers) as p:
            for m, msg in p.imap_unordered(do_month_wave if a.what == "wave" else do_month, ms):
                print(f"[{time.time() - t00:.0f}s] {m}: {msg}", flush=True)
