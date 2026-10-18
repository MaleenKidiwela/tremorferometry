#!/usr/bin/env python
"""s10_windows -- S10 pilot step 2 (PLAN §2): read, filter, despike and RMS-normalise every data and scrambled-control
window at every pilot station. Data rows in the cut store are read from it; everything else from the raw miniSEED day
files at the same offsets and over the same S-15..S+30 s span, so the processing is identical. No new cut stage.

    python s10_windows.py run --workers 12      # -> win/<NET.STA>.npz
    python s10_windows.py verify                # cut store vs raw, same windows
"""
from __future__ import annotations
import os, sys, glob, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S

CORD = {"Z": 0, "N": 1, "1": 1, "E": 2, "2": 2}
_SOS = {}


def sos_for(fs):
    from scipy.signal import butter
    if fs not in _SOS:
        _SOS[fs] = butter(4, [2.0 / (fs / 2), 8.0 / (fs / 2)], btype="band", output="sos")
    return _SOS[fs]


def process(x, fs, t0, tS):
    """x: native samples starting at t0 (epoch s) spanning ~tS-15..tS+30. -> (NW,) float on the tau grid, or None."""
    from scipy.signal import sosfiltfilt
    from scipy.signal.windows import tukey
    from scipy.interpolate import CubicSpline
    x = np.asarray(x, np.float64)
    if x.size < 100:
        return None
    tt = t0 + np.arange(x.size) / fs - tS
    if tt[0] > S.T0W - 0.05 or tt[-1] < S.T1W + 0.05:
        return None
    # gaps = runs of exact zeros (cut-store and merge convention) inside the stored span
    m = (tt >= S.T0W - 0.1) & (tt <= S.T1W + 0.1)
    z = x[m] == 0
    if z.any():
        run = np.convolve(z.astype(int), np.ones(5, int), "valid")
        if (run >= 5).any():
            return None
    good = x != 0
    x = np.where(good, x - x[good].mean(), 0.0) if good.any() else x
    x = x * tukey(x.size, 0.1)
    y = sosfiltfilt(sos_for(float(fs)), x)
    j0 = max(np.searchsorted(tt, S.T0W) - 10, 0); j1 = min(np.searchsorted(tt, S.T1W) + 10, tt.size)
    return CubicSpline(tt[j0:j1], y[j0:j1])(S.TAU)


class Raw:
    """lazy raw day-file reader for one station, with neighbouring days."""
    def __init__(self, netsta):
        self.ns = netsta; self.cache = {}

    def day(self, d):
        if d in self.cache:
            return self.cache[d]
        import obspy
        y, doy = d.year, d.timetuple().tm_yday
        out = None
        for base in ((S.RAW_BW, S.RAW_FW) if y <= 2013 else (S.RAW_FW, S.RAW_BW)):
            p = f"{base}/{self.ns}/{y}/{doy:03d}.mseed"
            if os.path.exists(p):
                try:
                    st = obspy.read(p)
                except Exception:
                    continue
                bands = {}
                for tr in st:
                    c = tr.stats.channel
                    if c[-1] in CORD:
                        bands.setdefault(c[:2], set()).add(CORD[c[-1]])
                full = [b for b, cs in bands.items() if 0 in cs]
                if not full:
                    continue
                band = sorted(full, key=lambda b: (-len(bands[b]), {"HH": 0, "EH": 1, "BH": 2}.get(b, 3)))[0]
                st = st.select(channel=f"{band}?")
                comps = {}
                for c in set(tr.stats.channel for tr in st):
                    s2 = st.select(channel=c).copy()
                    try:
                        s2.merge(method=0, fill_value=0)
                    except Exception:
                        continue
                    tr = s2[0]
                    comps[CORD[c[-1]]] = (tr.data.astype(np.float64), float(tr.stats.sampling_rate),
                                          float(tr.stats.starttime.timestamp))
                out = comps
                break
        self.cache[d] = out
        if len(self.cache) > 3:
            self.cache.pop(next(iter(self.cache)))
        return out

    def window(self, tS):
        """3 x (x, fs, t0) slices spanning tS-15..tS+30 s, or None per component."""
        a, b = tS - S.CUT_PRE, tS + S.CUT_POST
        res = [None, None, None]
        d0 = pd.Timestamp(a, unit="s").normalize(); d1 = pd.Timestamp(b, unit="s").normalize()
        for c in range(3):
            pieces = []
            for d in sorted({d0, d1, d0 - pd.Timedelta(days=1)}):
                D = self.day(d)
                if not D or c not in D:
                    continue
                x, fs, t0 = D[c]
                i0 = int(np.floor((a - t0) * fs)); i1 = int(np.ceil((b - t0) * fs)) + 1
                if i1 <= 0 or i0 >= x.size:
                    continue
                pieces.append((max(i0, 0), min(i1, x.size), x, fs, t0, i0 >= 0 and i1 <= x.size))
            full = [p for p in pieces if p[5]]
            if full:
                i0, i1, x, fs, t0, _ = full[0]
                res[c] = (x[i0:i1], fs, t0 + i0 / fs)
        return res


def cut_windows(month, det_ids, netsta):
    import h5py, pyarrow.parquet as pq
    p = f"{S.CUT}/index_{month}.parquet"
    if not os.path.exists(p) or not det_ids:
        return {}
    cols = ["det_id", "netsta", "cha", "fs", "t0", "offset_s", "n_samples", "n_gap", "dataset", "h5_start", "h5_file"]
    I = pq.read_table(p, columns=cols, filters=[("det_id", "in", list(det_ids)), ("netsta", "==", netsta)]).to_pandas()
    out = {}
    for h5f, g in I.groupby("h5_file"):
        with h5py.File(f"{S.CUT}/{h5f}", "r") as h:
            for r in g.itertuples():
                c = CORD.get(r.cha[-1])
                if c is None:
                    continue
                x = h[r.dataset][int(r.h5_start):int(r.h5_start) + int(r.n_samples)]
                out.setdefault(int(r.det_id), [None, None, None])[c] = (x, float(r.fs), float(r.t0), float(r.offset_s))
    return out


def task(arg):
    netsta, date, R = arg
    try:
        raw = Raw(netsta)
        n = len(R)
        W = np.full((n, 3, S.NW), np.nan, np.float32)
        src = np.zeros(n, np.int8)            # 1 cut store, 2 raw
        cutd = {}
        dids = [int(x) for x in R.det_id if x >= 0]
        if dids:
            cutd = cut_windows(date[:7], dids, netsta)
        for j, r in enumerate(R.itertuples()):
            cw = cutd.get(int(r.det_id)) if r.det_id >= 0 else None
            if cw is not None and any(v is not None for v in cw):
                for c in range(3):
                    if cw[c] is None:
                        continue
                    x, fs, t0, off = cw[c]
                    y = process(x, fs, t0, r.tS)
                    if y is not None:
                        W[j, c] = y
                src[j] = 1
            else:
                for c, v in enumerate(raw.window(r.tS)):
                    if v is None:
                        continue
                    y = process(v[0], v[1], v[2], r.tS)
                    if y is not None:
                        W[j, c] = y
                src[j] = 2
        # despike: clip at 8 x MAD (x1.4826), MAD = median over the station-day's data windows, per component
        rms_raw = np.full((n, 3), np.nan, np.float32)
        isdata = (R.series.to_numpy() == "data")
        for c in range(3):
            ok = np.isfinite(W[:, c, 0])
            if not ok.any():
                continue
            mads = 1.4826 * np.median(np.abs(W[ok, c] - np.median(W[ok, c], axis=1, keepdims=True)), axis=1)
            ref = np.median(mads[isdata[ok]]) if (isdata & ok).any() else np.median(mads)
            lim = 8.0 * ref
            W[ok, c] = np.clip(W[ok, c], -lim, lim)
            rr = np.sqrt(np.mean(W[ok, c].astype(np.float64) ** 2, axis=1))
            rms_raw[ok, c] = rr
            W[ok, c] = W[ok, c] / np.where(rr > 0, rr, np.nan)[:, None]
        return netsta, date, R.row_id.to_numpy(), W, src, rms_raw, ""
    except Exception as e:
        import traceback
        return netsta, date, R.row_id.to_numpy(), None, None, None, f"{e} {traceback.format_exc()[-400:]}"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("what"); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--stations", default="")
    ap.add_argument("--c1f", action="store_true", help="frame-time scrambled-control windows (Revision 1)")
    a = ap.parse_args()
    S.ensure()
    log = S.Log(f"{S.LOGD}/windows.log")
    RS = pd.read_parquet(f"{S.TB}/row_station.parquet")
    X = pd.read_parquet(f"{S.TB}/rows.parquet", columns=["row_id", "det_id", "series"])
    RS = RS.merge(X, on="row_id")
    sfx = ""
    if a.c1f:
        # Revision 1: scrambled-control row read at its detector frame time t + offset (no site shift) = a same-day
        # time unrelated to what this station contributed to the scrambled detection
        Xr = pd.read_parquet(f"{S.TB}/rows.parquet", columns=["row_id", "t", "fam"])
        O = pd.read_csv(f"{S.TB}/offsets_sp.csv").set_index(["fam", "netsta"]).offset
        RS = RS[RS.series == "c1"].drop(columns=["tS"]).merge(Xr, on="row_id")
        RS["tS"] = RS.t + np.array([O.get((f, ns), np.nan) for f, ns in zip(RS.fam, RS.netsta)])
        RS = RS[np.isfinite(RS.tS)].drop(columns=["t", "fam"]); RS["det_id"] = -1
        sfx = "__c1f"
    if a.what == "verify":
        V = RS[(RS.det_id >= 0) & (RS.date >= "2014")]
        V = V[V.netsta.isin(["CN.PGC", "PB.B011", "PB.B926", "CN.VGZ"])].groupby("netsta").head(25)
        res = []
        for (ns, d), g in V.groupby(["netsta", "date"]):
            _, _, rid, W1, _, _, e1 = task((ns, d, g))
            g2 = g.assign(det_id=-1)
            _, _, _, W2, _, _, e2 = task((ns, d, g2))
            for j in range(len(rid)):
                for c in range(3):
                    if np.isfinite(W1[j, c, 0]) and np.isfinite(W2[j, c, 0]):
                        res.append(dict(netsta=ns, row_id=rid[j], comp=c,
                                        maxdiff=float(np.max(np.abs(W1[j, c] - W2[j, c]))),
                                        cc=float(np.corrcoef(W1[j, c], W2[j, c])[0, 1])))
        Rv = pd.DataFrame(res); Rv.to_csv(f"{S.TB}/verify_cut_vs_raw.csv", index=False)
        log(f"verify cut vs raw: {len(Rv)} window-components, max |diff| (RMS units) {Rv.maxdiff.max():.2e}, min cc {Rv.cc.min():.6f}")
        return
    if a.stations:
        RS = RS[RS.netsta.isin(a.stations.split(","))]
    tasks = [(ns, d, g.reset_index(drop=True)) for (ns, d), g in RS.groupby(["netsta", "date"])]
    tasks.sort(key=lambda t: -len(t[2]))
    log(f"windows: {len(tasks)} station-days, {len(RS)} row-stations, {a.workers} workers")
    acc = {}
    import multiprocessing as mp
    t00 = time.time(); nerr = 0
    with mp.get_context("fork").Pool(a.workers) as p:
        for k, (ns, d, rid, W, src, rr, err) in enumerate(p.imap_unordered(task, tasks, chunksize=1)):
            if err:
                nerr += 1; log(f"ERROR {ns} {d}: {err}"); continue
            acc.setdefault(ns, []).append((rid, W, src, rr))
            if k % 200 == 0:
                log(f"  {k + 1}/{len(tasks)} {time.time() - t00:.0f}s")
    for ns, L in acc.items():
        rid = np.concatenate([x[0] for x in L]); W = np.concatenate([x[1] for x in L])
        src = np.concatenate([x[2] for x in L]); rr = np.concatenate([x[3] for x in L])
        np.savez(f"{S.WIN}/{ns}{sfx}.npz.tmp.npz", row_id=rid, W=W, src=src, rms_raw=rr, tau=S.TAU)
        os.replace(f"{S.WIN}/{ns}{sfx}.npz.tmp.npz", f"{S.WIN}/{ns}{sfx}.npz")
        okz = np.isfinite(W[:, 0, 0])
        log(f"{ns}: {len(rid)} windows, Z ok {okz.sum()}, cut {int((src == 1).sum())} raw {int((src == 2).sum())}")
    log(f"windows done, {nerr} station-day errors, {time.time() - t00:.0f}s")
    open(f"{S.OUT}/WINDOWS{sfx}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
