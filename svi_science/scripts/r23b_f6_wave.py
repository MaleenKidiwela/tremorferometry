#!/usr/bin/env python
"""r23b_f6_wave (2026-10-11) -- Finding 6b waveforms: S6 matched-filter amplitudes on the fixed CN station set
(PGC, SNB, VGZ, NLLB) for 2003-2013 ETS events, read from raw miniSEED with the S6 processing.
Rules: svi_science/results/r23b/PLAN.md §3b (pre-registered).

  python r23b_f6_wave.py targets            # tables/f6/wave_targets.parquet (backward majors + forward validation events)
  python r23b_f6_wave.py extract --workers 8
  python r23b_f6_wave.py stacks             # S6 family x station reference stacks for the 4 CN stations (S6 pass 1)
  python r23b_f6_wave.py project            # a_mf, cc, noise projections -> tables/f6/mf_cn_backward.parquet
  python r23b_f6_wave.py validate           # raw reader vs S6 store on the forward validation events
"""
from __future__ import annotations
import argparse, glob, json, os, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
TB = f"{ROOT}/svi_science/results/r23b/tables"
OUT = f"{TB}/f6"
WAV = f"{OUT}/wave"
S6 = f"{ROOT}/svi_science/results/s6"
CUT = f"{ROOT}/lfe_extension/results/cut/svi_ext"
STAS = ["CN.NLLB", "CN.PFB", "CN.PGC", "CN.SNB", "CN.VGZ", "PB.B001", "PB.B003", "PB.B004", "PB.B005", "PB.B006", "PB.B007",
        "PB.B009", "PB.B010", "PB.B011", "PB.B926", "UW.LRIV"]          # S6 index order
CN = ["CN.NLLB", "CN.PGC", "CN.SNB", "CN.VGZ"]
SIDX = {s: i for i, s in enumerate(STAS)}
WT0, WN = -14.5, 380
CORD = {"Z": 0, "N": 1, "1": 1, "E": 2, "2": 2}
S0_, L, LAG = 280, 80, 5
NWIN = [5, 115]
os.makedirs(WAV, exist_ok=True)


# ------------------------------------------------------------------------------------------- targets
def targets():
    W = pd.read_csv(f"{TB}/windows.csv", parse_dates=["w0", "w1"])
    import importlib.util
    spec = importlib.util.spec_from_file_location("s3b", f"{ROOT}/svi_science/scripts/s3_build.py")
    s3b = importlib.util.module_from_spec(spec); spec.loader.exec_module(s3b)
    V = pd.read_parquet(f"{ROOT}/svi_science/results/r23/tables/events_2003_2026.parquet",
                        columns=["series", "fam", "time_utc", "date", "seg", "P", "detector", "src"])
    V = V[(V.series == "data") & (V.src == "bw")].reset_index(drop=True)
    M = W[W.major]
    V["ep"] = s3b.assign(V.time_utc, M)
    B = V[V.ep.notna()].copy()
    B["ev_id"] = 10_000_000 + B.index.astype(np.int64)
    B["kind"] = "bw"
    # forward validation events: S6 headline targets in the 2015-12 and 2018-06 padded windows
    T = pd.read_parquet(f"{S6}/targets.parquet")
    T = T[(T.kind == "event") & (T.P_same >= 0.75) & T.keep_xfam8]
    T["ep"] = s3b.assign(T.time_utc, M[M.episode.isin(["2015-12", "2018-06"])])
    T = T[T.ep.notna()].copy()
    Fv = pd.DataFrame(dict(ev_id=T.det_id.astype(np.int64), fam=T.fam.astype(int), time_utc=T.time_utc, ep=T.ep, kind="fwval",
                           seg="FWD", P=T.P_same, detector="forward"))
    X = pd.concat([B[["ev_id", "fam", "time_utc", "ep", "kind", "seg", "P", "detector"]], Fv], ignore_index=True)
    X["time_utc"] = X.time_utc.astype("datetime64[ns]")
    X.to_parquet(f"{OUT}/wave_targets.parquet", index=False)
    print(X.groupby(["kind", "seg"]).size())


# ------------------------------------------------------------------------------------------- gains
_G = None


def gain_rows():
    global _G
    if _G is not None:
        return _G
    I = pd.read_csv(f"{ROOT}/svi_science/results/backward_scope/inventory_channels.csv")
    I = I[(I.centre == "NRCAN") & (I.net == "CN") & I.sta.isin([s.split(".")[1] for s in CN])]
    g = [dict(sta="CN." + r.sta, cha=r.cha, start=pd.Timestamp(r.start), end=pd.Timestamp(r.end), gain=float(r.scale), src="inventory")
         for r in I.itertuples()]
    E = pd.read_csv(f"{CUT}/channel_epochs.csv")
    E = E[E.netsta.isin(CN)]
    g += [dict(sta=r.netsta, cha=r.cha, start=pd.Timestamp(r.start), end=pd.Timestamp(r.end), gain=float(r.gain), src="cut_epochs")
          for r in E.itertuples()]
    _G = pd.DataFrame(g)
    return _G


def gain_of(G, sta, cha, t):
    t = pd.Timestamp(t)
    g = G[(G.sta == sta) & (G.cha == cha) & (G.start <= t) & (G.end > t)]
    if not len(g):
        return np.nan
    pref = g[g.src == ("cut_epochs" if t >= pd.Timestamp("2014-01-01") else "inventory")]
    return float((pref if len(pref) else g).gain.iloc[0])


# ------------------------------------------------------------------------------------------- extraction
def _job(args):
    sta, date, rows, kind = args
    from obspy import read, UTCDateTime, Stream
    from scipy.signal import butter, sosfiltfilt
    from scipy.interpolate import CubicSpline
    G = gain_rows()
    bases = ([f"{ROOT}/lfe_extension/traces_backward/{sta}"] if kind == "bw" else []) + [f"{ROOT}/lfe_extension/traces/{sta}"]
    d0 = pd.Timestamp(date)
    st = Stream()
    for dd in (d0 - pd.Timedelta(days=1), d0, d0 + pd.Timedelta(days=1)):
        f = next((f"{b}/{dd.year}/{dd.dayofyear:03d}.mseed" for b in bases if os.path.exists(f"{b}/{dd.year}/{dd.dayofyear:03d}.mseed")),
                 f"{bases[0]}/{dd.year}/{dd.dayofyear:03d}.mseed")
        lo = min(r[2] for r in rows) - 20; hi = max(r[2] for r in rows) + 35
        dstart = UTCDateTime(dd.strftime("%Y-%m-%d")); dend = dstart + 86400
        if os.path.exists(f) and (dend > lo and dstart < hi):
            try:
                st += read(f)
            except Exception:
                pass
    if not len(st):
        return sta, date, [], "nofile"
    # one band per station-day: prefer HH, else BH, else EH
    bands = sorted({tr.stats.channel[:2] for tr in st})
    band = "HH" if "HH" in bands else ("BH" if "BH" in bands else bands[0])
    st = st.select(channel=f"{band}?")
    try:
        st.merge(method=1, fill_value=None)
    except Exception:
        for tr in st:
            tr.data = tr.data.astype(np.float64)
        st.merge(method=1, fill_value=None)
    sos = {}
    out = []
    tg = WT0 + np.arange(WN) / 20.0
    for ev_id, fam, s_epoch in rows:
        Wv = np.full((3, WN), np.nan, np.float32); ok_any = False
        for tr in st:
            c = CORD.get(tr.stats.channel[-1])
            if c is None:
                continue
            ws = UTCDateTime(s_epoch - 15.0); we = UTCDateTime(s_epoch + 30.0)
            t2 = tr.slice(ws, we)
            if t2.stats.npts < 50:
                continue
            dat = t2.data
            gapm = np.ma.getmaskarray(dat) if np.ma.isMaskedArray(dat) else np.zeros(t2.stats.npts, bool)
            x = np.ma.filled(dat.astype(np.float64), 0.0) if np.ma.isMaskedArray(dat) else dat.astype(np.float64)
            fs = float(t2.stats.sampling_rate)
            if fs not in sos:
                sos[fs] = butter(4, [2.0 / (fs / 2), 8.0 / (fs / 2)], btype="band", output="sos")
            n_gap = int(gapm.sum())
            good = (x != 0) if n_gap > 0 else np.ones(x.size, bool)
            if good.sum() < 50:
                continue
            x = np.where(good, x - x[good].mean(), 0.0)
            gain = gain_of(G, sta, t2.stats.channel, t2.stats.starttime.datetime)
            y = sosfiltfilt(sos[fs], x) / gain
            tt = float(t2.stats.starttime.timestamp) + np.arange(x.size) / fs - s_epoch
            if tg[0] < tt[0] or tg[-1] > tt[-1]:
                continue
            if n_gap > 0:
                i0, i1 = np.searchsorted(tt, tg[0]) - 2, np.searchsorted(tt, tg[-1]) + 2
                if not good[max(i0, 0):i1].all():
                    continue
            j0, j1 = max(np.searchsorted(tt, tg[0]) - 10, 0), min(np.searchsorted(tt, tg[-1]) + 10, tt.size)
            Wv[c] = CubicSpline(tt[j0:j1], y[j0:j1])(tg).astype(np.float32); ok_any = True
        if ok_any:
            out.append((ev_id, SIDX[sta], Wv))
    return sta, date, out, band


def extract(workers):
    import multiprocessing as mp
    X = pd.read_parquet(f"{OUT}/wave_targets.parquet")
    off = pd.read_csv(f"{TB}/cn_offsets.csv")
    off = {(int(r.fam), r.netsta): float(r.off) for r in off.itertuples()}
    X["t_epoch"] = X.time_utc.astype("int64") / 1e9
    jobs = []
    for sta in CN:
        x = X.copy()
        x["s_epoch"] = [t + off.get((f, sta), np.nan) for t, f in zip(x.t_epoch, x.fam)]
        x = x[np.isfinite(x.s_epoch)]
        x["sdate"] = pd.to_datetime(x.s_epoch, unit="s").dt.strftime("%Y-%m-%d")
        for (d, k), g in x.groupby(["sdate", "kind"]):
            jobs.append((sta, d, list(zip(g.ev_id.astype(np.int64), g.fam.astype(int), g.s_epoch.astype(float))), "bw" if k == "bw" else "fw"))
    print(len(jobs), "station-day jobs", flush=True)
    t0 = time.time(); acc = {}; stat = []
    with mp.Pool(workers) as p:
        for i, (sta, d, out, band) in enumerate(p.imap_unordered(_job, jobs, chunksize=2)):
            yr = d[:4]
            acc.setdefault(yr, []).extend(out)
            stat.append(dict(sta=sta, date=d, n_ok=len(out), band=band))
            if i % 200 == 0:
                print(f"[{time.time() - t0:.0f}s] {i}/{len(jobs)}", flush=True)
    for yr, L_ in acc.items():
        if not L_:
            continue
        np.savez_compressed(f"{WAV}/wave_{yr}.npz", ev_id=np.array([a[0] for a in L_], np.int64),
                            sta=np.array([a[1] for a in L_], np.int8), W=np.stack([a[2] for a in L_]))
    pd.DataFrame(stat).to_csv(f"{OUT}/extract_status.csv", index=False)
    print("done", time.time() - t0)


# ------------------------------------------------------------------------------------------- cut store (S6 path), 2010-08..2013
def _cut_month(args):
    m, T = args
    import h5py, pyarrow.parquet as pq
    from scipy.signal import butter, sosfiltfilt
    from scipy.interpolate import CubicSpline
    G = gain_rows()
    fi = f"{CUT}/index_{m}.parquet"
    if not os.path.exists(fi):
        return m, None
    I = pq.read_table(fi, columns=["det_id", "netsta", "cha", "fs", "t0", "offset_s", "n_samples", "n_gap", "dataset", "h5_start", "h5_file"]).to_pandas()
    I = I[I.det_id.isin(T.det_id) & I.netsta.isin(CN)]
    te = dict(zip(T.det_id, T.t_epoch)); evid = dict(zip(T.det_id, T.ev_id))
    keys = I.groupby(["det_id", "netsta"]).size().index.tolist(); kpos = {k: i for i, k in enumerate(keys)}
    Wm = np.full((len(keys), 3, WN), np.nan, np.float32); sos = {}
    tg = WT0 + np.arange(WN) / 20.0
    I = I.sort_values(["h5_file", "dataset", "h5_start"])
    for h5f, g in I.groupby("h5_file"):
        with h5py.File(f"{CUT}/{h5f}", "r") as h:
            for r in g.itertuples():
                c = CORD.get(r.cha[-1])
                if c is None:
                    continue
                k = kpos[(r.det_id, r.netsta)]
                x = h[r.dataset][int(r.h5_start):int(r.h5_start) + int(r.n_samples)].astype(np.float64)
                if x.size < 50:
                    continue
                fs = float(r.fs)
                if fs not in sos:
                    sos[fs] = butter(4, [2.0 / (fs / 2), 8.0 / (fs / 2)], btype="band", output="sos")
                good = x != 0 if r.n_gap > 0 else np.ones(x.size, bool)
                if good.sum() < 50:
                    continue
                x = np.where(good, x - x[good].mean(), 0.0)
                gain = gain_of(G, r.netsta, r.cha, pd.Timestamp(float(r.t0), unit="s"))
                y = sosfiltfilt(sos[fs], x) / gain
                tt = float(r.t0) + np.arange(x.size) / fs - (te[r.det_id] + float(r.offset_s))
                if tg[0] < tt[0] or tg[-1] > tt[-1]:
                    continue
                if r.n_gap > 0:
                    i0, i1 = np.searchsorted(tt, tg[0]) - 2, np.searchsorted(tt, tg[-1]) + 2
                    if not good[max(i0, 0):i1].all():
                        continue
                j0, j1 = max(np.searchsorted(tt, tg[0]) - 10, 0), min(np.searchsorted(tt, tg[-1]) + 10, tt.size)
                Wm[k, c] = CubicSpline(tt[j0:j1], y[j0:j1])(tg).astype(np.float32)
    keep = np.isfinite(Wm).any(axis=(1, 2))
    return m, (np.array([evid[k[0]] for k in keys], np.int64)[keep], np.array([SIDX[k[1]] for k in keys], np.int8)[keep], Wm[keep])


def cutstore(workers):
    import multiprocessing as mp
    X = pd.read_parquet(f"{OUT}/wave_targets.parquet")
    X = X[(X.kind == "bw") & (X.time_utc >= "2010-08-01") & (X.time_utc < "2014-01-01")].sort_values("time_utc")
    D = pd.read_parquet(f"{CUT}/detections.parquet", columns=["det_id", "fam", "t", "t_epoch"])
    D["fam"] = D.fam.astype(int); D["t"] = D.t.astype("datetime64[ns]"); D = D.sort_values("t")
    m = pd.merge_asof(X, D, left_on="time_utc", right_on="t", by="fam", tolerance=pd.Timedelta("0.1s"), direction="nearest")
    print("cut-store match", float(m.det_id.notna().mean()), len(m), flush=True)
    m = m.dropna(subset=["det_id"]); m["det_id"] = m.det_id.astype(np.int64)
    m["month"] = m.time_utc.dt.strftime("%Y-%m")
    jobs = [(mo, g[["det_id", "ev_id", "t_epoch"]]) for mo, g in m.groupby("month")]
    with mp.Pool(workers) as p:
        for mo, res in p.imap_unordered(_cut_month, jobs):
            if res is None or not len(res[0]):
                print(mo, "none", flush=True); continue
            np.savez_compressed(f"{WAV}/cut_{mo}.npz", ev_id=res[0], sta=res[1], W=res[2])
            print(mo, len(res[0]), flush=True)


# ------------------------------------------------------------------------------------------- stacks (S6 pass 1, CN only)
def stacks():
    T = pd.read_parquet(f"{S6}/targets.parquet")
    T["head"] = (T.P_same >= 0.75) & T.keep_xfam8 & (T.kind == "event")
    fam = dict(zip(T.det_id, T.fam)); head = set(T.det_id[T["head"]]); ev05 = set(T.det_id[(T.kind == "event")])
    cn = {SIDX[s] for s in CN}
    acc, accb = {}, {}
    for f in sorted(glob.glob(f"{S6}/env/wave_*.npz")):
        z = np.load(f); W = z["W"][:, :, S0_ - LAG:S0_ + L + LAG].astype(np.float64)
        for i, (d, s) in enumerate(zip(z["det_id"], z["sta"])):
            if int(s) not in cn or d not in ev05:
                continue
            k = (fam[d], int(s)); x = W[i]; ok = np.isfinite(x).all(axis=1)
            if not ok.any():
                continue
            for D_, use in ((accb, True), (acc, d in head)):
                if not use:
                    continue
                a = D_.setdefault(k, [np.zeros_like(x), np.zeros(3)])
                a[0][ok] += x[ok]; a[1][ok] += 1
    keys, Ts, src, nst = [], [], [], []
    for k, (sm, n) in acc.items():
        if n.max() >= 20:
            keys.append(k); Ts.append(sm[:, LAG:LAG + L] / np.maximum(n, 1)[:, None]); src.append("headline"); nst.append(n.max())
    for k, (sm, n) in accb.items():
        if k not in set(keys) and n.max() >= 5:
            keys.append(k); Ts.append(sm[:, LAG:LAG + L] / np.maximum(n, 1)[:, None]); src.append("P050_fallback"); nst.append(n.max())
    np.savez(f"{OUT}/cn_stacks.npz", fam=np.array([k[0] for k in keys]), sta=np.array([k[1] for k in keys]), T=np.stack(Ts),
             src=np.array(src), n=np.array(nst))
    chk = pd.read_csv(f"{S6}/family_station_stacks.csv")
    mine = pd.DataFrame(dict(fam=[k[0] for k in keys], sta=[STAS[k[1]] for k in keys], rms=[float(np.sqrt(np.mean(t ** 2))) for t in Ts]))
    m = mine.merge(chk, on=["fam", "sta"])
    print("stacks", len(keys), "vs S6 rms max rel diff", float((m.rms / m.stack_rms - 1).abs().max()), "n matched", len(m))


def project():
    z = np.load(f"{OUT}/cn_stacks.npz", allow_pickle=True)
    stk = {(int(f), int(s)): (T, n, src) for f, s, T, n, src in zip(z["fam"], z["sta"], z["T"], z["n"], z["src"])}
    X = pd.read_parquet(f"{OUT}/wave_targets.parquet"); fam = dict(zip(X.ev_id, X.fam))
    rows = []
    for f in sorted(glob.glob(f"{WAV}/wave_*.npz")) + sorted(glob.glob(f"{WAV}/cut_*.npz")):
        reader = "cut" if os.path.basename(f).startswith("cut_") else "raw"
        q = np.load(f); Wall = q["W"].astype(np.float64)
        for i, (d, s) in enumerate(zip(q["ev_id"], q["sta"])):
            k = (fam[int(d)], int(s))
            if k not in stk:
                continue
            Tk, nst, src = stk[k]; x = Wall[i]
            ok = np.isfinite(x).all(axis=1) & (np.abs(Tk).sum(axis=1) > 0)
            if ok.sum() < 2:
                continue
            Tc = Tk[ok]; xc = x[ok]; tt = (Tc * Tc).sum()
            SW = np.lib.stride_tricks.sliding_window_view(xc, L, axis=1)

            def proj(st):
                seg = SW[:, st - LAG:st + LAG + 1, :]
                num = np.einsum("clk,ck->l", seg, Tc); j = int(np.argmax(num))
                return num[j] / tt, num[j] / np.sqrt((seg[:, j] ** 2).sum() * tt + 1e-60)
            a, cc = proj(S0_); n1, _ = proj(NWIN[0]); n2, _ = proj(NWIN[1])
            # 3-C envelope for the S-window check (0.1 s bins over -14.5..+4.5 s; normalised later)
            e = np.nansum(x ** 2, axis=0).reshape(-1, 2).mean(1)
            rows.append((int(d), STAS[int(s)], a, cc, n1, n2, int(ok.sum()), int(nst), str(src), reader, e.astype(np.float32)))
    M = pd.DataFrame(rows, columns=["ev_id", "sta", "a_mf", "cc_mf", "an1", "an2", "ncomp_mf", "n_stack", "stack_src", "reader", "env"])
    M["sig_noise"] = np.sqrt((M.an1 ** 2 + M.an2 ** 2) / 2)
    env = np.stack(M.pop("env").to_numpy())
    np.save(f"{OUT}/env_cn_backward.npy", env)
    M.to_parquet(f"{OUT}/mf_cn_backward.parquet")
    print(len(M), M.groupby("sta").size().to_dict())


def validate():
    M = pd.read_parquet(f"{OUT}/mf_cn_backward.parquet")
    X = pd.read_parquet(f"{OUT}/wave_targets.parquet")
    v = M[M.ev_id.isin(X[X.kind == "fwval"].ev_id)]
    S = pd.read_parquet(f"{S6}/mf_station.parquet")
    m = v.merge(S[["det_id", "sta", "a_mf", "sig_noise"]].rename(columns={"det_id": "ev_id"}), on=["ev_id", "sta"], suffixes=("", "_s6"))
    m = m[(m.a_mf > 0) & (m.a_mf_s6 > 0)]
    d = np.abs(np.log10(m.a_mf / m.a_mf_s6)); dn = np.abs(np.log10(m.sig_noise / m.sig_noise_s6))
    out = dict(n=int(len(m)), median_abs_dlog=float(d.median()), q90_abs_dlog=float(d.quantile(0.9)),
               frac_le_0p05=float((d <= 0.05).mean()), noise_median_abs_dlog=float(dn.median()),
               by_sta={s: float(g.median()) for s, g in d.groupby(m.sta)},
               pass_=bool(d.median() <= 0.02 and (d <= 0.05).mean() >= 0.9))
    bw = M[M.ev_id >= 10_000_000]
    rc = bw[bw.reader == "raw"].merge(bw[bw.reader == "cut"], on=["ev_id", "sta"], suffixes=("_raw", "_cut"))
    rc = rc[(rc.a_mf_raw > 0) & (rc.a_mf_cut > 0)]
    d2 = np.abs(np.log10(rc.a_mf_raw / rc.a_mf_cut))
    out["raw_vs_cut_2010_13"] = dict(n=int(len(rc)), median_abs_dlog=float(d2.median()) if len(rc) else None,
                                      q90_abs_dlog=float(d2.quantile(0.9)) if len(rc) else None,
                                      frac_le_0p05=float((d2 <= 0.05).mean()) if len(rc) else None,
                                      by_sta={s_: float(g.median()) for s_, g in d2.groupby(rc.sta)},
                                      pass_=bool(len(rc) and d2.median() <= 0.02 and (d2 <= 0.05).mean() >= 0.9))
    json.dump(out, open(f"{OUT}/reader_validation.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what"); ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    {"targets": targets, "stacks": stacks, "project": project, "validate": validate,
     "cutstore": lambda: cutstore(a.workers)}.get(a.what, lambda: extract(a.workers))()
