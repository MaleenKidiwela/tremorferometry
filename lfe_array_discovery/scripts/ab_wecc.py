#!/usr/bin/env python
"""Tremor DAYS before 2009-08-06 (no public location catalogue) by envelope cross-correlation,
after Wech & Creager (2008, GRL 35, L20302: "Automated detection and location of Cascadia
tremor") -- the method behind the PNSN catalogue that defines our tremor days from 2009-08-06.

WECC as implemented here (the published recipe, simplified where noted):
  * each station-day: 1.5-5.5 Hz band-pass (Wech & Creager 2008), envelope (|analytic signal|,
    root-sum-square over the available components), low-pass 0.1 Hz, decimated to 1 sps, log-free;
  * 5-minute windows (288 per day); per window each station envelope is demeaned/normalised;
  * all station-pair envelope cross-correlations over +-MAXLAG s;
  * grid search over epicentres (0.05 deg, subarray region + 40 km) at a fixed 40 km source
    depth: coherence(node) = mean over pairs of cc_ij at the predicted S-wave lag, S times from
    the GSC 1-D model (Vp / 1.73) -- the coherence-maximising epicentre is the window's location;
    (simplification: no bootstrap error ellipse; the coherence value is kept as a column)
  * a window is TREMOR if its coherence >= C_THR and its epicentre is within 25 km of the
    subarray's pool (the same in-region test as the PNSN day rule);
  * a day is a TREMOR DAY if it has >= N_WIN tremor windows (12 x 5 min = 1 hour, the same 1 h
    as the PNSN rule of >= 24 x 2.5-min epicentres).
C_THR is CALIBRATED against the PNSN day labels (same rule) on 2009-08-06..2011-12-31 (NW), the
agreement (accuracy, Cohen's kappa, TPR, FPR) is reported, and the calibrated detector is then
applied to 2007-01-01..2009-08-05 (NW) and 2007-01-01..2008-12-31 (SW).

Stations: the subarray's 7 sec.3.1 autocorrelation stations (continuous CAFE/PB/UW recording).

  python scripts/ab_wecc.py fetch            # raw station-days (ab_fetch)
  python scripts/ab_wecc.py env --workers 24 # 1-sps envelopes -> results/ab_royer/tremor/env/
  python scripts/ab_wecc.py detect           # window coherence -> tremor/wecc_windows_<SUB>.parquet
  python scripts/ab_wecc.py calibrate        # vs PNSN -> tremor/wecc_calibration.csv, wecc_days_*.csv
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

OUT = f"{A.RES}/tremor"
ENV = f"{OUT}/env"
SPAN = {"NW": ("2007-01-01", "2011-12-31"), "SW": ("2007-01-01", "2008-12-31")}
VALID = ("2009-08-06", "2011-12-31")          # PNSN catalogue era, NW
BAND = (1.5, 5.5)
WIN = 300
MAXLAG = 20
DEPTH = 40.0
N_WIN = 12
RADIUS = 25.0
VS_MODEL = ([0.0, 1.0, 6.0, 30.0, 46.0], [5.01, 6.01, 6.71, 7.11, 7.76])   # GSC onshore, /1.73


def stations(sub):
    return list(A.SUBARRAY[sub]["stas"])


def days(sub):
    a, b = SPAN[sub]
    return list(pd.date_range(a, b).strftime("%Y-%m-%d"))


def cmd_fetch(workers=8):
    from concurrent.futures import ThreadPoolExecutor
    import ab_fetch as F
    jobs = []
    for sub in ("NW", "SW"):
        for d in days(sub):
            for s in stations(sub):
                p = A.raw_path(s, d)
                if not os.path.exists(p) and not os.path.exists(p + ".nodata"):
                    jobs.append((s, d, p))
    jobs = sorted(set(jobs))
    A.log(f"[wecc fetch] {len(jobs)} station-days")
    rows = []
    with ThreadPoolExecutor(max_workers=min(workers, 8)) as ex:
        for i, r in enumerate(ex.map(F.fetch_one, jobs), 1):
            rows.append(r)
            if r["status"] != "data":
                p = A.raw_path(r["sta"], r["date"]) + ".nodata"
                os.makedirs(os.path.dirname(p), exist_ok=True)
                open(p, "w").write(r["err"][:200])
            if i % 200 == 0:
                A.log(f"  {i}/{len(jobs)}")
    pd.DataFrame(rows).to_csv(f"{OUT}/wecc_fetch_manifest.csv", index=False)


def env_one(job):
    sta, date = job
    dst = f"{ENV}/{A.netsta(sta)}_{date}.npy"
    if os.path.exists(dst):
        return dst
    src = A.raw_path(sta, date)
    if not os.path.exists(src):
        return None
    try:
        warnings.filterwarnings("ignore")
        from obspy import read
        from scipy.signal import butter, hilbert, sosfiltfilt
        from scipy.fft import next_fast_len
        st, _ = A.pick_band(read(src))
        # one trace per channel: many-fragment days otherwise cost one full-day transform PER fragment
        st.merge(method=0, fill_value=None)
        t_day = pd.Timestamp(date).timestamp()
        tot = None
        for tr in st:
            fs = float(tr.stats.sampling_rate)
            n = int(round(86400 * fs))
            x = np.zeros(n)
            v = np.zeros(n, bool)
            i0 = int(round((tr.stats.starttime.timestamp - t_day) * fs))
            a, b = max(i0, 0), min(i0 + tr.stats.npts, n)
            if b <= a:
                continue
            seg = tr.data[a - i0:b - i0]
            m = np.ma.getmaskarray(seg)
            x[a:b] = np.ma.filled(seg.astype(float), 0.0)
            v[a:b] = ~m
            x[v] -= x[v].mean()
            sos = butter(4, [BAND[0] / (fs / 2), BAND[1] / (fs / 2)], "band", output="sos")
            e = np.abs(hilbert(sosfiltfilt(sos, x), N=next_fast_len(n))[:n])
            k = int(round(fs))
            e = e[: (n // k) * k].reshape(-1, k).mean(1)          # 1 sps
            vv = v[: (n // k) * k].reshape(-1, k).all(1)
            e = np.where(vv, e, np.nan)
            tot = e ** 2 if tot is None else np.nansum([tot, e ** 2], axis=0)
        if tot is None:
            return None
        env = np.sqrt(tot)
        lp = butter(2, 0.1 / 0.5, "low", output="sos")
        good = np.isfinite(env)
        f = np.where(good, env, np.nanmedian(env) if good.any() else 0.0)
        env = np.where(good, sosfiltfilt(lp, f), np.nan).astype(np.float32)
        os.makedirs(ENV, exist_ok=True)
        np.save(dst, env)
        return dst
    except Exception:
        return None


def cmd_env(workers=24):
    import multiprocessing as mp
    jobs = sorted({(s, d) for sub in ("NW", "SW") for d in days(sub) for s in stations(sub)})
    A.log(f"[wecc env] {len(jobs)} station-days")
    with mp.Pool(min(workers, 30), maxtasksperchild=50) as pool:
        n = sum(1 for r in pool.imap_unordered(env_one, jobs, chunksize=4) if r)
    A.log(f"[wecc env] {n} envelopes")


def tt_table():
    sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
    import p_vel1d as V
    V.ZTOP = np.array(VS_MODEL[0])
    V.VP = np.array(VS_MODEL[1])
    xs = np.arange(0, 200.5, 1.0)
    ts = np.array([V.tt(float(x), DEPTH, "S", vpvs=1.73) for x in xs])
    return xs, ts


def cmd_detect():
    import torch
    dev = A.torch_dev()
    xs, ts = tt_table()
    for sub in ("NW", "SW"):
        stas = stations(sub)
        la0 = np.array([A.COORD[s][0] for s in A.POOL[sub]])
        lo0 = np.array([A.COORD[s][1] for s in A.POOL[sub]])
        glat = np.arange(la0.min() - 0.4, la0.max() + 0.4, 0.05)
        glon = np.arange(lo0.min() - 0.55, lo0.max() + 0.55, 0.05)
        G = np.array([(a, b) for a in glat for b in glon])
        inreg = A.epi_km(G[:, 0:1], G[:, 1:2], la0[None], lo0[None]).min(1) <= RADIUS
        T = np.stack([np.interp(A.epi_km(G[:, 0], G[:, 1], A.COORD[s][0], A.COORD[s][1]), xs, ts)
                      for s in stas], 1)                                    # (nG, nS)
        pairs = [(i, j) for i in range(len(stas)) for j in range(i + 1, len(stas))]
        rows = []
        for d in days(sub):
            E = []
            for s in stas:
                p = f"{ENV}/{A.netsta(s)}_{d}.npy"
                E.append(np.load(p) if os.path.exists(p) else np.full(86400, np.nan, np.float32))
            E = np.stack([e[:86400] if e.size >= 86400 else
                          np.r_[e, np.full(86400 - e.size, np.nan)] for e in E])
            nw = 86400 // WIN
            W = E[:, : nw * WIN].reshape(len(stas), nw, WIN)                 # (S, nw, WIN)
            ok = np.isfinite(W).all(2)                                      # (S, nw)
            Wn = np.nan_to_num(W - np.nanmean(W, 2, keepdims=True))
            Wn = Wn / np.maximum(np.linalg.norm(Wn, axis=2, keepdims=True), 1e-12)
            Wt = torch.from_numpy(Wn).to(dev).double()
            F = torch.fft.rfft(Wt, n=2 * WIN, dim=2)
            best = np.full(nw, np.nan)
            loc = np.full(nw, -1)
            npair = np.zeros(nw)
            acc = torch.zeros((nw, len(G)), device=dev, dtype=torch.float64)
            cnt = torch.zeros((nw, 1), device=dev, dtype=torch.float64)
            for i, j in pairs:
                cc = torch.fft.irfft(F[i] * F[j].conj(), n=2 * WIN, dim=1)  # (nw, 2WIN)
                lag = np.round(T[:, j] - T[:, i]).astype(int)               # t_j - t_i
                lag = np.clip(lag, -MAXLAG, MAXLAG)
                idx = torch.from_numpy((-lag) % (2 * WIN)).to(dev)          # cc at lag (j later)
                v = cc[:, idx]                                              # (nw, nG)
                okp = torch.from_numpy(ok[i] & ok[j]).to(dev).double()[:, None]
                acc += v * okp
                cnt += okp
                npair += ok[i] & ok[j]
            coh = (acc / cnt.clamp_min(1)).cpu().numpy()
            best = coh.max(1)
            loc = coh.argmax(1)
            for w in range(nw):
                rows.append(dict(date=d, win=w, coh=float(best[w]), n_pair=int(npair[w]),
                                 lat=float(G[loc[w], 0]), lon=float(G[loc[w], 1]),
                                 in_region=bool(inreg[loc[w]])))
        D = pd.DataFrame(rows, columns=["date", "win", "coh", "n_pair", "lat", "lon",
                                        "in_region"])
        D.loc[D.n_pair < 6, "coh"] = np.nan                                 # < 4 stations
        D.to_parquet(f"{OUT}/wecc_windows_{sub}.parquet", index=False)
        A.log(f"[wecc detect {sub}] {D.date.nunique()} days")


def pnsn_labels():
    import ab_tremor_days as TD
    T = TD.catalogue()
    out = {}
    for sub in ("NW",):
        la = np.array([A.COORD[s][0] for s in A.POOL[sub]])
        lo = np.array([A.COORD[s][1] for s in A.POOL[sub]])
        C = T[(T.time >= VALID[0]) & (T.time < pd.Timestamp(VALID[1]) + pd.Timedelta("1D"))]
        dk = A.epi_km(C.lat.to_numpy()[:, None], C.lon.to_numpy()[:, None], la, lo).min(1)
        n = C[dk <= TD.RADIUS_KM].groupby(C.time[dk <= TD.RADIUS_KM].dt.strftime("%Y-%m-%d")).size()
        dd = pd.date_range(*VALID).strftime("%Y-%m-%d")
        out[sub] = pd.Series([int(n.get(d, 0)) for d in dd], index=dd)
    return out


def cmd_calibrate():
    lab = pnsn_labels()["NW"]
    D = pd.read_parquet(f"{OUT}/wecc_windows_NW.parquet")
    V = D[D.date.isin(lab.index)]
    rows = []
    for thr in np.round(np.arange(0.30, 0.96, 0.01), 2):
        k = V[(V.coh >= thr) & V.in_region].groupby("date").size().reindex(lab.index).fillna(0)
        cov = V.groupby("date").coh.apply(lambda c: c.notna().mean()).reindex(lab.index)
        m = cov >= 0.5                                   # days the 7 stations actually recorded
        y, p = (lab[m] >= 24).to_numpy(), (k[m] >= N_WIN).to_numpy()
        tp, fp = int((y & p).sum()), int((~y & p).sum())
        fn, tn = int((y & ~p).sum()), int((~y & ~p).sum())
        n = tp + fp + fn + tn
        po = (tp + tn) / n
        pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n ** 2
        rows.append(dict(c_thr=thr, n_days=n, pnsn_tremor_days=int(y.sum()),
                         wecc_tremor_days=int(p.sum()), tp=tp, fp=fp, fn=fn, tn=tn,
                         accuracy=po, kappa=(po - pe) / (1 - pe) if pe < 1 else np.nan,
                         tpr=tp / max(tp + fn, 1), fpr=fp / max(fp + tn, 1)))
    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/wecc_calibration.csv", index=False)
    b = R.loc[R.kappa.idxmax()]
    A.log(f"[wecc calibrate] best kappa {b.kappa:.3f} at c_thr {b.c_thr}: accuracy {b.accuracy:.3f}"
          f" TPR {b.tpr:.3f} FPR {b.fpr:.3f} ({int(b.tp)} TP, {int(b.fp)} FP, {int(b.fn)} FN, "
          f"{int(b.tn)} TN over {int(b.n_days)} days)")
    for sub in ("NW", "SW"):
        D = pd.read_parquet(f"{OUT}/wecc_windows_{sub}.parquet")
        k = D[(D.coh >= b.c_thr) & D.in_region].groupby("date").size()
        cov = D.groupby("date").coh.apply(lambda c: c.notna().mean())
        out = pd.DataFrame(dict(n_tremor_windows=k.reindex(cov.index).fillna(0).astype(int),
                                coverage=cov))
        out["tremor_day"] = (out.n_tremor_windows >= N_WIN) & (out.coverage >= 0.5)
        out["c_thr"] = b.c_thr
        out.index.name = "date"
        out.to_csv(f"{OUT}/wecc_days_{sub}.csv")
        A.log(f"[wecc {sub}] tremor days by year: "
              + str(out[out.tremor_day].groupby(out[out.tremor_day].index.str[:4]).size()
                    .to_dict()))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    c = sys.argv[1]
    {"fetch": cmd_fetch, "env": cmd_env, "detect": cmd_detect,
     "calibrate": cmd_calibrate}[c]()
