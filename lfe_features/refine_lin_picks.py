"""Precise LFE arrival times by cross-correlation alignment (fixes the OT+11s label jitter).

Cut Z windows around each curated Lin LFE (search +/-6 s about OT+ANCHOR), iteratively align to
a reference stack by cross-correlation, take each event's best-lag as its REFINED arrival. Output
data/lin_precise_<sta>.csv (refined arrival epoch per event) for the per-sample picker labels.

Usage: PYTHONPATH=src python lfe_features/refine_lin_picks.py --net PB --sta B011 --slat 48.65 --slon -123.448
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy import signal
sys.path.insert(0, os.path.dirname(__file__))

ANCHOR = 11.0; SR = 50.0; SEARCH = 6.0      # +/- s search; window length = 2*SEARCH
HALF = int(SEARCH * SR); WL = 2 * HALF
HERE = os.path.dirname(__file__)


def cut_day(payload):
    net, sta, year, julday, items = payload   # items: (idx, center_epoch)
    from obspy import read, UTCDateTime
    path = f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"
    if not os.path.exists(path):
        return []
    try:
        st = read(path).select(channel="*Z"); st = st.merge(method=1, fill_value=0)
    except Exception:
        return []
    if len(st) == 0:
        return []
    tr = st[0]; fs = tr.stats.sampling_rate
    sos = signal.butter(4, [2, 8], "band", fs=fs, output="sos")
    day0 = UTCDateTime(year=year, julday=julday).timestamp
    out = []
    pad = 1.0
    for idx, cen in items:
        i0 = int(round((cen - day0 - SEARCH - pad) * fs)); n = int(round((2 * SEARCH + 2 * pad) * fs))
        if i0 < 0 or i0 + n > tr.stats.npts:
            continue
        x = signal.sosfiltfilt(sos, signal.detrend(tr.data[i0:i0 + n].astype(float)))
        if abs(fs - SR) > 0.1:
            g = math.gcd(int(round(fs)), int(SR)); x = signal.resample_poly(x, int(SR) // g, int(round(fs)) // g)
        s = int(round(pad * SR)); x = x[s:s + WL]
        if len(x) < WL or np.std(x) == 0:
            continue
        out.append((idx, (x / (np.std(x) + 1e-9)).astype(np.float32)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--n", type=int, default=6000); ap.add_argument("--workers", type=int, default=14)
    a = ap.parse_args(); s = a.sta.lower()
    files = glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed")
    have = set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0])) for f in files)
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dk = np.sqrt(((lin.lat - a.slat) * 111) ** 2 + ((lin.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    cur = lin[(lin.N >= 6) & (lin.residual < 1.2) & (dk < 35)].copy().reset_index(drop=True)
    cur["t"] = pd.to_datetime(cur.OT, utc=True); cur["cen"] = cur.t.values.astype("datetime64[ns]").astype("int64") / 1e9 + ANCHOR
    cur["y"] = cur.t.dt.year; cur["j"] = cur.t.dt.dayofyear
    cur = cur[[ (r.y, r.j) in have for r in cur.itertuples(index=False)]]
    if len(cur) > a.n:
        cur = cur.sample(a.n, random_state=1)
    cur = cur.reset_index(drop=True)
    jobs = defaultdict(list)
    for i, r in enumerate(cur.itertuples(index=False)):
        jobs[(int(r.y), int(r.j))].append((i, float(r.cen)))
    payloads = [(a.net, a.sta, k[0], k[1], v) for k, v in jobs.items()]
    print(f"[{a.sta}] cutting {len(cur)} Lin windows for alignment...", flush=True)
    idxs = []; W = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(cut_day, p) for p in payloads]):
            for idx, x in f.result():
                idxs.append(idx); W.append(x)
    W = np.stack(W); idxs = np.array(idxs)
    print(f"[{a.sta}] {len(W)} windows; iterative CC alignment...", flush=True)

    # iterative CC alignment to mean reference
    lags = np.zeros(len(W), int)
    ref = W.mean(0)
    for it in range(4):
        new = np.empty_like(W); newlags = np.zeros(len(W), int)
        rc = (ref - ref.mean()) / (np.std(ref) + 1e-9)
        for k in range(len(W)):
            c = signal.correlate(W[k], rc, mode="same")
            lag = c.argmax() - WL // 2
            lag = int(np.clip(lag, -HALF // 2, HALF // 2))
            newlags[k] = lag
            new[k] = np.roll(W[k], -lag)
        ref = new.mean(0); lags = newlags
        # alignment quality = mean CC of aligned windows with ref
        rc = (ref - ref.mean()) / (np.std(ref) + 1e-9)
        q = np.mean([np.corrcoef(np.roll(W[k], -lags[k]), rc)[0, 1] for k in range(0, len(W), 20)])
        print(f"  iter {it+1}: mean aligned CC vs ref = {q:.3f}", flush=True)

    refined = cur.copy()
    refined["lag_s"] = np.nan
    refined.loc[idxs, "lag_s"] = lags / SR
    refined["arrival"] = refined.cen + refined["lag_s"].fillna(0)
    refined = refined.dropna(subset=["lag_s"])
    refined[["OT", "lat", "lon", "arrival", "lag_s"]].to_csv(f"{HERE}/data/lin_precise_{s}.csv", index=False)
    np.save(f"{HERE}/data/lin_reference_{s}.npy", ref)
    print(f"[{a.sta}] saved data/lin_precise_{s}.csv ({len(refined)} refined picks; lag std {refined.lag_s.std():.2f}s)")
    print("REFINE DONE")


if __name__ == "__main__":
    main()
