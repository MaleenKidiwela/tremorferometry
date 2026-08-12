"""Score discovery candidates with the tremor-picker -> P(LFE) column, then write a
BASELINE parquet (all scorable candidates) and a FILTERED parquet (P(LFE) >= thr), so
discover_gpu can be run on each to compare the families that form.

Usage: PYTHONPATH=src python lfe_features/score_candidates.py --net PB --sta B011 \
         --cand data/b011_pnsn_candidates_100km.parquet --y0 2010 --y1 2013 --thr 0.4 --workers 14
"""
import argparse, os, sys
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all, FEATURE_ORDER
import joblib
HERE = os.path.dirname(__file__); PRE, POST = 10.0, 30.0


def cut_day(payload):
    net, sta, year, julday, items, cols, target_fs = payload      # items: (idx, epoch)
    from obspy import read, UTCDateTime
    path = f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"
    if not os.path.exists(path):
        return []
    try:
        st = read(path)
    except Exception:
        return []
    # decimate the full day to target_fs BEFORE slicing so features match the picker's training fs (40 Hz for
    # the broadband picker; native 100 Hz HHZ would give a different Nyquist -> wrong spectral features).
    if target_fs:
        from fractions import Fraction
        from scipy.signal import resample_poly
        for tr in st:
            nat = float(tr.stats.sampling_rate)
            if abs(nat - target_fs) > 1e-6:
                fr = Fraction(int(round(target_fs)), int(round(nat))).limit_denominator(1000)
                try:
                    tr.data = resample_poly(tr.data.astype(np.float64), fr.numerator, fr.denominator).astype(np.float32)
                    tr.stats.sampling_rate = target_fs
                except Exception:
                    pass

    def get(*ch):
        for c in ch:
            s = st.select(channel=f"*{c}")
            if len(s):
                try: s = s.merge(method=1, fill_value=0)
                except Exception: pass
                return s[0]
        return None
    trZ, trH1, trH2 = get("Z"), get("1", "N"), get("2", "E")
    if trZ is None:
        return []
    fs = float(trZ.stats.sampling_rate); exp = int(round((PRE + POST) * fs))
    out = []
    for idx, ep in items:
        a = UTCDateTime(ep)
        try:
            z = trZ.slice(a - PRE, a + POST).data.astype(float)
        except Exception:
            continue
        if len(z) < 0.9 * exp or np.std(z) == 0 or not np.all(np.isfinite(z)):
            continue
        z = z[:exp]
        h1 = h2 = None
        if trH1 is not None and trH2 is not None:
            try:
                aa = trH1.slice(a - PRE, a + POST).data.astype(float); bb = trH2.slice(a - PRE, a + POST).data.astype(float)
                if len(aa) >= 0.9 * exp and len(bb) >= 0.9 * exp:
                    n = min(len(z), len(aa), len(bb)); h1, h2, z = aa[:n], bb[:n], z[:n]
            except Exception:
                h1 = h2 = None
        try:
            f = compute_all(z, fs, PRE, h1, h2)
        except Exception:
            continue
        if any(c not in f or not np.isfinite(f[c]) for c in cols):
            continue
        out.append((idx, [f[c] for c in cols]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--cand", required=True); ap.add_argument("--y0", type=int, default=2010); ap.add_argument("--y1", type=int, default=2013)
    ap.add_argument("--thr", type=float, default=0.4); ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--target-fs", type=float, default=None, help="resample each day to this fs before feature "
                    "extraction (set 40 for the broadband picker; omit to keep native, e.g. boreholes)")
    a = ap.parse_args(); s = a.sta.lower()
    B = joblib.load(f"{HERE}/models/tremor_picker_{s}.joblib")
    model, scaler, cols, classes = B["model"], B["scaler"], B["cols"], B["classes"]
    iL = classes.index("LFE")

    c = pd.read_parquet(a.cand); c["t"] = pd.to_datetime(c["time"]);
    c = c[(c.t.dt.year >= a.y0) & (c.t.dt.year <= a.y1)].reset_index(drop=True)
    c["epoch"] = c.t.values.astype("datetime64[ns]").astype("int64") / 1e9
    c["yr"] = c.t.dt.year; c["jd"] = c.t.dt.dayofyear
    print(f"[{a.sta}] scoring {len(c)} candidates {a.y0}-{a.y1}...", flush=True)
    jobs = defaultdict(list)
    for i, r in enumerate(c.itertuples(index=False)):
        jobs[(int(r.yr), int(r.jd))].append((i, float(r.epoch)))
    payloads = [(a.net, a.sta, k[0], k[1], v, cols, a.target_fs) for k, v in jobs.items()]
    idxs = []; feats = []; done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(cut_day, p) for p in payloads]):
            for idx, ff in f.result():
                idxs.append(idx); feats.append(ff)
            done += 1
            if done % 300 == 0:
                print(f"  {done}/{len(payloads)} days, {len(idxs):,} scored", flush=True)
    idxs = np.array(idxs); F = np.array(feats)
    pL = model.predict_proba(scaler.transform(F))[:, iL]
    c["p_lfe"] = np.nan; c.loc[idxs, "p_lfe"] = pL
    scored = c.dropna(subset=["p_lfe"]).copy()
    keep_cols = ["time", "lat", "lon", "snr", "station", "p_lfe"]
    base = scored[keep_cols]
    filt = scored[scored.p_lfe >= a.thr][keep_cols]
    base.to_parquet(f"data/{s}_cand_baseline.parquet", index=False)
    filt.to_parquet(f"data/{s}_cand_filtered.parquet", index=False)
    print(f"[{a.sta}] scored {len(scored)} | BASELINE {len(base)} -> FILTERED {len(filt)} "
          f"(kept {100*len(filt)/len(base):.0f}% at P(LFE)>={a.thr})")
    print(f"  wrote data/{s}_cand_baseline.parquet + data/{s}_cand_filtered.parquet")
    print("SCORE CANDIDATES DONE")


if __name__ == "__main__":
    main()
