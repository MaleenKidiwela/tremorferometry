"""PROOF: does the tremor-picker make discovery cleaner?

Replicates the discovery candidate step — envelope-peak detection INSIDE PNSN tremor windows —
then clusters candidates into families (cc>=0.8, complete-linkage, >=3 members, >=2 years) TWO ways:
  BASELINE  = all envelope-peak candidates (current discovery)
  FILTERED  = only candidates with tremor-picker P(LFE) >= thr
Each resulting family is classified by its members' ground truth (LFE / EQ / BLAST / UNKNOWN) using
Lin LFE times + ANSS earthquake/explosion times. Shows the filter removes spurious EQ/blast families.

Usage: PYTHONPATH=src python lfe_features/discovery_filter_demo.py --net PB --sta B011 --slat 48.65 --slon -123.448
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy import signal
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import joblib
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all, FEATURE_ORDER

HERE = os.path.dirname(__file__)
FPRE, FPOST = 10.0, 30.0     # feature window
TPRE, TPOST = 1.0, 1.0       # 2-s template window for clustering
WINLEN = 300.0               # tremor-window length (s)
POL = ["hv_ratio", "rectilinearity", "planarity"]


def detect_day(payload):
    """Envelope-peak candidates inside the day's tremor windows; return (time, template2s, features)."""
    net, sta, year, julday, tremor_epochs, cols = payload
    from obspy import read, UTCDateTime
    path = f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"
    if not os.path.exists(path):
        return []
    try:
        st = read(path)
    except Exception:
        return []

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
    fs = float(trZ.stats.sampling_rate); day0 = UTCDateTime(year=year, julday=julday).timestamp
    Z = trZ.data.astype(float); L = len(Z)
    sos = signal.butter(4, [2, 8], "band", fs=fs, output="sos")
    H1 = trH1.data.astype(float) if trH1 is not None else None
    H2 = trH2.data.astype(float) if trH2 is not None else None
    out = []
    fexp = int((FPRE + FPOST) * fs); texp = int((TPRE + TPOST) * fs)
    for te in tremor_epochs:
        i0 = int((te - day0) * fs); i1 = i0 + int(WINLEN * fs)
        if i0 < 0 or i1 > L:
            continue
        seg = Z[i0:i1]
        if np.std(seg) == 0:
            continue
        env = np.abs(signal.hilbert(signal.sosfiltfilt(sos, signal.detrend(seg))))
        bg = np.median(env) + 1e-9
        pk, _ = signal.find_peaks(env, height=3 * bg, distance=int(1.0 * fs))
        for p in pk:
            gp = i0 + p; tt = day0 + gp / fs
            ti0 = gp - int(TPRE * fs); fi0 = gp - int(FPRE * fs)
            if ti0 < 0 or ti0 + texp > L or fi0 < 0 or fi0 + fexp > L:
                continue
            tmpl = Z[ti0:ti0 + texp].copy()
            if np.std(tmpl) == 0:
                continue
            tmpl = (tmpl - tmpl.mean()) / (np.linalg.norm(tmpl) + 1e-9)
            zf = Z[fi0:fi0 + fexp]
            h1 = H1[fi0:fi0 + fexp] if H1 is not None else None
            h2 = H2[fi0:fi0 + fexp] if H2 is not None else None
            try:
                feats = compute_all(zf, fs, FPRE, h1, h2)
            except Exception:
                continue
            if any(c not in feats or not np.isfinite(feats[c]) for c in cols):
                continue
            out.append((tt, tmpl.astype(np.float32), [feats[c] for c in cols]))
    return out


def cc_matrix(T, max_shift=20):
    """max-shift normalized cross-correlation distance matrix (condensed)."""
    n = len(T)
    best = np.full((n, n), -1.0)
    for sh in range(-max_shift, max_shift + 1, 2):
        Ts = np.roll(T, sh, axis=1)
        c = Ts @ T.T
        best = np.maximum(best, c)
    np.fill_diagonal(best, 1.0)
    best = np.clip(best, -1, 1)
    return 1.0 - best


def cluster(idx, T, times, min_members=3, min_years=2, cc=0.8):
    if len(idx) < min_members:
        return []
    D = cc_matrix(T[idx])
    Z = linkage(squareform(D, checks=False), method="complete")
    lab = fcluster(Z, t=1 - cc, criterion="distance")
    fams = []
    yr = np.array([pd.Timestamp(t, unit="s", tz="UTC").year for t in times[idx]])
    for c in np.unique(lab):
        m = np.where(lab == c)[0]
        if len(m) >= min_members and len(set(yr[m])) >= min_years:
            fams.append(idx[m])
    return fams


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--n-windows", type=int, default=900); ap.add_argument("--max-cand", type=int, default=3000)
    ap.add_argument("--thr", type=float, default=0.4); ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--min-members", type=int, default=3); ap.add_argument("--min-years", type=int, default=2)
    a = ap.parse_args(); s = a.sta.lower(); rng = np.random.RandomState(1)
    B = joblib.load(f"{HERE}/models/tremor_picker_{s}.joblib")
    model, scaler, cols, classes = B["model"], B["scaler"], B["cols"], B["classes"]
    iL = classes.index("LFE")

    have = set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0]))
               for f in glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed"))
    tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    dk = np.sqrt(((tr.lat - a.slat) * 111) ** 2 + ((tr.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    tr = tr[dk < 30].copy(); tr["t"] = pd.to_datetime(tr.time, utc=True)
    tr["e"] = tr.t.values.astype("datetime64[ns]").astype("int64") / 1e9
    tr["y"] = tr.t.dt.year; tr["j"] = tr.t.dt.dayofyear
    tr = tr[[ (r.y, r.j) in have for r in tr.itertuples(index=False)]]
    tr = tr.sample(min(a.n_windows, len(tr)), random_state=1)       # spread across episodes/years
    jobs = defaultdict(list)
    for r in tr.itertuples(index=False):
        jobs[(int(r.y), int(r.j))].append(float(r.e))
    payloads = [(a.net, a.sta, k[0], k[1], v, cols) for k, v in jobs.items()]
    print(f"[{a.sta}] running envelope-peak detection in {len(tr)} tremor windows over {len(payloads)} days...", flush=True)

    times = []; T = []; F = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(detect_day, p) for p in payloads]):
            for tt, tmpl, feats in f.result():
                times.append(tt); T.append(tmpl); F.append(feats)
    times = np.array(times); T = np.stack(T); F = np.array(F)
    if len(times) > a.max_cand:
        sel = rng.choice(len(times), a.max_cand, replace=False)
        times, T, F = times[sel], T[sel], F[sel]
    print(f"[{a.sta}] {len(times)} candidates (envelope-peaks in tremor windows)", flush=True)

    # ground truth labels
    def load_t(path, col, anchor=0.0, eqstyle=False):
        d = pd.read_csv(path)
        tt = pd.to_datetime(d[col], utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9
        return np.sort(tt + anchor)
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dl = np.sqrt(((lin.lat - a.slat) * 111) ** 2 + ((lin.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    lin_t = np.sort(pd.to_datetime(lin[(lin.N >= 6) & (lin.residual < 1.2) & (dl < 35)].OT, utc=True)
                    .values.astype("datetime64[ns]").astype("int64") / 1e9 + 11.0)
    def eq_arr(csv):
        d = pd.read_csv(csv); ep = np.sqrt(((d.lat - a.slat) * 111) ** 2 + ((d.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
        return np.sort(pd.to_datetime(d.ot, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9 +
                       np.sqrt(ep ** 2 + np.clip(d.depth, 0, 60) ** 2) / 3.6)
    eq_t = eq_arr(f"data/eq_catalog_{s}.csv"); bl_t = eq_arr(f"data/blast_catalog_{s}.csv")

    def nearest(arr, t, tol=6.0):
        if len(arr) == 0: return False
        k = np.searchsorted(arr, t); k = min(max(k, 1), len(arr) - 1)
        return min(abs(t - arr[k]), abs(t - arr[k - 1])) <= tol
    gt = np.array(["LFE" if nearest(lin_t, t) else "EQ" if nearest(eq_t, t) else "BLAST" if nearest(bl_t, t)
                   else "UNK" for t in times])
    print("  candidate ground-truth mix:", dict(Counter(gt)), flush=True)

    pL = model.predict_proba(scaler.transform(F))[:, iL]

    def famtype(members):
        c = Counter(gt[members])
        for k in ("LFE", "EQ", "BLAST"):
            if c[k] >= max(3, 0.3 * len(members)):  # dominant known label
                return k
        return "UNK"

    for tag, idx in [("BASELINE (envelope-only)", np.arange(len(times))),
                     (f"FILTERED (P(LFE)>={a.thr})", np.where(pL >= a.thr)[0])]:
        fams = cluster(idx, T, times, min_members=a.min_members, min_years=a.min_years)
        types = Counter(famtype(m) for m in fams)
        memb = sum(len(m) for m in fams)
        print(f"\n=== {tag}: {len(idx)} candidates -> {len(fams)} families ({memb} clustered members) ===")
        print("   family types:", dict(types))
    print("\nDISCOVERY FILTER DEMO DONE")


if __name__ == "__main__":
    main()
