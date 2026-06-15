"""Stage-3 prep: cut windows for LIN / EQ / RAND / FAM and save fixed-size log-spectrograms
for the autoencoder. Frequency axis fixed in Hz (0-25 Hz) so it is instrument-comparable.

Output: data/spec_<sta>.npz  with X (N, F, T) log-spectrograms, y (labels), fam (family id or '').

Usage: PYTHONPATH=src python lfe_features/extract_spectrograms.py --net PB --sta B011 \
         --slat 48.65 --slon -123.448 --per 2500 --workers 14
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy import signal
sys.path.insert(0, os.path.dirname(__file__))

PRE, POST = 10.0, 30.0
FMAX, NF, NT = 25.0, 40, 64        # fixed spectrogram grid: 40 freq bins (0-25Hz) x 64 time frames
HERE = os.path.dirname(__file__)


def logspec(z, fs):
    f, t, S = signal.spectrogram(signal.detrend(z), fs=fs, nperseg=int(fs * 1.5),
                                 noverlap=int(fs * 1.5 * 0.8))
    S = np.log10(S + 1e-12)
    fm = f <= FMAX
    S = S[fm]                                   # crop to 0-25 Hz
    # resample to NF x NT via linear interp
    fi = np.linspace(0, S.shape[0] - 1, NF); ti = np.linspace(0, S.shape[1] - 1, NT)
    from scipy.ndimage import map_coordinates
    gy, gx = np.meshgrid(fi, ti, indexing="ij")
    R = map_coordinates(S, [gy, gx], order=1, mode="nearest")
    R = (R - R.mean()) / (R.std() + 1e-9)       # per-window standardize
    return R.astype(np.float32)


def cut_day(payload):
    net, sta, year, julday, items = payload
    from obspy import read, UTCDateTime
    path = f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"
    if not os.path.exists(path):
        return []
    try:
        st = read(path).select(channel="*Z")
        st = st.merge(method=1, fill_value=0)
    except Exception:
        return []
    if len(st) == 0:
        return []
    tr = st[0]; fs = float(tr.stats.sampling_rate); exp = int(round((PRE + POST) * fs))
    out = []
    for label, fam, anchor in items:
        try:
            z = tr.slice(UTCDateTime(anchor) - PRE, UTCDateTime(anchor) + POST).data.astype(float)
        except Exception:
            continue
        if len(z) < exp or not np.all(np.isfinite(z)) or np.std(z) == 0:
            continue
        try:
            out.append((label, fam, logspec(z[:exp], fs)))
        except Exception:
            continue
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--per", type=int, default=2500); ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(); s = a.sta.lower(); rng = np.random.RandomState(a.seed)
    ANCHOR = 11.0
    files = glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed")
    have = set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0])) for f in files)
    daylist = sorted(have)
    jobs = defaultdict(list)

    # LIN
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dkm = np.sqrt(((lin.lat - a.slat) * 111) ** 2 + ((lin.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    cur = lin[(lin.N >= 6) & (lin.residual < 1.2) & (dkm < 35)].copy()
    cur["t"] = pd.to_datetime(cur.OT, utc=True); cur["y"] = cur.t.dt.year; cur["j"] = cur.t.dt.dayofyear
    cur = cur[[ (r.y, r.j) in have for r in cur.itertuples(index=False)]].sample(min(a.per, len(cur)), random_state=1)
    for r in cur.itertuples(index=False):
        jobs[(int(r.y), int(r.j))].append(("LIN", "", r.t.value / 1e9 + ANCHOR))
    # EQ
    if os.path.exists(f"data/eq_catalog_{s}.csv"):
        ev = pd.read_csv(f"data/eq_catalog_{s}.csv"); ev["t"] = pd.to_datetime(ev.ot, utc=True)
        epi = np.sqrt(((ev.lat - a.slat) * 111) ** 2 + ((ev.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
        ev["anch"] = ev.t.values.astype("datetime64[ns]").astype("int64") / 1e9 + np.sqrt(epi ** 2 + np.clip(ev.depth, 0, 60) ** 2) / 3.6
        ev["adt"] = pd.to_datetime(ev.anch, unit="s", utc=True); ev["y"] = ev.adt.dt.year; ev["j"] = ev.adt.dt.dayofyear
        ev = ev[[ (r.y, r.j) in have for r in ev.itertuples(index=False)]].sample(min(a.per, len(ev)), random_state=1)
        for r in ev.itertuples(index=False):
            jobs[(int(r.y), int(r.j))].append(("EQ", "", float(r.anch)))
    # RAND
    from datetime import datetime
    for _ in range(a.per):
        y, jd = daylist[rng.randint(len(daylist))]
        base = datetime(y, 1, 1).timestamp() + (jd - 1) * 86400 + rng.uniform(60, 86340)
        jobs[(y, jd)].append(("RAND", "", base + ANCHOR))
    # FAM (sampled densified detections, with grade in fam string)
    A = pd.read_csv("results/family_verdicts/ALL_families_labeled.csv"); A = A[A.station == a.sta]
    grade = dict(zip(A.fam, A.grade))
    m = pd.read_csv(f"data/mf_{s}_all.csv", usecols=["template", "time"])
    m["t"] = pd.to_datetime(m.time, utc=True); m["y"] = m.t.dt.year; m["j"] = m.t.dt.dayofyear
    m = m[[ (r.y, r.j) in have for r in m.itertuples(index=False)]]
    parts = [g.sample(min(40, len(g)), random_state=1) for _, g in m.groupby("template")]
    ms = pd.concat(parts, ignore_index=True)
    for r in ms.itertuples(index=False):
        jobs[(int(r.y), int(r.j))].append((f"FAM_{grade.get(r.template,'NA')}", r.template, r.t.value / 1e9))

    payloads = [(a.net, a.sta, k[0], k[1], v) for k, v in jobs.items()]
    print(f"[{a.sta}] spectrograms over {len(payloads)} days...", flush=True)
    X = []; Y = []; F = []; done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(cut_day, p) for p in payloads]):
            for label, fam, R in f.result():
                X.append(R); Y.append(label); F.append(fam)
            done += 1
            if done % 800 == 0:
                print(f"  {done}/{len(payloads)} days, {len(X):,} specs", flush=True)
    X = np.stack(X); Y = np.array(Y); F = np.array(F)
    np.savez_compressed(f"{HERE}/data/spec_{s}.npz", X=X, y=Y, fam=F)
    import collections
    print(f"[{a.sta}] saved data/spec_{s}.npz  X={X.shape}  classes={dict(collections.Counter(Y))}")
    print("SPEC EXTRACT DONE")


if __name__ == "__main__":
    main()
