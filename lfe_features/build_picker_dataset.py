"""Stage-4 prep: build labeled continuous SEGMENTS for a PhaseNet/EQTransformer-style
per-sample LFE picker.

Each segment = 3-comp, ~41 s, 50 Hz (2048 samples), bandpass 1-15 Hz, per-trace normalized.
Per-sample target = 3 softmax channels [LFE, EQ, NOISE]: Gaussian bumps (sigma 0.4 s) at every
curated Lin LFE S-arrival and every ANSS EQ S-arrival inside the segment; NOISE = 1 - LFE - EQ.

Segments are drawn LFE-centered / EQ-centered / noise so all three classes are represented, and
ALL catalog events in each window are labeled (so overlapping events aren't mislabeled noise).

Output: data/picker_ds_<sta>.npz  X (N,3,2048) float32, Y (N,3,2048) float32.

Usage: PYTHONPATH=src python lfe_features/build_picker_dataset.py --net PB --sta B011 \
         --slat 48.65 --slon -123.448 --n 5000 --workers 14
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy import signal
sys.path.insert(0, os.path.dirname(__file__))

FS = 50.0; SEG = 2048; T = SEG / FS          # ~40.96 s
SIG_LFE = 1.0     # LFE label width (s) — narrowed for CC-refined picks (was 1.5 for raw OT+11 anchor)
SIG_EQ = 0.6      # EQ label width (s) — impulsive, better located
ANCHOR = 11.0; VS = 3.6
HERE = os.path.dirname(__file__)


PAD = 2.0   # filter edge pad (s)


def _prep_slice(x, fs_raw):
    """Filter+resample a SHORT raw slice (fast); caller trims the pad."""
    x = signal.detrend(x)
    sos = signal.butter(4, [1, 15], "band", fs=fs_raw, output="sos")
    x = signal.sosfiltfilt(sos, x)
    if abs(fs_raw - FS) > 0.1:
        g = math.gcd(int(round(fs_raw)), int(FS))
        x = signal.resample_poly(x, int(FS) // g, int(round(fs_raw)) // g)
    return x


def build_day(payload):
    net, sta, year, julday, segs, lin_t, eq_t = payload
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
    if trZ is None or trH1 is None or trH2 is None:
        return []
    fs_raw = float(trZ.stats.sampling_rate)
    raws = [trZ.data.astype(float), trH1.data.astype(float), trH2.data.astype(float)]
    L = min(len(r) for r in raws)
    day0 = UTCDateTime(year=year, julday=julday).timestamp
    n_pad = int(round(PAD * fs_raw)); n_raw = int(round((T + 2 * PAD) * fs_raw))
    start = int(round(PAD * FS))
    out = []
    tt = np.arange(SEG) / FS
    for t0 in segs:
        i0 = int(round((t0 - day0 - PAD) * fs_raw))
        if i0 < 0 or i0 + n_raw > L:
            continue
        chans = []
        ok = True
        for r in raws:
            sl = r[i0:i0 + n_raw]
            if np.std(sl) == 0 or not np.all(np.isfinite(sl)):
                ok = False; break
            p = _prep_slice(sl, fs_raw)[start:start + SEG]
            if len(p) < SEG:
                ok = False; break
            chans.append(p)
        if not ok:
            continue
        X = np.stack(chans)
        X = X / (np.abs(X).max() + 1e-9)
        Y = np.zeros((3, SEG), np.float32)
        for arr, ch, sg in ((lin_t, 0, SIG_LFE), (eq_t, 1, SIG_EQ)):
            lo, hi = np.searchsorted(arr, [t0, t0 + T])
            for et in arr[lo:hi]:
                Y[ch] += np.exp(-0.5 * ((tt - (et - t0)) / sg) ** 2)
        Y[:2] = np.clip(Y[:2], 0, 1)
        Y[2] = np.clip(1 - Y[0] - Y[1], 0, 1)
        out.append((X.astype(np.float32), Y))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--n", type=int, default=5000); ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(); s = a.sta.lower(); rng = np.random.RandomState(a.seed)
    files = glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed")
    have = sorted(set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0])) for f in files))
    day0 = {(y, j): pd.Timestamp(f"{y}-01-01", tz="UTC").timestamp() + (j - 1) * 86400 for y, j in have}

    # LFE S-arrival times — use cross-correlation-REFINED picks if available (else OT+ANCHOR)
    pp = f"{HERE}/data/lin_precise_{s}.csv"
    if os.path.exists(pp):
        lin_t = np.sort(pd.read_csv(pp)["arrival"].values.astype(float))
        print(f"[{a.sta}] using {len(lin_t)} cross-correlation-refined LFE picks (SIG_LFE={SIG_LFE})", flush=True)
    else:
        lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
        dk = np.sqrt(((lin.lat - a.slat) * 111) ** 2 + ((lin.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
        cur = lin[(lin.N >= 6) & (lin.residual < 1.2) & (dk < 35)]
        lin_t = np.sort(pd.to_datetime(cur.OT, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9 + ANCHOR)
    # EQ S-arrival times
    ev = pd.read_csv(f"data/eq_catalog_{s}.csv")
    epi = np.sqrt(((ev.lat - a.slat) * 111) ** 2 + ((ev.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    eq_t = np.sort(pd.to_datetime(ev.ot, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9 +
                   np.sqrt(epi ** 2 + np.clip(ev.depth, 0, 60) ** 2) / VS)

    def day_ok(epoch):
        from datetime import datetime, timezone
        d = datetime.fromtimestamp(epoch, tz=timezone.utc); return (d.year, d.timetuple().tm_yday)

    # draw segment starts: 45% LFE-centered, 25% EQ-centered, 30% noise
    segs_by_day = defaultdict(list)
    def place(center_times, n):
        if len(center_times) == 0: return
        for c in rng.choice(center_times, min(n, len(center_times) * 3), replace=True):
            t0 = c - rng.uniform(8, T - 8)
            k = day_ok(t0)
            if k in day0: segs_by_day[k].append(t0)
    place(lin_t, int(a.n * 0.45)); place(eq_t, int(a.n * 0.25))
    for _ in range(int(a.n * 0.30)):
        y, j = have[rng.randint(len(have))]
        segs_by_day[(y, j)].append(day0[(y, j)] + rng.uniform(60, 86340 - T))

    payloads = [(a.net, a.sta, k[0], k[1], v, lin_t, eq_t) for k, v in segs_by_day.items()]
    print(f"[{a.sta}] building ~{sum(len(v) for v in segs_by_day.values())} segments over {len(payloads)} days "
          f"(LFE arrivals {len(lin_t)}, EQ {len(eq_t)})", flush=True)
    XS, YS = [], []; done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(build_day, p) for p in payloads]):
            for X, Y in f.result():
                XS.append(X); YS.append(Y)
            done += 1
            if done % 400 == 0:
                print(f"  {done}/{len(payloads)} days, {len(XS):,} segs", flush=True)
    X = np.stack(XS); Y = np.stack(YS)
    np.savez_compressed(f"{HERE}/data/picker_ds_{s}.npz", X=X, Y=Y)
    frac_lfe = (Y[:, 0].max(1) > 0.5).mean(); frac_eq = (Y[:, 1].max(1) > 0.5).mean()
    print(f"[{a.sta}] saved data/picker_ds_{s}.npz X={X.shape}  segs-with-LFE {frac_lfe:.2f}  with-EQ {frac_eq:.2f}")
    print("PICKER DATASET DONE")


if __name__ == "__main__":
    main()
