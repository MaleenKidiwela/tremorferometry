"""Extract the SAME Stage-1 features at ground-truth windows for a Lin-region station:

  LIN  = curated Lin (2023) LFE origin times (N>=6, residual<0.8, depth 25-45 km,
         within RADIUS km of the station) -> the POSITIVE reference (real LFEs)
  RAND = random times in the station's record -> NEGATIVE baseline (mostly noise)

Window anchored like the densified windows: S ~ OT + ANCHOR, window [anchor-PRE, anchor+POST]
so feature geometry matches feat_<sta>.parquet. Output: feat_ref_<sta>.parquet (label LIN/RAND).

Usage: PYTHONPATH=src python lfe_features/extract_lin_positives.py --net CN --sta PGC \
         --slat 48.6498 --slon -123.4521 --n-lin 8000 --n-rand 8000 --workers 16
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all, FEATURE_ORDER

PRE, POST, ANCHOR = 10.0, 30.0, 11.0      # S ~ OT+ANCHOR; window [anchor-PRE, anchor+POST]
HERE = os.path.dirname(__file__)


def day_file(net, sta, year, julday):
    return f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"


def process_day(payload):
    net, sta, year, julday, items = payload      # items: list of (label, epoch_seconds)
    from obspy import read, UTCDateTime
    path = day_file(net, sta, year, julday)
    if not os.path.exists(path):
        return []
    try:
        st = read(path)
    except Exception:
        return []

    def get(*chars):
        for ch in chars:
            s = st.select(channel=f"*{ch}")
            if len(s):
                try:
                    s = s.merge(method=1, fill_value=0)
                except Exception:
                    pass
                return s[0]
        return None

    trZ = get("Z")
    if trZ is None:
        return []
    fs = float(trZ.stats.sampling_rate)
    trH1 = get("1", "N"); trH2 = get("2", "E")
    exp = int(round((PRE + POST) * fs))
    rows = []
    for label, t_epoch in items:
        # window anchored at S ~ t_epoch + ANCHOR (Lin OT) ; RAND uses t_epoch directly + ANCHOR too (irrelevant)
        anchor = UTCDateTime(t_epoch) + ANCHOR
        try:
            z = trZ.slice(anchor - PRE, anchor + POST).data.astype(float)
        except Exception:
            continue
        if len(z) < 0.9 * exp or not np.all(np.isfinite(z)) or np.std(z) == 0:
            continue
        z = z[:exp]
        h1 = h2 = None
        if trH1 is not None and trH2 is not None:
            try:
                a = trH1.slice(anchor - PRE, anchor + POST).data.astype(float)
                b = trH2.slice(anchor - PRE, anchor + POST).data.astype(float)
                if len(a) >= 0.9 * exp and len(b) >= 0.9 * exp:
                    n = min(len(z), len(a), len(b)); h1, h2, z = a[:n], b[:n], z[:n]
            except Exception:
                h1 = h2 = None
        try:
            feats = compute_all(z, fs, PRE, h1, h2)
        except Exception:
            continue
        feats.update(label=label, time=float(t_epoch))
        rows.append(feats)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--radius", type=float, default=35.0)
    ap.add_argument("--min-n", type=int, default=6, help="min Lin station count (quality)")
    ap.add_argument("--max-residual", type=float, default=1.2, help="max Lin location residual (s)")
    ap.add_argument("--depth-min", type=float, default=None, help="optional Lin depth lower bound (negative km)")
    ap.add_argument("--depth-max", type=float, default=None, help="optional Lin depth upper bound (negative km)")
    ap.add_argument("--n-lin", type=int, default=8000); ap.add_argument("--n-rand", type=int, default=8000)
    ap.add_argument("--workers", type=int, default=16); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    s = a.sta.lower(); rng = np.random.RandomState(a.seed)

    # which day-files exist for this station (defines the usable time span)
    files = glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed")
    have_days = set()
    for f in files:
        yr = int(f.split("/")[-2]); jd = int(os.path.basename(f).split(".")[0]); have_days.add((yr, jd))
    if not have_days:
        sys.exit(f"no waveforms for {a.net}.{a.sta}")
    print(f"[{a.sta}] {len(have_days)} day-files on disk", flush=True)

    # ---- curated Lin positives within radius ----
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dkm = np.sqrt(((lin.lat - a.slat) * 111.0) ** 2 +
                  ((lin.lon - a.slon) * 111.0 * math.cos(math.radians(a.slat))) ** 2)
    q = (lin.N >= a.min_n) & (lin.residual < a.max_residual) & (dkm < a.radius)
    if a.depth_min is not None:
        q &= (lin.depth >= a.depth_min)
    if a.depth_max is not None:
        q &= (lin.depth <= a.depth_max)
    cur = lin[q].copy()
    cur["t"] = pd.to_datetime(cur.OT, utc=True)
    cur["epoch"] = cur.t.values.astype("datetime64[ns]").astype("int64") / 1e9
    cur["year"] = cur.t.dt.year; cur["julday"] = cur.t.dt.dayofyear
    cur = cur[[ (r.year, r.julday) in have_days for r in cur.itertuples(index=False)]]
    if len(cur) > a.n_lin:
        cur = cur.sample(a.n_lin, random_state=a.seed)
    print(f"[{a.sta}] curated Lin positives in-window: {len(cur):,} "
          f"(of {((lin.N>=6)&(lin.residual<0.8)&(dkm<a.radius)).sum():,} curated within radius)", flush=True)

    # ---- random-time negatives on existing days ----
    daylist = sorted(have_days)
    rand_items = []
    for _ in range(a.n_rand):
        yr, jd = daylist[rng.randint(len(daylist))]
        from datetime import datetime
        base = datetime(yr, 1, 1).timestamp() + (jd - 1) * 86400 + rng.uniform(60, 86340)
        rand_items.append((yr, jd, base))

    # ---- build per-day jobs ----
    jobs = defaultdict(list)
    for r in cur.itertuples(index=False):
        jobs[(int(r.year), int(r.julday))].append(("LIN", float(r.epoch)))
    for yr, jd, base in rand_items:
        jobs[(yr, jd)].append(("RAND", float(base)))
    payloads = [(a.net, a.sta, k[0], k[1], v) for k, v in jobs.items()]
    print(f"[{a.sta}] extracting {sum(len(v) for v in jobs.values()):,} ref windows over {len(payloads)} days...", flush=True)

    rows = []; done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(process_day, p) for p in payloads]
        for f in as_completed(futs):
            rows.extend(f.result()); done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(payloads)} days, {len(rows):,} windows", flush=True)
    df = pd.DataFrame(rows)
    os.makedirs(f"{HERE}/data", exist_ok=True)
    out = f"{HERE}/data/feat_ref_{s}.parquet"
    df.to_parquet(out, index=False)
    print(f"[{a.sta}] wrote {out}: " + ", ".join(f"{k}={v}" for k, v in df.label.value_counts().items()) +
          f"  3-comp={'yes' if 'hv_ratio' in df.columns else 'no'}", flush=True)
    print("REF EXTRACT DONE", flush=True)


if __name__ == "__main__":
    main()
