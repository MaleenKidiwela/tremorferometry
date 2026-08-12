"""Extract the 'in-tremor NOISE/cultural' negative class: random times INSIDE PNSN tremor
windows that are NOT near a curated Lin LFE. These represent the ambient/cultural energy the
discovery envelope-peak detector grabs DURING tremor (the real confuser, not quiet-time noise).

Output: feat_tnoise_<sta>.parquet (label TNOISE).

Usage: PYTHONPATH=src python lfe_features/extract_tremor_noise.py --net PB --sta B011 \
         --slat 48.65 --slon -123.448 --n 8000 --workers 14
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(__file__))
from extract_eq_negatives import process_day            # reuse the windowed feature extractor
HERE = os.path.dirname(__file__); ANCHOR = 11.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--radius", type=float, default=30.0); ap.add_argument("--winlen", type=float, default=300.0)
    ap.add_argument("--n", type=int, default=8000); ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(); s = a.sta.lower(); rng = np.random.RandomState(a.seed)
    have = set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0]))
               for f in glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed"))

    # PNSN tremor windows near station, on days we have
    tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    dk = np.sqrt(((tr.lat - a.slat) * 111) ** 2 + ((tr.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    tr = tr[dk < a.radius].copy(); tr["t"] = pd.to_datetime(tr.time, utc=True)
    tr["epoch"] = tr.t.values.astype("datetime64[ns]").astype("int64") / 1e9
    tr["y"] = tr.t.dt.year; tr["j"] = tr.t.dt.dayofyear
    tr = tr[[ (r.y, r.j) in have for r in tr.itertuples(index=False)]]
    print(f"[{a.sta}] {len(tr)} PNSN tremor windows on-disk", flush=True)

    # curated Lin LFE arrival times (to AVOID) — sorted
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dl = np.sqrt(((lin.lat - a.slat) * 111) ** 2 + ((lin.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    cur = lin[(lin.N >= 6) & (lin.residual < 1.2) & (dl < 35)]
    lfe_t = np.sort(pd.to_datetime(cur.OT, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9 + ANCHOR)

    # draw random times inside tremor windows, >20s from any LFE
    trv = tr.epoch.values
    jobs = defaultdict(list); ntries = 0
    while sum(len(v) for v in jobs.values()) < a.n and ntries < a.n * 20:
        ntries += 1
        e = trv[rng.randint(len(trv))] + rng.uniform(0, a.winlen)
        k = np.searchsorted(lfe_t, e)
        near = min(abs(e - lfe_t[min(k, len(lfe_t) - 1)]), abs(e - lfe_t[max(k - 1, 0)]))
        if near < 20:
            continue
        from datetime import datetime, timezone
        d = datetime.fromtimestamp(e, tz=timezone.utc)
        if (d.year, d.timetuple().tm_yday) in have:
            jobs[(d.year, d.timetuple().tm_yday)].append(("TNOISE", e))
    payloads = [(a.net, a.sta, k[0], k[1], v) for k, v in jobs.items()]
    print(f"[{a.sta}] extracting {sum(len(v) for v in jobs.values())} in-tremor noise windows...", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(process_day, p) for p in payloads]):
            rows.extend(f.result())
    df = pd.DataFrame(rows)
    df.to_parquet(f"{HERE}/data/feat_tnoise_{s}.parquet", index=False)
    print(f"[{a.sta}] wrote feat_tnoise_{s}.parquet: TNOISE={len(df)}  3-comp={'yes' if 'hv_ratio' in df.columns else 'no'}")
    print("TNOISE EXTRACT DONE")


if __name__ == "__main__":
    main()
