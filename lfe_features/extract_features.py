"""Stage-1 feature extraction over densified detection windows for one station.

For each family in mf_<sta>_all.csv, sample detections, cut a long window
[t_d-PRE, t_d+POST] from the daily mseed, compute hand-crafted features
(feature_defs.compute_all), and write per-detection + per-family tables.

Outputs:
  lfe_features/data/feat_<sta>.parquet         per-detection features (+ fam, grade, cc, time)
  lfe_features/data/feat_<sta>_family.csv       per-family median/std fingerprint

Usage:
  PYTHONPATH=src python lfe_features/extract_features.py --net PB --sta B927 \
      --per-fam 150 --workers 24
"""
import argparse, os, sys, glob
import numpy as np
import pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all, FEATURE_ORDER

PRE, POST = 10.0, 30.0
HERE = os.path.dirname(__file__)


def day_file(net, sta, year, julday):
    return f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"


def process_day(payload):
    net, sta, year, julday, dets = payload   # dets: list of (fam, epoch_seconds, cc)
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
    trH1 = get("1", "N")
    trH2 = get("2", "E")
    exp = int(round((PRE + POST) * fs))
    rows = []
    for fam, t_epoch, cc in dets:
        tt = UTCDateTime(t_epoch)
        try:
            z = trZ.slice(tt - PRE, tt + POST).data.astype(float)
        except Exception:
            continue
        if len(z) < 0.9 * exp or not np.all(np.isfinite(z)) or np.std(z) == 0:
            continue
        z = z[:exp] if len(z) > exp else z
        h1 = h2 = None
        if trH1 is not None and trH2 is not None:
            try:
                a = trH1.slice(tt - PRE, tt + POST).data.astype(float)
                b = trH2.slice(tt - PRE, tt + POST).data.astype(float)
                if len(a) >= 0.9 * exp and len(b) >= 0.9 * exp:
                    n = min(len(z), len(a), len(b))
                    h1, h2 = a[:n], b[:n]
                    z = z[:n]
            except Exception:
                h1 = h2 = None
        try:
            feats = compute_all(z, fs, PRE, h1, h2)
        except Exception:
            continue
        feats.update(fam=fam, time=float(t_epoch), cc=float(cc))
        rows.append(feats)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True)
    ap.add_argument("--sta", required=True)
    ap.add_argument("--per-fam", type=int, default=150)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--max-days", type=int, default=0, help="0 = no cap (demo speed control)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    s = a.sta.lower()

    print(f"[{a.sta}] loading detections + grades...", flush=True)
    m = pd.read_csv(f"data/mf_{s}_all.csv", usecols=["template", "time", "cc"])
    m["t"] = pd.to_datetime(m.time, utc=True)
    A = pd.read_csv("results/family_verdicts/ALL_families_labeled.csv")
    grade = dict(zip(A[A.station == a.sta].fam, A[A.station == a.sta].grade))

    # sample per family: half highest-cc, half random (capture range)
    samp = []
    rng = np.random.RandomState(a.seed)
    for fam, g in m.groupby("template"):
        n = min(a.per_fam, len(g))
        top = g.nlargest(n // 2, "cc")
        rest = g.drop(top.index)
        rnd = rest.sample(min(n - len(top), len(rest)), random_state=a.seed) if len(rest) else rest
        samp.append(pd.concat([top, rnd]))
    samp = pd.concat(samp, ignore_index=True)
    samp["epoch"] = samp.t.values.astype("datetime64[ns]").astype("int64") / 1e9
    print(f"[{a.sta}] {len(samp):,} sampled windows across {samp.template.nunique()} families", flush=True)

    # group by (year, julday)
    samp["year"] = samp.t.dt.year
    samp["julday"] = samp.t.dt.dayofyear
    jobs = defaultdict(list)
    for r in samp.itertuples(index=False):
        jobs[(r.year, r.julday)].append((r.template, r.epoch, r.cc))
    keys = list(jobs.keys())
    if a.max_days and len(keys) > a.max_days:
        keys = list(rng.permutation(keys))[: a.max_days]
        kept = sum(len(jobs[k]) for k in keys)
        print(f"[{a.sta}] max-days cap: {len(keys)} days / {kept:,} windows", flush=True)
    payloads = [(a.net, a.sta, int(k[0]), int(k[1]), jobs[k]) for k in keys]

    print(f"[{a.sta}] extracting over {len(payloads)} day-files (workers={a.workers})...", flush=True)
    rows = []
    done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(process_day, p) for p in payloads]
        for f in as_completed(futs):
            rows.extend(f.result())
            done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(payloads)} days, {len(rows):,} windows", flush=True)

    if not rows:
        sys.exit(f"[{a.sta}] no windows extracted (waveforms present?)")
    df = pd.DataFrame(rows)
    df["grade"] = df.fam.map(grade).fillna("NA")
    os.makedirs(f"{HERE}/data", exist_ok=True)
    out = f"{HERE}/data/feat_{s}.parquet"
    df.to_parquet(out, index=False)
    feat_cols = [c for c in FEATURE_ORDER if c in df.columns]
    fam_agg = df.groupby(["fam", "grade"])[feat_cols].agg(["median", "std", "count"])
    fam_agg.columns = ["_".join(c) for c in fam_agg.columns]
    fam_agg = fam_agg.reset_index()
    fam_agg.to_csv(f"{HERE}/data/feat_{s}_family.csv", index=False)
    print(f"[{a.sta}] wrote {out} ({len(df):,} windows) + feat_{s}_family.csv "
          f"({len(fam_agg)} families). 3-comp={'yes' if 'hv_ratio' in df.columns else 'no (Z only)'}", flush=True)
    print("EXTRACT DONE", flush=True)


if __name__ == "__main__":
    main()
