"""Extract the SAME features at regional EARTHQUAKE arrivals -> hard-negative class 'EQ'.

A real LFE picker must reject impulsive/broadband earthquakes, not just silence. We query
ANSS/USGS regional events near the station, anchor each window at the event's estimated
S-arrival (OT + hypocentral-distance/Vs) so geometry matches the LIN/family windows, and
extract identical features. Output: feat_eq_<sta>.parquet (label EQ). Caches the event
catalog to data/eq_catalog_<sta>.csv.

Usage: PYTHONPATH=src python lfe_features/extract_eq_negatives.py --net PB --sta B011 \
         --slat 48.6500 --slon -123.4480 --maxrad 1.3 --minmag 1.0 --n 8000 --workers 14
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all, FEATURE_ORDER

PRE, POST = 10.0, 30.0      # window [S_anchor-PRE, S_anchor+POST]; S at +10s
VS = 3.6                    # km/s crustal S
HERE = os.path.dirname(__file__)


def day_file(net, sta, year, julday):
    return f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"


def process_day(payload):
    net, sta, year, julday, items = payload     # items: (label, anchor_epoch)
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
                try: s = s.merge(method=1, fill_value=0)
                except Exception: pass
                return s[0]
        return None
    trZ = get("Z")
    if trZ is None:
        return []
    fs = float(trZ.stats.sampling_rate)
    trH1 = get("1", "N"); trH2 = get("2", "E")
    exp = int(round((PRE + POST) * fs))
    rows = []
    for label, anchor in items:
        a = UTCDateTime(anchor)
        try:
            z = trZ.slice(a - PRE, a + POST).data.astype(float)
        except Exception:
            continue
        if len(z) < 0.9 * exp or not np.all(np.isfinite(z)) or np.std(z) == 0:
            continue
        z = z[:exp]
        h1 = h2 = None
        if trH1 is not None and trH2 is not None:
            try:
                aa = trH1.slice(a - PRE, a + POST).data.astype(float)
                bb = trH2.slice(a - PRE, a + POST).data.astype(float)
                if len(aa) >= 0.9 * exp and len(bb) >= 0.9 * exp:
                    n = min(len(z), len(aa), len(bb)); h1, h2, z = aa[:n], bb[:n], z[:n]
            except Exception:
                h1 = h2 = None
        try:
            feats = compute_all(z, fs, PRE, h1, h2)
        except Exception:
            continue
        feats.update(label=label, time=float(anchor))
        rows.append(feats)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--maxrad", type=float, default=1.3); ap.add_argument("--minmag", type=float, default=1.0)
    ap.add_argument("--n", type=int, default=8000); ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--eventtype", default=None, help="USGS eventtype filter, e.g. 'explosion'")
    ap.add_argument("--label", default="EQ", help="class label written to the parquet")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(); s = a.sta.lower(); rng = np.random.RandomState(a.seed)
    lab = a.label.lower()

    files = glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed")
    have = set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0])) for f in files)
    yrs = sorted(set(y for y, _ in have))

    cat_csv = f"data/{lab}_catalog_{s}.csv"
    if os.path.exists(cat_csv):
        ev = pd.read_csv(cat_csv)
    else:
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        kw = dict(starttime=UTCDateTime(f"{yrs[0]}-01-01"), endtime=UTCDateTime(f"{yrs[-1]}-12-31"),
                  latitude=a.slat, longitude=a.slon, maxradius=a.maxrad, minmagnitude=a.minmag, limit=20000)
        if a.eventtype:
            kw["eventtype"] = a.eventtype
        cat = Client("USGS").get_events(**kw)
        recs = []
        for e in cat:
            o = e.origins[0]; m = e.magnitudes[0].mag if e.magnitudes else np.nan
            recs.append(dict(ot=str(o.time), lat=o.latitude, lon=o.longitude,
                             depth=(o.depth or 0) / 1000.0, mag=m))
        ev = pd.DataFrame(recs); ev.to_csv(cat_csv, index=False)
    ev["t"] = pd.to_datetime(ev.ot, utc=True)
    # anchor = OT + S-traveltime (hypocentral)
    epi = np.sqrt(((ev.lat - a.slat) * 111.0) ** 2 +
                  ((ev.lon - a.slon) * 111.0 * math.cos(math.radians(a.slat))) ** 2)
    hyp = np.sqrt(epi ** 2 + np.clip(ev.depth, 0, 60) ** 2)
    ev["anchor"] = ev.t.values.astype("datetime64[ns]").astype("int64") / 1e9 + hyp / VS
    ev["adt"] = pd.to_datetime(ev.anchor, unit="s", utc=True)
    ev["year"] = ev.adt.dt.year; ev["julday"] = ev.adt.dt.dayofyear
    ev = ev[[ (r.year, r.julday) in have for r in ev.itertuples(index=False)]]
    if len(ev) > a.n:
        ev = ev.sample(a.n, random_state=a.seed)
    print(f"[{a.sta}] {len(ev)} EQ windows on-disk (of catalog {os.path.getsize(cat_csv)//1024}KB)", flush=True)

    jobs = defaultdict(list)
    for r in ev.itertuples(index=False):
        jobs[(int(r.year), int(r.julday))].append((a.label, float(r.anchor)))
    payloads = [(a.net, a.sta, k[0], k[1], v) for k, v in jobs.items()]
    rows = []; done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(process_day, p) for p in payloads]):
            rows.extend(f.result()); done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(payloads)} days, {len(rows):,} windows", flush=True)
    df = pd.DataFrame(rows)
    out = f"{HERE}/data/feat_{lab}_{s}.parquet"; df.to_parquet(out, index=False)
    print(f"[{a.sta}] wrote {out}: {a.label}={len(df)}  3-comp={'yes' if 'hv_ratio' in df.columns else 'no'}", flush=True)
    print("EQ EXTRACT DONE", flush=True)


if __name__ == "__main__":
    main()
