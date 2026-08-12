"""Decisive follow-up to lin_coincidence_b011.py: a LAG HISTOGRAM.

Zero-lag coincidence was null, but (1) B011 mf detection time ~= S-arrival ~= Lin OT + ~10s
(travel time), so genuine matches sit at a NON-zero lag; (2) raw lists are noise-diluted.
Both are fixed by histogramming the signed offset to the nearest Lin event over +-60s and
looking for a PEAK above a flat noise background. Genuine shared events concentrate at the
true lag; pure noise stays flat. Peak/background ratio = strength of real coincidence.

Caches GOLD per-family detection times to data/b011_gold_times.npz on first run.
Read-only otherwise. Output: printed histogram + peak stats, by grade and location-matched.
"""
import pandas as pd, numpy as np, math, os

LAGWIN = 60.0          # +/- seconds searched
BIN = 0.5              # histogram bin (s)
edges = np.arange(-LAGWIN, LAGWIN + BIN, BIN)
centers = (edges[:-1] + edges[1:]) / 2

def to_sec(s):
    return pd.to_datetime(s, utc=True).values.astype('datetime64[ns]').astype('int64') / 1e9

lin = pd.read_csv('data/raw_lfe/lin2023_lfe.csv', usecols=['OT', 'lat', 'lon'])
lin_t = to_sec(lin.OT)
o = np.argsort(lin_t)
lin_t = lin_t[o]; lin_lat = lin.lat.to_numpy()[o]; lin_lon = lin.lon.to_numpy()[o]
TMIN, TMAX = lin_t.min(), lin_t.max()

A = pd.read_csv('results/family_verdicts/ALL_families_labeled.csv')
A = A[A.station == 'B011']
grade = dict(zip(A.fam, A.grade))
gold_fams = [f for f, g in grade.items() if g == 'GOLD']

CACHE = 'data/b011_gold_times.npz'
if os.path.exists(CACHE):
    z = np.load(CACHE, allow_pickle=True)
    gold_times = {k: z[k] for k in z.files}
    print(f"loaded cached GOLD times for {len(gold_times)} families", flush=True)
else:
    acc = {f: [] for f in gold_fams}
    n = 0
    for ch in pd.read_csv('data/mf_b011_all.csv', usecols=['template', 'time'], chunksize=5_000_000):
        n += len(ch)
        ch = ch[ch.template.isin(gold_fams)]
        if len(ch):
            t = to_sec(ch.time); fam = ch.template.to_numpy()
            inwin = (t >= TMIN) & (t <= TMAX); t = t[inwin]; fam = fam[inwin]
            for f in np.unique(fam):
                acc[f].append(t[fam == f])
        print(f"  ...{n:,} rows", flush=True)
    gold_times = {f: (np.concatenate(v) if v else np.array([])) for f, v in acc.items()}
    np.savez(CACHE, **gold_times)
    print(f"cached -> {CACHE}", flush=True)

def signed_nearest(times, ref):
    """signed (times - nearest ref) for each time; ref sorted."""
    if len(ref) == 0 or len(times) == 0:
        return np.array([])
    idx = np.clip(np.searchsorted(ref, times), 1, len(ref) - 1)
    left = times - ref[idx - 1]; right = times - ref[idx]
    return np.where(np.abs(left) < np.abs(right), left, right)

# ---- GLOBAL lag histogram (all GOLD dets vs all Lin) ----
allg = np.concatenate([v for v in gold_times.values() if len(v)])
d = signed_nearest(allg, lin_t)
h, _ = np.histogram(d, bins=edges)
bg = np.median(h[(np.abs(centers) > 20)])      # flat background away from any peak
pk = h.argmax()
print("\n=== GLOBAL lag histogram (GOLD dets vs all Lin OT) ===", flush=True)
print(f"n_det={len(allg):,}  peak at lag {centers[pk]:+.1f}s  height={h[pk]:,}  bg/bin~{bg:.0f}  peak/bg={h[pk]/max(bg,1):.1f}", flush=True)
# show the +-15s region
for c, cnt in zip(centers, h):
    if abs(c) <= 15 and cnt > 0:
        bar = '#' * int(40 * cnt / h[pk]) if h[pk] else ''
        print(f"  {c:+6.1f}s {cnt:8d} {bar}", flush=True)

# ---- LOCATION-MATCHED aggregate lag histogram ----
def parse_ll(fid):
    a = fid.split('__')[0].split('_'); return float(a[0]), float(a[1])
R = 8.0
hist_loc = np.zeros(len(centers))
ndet_loc = 0
for f, tt in gold_times.items():
    if len(tt) == 0: continue
    flat, flon = parse_ll(f)
    dkm = np.sqrt(((lin_lat - flat) * 111.0) ** 2 +
                  ((lin_lon - flon) * 111.0 * math.cos(math.radians(flat))) ** 2)
    near = np.sort(lin_t[dkm < R])
    if len(near) < 50: continue
    d = signed_nearest(tt, near)
    hh, _ = np.histogram(d, bins=edges)
    hist_loc += hh; ndet_loc += len(tt)
bg = np.median(hist_loc[(np.abs(centers) > 20)])
pk = hist_loc.argmax()
print("\n=== LOCATION-MATCHED aggregate lag histogram (R=8km) ===", flush=True)
print(f"n_det={ndet_loc:,}  peak at lag {centers[pk]:+.1f}s  height={hist_loc[pk]:,.0f}  bg/bin~{bg:.0f}  peak/bg={hist_loc[pk]/max(bg,1):.1f}", flush=True)
for c, cnt in zip(centers, hist_loc):
    if abs(c) <= 15 and cnt > 0:
        bar = '#' * int(40 * cnt / hist_loc[pk]) if hist_loc[pk] else ''
        print(f"  {c:+6.1f}s {cnt:8.0f} {bar}", flush=True)
print("\nLAG HISTOGRAM DONE", flush=True)
