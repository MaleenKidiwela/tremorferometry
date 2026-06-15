"""External validation: do B011's matched-filter detections coincide in TIME with
Lin (2023) LFE origin times? B011 is a PB borehole co-located with PGC in the
densest part of Lin's catalog, and it's already trust-graded (38 GOLD).

If GOLD-family detections cluster on independent Lin LFE times far above a
time-shifted null -> "GOLD" is externally confirmed as real LFE (first outside
ground truth the trust battery has had). Stratify by grade; also do a
location-matched per-GOLD-family pass to label each family Lin-confirmed or not.

Read-only. Outputs printed summary + data/lin_coinc_B011_byfamily.csv.
"""
import pandas as pd, numpy as np, math

DT = 4.0          # +/- seconds for an event-level time match
NULL_SHIFT = 123.4 * 86400.0   # null: shift B011 times by ~123 days

def to_sec(s):
    """Unix seconds, robust to pandas us/ns datetime resolution."""
    return pd.to_datetime(s, utc=True).values.astype('datetime64[ns]').astype('int64') / 1e9

# --- Lin LFE origin times (sorted unix seconds) + locations ---
lin = pd.read_csv('data/raw_lfe/lin2023_lfe.csv', usecols=['OT', 'lat', 'lon'])
lin_t = to_sec(lin.OT)
o = np.argsort(lin_t)
lin_t = lin_t[o]; lin_lat = lin.lat.to_numpy()[o]; lin_lon = lin.lon.to_numpy()[o]
TMIN, TMAX = lin_t.min(), lin_t.max()
print(f"Lin: {len(lin_t):,} LFEs  {pd.Timestamp(TMIN,unit='s')} .. {pd.Timestamp(TMAX,unit='s')}", flush=True)

# --- grades for B011 ---
A = pd.read_csv('results/family_verdicts/ALL_families_labeled.csv')
A = A[A.station == 'B011']
grade = dict(zip(A.fam, A.grade))
print(f"B011 graded families: {len(grade)}  ({(A.grade=='GOLD').sum()} GOLD)", flush=True)

def nearest_dt(times):
    """min |t - nearest Lin OT| for each time (seconds)."""
    idx = np.clip(np.searchsorted(lin_t, times), 1, len(lin_t) - 1)
    return np.minimum(np.abs(times - lin_t[idx - 1]), np.abs(times - lin_t[idx]))

# --- chunked pass over B011's 44.7M detections ---
times_by_grade = {g: [] for g in ['GOLD', 'TRUSTED', 'UNDETERMINED', 'FAIL']}
gold_times = {}     # family -> list of time-arrays (location-matched pass)
nrows = 0
for ch in pd.read_csv('data/mf_b011_all.csv', usecols=['template', 'time'], chunksize=5_000_000):
    nrows += len(ch)
    g = ch.template.map(grade)
    ch = ch[g.notna()].copy(); ch['g'] = g[g.notna()].to_numpy()
    t = to_sec(ch.time)
    inwin = (t >= TMIN) & (t <= TMAX)
    t = t[inwin]; gg = ch.g.to_numpy()[inwin]; fam = ch.template.to_numpy()[inwin]
    for grd in times_by_grade:
        times_by_grade[grd].append(t[gg == grd])
    for f in np.unique(fam[gg == 'GOLD']):
        gold_times.setdefault(f, []).append(t[fam == f])
    print(f"  ...{nrows:,} rows scanned", flush=True)

print("\n=== GRADE-STRATIFIED COINCIDENCE (B011 dets vs Lin OT, +/-%.0fs) ===" % DT, flush=True)
print(f"{'grade':14s}{'n_det':>12s}{'real%':>9s}{'null%':>9s}{'enrich':>8s}", flush=True)
for grd, arrs in times_by_grade.items():
    tt = np.concatenate(arrs) if arrs else np.array([])
    if len(tt) == 0:
        print(f"{grd:14s}{0:>12d}", flush=True); continue
    real = (nearest_dt(tt) < DT).mean()
    null = (nearest_dt(tt + NULL_SHIFT) < DT).mean()
    print(f"{grd:14s}{len(tt):>12,}{100*real:>9.2f}{100*null:>9.2f}{(real/null if null>0 else np.inf):>8.1f}", flush=True)

# --- per-GOLD-family, LOCATION-matched (Lin LFEs within R km of the family patch) ---
def parse_ll(fid):
    a = fid.split('__')[0].split('_'); return float(a[0]), float(a[1])
R = 8.0
rows = []
for f, arrs in gold_times.items():
    tt = np.concatenate(arrs)
    flat, flon = parse_ll(f)
    dkm = np.sqrt(((lin_lat - flat) * 111.0) ** 2 +
                  ((lin_lon - flon) * 111.0 * math.cos(math.radians(flat))) ** 2)
    near = np.sort(lin_t[dkm < R])
    if len(near) < 50:
        loc_real = np.nan
    else:
        idx = np.clip(np.searchsorted(near, tt), 1, len(near) - 1)
        nd = np.minimum(np.abs(tt - near[idx - 1]), np.abs(tt - near[idx]))
        loc_real = (nd < DT).mean()
    glob = (nearest_dt(tt) < DT).mean()
    rows.append(dict(fam=f, n_det=len(tt), lin_near=len(near),
                     coinc_global=round(glob, 4), coinc_locmatch=round(loc_real, 4)))
out = pd.DataFrame(rows).sort_values('coinc_locmatch', ascending=False)
out.to_csv('data/lin_coinc_B011_byfamily.csv', index=False)
print(f"\nPer-GOLD-family (location-matched, R={R}km) -> data/lin_coinc_B011_byfamily.csv", flush=True)
print(out.to_string(index=False), flush=True)
conf = (out.coinc_locmatch > 0.2).sum()
print(f"\nLin-CONFIRMED GOLD families (locmatch coincidence >20%): {conf}/{len(out)}", flush=True)
print("COINCIDENCE TEST DONE", flush=True)
