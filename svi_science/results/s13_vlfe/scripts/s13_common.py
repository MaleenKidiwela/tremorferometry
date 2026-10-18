"""S13 shared paths, catalogue loaders, origin times, waveform location (rules: ../PLAN.md)."""
from __future__ import annotations
import os, glob, datetime as dt
import numpy as np
import pandas as pd

ROOT = '/home/jovyan/tremorferometry'
OUT = f'{ROOT}/svi_science/results/s13_vlfe'
TB, FIG, LOGD, DATA = f'{OUT}/tables', f'{OUT}/figures', f'{OUT}/logs', f'{OUT}/data'
DISP = os.environ.get('S13_DISP', f'{DATA}/disp1sps')
RESP = f'{DATA}/resp'
TRACE_ROOTS = ['lfe_extension/traces', 'lfe_extension/traces_backward', 'lfe_extension_wa/traces',
               'lfe_extension_wa/traces_dl', 'lfe_extension_wa/traces_dl2', 'lfe_catalog_repair/traces',
               'lfe_array_discovery/traces', 'data/waveforms', f'svi_science/results/s13_vlfe/data/raw']
VS_KMS = 3.95
SEED = 20261013


def log(msg, f=None):
    s = f'{dt.datetime.now():%Y-%m-%d %H:%M:%S} {msg}'
    print(s, flush=True)
    if f:
        with open(f, 'a') as fh:
            fh.write(s + '\n')


def day_files(netsta, date):
    """All raw day files for NET.STA on date (YYYY-MM-DD), in root order."""
    d = pd.Timestamp(date)
    out = []
    for r in TRACE_ROOTS:
        p = f'{ROOT}/{r}/{netsta}/{d.year}/{d.dayofyear:03d}.mseed'
        if os.path.isfile(p) and os.path.getsize(p) > 0:
            out.append(p)
    return out


def svi_rows():
    """SVI LFE rows (data + c1), as s10_multi_inventory.load_rows (backward primary&keep_xfam8 + forward prod)."""
    RES = f'{ROOT}/svi_science/results'
    fc = ['run', 'date', 'series', 'k', 'fam', 'P_same', 'keep_xfam8', 'time_utc']
    A = pd.read_parquet(f'{RES}/s0/analysis_table.parquet', columns=fc)
    A = A[(A.run.astype(str) == 'prod') & A.keep_xfam8].copy()
    A = A.assign(series=A.series.astype(str), P=A.P_same, counted_once=True,
                 t=A.time_utc.astype('datetime64[ns]').astype('int64') / 1e9, src='fw')
    A = A[['date', 'series', 'fam', 'P', 't', 'counted_once', 'src']]
    B = pd.read_parquet(f'{ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet',
                        columns=['date', 'series', 'fam', 'P', 't', 'primary', 'keep_xfam8', 'counted_once'])
    B = B[B.primary & (B.keep_xfam8 == 1)].copy()
    B = B.assign(series=B.series.astype(str), src='bw')[['date', 'series', 'fam', 'P', 't', 'counted_once', 'src']]
    R = pd.concat([B, A], ignore_index=True)
    R['fam'] = R.fam.astype(int)
    R['date'] = R.date.astype(str)
    return R


def svi_families():
    F = pd.read_csv(f'{ROOT}/data/catalogue/families_v1_3.csv', low_memory=False)
    F = F[F.region == 'S_Vancouver_Is'].copy()
    F['fam'] = pd.to_numeric(F.family_id, errors='coerce').astype(int)
    return F[['fam', 'lat', 'lon', 'depth_km', 'location_quality']].set_index('fam')


def station_coords():
    """lat/lon/elev per NET.STA from the cached StationXML (first epoch)."""
    from obspy import read_inventory
    rows = {}
    for f in sorted(glob.glob(f'{RESP}/*.xml')):
        inv = read_inventory(f)
        for net in inv:
            for sta in net:
                rows.setdefault(f'{net.code}.{sta.code}', (sta.latitude, sta.longitude, sta.elevation))
    return pd.DataFrame(rows, index=['lat', 'lon', 'elev']).T


def origin_shift_svi(fams):
    """Per family: origin = t + shift, shift = median_sta(offset_S - r/VS)."""
    from s13_green import ned_offset
    O = pd.read_csv(f'{ROOT}/lfe_extension/results/cut/svi_ext/offsets_used.csv')
    F = svi_families(); C = station_coords()
    out = {}
    for f in fams:
        o = O[O.fam == f]
        v = []
        for r in o.itertuples():
            if r.netsta not in C.index or f not in F.index:
                continue
            n, e = ned_offset(F.loc[f, 'lat'], F.loc[f, 'lon'], C.loc[r.netsta, 'lat'], C.loc[r.netsta, 'lon'])
            rr = np.sqrt(n ** 2 + e ** 2 + (F.loc[f, 'depth_km'] * 1e3) ** 2) / 1e3
            v.append(r.offset_s - rr / VS_KMS)
        out[f] = float(np.median(v)) if v else np.nan
    return pd.Series(out, name='shift')


GUARD_FRAC, GUARD_CAP = 1e-3, 0.3          # PLAN.md Amendment 4 (2026-10-14 ~5:52 pm PT)


def guard_weights(sig, frac=GUARD_FRAC, cap=GUARD_CAP):
    """Amendment 4 channel guard. sig: dict sta -> (3,) control sigma (N, E, Z). Per component, a station whose sigma is
    <= frac x the median sigma of that component over the unit's stations (or non-finite) gets weight 0; the rest get
    1/sigma^2, then no single station-component holds more than `cap` of the total (water-filling; < 4 live components ->
    equal weights). Returns (wts dict sta -> (3,), share dict sta -> (3,), dropped dict sta -> (3,) bool, capped flag)."""
    stas = list(sig)
    S = np.array([np.asarray(sig[s], float) for s in stas]).reshape(len(stas), 3)
    med = np.array([np.median(S[np.isfinite(S[:, c]), c]) if np.isfinite(S[:, c]).any() else np.nan for c in range(3)])
    live = np.isfinite(S) & (S > frac * np.where(np.isfinite(med), med, np.inf)[None, :]) & (S > 0)
    W = np.where(live, 1.0 / np.where(live, S, 1.0) ** 2, 0.0)
    v = W[live]
    capped = False
    if v.size and v.size * cap < 1:
        W[live] = 1.0; capped = True
    elif v.size:
        vs = np.sort(v)[::-1]
        c = None
        for k in range(vs.size):
            ck = cap * vs[k:].sum() / (1 - cap * k)
            if vs[k] <= ck * (1 + 1e-12):
                c = ck if k > 0 else None
                break
        if c is not None:
            W = np.minimum(W, c); capped = True
    tot = W.sum()
    share = W / tot if tot > 0 else W
    return ({s: W[i] for i, s in enumerate(stas)}, {s: share[i] for i, s in enumerate(stas)},
            {s: ~live[i] for i, s in enumerate(stas)}, capped)


def guard_violation(sig, frac=GUARD_FRAC, cap=GUARD_CAP, tol=1e-6):
    """True if plain 1/sigma^2 weights (the pre-Amendment-4 fit) break the guard: a flat/dead component in use, or one
    component above `cap` of the total weight."""
    _, _, dropped, _ = guard_weights(sig, frac, cap)
    if any(d.any() for d in dropped.values()):
        return True
    W = np.array([1.0 / np.asarray(sig[s], float) ** 2 for s in sig])
    return bool(W.size and (W / W.sum()).max() > cap * (1 + tol))


# ---- Amendment 4 v2 (2026-10-14 ~6:25 pm PT, after Merlin's review; replaces the 0.3 cap rule above, NO weight cap) ----
GUARD_RULE = 'A4v2'


def guard2(sig, nused, min_used=5, f_sta=0.05, f_unit=0.02):
    """Drop a station-component if n_used < 5, sigma <= 0 (or non-finite), sigma < 0.05 x median sigma of that station's
    other two components, or sigma < 0.02 x median sigma over all channels in the unit. Kept channels get 1/sigma^2.
    sig, nused: dict sta -> (3,) (N, E, Z). Returns (wts dict, dropped dict sta -> (3,) bool, reasons dict sta -> list)."""
    stas = list(sig)
    S = np.array([np.asarray(sig[s], float) for s in stas]).reshape(len(stas), 3)
    N = np.array([np.asarray(nused[s], float) for s in stas]).reshape(len(stas), 3)
    pos = np.isfinite(S) & (S > 0)
    med_unit = np.median(S[pos]) if pos.any() else np.nan
    reasons = {s: ['', '', ''] for s in stas}
    drop = np.zeros_like(S, bool)
    for i, s in enumerate(stas):
        for c in range(3):
            r = []
            if not N[i, c] >= min_used:
                r.append('n_used<5')
            if not pos[i, c]:
                r.append('sigma<=0')
            else:
                o = S[i, [k for k in range(3) if k != c]]
                o = o[np.isfinite(o)]
                if o.size and S[i, c] < f_sta * np.median(o):
                    r.append('<0.05xsta')
                if np.isfinite(med_unit) and S[i, c] < f_unit * med_unit:
                    r.append('<0.02xunit')
            if r:
                drop[i, c] = True; reasons[s][c] = '+'.join(r)
    W = np.where(drop, 0.0, 1.0 / np.where(pos, S, 1.0) ** 2)
    return ({s: W[i] for i, s in enumerate(stas)}, {s: drop[i] for i, s in enumerate(stas)}, reasons)
