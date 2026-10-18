"""S15 shared loaders: the 540 closed-catalogue ETS VLFEs, the LFE pool used by cat_e2.lfe_status, family positions,
episode windows and front directions. Rules: ../PLAN.md."""
import sys, os, json
import numpy as np, pandas as pd

ROOT = '/home/jovyan/tremorferometry'
RES = f'{ROOT}/svi_science/results'
CAT = f'{RES}/s13_vlfe/catalogue'
TB = f'{CAT}/tables'
S15 = f'{RES}/s15_vlfe_science'
PMIN = float(os.environ.get('S15_PMIN', '0.75'))          # Amendment 3: P >= 0.5 concurrency column
OUT_T = f'{S15}/tables' if PMIN == 0.75 else f'{S15}/tables/p%03d' % round(PMIN * 100)
os.makedirs(OUT_T, exist_ok=True)
OUT_F = f'{S15}/figures'
sys.path[:0] = [f'{RES}/s13_vlfe/scripts', f'{RES}/s13_vlfe/hg_repro/scripts', f'{CAT}/scripts', f'{ROOT}/svi_science/scripts']

RUNS = ['prod2010', 'prod2011', 'prod2012', 'prodq2013', 'prod2014', 'prod2015', 'prodq2016', 'prodq2017', 'prodrq2018',
        'prodf2019', 'prodf2020', 'prodfy2021', 'prodfy2022', 'prodfy2023', 'prodfy2024', 'prodfy2025', 'prodfy2026']
BIG9 = ['2011-08', '2012-08', '2013-09', '2018-06', '2019-08', '2020-11', '2023-04', '2024-08', '2025-09b1']
DAY = 86400.0
KM = 111.19
LAT0, LON0 = 48.3, -123.6


def xy(lat, lon):
    return (np.asarray(lon) - LON0) * KM * np.cos(np.radians(LAT0)), (np.asarray(lat) - LAT0) * KM


def to_pt(ep):
    t = pd.to_datetime(np.asarray(ep, float), unit='s', utc=True).tz_convert('America/Los_Angeles')
    return t.strftime('%Y-%m-%d %H:%M:%S %Z')


def episodes():
    """union of episode windows (UTC epoch s; window = start 00:00 .. end+1 00:00 UTC)."""
    E = pd.concat([pd.read_csv(f'{TB}/{r}_episodes_v2.csv').assign(run=r) for r in RUNS], ignore_index=True)
    g = E.groupby('episode').agg(start=('start', 'min'), end=('end', 'max'), major=('major', 'max'),
                                 compl_38=('compl_3.8', 'mean')).reset_index()
    tz = 'UTC'   # table dates are UTC dates (an ETS-labelled event at 2022-11-30 03:37 UTC sits on the start date)
    g['t0'] = pd.to_datetime(g.start).dt.tz_localize(tz).dt.tz_convert('UTC').dt.as_unit('ns').astype('int64') / 1e9
    g['t1'] = (pd.to_datetime(g.end) + pd.Timedelta(days=1)).dt.tz_localize(tz).dt.tz_convert('UTC').dt.as_unit('ns').astype('int64') / 1e9
    return g


def vlfes():
    rows = []
    for r in RUNS:
        e = pd.read_csv(f'{TB}/{r}_ets_events_v2.csv')
        a = pd.read_csv(f'{TB}/{r}_events_v2.csv')[['epoch', 'episode', 'period', 'Mw_scan']]
        m = e.merge(a.drop_duplicates('epoch'), on='epoch', how='left')
        assert len(m) == len(e)
        m['run'] = r
        rows.append(m)
    V = pd.concat(rows, ignore_index=True)
    assert len(V) == 540, len(V)
    E = episodes().set_index('episode')
    V['t0_ep'] = V.episode.map(E.t0); V['t1_ep'] = V.episode.map(E.t1)
    V['flag_2016_17'] = V.run.isin(['prodq2016', 'prodq2017'])
    V['year'] = pd.to_datetime(V.epoch, unit='s').dt.year
    V['vx'], V['vy'] = xy(V.lat_best, V.lon_best)
    return V.sort_values('epoch').reset_index(drop=True)


def families():
    """all families of the pool: key, src, lat, lon, x, y, along (S0 frame, + = NW), svi_fid (nearest SVI family)."""
    F = pd.read_csv(f'{RES}/s0/families.csv')
    S = pd.DataFrame(dict(key=['S:%d' % f for f in F.fid], src='svi', fid=F.fid, lat=F.lat, lon=F.lon, along=F.along_strike_km))
    S['x'], S['y'] = xy(S.lat, S.lon)
    # the 15 Bostock families outside the production 115 (backward catalogue): positions from e5r_families_2012.csv
    B = pd.read_csv(f'{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv')
    B = B[~B.bkey.isin(F.fid)]
    SB = pd.DataFrame(dict(key=['S:%d' % f for f in B.bkey], src='svi_b', fid=B.bkey, lat=B.lat, lon=B.lon, along=np.nan))
    SB['x'], SB['y'] = xy(SB.lat, SB.lon)
    # linear map (x, y) -> along for the S0 frame
    X = np.c_[np.ones(len(S)), S.x, S.y]
    coef, *_ = np.linalg.lstsq(X, S.along.values, rcond=None)
    resid = S.along.values - X @ coef
    c = pd.read_parquet(f'{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet',
                        columns=['family', 'family_uid', 'region', 'selected_v2'])
    c = c[c.selected_v2 & (c.region == 'NW')].drop_duplicates('family')
    FV = pd.read_csv(f'{ROOT}/data/catalogue/families_v1_3.csv', low_memory=False)
    FV = FV[FV.region == 'Royer_NW_WA'].set_index('family_uid')
    lat = c.family_uid.map(FV.lat).fillna(c.family_uid.map(FV.pos_lat)).fillna(c.family_uid.map(FV.approx_lat))
    lon = c.family_uid.map(FV.lon).fillna(c.family_uid.map(FV.pos_lon)).fillna(c.family_uid.map(FV.approx_lon))
    W = pd.DataFrame(dict(key=['W:%d' % f for f in c.family], src='wa', fid=-1, lat=lat.values, lon=lon.values))
    W['x'], W['y'] = xy(W.lat, W.lon)
    W['along'] = coef[0] + coef[1] * W.x + coef[2] * W.y
    SB['along'] = coef[0] + coef[1] * SB.x + coef[2] * SB.y
    A = pd.concat([S, SB, W], ignore_index=True)
    sx, sy = S.x.values, S.y.values
    d = np.hypot(A.x.values[:, None] - sx[None], A.y.values[:, None] - sy[None])
    A['svi_fid'] = S.fid.values[np.nanargmin(np.where(np.isfinite(d), d, 1e9), 1)]
    A.loc[A.src == 'svi', 'svi_fid'] = A.loc[A.src == 'svi', 'fid']
    info = dict(along_fit_coef=coef.tolist(), along_fit_rms_km=float(np.sqrt(np.mean(resid ** 2))),
                n_wa=int(len(W)), n_wa_no_pos=int(W.lat.isna().sum()))
    return A, info


def lfes():
    """pool LFEs: t (UTC epoch s), key; SVI P>=0.75 data counted_once 2003+, WA NW selected_v2 (as hg_r3.lfe_lists)."""
    import s13_common as C
    Rw = C.svi_rows()
    Rw = Rw[(Rw.series == 'data') & Rw.counted_once & (Rw.P >= PMIN)]
    t0 = (pd.Timestamp('2003-01-01') - pd.Timestamp('1970-01-01')).total_seconds()
    Rw = Rw[Rw.t >= t0]
    S = pd.DataFrame(dict(t=Rw.t.values.astype(float), key=['S:%d' % f for f in Rw.fam.values]))
    c = pd.read_parquet(f'{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet',
                        columns=['time', 'family', 'region', 'selected_v2'])
    c = c[c.selected_v2 & (c.region == 'NW')]
    tw = c.time.astype('datetime64[ns]').astype('int64').values / 1e9
    W = pd.DataFrame(dict(t=tw, key=['W:%d' % f for f in c.family.values]))
    W = W[W.t >= t0]
    L = pd.concat([S, W], ignore_index=True).sort_values('t').reset_index(drop=True)
    return L


def directions(L, F, E):
    """front direction per episode: s3 table where present; else sign of Spearman(first-LFE time, along) over families
    with >= 3 LFEs in the window (+ = NW-ward)."""
    from scipy.stats import spearmanr
    k = pd.read_csv(f'{RES}/s3/episodes_kinematics.csv')
    k = k[k.region == 'SVI'].set_index('episode')
    al = F.set_index('key').along
    out = []
    for _, e in E.iterrows():
        m = (L.t >= e.t0) & (L.t < e.t1)
        g = L[m].groupby('key').t.agg(['min', 'size'])
        g = g[g['size'] >= 3]
        g['along'] = g.index.map(al)
        g = g.dropna()
        rho = spearmanr(g['min'], g['along']).statistic if len(g) >= 5 else np.nan
        own = 'NW-ward' if rho > 0 else ('SE-ward' if rho < 0 else None)
        s3 = k.direction.get(e.episode) if e.episode in k.index else None
        out.append(dict(episode=e.episode, n_fam_onset=len(g), rho_onset_along=rho, dir_own=own, dir_s3=s3,
                        direction=s3 if s3 is not None else own, dir_source='s3' if s3 is not None else 'own_onsets'))
    return pd.DataFrame(out)
