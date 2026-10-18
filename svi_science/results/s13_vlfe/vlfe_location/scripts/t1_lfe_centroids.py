#!/usr/bin/env python
"""Test 1 (PLAN.md): concurrent-LFE centroids for every final tier VLFE, offsets, scrambled control (within-episode
re-pairing), bootstrap over episodes. Read-only on ../catalogue.
-> tables/events_final.csv, tables/t1_offsets.csv, tables/t1_summary_year.csv, tables/t1_summary.json,
   tables/t1_scrambled.csv"""
import os, sys, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
S13 = os.path.dirname(OUT)
ROOT = '/home/jovyan/tremorferometry'
sys.path.insert(0, f'{S13}/scripts')
CTB = f'{S13}/catalogue/tables'
LAT0, LON0 = 48.37, -123.61
KX, KY = 111.19 * np.cos(np.radians(LAT0)), 111.19


def xy(lat, lon):
    return (np.asarray(lon) - LON0) * KX, (np.asarray(lat) - LAT0) * KY


def events():
    rows = []
    for y in range(2010, 2019):
        nm = 'prodr' if y >= 2016 else 'prod'
        E = pd.read_csv(f'{CTB}/{nm}{y}_events_v2.csv')
        E['year'] = y; E['run'] = nm
        rows.append(E)
    E = pd.concat(rows, ignore_index=True)
    E['episode'] = E.episode.fillna('')
    E['ep_key'] = np.where(E.episode == '', 'between_' + E.year.astype(str), E.episode)
    E['eid'] = np.arange(len(E))
    return E


def lfe_table():
    import s13_common as C
    R = C.svi_rows()
    R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.75)]
    F = pd.read_csv(f'{ROOT}/data/catalogue/families_v1_3.csv', low_memory=False)
    Fs = F[F.region == 'S_Vancouver_Is'].copy()
    Fs['fam'] = pd.to_numeric(Fs.family_id, errors='coerce')
    Fs = Fs.dropna(subset=['fam', 'lat', 'lon']).set_index(Fs.dropna(subset=['fam', 'lat', 'lon']).fam.astype(int))
    A = pd.DataFrame(dict(t=R.t.values, fam=R.fam.values.astype(int)))
    A['lat'] = Fs.lat.reindex(A.fam).values; A['lon'] = Fs.lon.reindex(A.fam).values
    A['famkey'] = 'SVI:' + A.fam.astype(str); A['src'] = 'svi'
    c = pd.read_parquet(f'{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet',
                        columns=['time', 'region', 'selected_v2', 'family_uid'])
    c = c[c.selected_v2 & (c.region == 'NW')]
    Fw = F[F.region == 'Royer_NW_WA'].dropna(subset=['lat', 'lon']).set_index('family_uid')
    B = pd.DataFrame(dict(t=c.time.astype('datetime64[ns]').astype('int64').values / 1e9, famkey='WA:' + c.family_uid.astype(str)))
    B['lat'] = Fw.lat.reindex(c.family_uid).values; B['lon'] = Fw.lon.reindex(c.family_uid).values
    B['src'] = 'wa'; B['fam'] = -1
    L = pd.concat([A, B], ignore_index=True).sort_values('t').reset_index(drop=True)
    return L


def centroid(L, t, w, src=None):
    a, b = np.searchsorted(L.t.values, t - w, 'left'), np.searchsorted(L.t.values, t + w, 'right')
    x = L.iloc[a:b]
    if src is not None:
        x = x[x.src == src]
    n_all = len(x)
    x = x[np.isfinite(x.lat)]
    if not len(x):
        return dict(n=n_all, n_loc=0, n_fam=0, lat=np.nan, lon=np.nan, med_lat=np.nan, med_lon=np.nan)
    fl = x.groupby('famkey')[['lat', 'lon']].first()
    return dict(n=n_all, n_loc=len(x), n_fam=len(fl), lat=x.lat.mean(), lon=x.lon.mean(),
                med_lat=fl.lat.median(), med_lon=fl.lon.median())


def main():
    E = events()
    L = lfe_table()
    rows = []
    for r in E.itertuples():
        o = dict(eid=r.eid)
        for tag, w, src in (('c10', 600, None), ('c2', 120, None), ('svi10', 600, 'svi')):
            c = centroid(L, r.epoch, w, src)
            for k, v in c.items():
                o[f'{tag}_{k}'] = v
        rows.append(o)
    E = E.merge(pd.DataFrame(rows), on='eid')
    ex, ey = xy(E.lat_best, E.lon_best)
    for tag in ('c10', 'c2', 'svi10'):
        cx, cy = xy(E[f'{tag}_lat'], E[f'{tag}_lon'])
        E[f'{tag}_dE'] = ex - cx; E[f'{tag}_dN'] = ey - cy; E[f'{tag}_dist'] = np.hypot(ex - cx, ey - cy)
    E.to_csv(f'{OUT}/tables/events_final.csv', index=False)
    D = E[np.isfinite(E.c10_lat)].copy()
    D.to_csv(f'{OUT}/tables/t1_offsets.csv', index=False)

    # per year / overall summaries
    def summ(g, tag='c10'):
        return dict(n=len(g), mean_dE=g[f'{tag}_dE'].mean(), mean_dN=g[f'{tag}_dN'].mean(),
                    mean_vec_km=float(np.hypot(g[f'{tag}_dE'].mean(), g[f'{tag}_dN'].mean())),
                    mean_vec_az_deg=float(np.degrees(np.arctan2(g[f'{tag}_dE'].mean(), g[f'{tag}_dN'].mean())) % 360),
                    sd_dE=g[f'{tag}_dE'].std(), sd_dN=g[f'{tag}_dN'].std(), median_dist=g[f'{tag}_dist'].median(),
                    q25_dist=g[f'{tag}_dist'].quantile(.25), q75_dist=g[f'{tag}_dist'].quantile(.75),
                    mean_dist=g[f'{tag}_dist'].mean())
    Y = []
    for y, g in D.groupby('year'):
        Y.append(dict(year=y, **summ(g)))
    for tag in ('c10', 'c2', 'svi10'):
        g = D[np.isfinite(D[f'{tag}_lat'])]
        Y.append(dict(year=f'all_{tag}', **summ(g, tag)))
    g = D[D.period == 'ETS']; Y.append(dict(year='all_ETS_c10', **summ(g)))
    g = D[D.year <= 2015]; Y.append(dict(year='2010-2015_c10', **summ(g)))
    Ys = pd.DataFrame(Y); Ys.to_csv(f'{OUT}/tables/t1_summary_year.csv', index=False)

    # scrambled control: within-episode re-pairing
    rng = np.random.default_rng(20261016)
    exD, eyD = xy(D.lat_best.values, D.lon_best.values)
    cxD, cyD = xy(D.c10_lat.values, D.c10_lon.values)
    groups = [np.where(D.ep_key.values == k)[0] for k in D.ep_key.unique()]
    groups = [g for g in groups if len(g) >= 2]
    usable = np.concatenate(groups)
    real_med = np.median(np.hypot(exD - cxD, eyD - cyD)[usable])
    real_mean_vec = (np.mean((exD - cxD)[usable]), np.mean((eyD - cyD)[usable]))
    sc = []
    for k in range(1000):
        perm = np.arange(len(D))
        for g in groups:
            p = rng.permutation(g)
            # derangement within group where possible
            for _ in range(20):
                if len(g) < 2 or not np.any(p == g):
                    break
                p = rng.permutation(g)
            perm[g] = p
        dE = exD[usable] - cxD[perm[usable]]; dN = eyD[usable] - cyD[perm[usable]]
        sc.append(dict(draw=k, median_dist=np.median(np.hypot(dE, dN)), mean_dE=dE.mean(), mean_dN=dN.mean(),
                       sd_dE=dE.std(ddof=1), sd_dN=dN.std(ddof=1)))
    SC = pd.DataFrame(sc); SC.to_csv(f'{OUT}/tables/t1_scrambled.csv', index=False)
    # one pooled scrambled offset sample (for histograms)
    perm = np.arange(len(D))
    for g in groups:
        perm[g] = rng.permutation(g)
    pd.DataFrame(dict(eid=D.eid.values[usable], dE=exD[usable] - cxD[perm[usable]], dN=eyD[usable] - cyD[perm[usable]])).to_csv(
        f'{OUT}/tables/t1_scrambled_one_draw.csv', index=False)

    # bootstrap over episodes for the mean offset vector
    eps = D.ep_key.unique()
    bs = []
    for k in range(2000):
        pick = rng.choice(eps, len(eps), replace=True)
        g = pd.concat([D[D.ep_key == e] for e in pick])
        bs.append((g.c10_dE.mean(), g.c10_dN.mean(), g.c10_dist.median()))
    bs = np.array(bs)
    J = dict(n_events_total=int(len(E)), n_with_located_lfe_10m=int(len(D)), n_with_located_lfe_2m=int(np.isfinite(D.c2_lat).sum()),
             n_in_repairing=int(len(usable)), n_groups=len(groups),
             real_median_dist=float(real_med), real_mean_dE=float(real_mean_vec[0]), real_mean_dN=float(real_mean_vec[1]),
             scr_median_dist_mean=float(SC.median_dist.mean()), scr_median_dist_p05=float(SC.median_dist.quantile(.05)),
             p_tracking=float((SC.median_dist <= real_med).mean()),
             scr_mean_dE=float(SC.mean_dE.mean()), scr_mean_dN=float(SC.mean_dN.mean()),
             scr_sd_dE=float(SC.sd_dE.mean()), scr_sd_dN=float(SC.sd_dN.mean()),
             real_sd_dE=float(np.std((exD - cxD)[usable], ddof=1)), real_sd_dN=float(np.std((eyD - cyD)[usable], ddof=1)),
             boot_dE_95=[float(np.percentile(bs[:, 0], 2.5)), float(np.percentile(bs[:, 0], 97.5))],
             boot_dN_95=[float(np.percentile(bs[:, 1], 2.5)), float(np.percentile(bs[:, 1], 97.5))],
             boot_median_dist_95=[float(np.percentile(bs[:, 2], 2.5)), float(np.percentile(bs[:, 2], 97.5))],
             n_episodes=int(len(eps)))
    json.dump(J, open(f'{OUT}/tables/t1_summary.json', 'w'), indent=1)
    pd.set_option('display.width', 250)
    print(Ys.round(2).to_string()); print(json.dumps(J, indent=1))


if __name__ == '__main__':
    main()
