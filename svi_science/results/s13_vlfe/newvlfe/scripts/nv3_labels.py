#!/usr/bin/env python
"""NV3 blind eye-label set (PLAN Amendment NV3 §5): every real 'new' level-8 row (all periods, all years) + the same number of
control level-8 'new' rows (random, matched by year where possible), mixed under random IDs. Each figure = raw 0.02-0.05 Hz
Z displacement record section at the row's node, every station on ITS OWN source day (so shift controls show what the
detector saw), -120..300 s; title = ID + scan numbers only (no time, no kind). Key + LFE / storm / template columns ->
data/nv3_label_key.csv (do not open before labelling); blank label sheet -> tables/nv3_labels.csv; figures/nv3_labels/<ID>.png"""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
S13 = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe'; NV = f'{S13}/newvlfe'
sys.path.insert(0, f'{S13}/scripts'); sys.path.insert(0, f'{S13}/catalogue/scripts'); sys.path.insert(0, f'{S13}/hg_repro/scripts')
import cat_common as K
FD = f'{NV}/figures/nv3_labels'; os.makedirs(FD, exist_ok=True)
FARS = ('IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR', 'UW.OMAK', 'TA.I02D')
MFDIR = {2011: 'prod2011', 2013: 'prodq2013', 2016: 'prodq2016', 2018: 'prodr2018', 2019: 'prodf2019'}


def sd(y):
    return f'{S13}/data/scan' if y == 2011 else f'{S13}/data/scan_{y}'


_cache = {}


def arr(f):
    if f not in _cache:
        _cache[f] = np.load(f, mmap_mode='r') if os.path.isfile(f) else None
    return _cache[f]


def fig(e, fid, cc, RT):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    y = int(e.year); rd = 'polBW' if e.run == 'realW' else e.run
    meta = json.load(open(f'{NV}/data/runs/{y}/{rd}/meta.json'))
    cfg = pd.read_parquet(f'{NV}/data/runs/{y}/{rd}/daycfg.parquet', filters=[('comp', '==', 'Z')])
    d, u = divmod(int(e.tsec), 86400); nday = 366 if y % 4 == 0 else 365
    on = dict(zip(cfg[cfg.day == d].sta, cfg[cfg.day == d].w_on))
    refbp = f'{sd(y)}/bp2_ref_v2' if y in (2011, 2012) else f'{sd(y)}/bp2_ref'
    tr = []
    for k, s in enumerate(meta['stations']):
        ts = int((d + meta['offsets'][k]) % nday) * 86400 + u
        tr.append((s, f'{meta["bp"]}/{s}.npy', ts, bool(on.get(s, False))))
    for r in RT.itertuples():
        if r.netsta in meta['stations']:
            continue
        kk = int(r.k_ctrl) + (30 if meta['shift'] else 0)
        ts = int((d + meta['shift'] * kk) % nday) * 86400 + u
        tr.append((r.netsta, f'{refbp}/{r.netsta}.npy', ts, True))
    rows = []
    for s, f, ts, act in tr:
        X = arr(f)
        if X is None or s not in cc.index or ts - 120 < 0 or ts + 300 > X.shape[1]:
            continue
        x = np.asarray(X[2, ts - 120:ts + 300], float) * 1e9
        if np.isfinite(x).sum() < 0.8 * len(x):
            continue
        rows.append((float(K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], e.lat, e.lon)), s, x, act))
    rows.sort()
    pk = {s: np.nanmax(np.abs(x[120:300])) for _, s, x, act in rows if act and s not in FARS}
    med = np.nanmedian(list(pk.values())) if pk else 1
    amp = np.nanmax([v for v in pk.values() if v <= 10 * med] or [1.0])
    fg, ax = plt.subplots(figsize=(9, 0.42 * len(rows) + 1.6))
    t = np.arange(-120, 300)
    for k, (dist, s, x, act) in enumerate(rows):
        col = '0.7' if not act else ('tab:red' if dist <= 40 else ('tab:blue' if dist <= 150 else 'k'))
        ax.plot(t, k + 0.45 * np.clip(x / amp, -2.5, 2.5), color=col, lw=0.8)
        ax.text(303, k, f'{s.split(".")[1]} {dist:.0f} km {np.nanmax(np.abs(x[120:300])):.0f} nm' + ('' if act else ' (off)'), fontsize=7, va='center')
    ax.axvspan(-60, 120, color='0.93', zorder=0); ax.axvline(0, color='0.5', lw=0.5)
    ax.set_yticks([]); ax.set_xlim(-120, 390); ax.set_xlabel('s after scan time (fit window shaded)')
    ax.set_title(f'ID {fid}   z {e.z_best:.1f}  z_loso {e.z_loso:.1f}  ratio_far {e.ratio_far:.2f}  emerg_near {e.emerg_med_near:.1f} '
                 f'far {e.emerg_far_med:.1f}  max/2nd {e.emerg_max_over_second:.1f}\nscale {amp:.0f} nm; red ≤ 40 km, blue ≤ 150 km, black far',
                 fontsize=8)
    fg.tight_layout(); fg.savefig(f'{FD}/{fid}.png', dpi=100); plt.close(fg)


def main():
    D = pd.read_parquet(f'{NV}/tables/nv3_rows.parquet')
    X = D[(D.level == 8) & (D['loc'] == 'new')]
    R = X[X.kind == 'real']; C = X[X.kind == 'control']
    rng = np.random.default_rng(20261016)
    NMAX = 120                                             # NV3 §5 amendment: stratified random sample if more real rows
    if len(R) > NMAX:
        nall = len(R)
        R = pd.concat([g.sample(min(len(g), max(1, int(round(NMAX * len(g) / nall)))), random_state=int(rng.integers(1e9)))
                       for _, g in R.groupby('year')])
    pick = []
    for y, n in R.year.value_counts().items():
        cy = C[C.year == y]
        pick.append(cy.sample(min(n, len(cy)), random_state=int(rng.integers(1e9))))
    Cs = pd.concat(pick) if pick else C.iloc[:0]
    short = len(R) - len(Cs)
    if short > 0:
        rest = C.drop(Cs.index)
        Cs = pd.concat([Cs, rest.sample(min(short, len(rest)), random_state=1)])
    S = pd.concat([R, Cs]).sample(frac=1, random_state=2).reset_index(drop=True)
    ids = set()
    while len(ids) < len(S):
        ids.add(''.join(rng.choice(list('ABCDEFGHJKLMNPQRSTUVWXYZ23456789'), 4)))
    S['ID'] = sorted(ids, key=lambda _: rng.random())
    import scan_refs as SR
    cc = SR.all_coords()
    for y in sorted(S.year.unique()):
        os.environ['SCAN_YEAR'] = str(y)
        RT = SR.load_table(); RT = RT[RT.lat.notna()]
        for e in S[S.year == y].itertuples():
            try:
                fig(e, e.ID, cc, RT)
            except Exception as ex:  # noqa
                print('fig fail', e.ID, ex)
    # columns (timeline-based; for controls they are the scrambled-control reference)
    import nv_survivors as SV
    Xl, LF, Rs, WA, F, P = SV.lfe_tables()
    out = []
    for e in S.itertuples():
        ep = float(e.epoch)
        n_any = int(Xl.cnt(LF['svi075'], [ep], 600)[0] + Xl.cnt(LF['wa_nw'], [ep], 600)[0])
        Fn = F[K.dkm(F.lat.values, F.lon.values, e.lat, e.lon) <= 20]
        sv = set(pd.to_numeric(Fn[Fn.region == 'S_Vancouver_Is'].family_id, errors='coerce').dropna().astype(int))
        wa = set(Fn[Fn.region == 'Royer_NW_WA'].family_id.astype(str))
        tt = np.sort(np.r_[Rs[Rs.fam.isin(sv)].t.values, WA[WA.fam.isin(wa)].t.values])
        obs = int(Xl.cnt(tt, [ep], 600)[0]); null = Xl.cnt(tt, np.floor(ep / 86400) * 86400 + rng.uniform(0, 86400, 1000), 600)
        dP = K.dkm(P.lat.values, P.lon.values, e.lat, e.lon)
        mf = np.nan
        try:
            M = pd.read_parquet(f'{S13}/catalogue/data/{MFDIR[int(e.year)]}/rows_scored_v2.parquet', columns=['template', 'source', 'epoch', 'zmad'])
            M = M[(M.source == 'real') & ~M.template.astype(str).str.startswith('N')]
            h = M[np.abs(M.epoch - ep) <= 30]; mf = float(h.zmad.max()) if len(h) else np.nan
        except Exception:  # noqa
            pass
        stm = SV.storm(int(e.year)); q = stm[stm.day == str(pd.Timestamp(ep, unit='s').date())]
        out.append(dict(ID=e.ID, n_lfe_any_10m=n_any, patch_lfe_10m=obs, patch_p=float(np.mean(null >= obs)) if obs else 1.0,
                        n_pnsn_100km_10m=int(((dP <= 100) & (np.abs(P.t.values - ep) <= 600)).sum()), mf_best_z=mf,
                        storm_pct=float(q.pct.iloc[0]) if len(q) else np.nan))
    S = S.merge(pd.DataFrame(out), on='ID')
    keep = ['ID', 'kind', 'ctype', 'run', 'variant', 'year', 'period', 'time_pt', 'tsec', 'epoch', 'lat', 'lon', 'dep', 'z_best', 'z_loso', 'ratio_far',
            'n_far', 'vr_thrust', 'share_excess', 'pw2_n', 'pw2_vr', 'emerg_med_near', 'n_near_emerg_ge2', 'n_near', 'emerg_max_over_second',
            'emerg_lead_sta', 'emerg_far_med', 't3_near_mid', 't3_far_over_nearmax', 't3_spearman', 'eq_wave', 'eq_far_ge_near', 'eq_m45_90min',
            'dist_tpl_km', 'nearest_tpl', 'nearest_family', 'dist_family_km', 'n_lfe_any_10m', 'patch_lfe_10m', 'patch_p', 'n_pnsn_100km_10m',
            'mf_best_z', 'storm_pct']
    S[[c for c in keep if c in S]].to_csv(f'{NV}/data/nv3_label_key.csv', index=False)
    pd.DataFrame(dict(ID=sorted(S.ID), label='', note='')).to_csv(f'{NV}/tables/nv3_labels.csv', index=False)
    print(f'label set: {len(S)} figures ({len(R)} real, {len(Cs)} control) in {FD}')


if __name__ == '__main__':
    main()
