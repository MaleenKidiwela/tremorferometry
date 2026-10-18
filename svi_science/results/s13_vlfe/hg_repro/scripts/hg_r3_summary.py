#!/usr/bin/env python
"""Amendment R3 summary: merge rows/events/aux, apply the pre-registered pass rules (R3.2) and classification (R3.4),
write tables/r3_events.csv (every candidate, all columns), tables/r3_pass_counts.csv, tables/r3_real_looking.csv,
figures/r3_overview.png, figures/r3_aggregate.png, figures/r3_wave_*.png (raw waveforms, near + far stations)."""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
import hg_common as H
import hg_r3 as X

FAMREG = ['S_Vancouver_Is', 'Royer_NW_WA', 'Royer_SW_WA', 'Royer_N_Vancouver_Is', 'Olympic_2010', 'Olympic_2010_armB']


def families():
    F = pd.read_csv(f'{H.ROOT}/data/catalogue/families_v1_3.csv', low_memory=False)
    F = F[F.region.isin(FAMREG)].copy()
    for a, b in (('pos_lat', 'pos_lon'), ('approx_lat', 'approx_lon')):
        m = F.lat.isna() & F[a].notna()
        F.loc[m, 'lat'] = F.loc[m, a]; F.loc[m, 'lon'] = F.loc[m, b]
    return F[F.lat.notna()][['region', 'family_id', 'lat', 'lon']].reset_index(drop=True)


def dkm(la1, lo1, la2, lo2):
    return np.hypot((la1 - la2) * 111.19, (lo1 - lo2) * 111.19 * np.cos(np.radians(48.3)))


def load():
    R = pd.read_parquet(f'{X.R3}/rows{X.SFX}.parquet')
    E = pd.concat([pd.read_parquet(f'{X.R3}/events{X.SFX}_{y}.parquet') for y in (2011, 2012)], ignore_index=True)
    A = pd.read_parquet(f'{X.R3}/aux{X.SFX}.parquet')
    D = R.merge(E, on='row_id', how='left').merge(A, on='row_id', how='left')
    F = families()
    P = pd.read_csv(f'{H.S13}/tables/scan_footprint_polygon.csv')
    poly = Path(np.c_[P.lon, P.lat])
    TT = H.templates().set_index('template')
    ok = D.lat_best.notna()
    D['dist_family_km'] = np.nan
    D.loc[ok, 'dist_family_km'] = [float(np.min(dkm(la, lo, F.lat.values, F.lon.values)))
                                   for la, lo in zip(D.lat_best[ok], D.lon_best[ok])]
    D['nearest_family'] = ''
    D.loc[ok, 'nearest_family'] = [f'{F.region.values[k]}:{F.family_id.values[k]}' for k in
                                   [int(np.argmin(dkm(la, lo, F.lat.values, F.lon.values))) for la, lo in
                                    zip(D.lat_best[ok], D.lon_best[ok])]]
    D['dist_template_km'] = dkm(D.lat_best, D.lon_best, D.template.map(TT.lat), D.template.map(TT.lon))
    D['in_footprint_poly'] = False
    D.loc[ok, 'in_footprint_poly'] = poly.contains_points(np.c_[D.lon_best[ok], D.lat_best[ok]])
    D['Mw_scan'] = (np.log10(D.m0_cols.abs()) - 9.1) / 1.5
    # ---- pre-registered pass rules (R3.2)
    D['a_testable'] = D.n6 >= 2
    D['pass_a'] = D.a_testable & (D.ratio6 <= 0.551)
    D['pass_b'] = ~D.travelling.fillna(False).astype(bool)
    D['locatable'] = D.z_best >= 5
    D['pass_c'] = D.locatable & ~D.grid_edge.fillna(True).astype(bool) & (D.fp_interior.fillna(False).astype(bool)
                                                                          | (D.dist_family_km <= 20))
    D['pass_e'] = ~(D.in_head2.fillna(True) | D.in_bw5.fillna(True) | D.in_reg.fillna(True) | D.in_reg3.fillna(True)).astype(bool)
    D['pass_f'] = (D.z_loso >= 3) & (D.z_loso >= 0.5 * D.z_best)
    D['real_looking'] = D.pass_a & D.pass_b & D.pass_c & D.pass_e & D.pass_f

    def cls(r):
        if not r.pass_e:
            return 'catalogued earthquake'
        if not r.a_testable:
            return 'far/near untestable'
        if not (r.pass_a and r.pass_b):
            return 'travelling wave'
        if r.locatable and not r.pass_f:
            return 'single-station'
        if not r.locatable:
            return 'not locatable (z<5)'
        if not r.pass_c:
            return 'outside footprint/family zone or grid edge'
        return 'real-looking local VLFE'
    D['class'] = D.apply(cls, axis=1)
    lfe = (D.n_svi075_10m.fillna(0) + D.n_wa_nw_10m.fillna(0)) >= 1
    cand = (D.n_c010_075_10m.fillna(0) + D.n_pnsn100_10m.fillna(0)) >= 1
    D['lfe_status'] = np.where(lfe, 'with LFEs', np.where(cand, 'candidates/tremor only', 'quiet'))
    return D


def wave_fig(r, fn):
    """raw 0.02-0.05 Hz displacement, vertical: scan near network, H&G network, far references; t0-300 .. t0+400 s."""
    import hg_r2 as R2
    import scan_refs as SR
    cc = SR.all_coords()
    ep = float(r.epoch); y = pd.Timestamp(ep, unit='s').year
    sd = f'{X.S13D}/scan' if y == 2011 else f'{X.S13D}/scan_2012'
    y0 = (pd.Timestamp(f'{y}-01-01') - pd.Timestamp('1970-01-01')).total_seconds()
    a = int(round(ep - y0)) - 300; n = 700
    la0, lo0 = (r.lat_best, r.lon_best) if np.isfinite(r.lat_best) else (48.3, -123.5)
    tr = []
    for s in ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV']:
        f = f'{sd}/bp_{y}g/{s}.npy'
        if os.path.isfile(f):
            Xa = np.load(f, mmap_mode='r'); tr.append(('near (scan)', s, np.asarray(Xa[2, a:a + n], float) * 1e9))
    i = int(r.i) - 300
    for s in H.STAS:
        Xc = np.load(f'{H.MF}/cont_{s}.npy', mmap_mode='r')
        tr.append(("H&G network", s, np.asarray(Xc[2, i:i + n], float) * 1e9))
    for s in ['UW.GNW', 'IU.COR', 'US.NEW', 'US.HAWA', 'UW.OMAK', 'TA.I02D']:
        f = f'{sd}/bp2_ref/{s}.npy'
        if os.path.isfile(f):
            Xa = np.load(f, mmap_mode='r'); tr.append(('far reference', s, np.asarray(Xa[2, a:a + n], float) * 1e9))
    dist = []
    for g, s, x in tr:
        try:
            dist.append(dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], la0, lo0))
        except KeyError:
            dist.append(np.nan)
    o = np.lexsort((np.nan_to_num(dist, nan=9e9), [{'near (scan)': 0, 'H&G network': 1, 'far reference': 2}[g] for g, _, _ in tr]))
    t = np.arange(n) - 300
    fig, ax = plt.subplots(figsize=(10, 0.45 * len(tr) + 1.5))
    col = {'near (scan)': 'k', 'H&G network': 'tab:blue', 'far reference': 'tab:red'}
    # template overlay on H&G Z channels (displacement, LS-scaled)
    TT = H.templates().set_index('template'); io = H.tidx(TT.loc[int(r.template)].time)
    for k, j in enumerate(o):
        g, s, x = tr[j]
        if not np.isfinite(x).any():
            ax.text(-300, -k, f'{s} no data', fontsize=7, va='center'); continue
        x = x - np.nanmean(x)
        pk = np.nanmax(np.abs(x))
        ax.plot(t, x / (2.2 * pk) - k, color=col[g], lw=0.7)
        if g == 'H&G network':
            Xc = np.load(f'{H.MF}/cont_{s}.npy', mmap_mode='r')
            w = np.asarray(Xc[2, io:io + H.L], float) * 1e9
            seg = x[300:300 + H.L]
            if np.isfinite(w).all() and np.isfinite(seg).all():
                w = w - w.mean(); sc = (seg @ w) / (w @ w)
                ax.plot(np.arange(H.L), sc * w / (2.2 * pk) - k, color='orange', lw=1.0, alpha=0.9)
        ax.text(402, -k, f'{s}  {dist[j]:.0f} km  peak {pk:.0f} nm', fontsize=6.5, va='center', color=col[g])
    ax.axvspan(0, H.L, color='0.9', zorder=0)
    ax.set_xlim(-300, 560); ax.set_yticks([]); ax.set_xlabel('s after t0 (origin of our matched detection)')
    ax.set_title(f'cand {int(r.cand_id)}  T{int(r.template)}  {r.time_utc} UTC = {r.time_pt}\n'
                 f'z_best {r.z_best:.1f}, ratio6 {r.ratio6:.2f}, class: {r["class"]}; Z, 0.02-0.05 Hz displacement, '
                 f'each trace scaled to its own peak (orange = template, LS-scaled)', fontsize=8)
    fig.tight_layout(); fig.savefig(fn, dpi=120); plt.close(fig)


def main():
    D = load()
    C = pd.read_csv(f'{H.TB}/r3{X.SFX}_candidates.csv')
    cand = D[D.kind == 'cand'].merge(C[['cand_id', 'i_S1', 't_S1', 'S', 'zmad', 'd_ours_minus_theirs', 'matched_raw30',
                                        'rank_in_ours', 'n_ours', 'top30', 'time_utc', 'time_pt', 'in_SVI_episode_pm2d',
                                        'in_NW_episode_pm2d']], on='cand_id').sort_values('zmad', ascending=False)
    cols_first = ['cand_id', 'template', 'time_utc', 'time_pt', 't_S1', 'zmad', 'top30', 'class', 'real_looking', 'pass_a',
                  'pass_b', 'pass_c', 'pass_e', 'pass_f', 'lfe_status', 'z_best', 'z_loso', 'ratio6', 'n6', 'ratio_far',
                  'ref_ratio', 'ratio_GNW', 'travelling', 'pw2_n', 'pw2_nfar250', 'pw2_vr', 'pw2_vapp', 'pw2_baz',
                  'lat_best', 'lon_best', 'dep_best', 'fp_interior', 'grid_edge', 'dist_family_km', 'nearest_family',
                  'dist_template_km', 'z_fixnode', 'vr_thrust', 'vr_free', 'Mw_tpl', 'Mw_scan', 'g_cc', 'g_pct_random',
                  'hf1_60s', 'hf1_60s_pct', 'hf2_60s', 'hf1_10m', 'n_svi075_60s', 'n_svi075_10m', 'p_svi075_10m',
                  'n_wa_nw_10m', 'n_c010_075_10m', 'n_pnsn100_10m', 'p_pnsn100_10m', 'in_head2', 'in_bw5', 'in_reg',
                  'in_reg3', 'in_m45', 'eq_names', 'hg_top_sta', 'hg_top_share', 'hg_zmad_without_top']
    rest = [c for c in cand.columns if c not in cols_first]
    cand[cols_first + rest].to_csv(f'{H.TB}/r3{X.SFX}_events.csv', index=False)
    # pass counts, candidates vs scrambled control vs selection-matched noise
    rows = []
    for kind, lab in (('cand', 'candidates'), ('scr', 'scrambled control (same-day random times)'),
                      ('noise', 'selection-matched noise-template detections (10 noise templates)')):
        for sub in ('ALL', 'TOP30'):
            d = D[D.kind == kind]
            if sub == 'TOP30':
                d = d[d.cand_id.isin(C[C.top30].cand_id)]
            r = dict(set=lab, subset=sub, n=len(d))
            for p in ['a_testable', 'pass_a', 'pass_b', 'locatable', 'pass_c', 'pass_e', 'pass_f', 'real_looking']:
                r[p] = int(d[p].sum()); r[p + '_frac'] = float(d[p].mean()) if len(d) else np.nan
            for c_ in d['class'].value_counts().items():
                r['class: ' + c_[0]] = int(c_[1])
            rows.append(r)
    PC = pd.DataFrame(rows)
    PC.to_csv(f'{H.TB}/r3{X.SFX}_pass_counts.csv', index=False)
    RL = cand[cand.real_looking]
    RL[cols_first].to_csv(f'{H.TB}/r3{X.SFX}_real_looking.csv', index=False)
    # base rate of the travelling flag at random between-ETS times
    pb = D[D.kind == 'pwbase']
    base = dict(n=len(pb), travelling_rate=float(pb.travelling.mean()), pw2_n_median=float(pb.pw2_n.median()))
    pd.Series(base).to_csv(f'{H.TB}/r3{X.SFX}_pw_baserate.csv')
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
    print(PC.T.to_string())
    print('pw base', base)
    print(cand[cols_first[:30]].head(40).round(2).to_string(index=False))
    # ---- figures
    fig, ax = plt.subplots(2, 2, figsize=(13, 10))
    s = D[D.kind == 'scr']; nz = D[D.kind == 'noise']
    a0 = ax[0, 0]
    a0.scatter(s.z_best, s.ratio6, s=6, c='0.75', label='scrambled control')
    a0.scatter(nz.z_best, nz.ratio6, s=8, facecolors='none', edgecolors='tab:cyan', lw=0.5, label='selection-matched noise')
    sc = a0.scatter(cand.z_best, cand.ratio6, s=28, c=cand.zmad, cmap='viridis', edgecolors='k', lw=0.4, label='candidates')
    plt.colorbar(sc, ax=a0, label='our CC at the candidate (MAD units)')
    a0.axhline(0.551, color='r', ls='--', lw=0.8); a0.axvline(5, color='r', ls='--', lw=0.8)
    a0.set_xscale('log'); a0.set_yscale('log'); a0.set_xlabel('scan z_best (t0 ± 30 s)'); a0.set_ylabel('ratio6 (GNW + far LHZ / near)')
    a0.legend(fontsize=7)
    a1 = ax[0, 1]
    P = pd.read_csv(f'{H.S13}/tables/scan_footprint_polygon.csv'); F = families()
    a1.plot(np.r_[P.lon, P.lon[:1]], np.r_[P.lat, P.lat[:1]], 'g-', lw=1, label='scan LFE footprint')
    a1.scatter(F.lon, F.lat, s=4, c='0.6', label='LFE families')
    sc = a1.scatter(cand.lon_best, cand.lat_best, s=20 + 4 * np.clip(cand.z_best.fillna(0), 0, 20), c=cand.z_best, cmap='magma_r',
                    edgecolors='k', lw=0.4, label='candidate best node (size/colour = z)')
    plt.colorbar(sc, ax=a1, label='z_best')
    TT = H.templates()
    a1.scatter(TT.lon, TT.lat, marker='*', s=150, c='tab:blue', edgecolors='k', label='H&G templates')
    a1.set_xlabel('lon'); a1.set_ylabel('lat'); a1.legend(fontsize=7)
    a2 = ax[1, 0]
    sc = a2.scatter(cand.zmad, cand.g_cc, c=cand.z_best, cmap='magma_r', s=24, edgecolors='k', lw=0.4)
    plt.colorbar(sc, ax=a2, label='z_best')
    a2.axhline(0, color='0.5', lw=0.6)
    a2.set_xlabel('our CC at the candidate (MAD units)'); a2.set_ylabel('held-out-station CC with the template event (PGC SNB PFB YOUB)')
    a3 = ax[1, 1]
    keys = ['a_testable', 'pass_a', 'pass_b', 'locatable', 'pass_c', 'pass_e', 'pass_f', 'real_looking']
    xx = np.arange(len(keys))
    for k, (lab, c_) in enumerate([('candidates', 'k'), ('scrambled control (same-day random times)', '0.6'),
                                   ('selection-matched noise-template detections (10 noise templates)', 'tab:cyan')]):
        v = PC[(PC.set == lab) & (PC.subset == 'ALL')].iloc[0]
        a3.bar(xx + (k - 1) * 0.27, [v[p + '_frac'] for p in keys], 0.27, color=c_, label=f'{lab.split(" (")[0]} (n={v.n})')
    a3.set_xticks(xx); a3.set_xticklabels(keys, rotation=30, fontsize=8); a3.set_ylabel('fraction passing'); a3.legend(fontsize=7)
    fig.suptitle(f'R3{X.SFX}: Table-S1 {"ETS-window (post-hoc positive control)" if X.SFX else "between-ETS"} candidates matched by our run — per-event tests (raw; no error bars)')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r3{X.SFX}_overview.png', dpi=120); plt.close(fig)
    # aggregate figure
    G = pd.read_parquet(f'{X.R3}/agg_all{X.SFX}.parquet')
    fig, ax = plt.subplots(2, 3, figsize=(14, 7))
    for i_, st in enumerate(['TOP30', 'ALL']):
        for j_, m in enumerate(['G', 'G_stack', 'H']):
            a_ = ax[i_, j_]; g = G[G.set == st]
            for src, c_ in (('random same-day', '0.6'), ('selection-matched noise template', 'tab:cyan')):
                v = g[g.source == src][m].dropna()
                a_.hist(v, bins=25, color=c_, alpha=0.7, label=f'{src} ({len(v)})')
            a_.axvline(g[g.source == 'candidates'][m].iloc[0], color='k', lw=2, label='candidates')
            a_.set_title(f'{st}: {m}', fontsize=9); a_.legend(fontsize=6)
    ax[0, 0].set_xlabel('G = mean held-out-station CC'); ax[0, 2].set_xlabel('H = H&G-network stack CC (selection-dependent)')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r3{X.SFX}_aggregate.png', dpi=120); plt.close(fig)
    # waveform figures: every real-looking candidate, plus the top 3 by rank as raw examples
    show = pd.concat([RL, cand.head(3)]).drop_duplicates('cand_id')
    for r in show.itertuples():
        try:
            wave_fig(cand.loc[cand.cand_id == r.cand_id].iloc[0], f'{H.FIG}/r3{X.SFX}_wave_cand{int(r.cand_id):03d}.png')
        except Exception as e:  # noqa
            print('wave fig failed', r.cand_id, e)


if __name__ == '__main__':
    main()
