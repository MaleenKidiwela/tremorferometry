#!/usr/bin/env python
"""Tests 2 + 3 summaries and the raw figures for Test 1/2 (PLAN.md).
-> tables/t2_zsurface_all.csv, tables/t2_summary.json, tables/t3_template_copy.csv, tables/t3_summary.json,
   figures/t1_map_offsets.png, figures/t1_offset_hist.png, figures/t2_zsurface_examples.png, figures/t2_width_vs_offset.png"""
import os, sys, glob, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S13 = os.path.dirname(OUT)
CTB = f'{S13}/catalogue/tables'
ROOT = '/home/jovyan/tremorferometry'
LAT0, LON0 = 48.37, -123.61
KX, KY = 111.19 * np.cos(np.radians(LAT0)), 111.19
BLUE, ORANGE, INK, MUTED = '#2a78d6', '#eb6834', '#0b0b0b', '#8a8984'


def xy(lat, lon):
    return (np.asarray(lon) - LON0) * KX, (np.asarray(lat) - LAT0) * KY


def template_locs():
    T = pd.read_csv(f'{S13}/tables/hg2018_table1.csv', comment='#')
    L = {f'T{r.template}': (r.lat, r.lon) for r in T[T.template <= 5].itertuples()}
    B = pd.read_csv(f'{CTB}/e3_bank.csv')
    for r in B[B.primary & B.name.str.match(r'^C\d$')].itertuples():
        L[r.name] = (r.lat, r.lon)
    R = pd.read_csv(f'{CTB}/prod2018_relay_events.csv')
    for c, g in R.groupby('cluster'):
        L[f'RC{c}'] = (g.lat_best.mean(), g.lon_best.mean())
    return L


def main():
    E = pd.read_csv(f'{OUT}/tables/events_final.csv')
    T = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f'{OUT}/tables/t2_zsurface_2*.csv'))])
    E = E.merge(T, on=['eid', 'year'], how='left')
    E.to_csv(f'{OUT}/tables/t2_zsurface_all.csv', index=False)
    ets = E[E.period == 'ETS']; lf = ets[np.isfinite(ets.c10_lat)]
    J = dict(n_all=len(E), n_ets=len(ets), node_reproduced=float((E.node_stored == E.node_local).mean()))
    for nm, g in (('all', E), ('ets', ets), ('ets_lfe', lf)):
        for c in ('n_dz1', 'req_dz1_km', 'maxd_dz1_km', 'ext_E_dz1_km', 'ext_N_dz1_km', 'n_dz2', 'req_dz2_km', 'maxd_dz2_km'):
            J[f'{nm}_{c}_median'] = float(g[c].median())
        J[f'{nm}_frac_single_node_dz1'] = float((g.n_dz1 == 1).mean())
    J['ets_lfe_lfe_node_in_dz1'] = float(lf.lfe_node_in_dz1.mean()); J['ets_lfe_lfe_node_in_dz2'] = float(lf.lfe_node_in_dz2.mean())
    J['ets_lfe_dz_lfe_node_median'] = float(lf.dz_lfe_node.median())
    J['ets_lfe_offset_median_km'] = float(lf.c10_dist.median())
    # width vs z
    for lo_, hi_ in ((5, 7), (7, 10), (10, 100)):
        g = lf[(lf.z_stored >= lo_) & (lf.z_stored < hi_)]
        J[f'z{lo_}-{hi_}'] = dict(n=len(g), req_dz1_median=float(g.req_dz1_km.median()), maxd_dz1_median=float(g.maxd_dz1_km.median()),
                                  offset_median=float(g.c10_dist.median()), lfe_node_in_dz1=float(g.lfe_node_in_dz1.mean()))
    json.dump(J, open(f'{OUT}/tables/t2_summary.json', 'w'), indent=1)
    print(json.dumps(J, indent=1))

    # -------- Test 3 template dependence
    TL = template_locs()
    Z = np.load(f'{OUT}/data/zsurf_2011.npz')
    nla, nlo = Z['lat'], Z['lon']; nx, ny = xy(nla, nlo)

    def node_of(lat, lon):
        x, y = xy(lat, lon); return int(np.argmin(np.hypot(nx - x, ny - y)))
    TN = {k: node_of(*v) for k, v in TL.items()}
    rng = np.random.default_rng(7)
    rows = []
    for r in ets.itertuples():
        bank = [k for k in TN if not k.startswith('RC')] + ([k for k in TN if k.startswith('RC')] if r.run == 'prodr' else [])
        ex, ey = xy(r.lat_best, r.lon_best)
        def d_t(k):
            j = TN[k]; return float(np.hypot(nx[j] - ex, ny[j] - ey))
        o = dict(eid=r.eid, year=r.year, template=r.template, agree=r.agree, d_template=d_t(r.template),
                 d_agree=d_t(r.agree) if isinstance(r.agree, str) and r.agree in TN else np.nan, c10_dist=r.c10_dist)
        o['d_rand_mean'] = float(np.mean([d_t(k) for k in bank]))
        o['copy10'] = o['d_template'] <= 10.5; o['copy14'] = o['d_template'] <= 14.2
        o['copy10_either'] = o['copy10'] or (np.isfinite(o['d_agree']) and o['d_agree'] <= 10.5)
        o['rand_copy10'] = float(np.mean([d_t(k) <= 10.5 for k in bank]))
        o['rand_copy10_either'] = float(np.mean([(d_t(a) <= 10.5) or (d_t(b) <= 10.5) for a in bank for b in bank if a != b]))
        o['nearer_template_than_lfe'] = (o['d_template'] < r.c10_dist) if np.isfinite(r.c10_dist) else np.nan
        tx, ty = xy(*TL[r.template])
        o['tpl_dist_to_lfe'] = float(np.hypot(tx - xy(r.c10_lat, r.c10_lon)[0], ty - xy(r.c10_lat, r.c10_lon)[1])) if np.isfinite(r.c10_lat) else np.nan
        rows.append(o)
    C = pd.DataFrame(rows)
    RL = pd.read_csv(f'{CTB}/prod2018_relay_events.csv')
    ep_ = ets.set_index('eid').epoch
    C['relay_self'] = [bool(np.any(np.abs(RL.epoch.values - ep_[e]) <= 30)) and y == 2018 for e, y in zip(C.eid, C.year)]
    C.to_csv(f'{OUT}/tables/t3_template_copy.csv', index=False)
    Cl = C[np.isfinite(C.c10_dist)]
    J3 = dict(n=len(C), copy10=float(C.copy10.mean()), copy14=float(C.copy14.mean()), copy10_either=float(C.copy10_either.mean()),
              scrambled_copy10=float(C.rand_copy10.mean()), scrambled_copy10_either=float(C.rand_copy10_either.mean()),
              median_d_template=float(C.d_template.median()), median_d_random_template=float(C.d_rand_mean.median()),
              n_lfe=len(Cl), nearer_template_than_lfe=float(Cl.nearer_template_than_lfe.astype(float).mean()),
              median_tpl_dist_to_lfe=float(Cl.tpl_dist_to_lfe.median()), median_vlfe_dist_to_lfe=float(Cl.c10_dist.median()),
              n_relay_self=int(C.relay_self.sum()), copy10_no_relay_self=float(C[~C.relay_self].copy10.mean()),
              scrambled_copy10_no_relay_self=float(C[~C.relay_self].rand_copy10.mean()),
              by_template={k: dict(n=int(len(g)), copy10=float(g.copy10.mean())) for k, g in C.groupby('template')},
              by_run={k: float(g.copy10.mean()) for k, g in C.merge(ets[['eid', 'run']], on='eid').groupby('run')})
    json.dump(J3, open(f'{OUT}/tables/t3_summary.json', 'w'), indent=1)
    print(json.dumps(J3, indent=1))

    # -------- figures
    F = pd.read_csv(f'{ROOT}/data/catalogue/families_v1_3.csv', low_memory=False)
    W = np.load(f'{OUT}/data/win_2012.npz')
    stas, sll = W['stas'], W['sll']
    D = lf
    yrs = sorted(D.year.unique())
    cm = plt.get_cmap('viridis', len(yrs)); col = {y: cm(i) for i, y in enumerate(yrs)}
    fig, ax = plt.subplots(1, 2, figsize=(14, 6.4), gridspec_kw=dict(width_ratios=[1.25, 1]))
    a = ax[0]
    a.scatter(nlo, nla, s=4, c='#cfcfcb', zorder=0, label='scan grid nodes (10 km, Slab2 depth)')
    for reg, mk in (('S_Vancouver_Is', 'o'), ('Royer_NW_WA', 's')):
        f = F[F.region == reg]
        a.scatter(f.lon, f.lat, s=10, marker=mk, facecolor='none', edgecolor=MUTED, lw=0.6, zorder=1,
                  label=f'LFE families {reg}')
    for y in yrs:
        g = D[D.year == y]
        for r in g.itertuples():
            a.plot([r.lon_best, r.c10_lon], [r.lat_best, r.c10_lat], '-', color=col[y], lw=0.8, alpha=0.7, zorder=2)
        a.scatter(g.c10_lon, g.c10_lat, s=14, color=col[y], marker='x', lw=1, zorder=3)
        jit = np.random.default_rng(y).normal(0, 0.006, (2, len(g)))
        a.scatter(g.lon_best + jit[0], g.lat_best + jit[1], s=26, color=col[y], edgecolor='white', lw=0.5, zorder=4,
                  label=f'{y} (n={len(g)})')
    for s, (la, lo) in zip(stas, sll):
        if 47.5 < la < 49.6 and -125.3 < lo < -122.2:
            a.plot(lo, la, '^', color=INK, ms=6, zorder=5); a.annotate(s.split('.')[1], (lo, la), fontsize=7, xytext=(3, 3), textcoords='offset points')
    for k, (la, lo) in TL.items():
        if not k.startswith('RC'):
            a.plot(lo, la, 'P', color='#e87ba4', ms=7, mec=INK, mew=0.5, zorder=5)
            a.annotate(k, (lo, la), fontsize=6, xytext=(3, -8), textcoords='offset points', color=INK)
    a.set_xlim(-125.0, -122.4); a.set_ylim(47.4, 49.3); a.set_aspect(1 / np.cos(np.radians(48.3)))
    a.set_xlabel('lon'); a.set_ylabel('lat')
    a.set_title(f'{len(D)} tier VLFEs (dots, jittered 0.6 km) → LFE centroid ±10 min (x), coloured by year\n'
                'pink + = H&G T1–T5 and E3 C0–C5 template locations; triangles = scan stations', fontsize=9)
    a.legend(fontsize=6.5, loc='lower left', ncol=2, frameon=False)
    a = ax[1]
    for y in yrs:
        g = D[D.year == y]
        a.scatter(g.c10_dE, g.c10_dN, s=22, color=col[y], edgecolor='white', lw=0.5, label=str(y))
    a.axhline(0, color=MUTED, lw=0.8); a.axvline(0, color=MUTED, lw=0.8)
    m = (D.c10_dE.mean(), D.c10_dN.mean())
    a.annotate('', xy=m, xytext=(0, 0), arrowprops=dict(arrowstyle='->', color=INK, lw=2))
    a.text(m[0], m[1] + 3, f'mean {np.hypot(*m):.1f} km', fontsize=8, ha='center')
    a.set_xlabel('VLFE − LFE centroid, east (km)'); a.set_ylabel('north (km)'); a.set_aspect('equal')
    a.set_xlim(-80, 80); a.set_ylim(-80, 80)
    a.set_title('offset vectors (each dot one event)', fontsize=9); a.legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(f'{OUT}/figures/t1_map_offsets.png', dpi=140); plt.close(fig)

    S1 = pd.read_csv(f'{OUT}/tables/t1_scrambled_one_draw.csv')
    SC = pd.read_csv(f'{OUT}/tables/t1_scrambled.csv')
    fig, ax = plt.subplots(1, 4, figsize=(15, 3.6))
    bins = np.arange(-90, 91, 5)
    for a, c, s, lab in ((ax[0], 'c10_dE', 'dE', 'east offset (km)'), (ax[1], 'c10_dN', 'dN', 'north offset (km)')):
        a.hist(S1[s], bins, color=ORANGE, alpha=0.55, label='scrambled control (one re-pairing)')
        a.hist(D[c], bins, histtype='step', color=BLUE, lw=2, label='real pairing')
        a.axvline(0, color=MUTED, lw=0.8); a.set_xlabel(lab); a.set_ylabel('events')
    ax[0].legend(fontsize=7, frameon=False)
    bd = np.arange(0, 101, 4)
    ax[2].hist(np.hypot(S1.dE, S1.dN), bd, color=ORANGE, alpha=0.55, label='scrambled (one re-pairing)')
    ax[2].hist(D.c10_dist, bd, histtype='step', color=BLUE, lw=2, label='real')
    ax[2].set_xlabel('distance VLFE – LFE centroid (km)'); ax[2].legend(fontsize=7, frameon=False)
    ax[3].hist(SC.median_dist, 30, color=ORANGE, alpha=0.7, label='scrambled medians (1000)')
    ax[3].axvline(np.median(D.c10_dist), color=BLUE, lw=2, label='real median')
    ax[3].set_xlabel('median distance (km)'); ax[3].legend(fontsize=7, frameon=False)
    fig.suptitle('VLFE minus concurrent-LFE centroid (±10 min, P ≥ 0.75 SVI + WA-NW selected)', fontsize=10)
    fig.tight_layout(); fig.savefig(f'{OUT}/figures/t1_offset_hist.png', dpi=140); plt.close(fig)

    # z-surface examples: 3 strongest, 3 median-z events with LFEs
    Zs = {y: np.load(f'{OUT}/data/zsurf_{y}.npz') for y in yrs}
    lfs = lf.sort_values('z_stored')
    pick = list(lfs.tail(3).itertuples()) + list(lfs.iloc[len(lfs) // 2 - 1:len(lfs) // 2 + 1].itertuples()) + list(lfs.head(2).itertuples())
    fig, ax = plt.subplots(2, 4, figsize=(15, 7.6)); ax = ax.ravel()
    for a, r in zip(ax[:-1], pick):
        z = Zs[r.year]; k = list(z['eid']).index(r.eid); zz = z['z'][k]
        sc = a.scatter(z['lon'], z['lat'], c=zz[:len(z['lon'])] - np.nanmax(zz), s=34, marker='s', cmap='magma', vmin=-4, vmax=0)
        m1 = zz >= np.nanmax(zz) - 1
        a.scatter(z['lon'][m1], z['lat'][m1], s=34, marker='s', facecolor='none', edgecolor='white', lw=0.8)
        a.plot(r.lon_best, r.lat_best, 'o', color=BLUE, mec='white', ms=8)
        a.plot(r.c10_lon, r.c10_lat, 'X', color='#1baf7a', mec='white', ms=9)
        a.set_title(f'{r.time_pt[:16]} z={r.z_stored:.1f}\nΔz≤1: {r.n_dz1} nodes; offset {r.c10_dist:.0f} km', fontsize=8)
        a.set_xlim(-124.9, -122.5); a.set_ylim(47.75, 49.35); a.set_aspect(1 / np.cos(np.radians(48.3))); a.tick_params(labelsize=7)
    a = ax[-1]
    a.scatter(lf.req_dz1_km, lf.c10_dist, s=16, c=lf.z_stored, cmap='viridis', vmin=5, vmax=15)
    a.plot([0, 60], [0, 60], color=MUTED, lw=0.8); a.set_xlabel('Δz≤1 equivalent radius (km)'); a.set_ylabel('offset to LFE centroid (km)')
    a.set_title('width vs offset (colour = scan z 5–15)', fontsize=8)
    fig.subplots_adjust(right=0.9, wspace=0.3); cax = fig.add_axes([0.93, 0.56, 0.012, 0.3]); fig.colorbar(sc, cax=cax, label='z − z_best (white: Δz ≤ 1)')
    fig.suptitle('Production scan z surface at the best time (blue o = catalogue location, green X = LFE centroid ±10 min)', fontsize=10)
    fig.savefig(f'{OUT}/figures/t2_zsurface_examples.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
