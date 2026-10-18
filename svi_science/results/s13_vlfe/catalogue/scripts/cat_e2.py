#!/usr/bin/env python
"""Catalogue E2 (PLAN_CAT §3-§6; Merlin step 2): the R3 battery on every V4 T1/T2/T3 detection >= 5.5 MAD,
2011-06-01 - 2013-01-01, on the 8 + 8 scrambled-control catalogues and on 10 noise templates (same absolute summed CC
as the real template), then counts vs MF level (5.5, 6, 7, 8, 10 MAD), (e) / e', +- S5, single / two-template tiers.
Modes
  rows            -> data/e2/rows.parquet (all T1-T5 detections >= 5.5 MAD; T4/T5 used for agreement only)
  battery YEAR    -> data/e2/bat_YEAR.parquet (run with SCAN_YEAR=YEAR; cheap columns for every unique time, expensive
                     columns for times passing (c) and (e) or the e' time part)
  summary         -> tables/e2_*.csv, figures/e2_*.png"""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

E2D = f'{K.DATA}/e2'
os.makedirs(E2D, exist_ok=True)
LEVELS = [5.5, 6, 7, 8, 10]
NNOISE = 10
F0, F1 = '2011-06-01', '2013-01-01'
NW = int(os.environ.get('CAT_NW', '8'))


def H():
    import hg_common
    return hg_common


def run_rows():
    Hc = H()
    lo, hi = Hc.tidx(F0), Hc.tidx(F1)
    out, thr = [], {}
    for k in range(1, 6):
        z = np.load(f'{K.HGR}/data/r2/real_vel_0_sum_all_90_T{k}.npz')
        med, mad = float(z['real_med']), float(z['real_mad'])
        thr[k] = (med, mad)
        for nm in ['real'] + [f'A1_{r}' for r in range(8)] + [f'A2_{r}' for r in range(8)]:
            pi, pv, pn = z[nm + '_pi'], z[nm + '_pv'], z[nm + '_pnl']
            m = (pi >= lo) & (pi < hi) & (pv >= med + LEVELS[0] * mad)
            src = 'real' if nm == 'real' else nm[:2]
            draw = -1 if nm == 'real' else int(nm[3:])
            out.append(pd.DataFrame(dict(source=src, draw=draw, template=k, i=pi[m], S=pv[m], n_live=pn[m],
                                         nch=len(z['chans']), zmad=(pv[m] - med) / mad)))
    for j in range(NNOISE):
        z = np.load(f'{K.HGR}/data/r2/noise_vel_0_sum_all_90_{j:03d}.npz')
        i0 = int(z['i0'])
        for q in range(int(z['ncfg'])):
            pi, pv, pn = z[f'c{q}_pi'], z[f'c{q}_pv'], z[f'c{q}_pnl']
            for k in z[f'c{q}_templates']:
                k = int(k); med, mad = thr[k]
                m = (pi >= lo) & (pi < hi) & (pv >= med + LEVELS[0] * mad) & (np.abs(pi - i0) > 90)
                out.append(pd.DataFrame(dict(source='noise', draw=j, template=k, i=pi[m], S=pv[m], n_live=pn[m],
                                             nch=int(z[f'c{q}_n']), zmad=(pv[m] - med) / mad)))
    R = pd.concat(out, ignore_index=True)
    EP0 = (Hc.T0 - pd.Timestamp('1970-01-01')).total_seconds()
    R['epoch'] = EP0 + R.i
    R['year'] = pd.to_datetime(R.epoch, unit='s').dt.year
    R['period'] = K.period_of(R.epoch.values)
    R.to_parquet(f'{E2D}/rows.parquet')
    K.log(f'rows: {len(R)}; by source {R.source.value_counts().to_dict()}; battery rows (T1-T3) '
          f'{int((R.template <= 3).sum())}; unique battery times {R[R.template <= 3].i.nunique()}', 'e2')


# ------------------------------------------------------------------------------ battery (scan side)
def cheap(year, I, EPS=None):
    """I: unique H-frame indices in calendar YEAR (EPS: their epochs, for other frames) -> DataFrame of cheap columns."""
    import hg_r3 as X
    X.ev_init()
    EV = X.EV
    Hc = H()
    EP0 = (Hc.T0 - pd.Timestamp('1970-01-01')).total_seconds()
    zmax, amax = EV['zmax'], EV['argmax']
    n = len(zmax)
    rows = []
    masks = EV['masks']
    for k, i in enumerate(I):
        ep = EP0 + i if EPS is None else float(EPS[k])
        tsec = int(round(ep - EV['y0']))
        lo, hi = max(tsec - 30, 0), min(tsec + 31, n)
        zw = np.asarray(zmax[lo:hi], np.float64)
        if np.isfinite(zw).any():
            k = int(np.nanargmax(zw)); ts = lo + k; node = int(amax[ts]); zb = float(zw[k])
        else:
            ts, node, zb = tsec, -1, np.nan
        a0, a1 = max(tsec - 60, 0), min(tsec + 121, n)
        o = dict(i=int(i), epoch=ep, tsec_t0=tsec, tsec_best=ts, node_best=node, z_best=zb)
        for t in ('head2', 'bw5', 'reg', 'reg3', 'm45', 'v2'):
            o[f'in_{t}'] = bool(masks[t][a0:a1].any())
        rows.append(o)
    return pd.DataFrame(rows)


def exp_init():
    import hg_r3 as X
    X.ev_init()


def exp_row(r):
    import hg_r3 as X
    SC = X.EV['SC']
    ts, node = int(r['tsec_best']), int(r['node_best'])
    o = dict(i=int(r['i']))
    c1 = SC.cols(dict(variant=0, tsec=ts, node=node, m0=1e15))
    m0 = c1.get('m0_cols', np.nan)
    c2 = SC.cols(dict(variant=0, tsec=ts, node=node, m0=m0)) if np.isfinite(m0) else c1
    keep = ['n_sta', 'n_sc', 'fp_interior', 'grid_edge', 'm0_cols', 'vr_thrust', 'cc', 'max_share', 'max_share_ch',
            'share_excess', 'vr_free', 'vr_strikeslip', 'vr_steepdip', 'n_Z_agree', 'z_chk', 'z_loso', 'loso_dropped',
            'isolation', 'hf1', 'hf2', 'hf_loud', 'A_near_nm', 'n_ref', 'ref_excess', 'ref_ratio']
    for kk in keep:
        o[kk] = c2.get(kk, np.nan)
    for tg in K.REF_ALL:
        o[f'ratio_{tg}'] = c2.get(f'ratio_{tg}', np.nan); o[f'ex_{tg}'] = c2.get(f'ex_{tg}', np.nan)
    r6 = [o[f'ratio_{tg}'] for tg in K.REF_ALL if np.isfinite(o[f'ratio_{tg}'])]
    rf = [o[f'ratio_{tg}'] for tg in K.FAR5 if np.isfinite(o[f'ratio_{tg}'])]
    o['n6'] = len(r6); o['ratio6'] = float(np.median(r6)) if r6 else np.nan
    o['n_far'] = len(rf); o['ratio_far'] = float(np.median(rf)) if rf else np.nan
    o.update(X.plane_wave2(ts, node))
    return o


def exp_chunk(rows):
    out = []
    for r in rows:
        try:
            out.append(exp_row(r))
        except Exception as e:  # noqa
            out.append(dict(i=int(r['i']), ev_error=f'{type(e).__name__}: {str(e)[:80]}'))
    return out


def run_battery(year, rows_file=f'{E2D}/rows.parquet', out_file=None, tag='e2', gate_all=False):
    assert os.environ.get('SCAN_YEAR') == str(year), 'set SCAN_YEAR'
    R = pd.read_parquet(rows_file)
    if 'template' in R and tag == 'e2':
        R = R[R.template <= 3]
    R = R[R.year == year]
    U = R.drop_duplicates('i').sort_values('i')
    I = U.i.values
    C = cheap(year, I, U.epoch.values if tag != 'e2' else None)
    N = K.node_table()
    ok = C.node_best >= 0
    nb = C.node_best.clip(lower=0)
    C['grid_edge_node'] = np.where(ok, N.grid_edge.values[nb], True)
    C['in_footprint_poly'] = np.where(ok, N.fp_interior.values[nb], False)
    C['dist_family_km'] = np.where(ok, N.dist_family_km.values[nb], np.nan)
    C['nearest_family'] = np.where(ok, N.nearest_family.values[nb], '')
    C['lat_best'] = np.where(ok, N.lat.values[nb], np.nan); C['lon_best'] = np.where(ok, N.lon.values[nb], np.nan)
    C['dep_best'] = np.where(ok, N.dep.values[nb], np.nan)
    C['head2_maxmag'] = K.head2_maxmag(C.epoch.values, np.full(len(C), year))
    pc = (C.z_best >= 5) & ~C.grid_edge_node & (C.in_footprint_poly | (C.dist_family_km <= 20))
    e_ok = ~(C.in_head2 | C.in_bw5 | C.in_reg | C.in_reg3)
    e2t = ~(C.in_bw5 | C.in_reg | C.in_reg3 | (C.in_head2 & (C.head2_maxmag >= 5.5).fillna(False)))
    C['gate'] = (pc & (e_ok | e2t)) | gate_all            # gate_all: Amendment Y8b computes every row
    K.log(f'{tag} battery {year}: {len(C)} unique times; gated for expensive columns {int(C.gate.sum())}', tag)
    G = C[C.gate][['i', 'tsec_best', 'node_best']].to_dict('records')
    out = []
    if G:
        from multiprocessing import Pool
        nch = max(1, min(len(G), NW * 12))
        chunks = [G[k::nch] for k in range(nch)]
        with Pool(NW, initializer=exp_init) as p:
            for k, r in enumerate(p.imap_unordered(exp_chunk, chunks)):
                out += r
                if k % 24 == 0:
                    K.log(f'{tag} battery {year}: chunk {k + 1}/{nch}', tag)
    Xp = pd.DataFrame(out)
    D = C.merge(Xp, on='i', how='left') if len(Xp) else C
    D.to_parquet(out_file or f'{E2D}/bat_{year}.parquet')
    K.log(f'{tag} battery {year}: done ({len(D)} rows)', tag)


# ------------------------------------------------------------------------------ summary
def lfe_status(ep):
    import hg_r3 as X
    LF = X.lfe_lists()
    ep = np.asarray(ep, float)
    n75 = X.cnt(LF['svi075'], ep, 600) + X.cnt(LF['wa_nw'], ep, 600)
    nc = X.cnt(LF['svi_c010_075'], ep, 600)
    # PNSN tremor within 100 km of the scan centre stands in for the template epicentre (T2 lies 10 km away)
    nt = X.cnt(LF['pnsn100_T2'], ep, 600)
    st = np.where(n75 >= 1, 'with LFEs', np.where(nc + nt >= 1, 'candidates/tremor only', 'quiet'))
    return st, n75


def events_dedupe(d, by):
    """count distinct events: rows within +-30 s of each other (same group) are one event."""
    n = 0
    for _, g in d.groupby(by):
        t = np.sort(g.i.values)
        n += int(1 + np.sum(np.diff(t) > 30)) if len(t) else 0
    return n


def agreement(R):
    """for every row: highest MF level (MAD) at which another template of the same source/draw has a detection within
    +-30 s (real / A1_r / A2_r), or noise j+1 (mod NNOISE) under the same template's configuration."""
    R = R.copy()
    R['agree_z'] = -np.inf
    keyed = {}
    for (src, dr, k), g in R.groupby(['source', 'draw', 'template']):
        o = np.argsort(g.i.values)
        keyed[(src, dr, k)] = (g.i.values[o], g.zmad.values[o])
    for (src, dr, k), g in R.groupby(['source', 'draw', 'template']):
        if k > 3:
            continue
        others = ([(src, dr, kk) for kk in range(1, 6) if kk != k] if src != 'noise'
                  else [(src, (dr + 1) % NNOISE, k)])
        best = np.full(len(g), -np.inf)
        for key in others:
            if key not in keyed:
                continue
            ti, tz = keyed[key]
            for j, i in enumerate(g.i.values):
                a, b = np.searchsorted(ti, i - 30), np.searchsorted(ti, i + 31)
                if b > a:
                    best[j] = max(best[j], tz[a:b].max())
        R.loc[g.index, 'agree_z'] = best
    return R


def run_summary():
    R = pd.read_parquet(f'{E2D}/rows.parquet')
    R = agreement(R)
    B = pd.concat([pd.read_parquet(f'{E2D}/bat_{y}.parquet') for y in (2011, 2012)], ignore_index=True)
    B = B.drop(columns=['epoch'])
    if 'grid_edge' in B:
        B['grid_edge'] = B.grid_edge.where(B.grid_edge.notna(), B.grid_edge_node)
    else:
        B['grid_edge'] = B.grid_edge_node
    D = R[R.template <= 3].merge(B, on='i', how='left')
    for c in ('n_far', 'n6'):
        D[c] = D[c].fillna(0)
    D = K.apply_rules(D)
    D['lfe_status'], D['n_lfe_10m'] = lfe_status(D.epoch.values)
    D.to_parquet(f'{E2D}/e2_rows_scored.parquet')
    yr_btw = K.BETWEEN_YR
    rows = []
    for L in LEVELS:
        for ecol in ('battery', 'battery_e2'):
            for s5 in (False, True):
                for agree in (False, True):
                    m = (D.zmad >= L) & D[ecol]
                    if s5:
                        m &= D.pass_S5
                    if agree:
                        m &= D.agree_z >= L
                    X = D[m]
                    r = dict(level_MAD=L, earthquake_test='e' if ecol == 'battery' else "e'", S5=s5,
                             two_template=agree)
                    for per in ('ETS', 'between'):
                        x = X[X.period == per]
                        r[f'real_{per}'] = events_dedupe(x[x.source == 'real'], ['source'])
                        cc = [events_dedupe(x[(x.source == s) & (x.draw == d)], ['source']) for s in ('A1', 'A2') for d in range(8)]
                        r[f'ctrl_mean_{per}'] = float(np.mean(cc)); r[f'ctrl_max_{per}'] = int(np.max(cc))
                        nn = [events_dedupe(x[(x.source == 'noise') & (x.draw == j)], ['source']) for j in range(NNOISE)]
                        r[f'noise_mean_{per}'] = float(np.mean(nn))
                    r['ctrl_rate_between_per_yr'] = r['ctrl_mean_between'] / yr_btw
                    r['noise_rate_between_per_tpl_yr'] = r['noise_mean_between'] / yr_btw
                    r['real_between_per_yr'] = r['real_between'] / yr_btw
                    r['accept'] = bool(r['real_ETS'] >= 8 and r['real_between'] <= r['ctrl_max_between'] + 2
                                       and r['ctrl_rate_between_per_yr'] <= 1)
                    rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(f'{K.TB}/e2_separation.csv', index=False)
    # detection counts (no battery) per level, for context
    cnt = []
    for L in LEVELS:
        x = D[D.zmad >= L]
        for per in ('ETS', 'between'):
            y = x[x.period == per]
            cnt.append(dict(level_MAD=L, period=per, real=events_dedupe(y[y.source == 'real'], ['source']),
                            ctrl_mean=np.mean([events_dedupe(y[(y.source == s) & (y.draw == d)], ['source'])
                                               for s in ('A1', 'A2') for d in range(8)]),
                            noise_mean=np.mean([events_dedupe(y[(y.source == 'noise') & (y.draw == j)], ['source'])
                                                for j in range(NNOISE)])))
    pd.DataFrame(cnt).to_csv(f'{K.TB}/e2_detection_counts.csv', index=False)
    # per-test pass fractions at 5.5 MAD
    pf = []
    for (src, per), g in D[D.period.isin(['ETS', 'between'])].groupby(['source', 'period']):
        pf.append(dict(source=src, period=per, n=len(g), **{c: float(g[c].mean()) for c in
                  ('pass_a', 'pass_b', 'pass_c', 'pass_e', 'pass_e2', 'pass_f', 'pass_S5', 'battery', 'battery_e2')}))
    pd.DataFrame(pf).to_csv(f'{K.TB}/e2_pass_fractions.csv', index=False)
    # the real events passing battery' at 5.5 MAD (any period)
    ev = D[(D.source == 'real') & D.battery_e2].sort_values('epoch')
    ev = ev.assign(time_pt=pd.to_datetime(ev.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles')
                   .dt.strftime('%Y-%m-%d %H:%M:%S %Z'))
    cols = ['time_pt', 'period', 'template', 'zmad', 'agree_z', 'battery', 'battery_e2', 'pass_S5', 'z_best', 'z_loso',
            'ratio_far', 'n_far', 'ratio_GNW', 'lat_best', 'lon_best', 'dep_best', 'nearest_family', 'dist_family_km',
            'head2_maxmag', 'n_far_excited', 'lfe_status', 'n_lfe_10m', 'm0_cols', 'i', 'epoch']
    ev[cols].to_csv(f'{K.TB}/e2_real_battery_events.csv', index=False)
    pd.set_option('display.width', 250)
    print(T.to_string(index=False))
    figs(D, T)


def figs(D, T):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # raw: MF level vs scan z for every T1-T3 detection, real vs controls, ETS / between, battery' ringed
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5), sharex=True, sharey=True)
    for a, per in enumerate(('ETS', 'between')):
        for b, (src, nm) in enumerate((('real', 'real'), ('A', 'scrambled controls (16)'), ('noise', 'noise templates (10)'))):
            ax = axs[a, b]
            x = D[(D.period == per) & (D.source.str.startswith(src))]
            sc = ax.scatter(x.zmad, x.z_best.clip(0.5, 100), c=x.ratio_far.clip(0.05, 3), cmap='viridis_r', s=5,
                            norm=matplotlib.colors.LogNorm(0.05, 3), lw=0)
            y = x[x.battery_e2]
            ax.scatter(y.zmad, y.z_best, s=40, facecolors='none', edgecolors='r', lw=1)
            ax.set_yscale('log'); ax.axhline(5, color='k', lw=0.5)
            ax.set_title(f'{per}: {nm} — {len(x)} dets, red ring = battery′ ({len(y)})', fontsize=9)
            if a == 1:
                ax.set_xlabel('MF level (MAD of the real template)')
            if b == 0:
                ax.set_ylabel('scan z_best (±30 s)')
    fig.colorbar(sc, ax=axs, label='ratio_far')
    fig.savefig(f'{K.FIG}/e2_mf_vs_scanz.png', dpi=120)
    # separation curves
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, ecol in zip(axs, ('e', "e'")):
        for agree, ls in ((False, '-'), (True, '--')):
            t = T[(T.earthquake_test == ecol) & (~T.S5) & (T.two_template == agree)]
            ax.plot(t.level_MAD, t.real_ETS, 'r' + ls, marker='o', label=f'real ETS{" 2-tpl" if agree else ""}')
            ax.plot(t.level_MAD, t.real_between, 'k' + ls, marker='o', label=f'real between{" 2-tpl" if agree else ""}')
            ax.plot(t.level_MAD, t.ctrl_mean_between, color='tab:purple', ls=ls, marker='s', label=f'control mean between{" 2-tpl" if agree else ""}')
            ax.plot(t.level_MAD, t.ctrl_max_between, color='tab:purple', ls=ls, marker='^', alpha=0.5, label=f'control max between{" 2-tpl" if agree else ""}')
        ax.set_xlabel('MF level (MAD)'); ax.set_ylabel('events passing battery'); ax.set_title(f'battery with ({ecol}), no S5')
        ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/e2_separation.png', dpi=130)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'rows':
        run_rows()
    elif mode == 'battery':
        run_battery(int(sys.argv[2]))
    elif mode == 'summary':
        run_summary()
