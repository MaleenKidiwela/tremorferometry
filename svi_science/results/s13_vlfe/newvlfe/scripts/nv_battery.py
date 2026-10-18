#!/usr/bin/env python
"""newvlfe stage 2 (PLAN.md §1-2): battery columns for every positive scan-v2 candidate (z >= 5) of one run.

usage: SCAN_YEAR=Y nv_battery.py Y RUN [--workers 2]
  RUN = real   -> data/scan_Y/runs/real_v2 (read only), variant 0
        polB   -> newvlfe/data/runs/Y/polB, variants 1..4 (polarity-scrambled)
        shA<s> -> newvlfe/data/runs/Y/shA<s>, variant 0 (station-day shift, all stations; references k + 30)
Stages: cheap (node table, time masks, gate = c AND e' time part), light (z_loso, f), full (scan_columns2.cols +
hg_r3.plane_wave2) for rows passing gate and f. Output newvlfe/data/bat_Y_RUN.parquet and live time
newvlfe/data/live_Y_RUN.json. Imports catalogue/hg_repro/scan code read-only; writes only under newvlfe/."""
import os, sys, json, time, argparse
os.environ['OMP_NUM_THREADS'] = '1'
S13 = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe'
NV = f'{S13}/newvlfe'
for p in (f'{S13}/scripts', f'{S13}/hg_repro/scripts', f'{S13}/catalogue/scripts'):
    sys.path.insert(0, p)
import numpy as np, pandas as pd

Y = int(os.environ['SCAN_YEAR'])
_SD = f'{S13}/data/scan' if Y == 2011 else f'{S13}/data/scan_{Y}'
os.environ.setdefault('SCAN_GRID', f'{_SD}/templates_grid_v2.npz')
os.environ['SCAN_GAIN_FALLBACK'] = '1'
import scan_common as S
import scan_columns2 as SC
import scan_masks2 as M2
import cat_common as K

BP = f'{S.SD}/bp_{Y}v2'
REFBP = f'{S.SD}/bp2_ref_v2' if Y in (2011, 2012) else f'{S.SD}/bp2_ref'     # production rule for 2011-12
LOG = f'{NV}/logs/nv_battery.log'
RUNDIR = {'realW': 'polBW'}
FARS = ('IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR', 'UW.OMAK', 'TA.I02D')                          # W2b re-run: real = variant 0 of the polBW pass


def log(m):
    with open(LOG, 'a') as f:
        f.write(f"{pd.Timestamp.now(tz='America/Los_Angeles').strftime('%F %H:%M:%S PT')} {Y} {m}\n")


def run_path(run):
    """run name as seen from S.SD/runs (scan_columns2.init reads f'{S.SD}/runs/{run}')."""
    if run == 'real':
        return 'real_v2'
    return f'../../../newvlfe/data/runs/{Y}/{RUNDIR.get(run, run)}'


def winit(run):
    SC.init(run_path(run), BP, REFBP)
    if SC.W['shift'] != 0:
        for r in SC.W['refs']:
            r['k'] = int(r['k']) + 30                     # PLAN §1: no reference shares a station's offset
    import hg_r3 as X
    X.EV['S'], X.EV['SC'] = S, SC


def light(v, tsec, node):
    """z_chk / z_loso at the node (first part of scan_columns2.cols, unchanged arithmetic)."""
    W = SC.W
    stas, offs = W['stas'], W['offs']
    d, u = divmod(int(tsec), 86400)
    sv = W['signs'][v]; g = W['g'][d]
    srcday = (d + offs) % S.NDAY
    nums = np.zeros(len(stas)); dens = np.zeros(len(stas))
    T = W['T']['thrust'][node]
    for i, s in enumerate(stas):
        ts = int(srcday[i]) * 86400 + u
        x = SC.cleaned_window(i, s, ts, tsec, S.FIT0, S.FIT1 + 1)
        sig = W['DS'][s]['sig'][srcday[i]] * 1e9
        for c in range(3):
            if W['won'][d, i, c] and np.isfinite(x[c]).all() and np.isfinite(sig[c]) and sig[c] > 0:
                w = 1.0 / sig[c] ** 2; Tt = T[i, c] * g[i, c]
                nums[i] += w * (Tt * x[c] * sv[i, c]).sum(); dens[i] += w * (Tt * Tt).sum()
    mad = W['mad'][v, node, d]
    if dens.sum() <= 0:
        return np.nan, np.nan
    kd = int(np.argmax(nums)); keep = np.arange(len(stas)) != kd
    zc = nums.sum() / np.sqrt(dens.sum()) / mad
    zl = nums[keep].sum() / np.sqrt(dens[keep].sum()) / mad if dens[keep].sum() > 0 else np.nan
    return zc, zl


def full(r):
    import hg_r3 as X
    v, ts, node = int(r['variant']), int(r['tsec']), int(r['node'])
    c1 = SC.cols(dict(variant=v, tsec=ts, node=node, m0=1e15))
    m0 = c1.get('m0_cols', np.nan)
    c2 = c1
    if r['need_true_m0'] and np.isfinite(m0):
        c2 = SC.cols(dict(variant=v, tsec=ts, node=node, m0=m0))
    o = {k: c2[k] for k in c2 if not k.startswith(('pred_', 'lag_'))}
    o['ex_at_true_m0'] = bool(r['need_true_m0'] and np.isfinite(m0))
    rf = [o[f'ratio_{t}'] for t in K.FAR5 if np.isfinite(o.get(f'ratio_{t}', np.nan))]
    r6 = [o[f'ratio_{t}'] for t in K.REF_ALL if np.isfinite(o.get(f'ratio_{t}', np.nan))]
    o['n_far'] = len(rf); o['ratio_far'] = float(np.median(rf)) if rf else np.nan
    o['n6'] = len(r6); o['ratio6'] = float(np.median(r6)) if r6 else np.nan
    o.update(X.plane_wave2(ts, node))
    o.update(emergence(ts, node))
    return o


def emergence(tsec, node):
    """Amendment NV3: catalogue Y9 emergence / amplitude-distance columns, each station on its own source day."""
    from scipy.signal import hilbert
    from scipy.stats import spearmanr
    W = SC.W
    d, u = divmod(int(tsec), 86400)
    nla, nlo = W['nodes'].lat[node], W['nodes'].lon[node]
    tr = []
    for i, s in enumerate(W['stas']):
        if not W['won'][d, i, 2]:
            continue
        ts = int((d + W['offs'][i]) % S.NDAY) * 86400 + u
        X = W['BP'][s]
        if ts - 1800 < 0 or ts + 1800 > X.shape[1]:
            continue
        tr.append((s, np.asarray(X[2, ts - 1800:ts + 1800], float) * 1e9, float(S.dist_km(nla, nlo, *W['sll'][i])), False))
    for ref in W['refs']:
        if ref['sta'] in W['stas']:
            continue
        dr = int((d + W['shift'] * ref['k']) % S.NDAY); ts = dr * 86400 + u
        X = ref['X']
        if ts - 1800 < 0 or ts + 1800 > X.shape[1]:
            continue
        tr.append((ref['sta'], np.asarray(X[2, ts - 1800:ts + 1800], float) * 1e9, float(S.dist_km(nla, nlo, ref['lat'], ref['lon'])),
                   ref['sta'] in FARS))
    em, amp, dist, isfar = {}, {}, {}, {}
    for s, x, dd, fa in tr:
        if np.isfinite(x).sum() < 0.9 * len(x) or not np.isfinite(x[1800:1980]).all():
            continue
        x = np.nan_to_num(x - np.nanmean(x)); env = np.abs(hilbert(x)); med = float(np.median(env))
        pk = float(np.max(np.abs(x[1800:1980])))
        amp[s] = pk; em[s] = pk / med if med > 0 else np.nan; dist[s] = dd; isfar[s] = fa
    nf = [s for s in em if not isfar[s]]
    netw = [s for s in nf if dist[s] <= 150]
    ev = np.array([em[s] for s in netw]); allnf = sorted([em[s] for s in nf if np.isfinite(em[s])], reverse=True)
    far = [s for s in em if isfar[s]]
    nr = [s for s in nf if dist[s] <= 40]
    if len(nr) < 2:
        nr = sorted(nf, key=dist.get)[:2]
    mid = [s for s in nf if 40 < dist[s] <= 150 and s not in nr]
    r1 = np.median([amp[s] for s in nr]) / np.median([amp[s] for s in mid]) if nr and mid else np.nan
    r2 = np.median([amp[s] for s in far]) / max(amp[s] for s in nr) if far and nr else np.nan
    allS = nr + mid + far
    rho = spearmanr([dist[s] for s in allS], [amp[s] for s in allS]).correlation if len(allS) >= 4 else np.nan
    return dict(emerg_med_near=float(np.nanmedian(ev)) if len(ev) else np.nan, n_near_emerg_ge2=int((ev >= 2).sum()), n_near=len(ev),
                emerg_max_over_second=allnf[0] / allnf[1] if len(allnf) > 1 and allnf[1] > 0 else np.nan,
                emerg_lead_sta=max(nf, key=lambda s: em[s]) if nf else '',
                emerg_far_med=float(np.nanmedian([em[s] for s in far])) if far else np.nan,
                t3_near_mid=r1, t3_far_over_nearmax=r2, t3_spearman=rho)


def chunk_light(rows):
    return [(i, *light(v, t, n)) for i, v, t, n in rows]


def chunk_full(rows):
    out = []
    for i, r in rows:
        try:
            out.append(dict(idx=i, **full(r)))
        except Exception as e:  # noqa
            out.append(dict(idx=i, col_error=f'{type(e).__name__}: {str(e)[:80]}'))
    return out


def anyw(mk, t):
    c = np.concatenate([[0], np.cumsum(mk)]); n = len(mk)
    a = np.clip(t - 60, 0, n); b = np.clip(t + 121, 0, n)
    return (c[b] - c[a]) > 0


def head2_m55_mask():
    A = pd.read_csv(S.tab('scan2_eq_events.csv'))
    A = A[A.head2_t0.notna() & (A.mag >= 5.5)]
    y0 = (S.Y0 - pd.Timestamp('1970-01-01')).total_seconds()
    n = S.NDAY * 86400
    m = np.zeros(n, bool)
    for a, b in zip(A.head2_t0.values - y0, A.head2_t1.values - y0):
        a, b = int(max(a, 0)), int(min(b, n))
        if b > a:
            m[a:b] = True
    return m


def period_vec(tsec):
    """ETS = SVI episode (start-2 d, end+3 d); other = NW-only episode; else between (production rule)."""
    y0 = S.Y0
    t = y0 + pd.to_timedelta(np.asarray(tsec), unit='s')
    out = np.full(len(t), 'between', dtype=object)
    E = pd.read_csv(f'{K.ROOT}/svi_science/results/b_ets_record/episodes_b.csv')
    for reg, lab in (('NW', 'other'), ('SVI', 'ETS')):
        for r in E[E.region == reg].itertuples():
            a = pd.Timestamp(r.start) - pd.Timedelta(days=2); b = pd.Timestamp(r.end) + pd.Timedelta(days=3)
            out[(t >= a) & (t < b)] = lab
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('year', type=int); ap.add_argument('run')
    ap.add_argument('--workers', type=int, default=2); ap.add_argument('--limit', type=int, default=0)
    a = ap.parse_args()
    assert a.year == Y
    t0 = time.time()
    rdir = f'{S.SD}/runs/{run_path(a.run)}'
    Cd = pd.read_parquet(f'{rdir}/cand.parquet')
    Cd = Cd[(Cd.sign > 0) & (Cd.z >= 5)]
    Cd = Cd[Cd.variant >= 1] if a.run.startswith('polB') else Cd[Cd.variant == 0]
    Cd = Cd.reset_index(drop=True)
    if a.limit:
        Cd = Cd.sample(a.limit, random_state=0).sort_values('tsec').reset_index(drop=True)
    N = K.node_table()
    nb = Cd.node.values
    Cd['lat'] = N.lat.values[nb]; Cd['lon'] = N.lon.values[nb]; Cd['dep'] = N.dep.values[nb]
    Cd['grid_edge_node'] = N.grid_edge.values[nb]; Cd['in_footprint_poly'] = N.fp_interior.values[nb]
    Cd['dist_family_km'] = N.dist_family_km.values[nb]; Cd['nearest_family'] = N.nearest_family.values[nb]
    t = Cd.tsec.values.astype(int)
    for nm in ('head2', 'bw5', 'reg', 'reg3', 'm45', 'v2'):
        Cd[f'in_{nm}'] = anyw(M2.mask(nm), t)
    y0e = (S.Y0 - pd.Timestamp('1970-01-01')).total_seconds()
    Cd['epoch'] = y0e + t
    Cd['head2_maxmag'] = K.head2_maxmag(Cd.epoch.values, np.full(len(Cd), Y))
    Cd['z_best'] = Cd.z
    pc = (Cd.z_best >= 5) & ~Cd.grid_edge_node & (Cd.in_footprint_poly | (Cd.dist_family_km <= 20))
    e2t = ~(Cd.in_bw5 | Cd.in_reg | Cd.in_reg3 | (Cd.in_head2 & (Cd.head2_maxmag >= 5.5).fillna(False)))
    Cd['gate'] = pc & e2t
    Cd['period'] = period_vec(t)
    Cd['need_true_m0'] = Cd.in_head2
    log(f'{a.run}: {len(Cd)} rows z>=5, gate {int(Cd.gate.sum())}')
    from multiprocessing import Pool
    G = Cd[Cd.gate]
    rows = list(zip(G.index, G.variant, G.tsec, G.node))
    res = []
    if rows:
        nch = max(1, min(len(rows), a.workers * 16)); ch = [rows[k::nch] for k in range(nch)]
        with Pool(a.workers, initializer=winit, initargs=(a.run,)) as p:
            for r in p.imap_unordered(chunk_light, ch):
                res += r
    L = pd.DataFrame(res, columns=['idx', 'z_chk_light', 'z_loso_light']).set_index('idx')
    Cd = Cd.join(L)
    Cd['pass_f_light'] = (Cd.z_loso_light >= 3) & (Cd.z_loso_light >= 0.5 * Cd.z_best)
    log(f'{a.run}: light done, f-pass {int((Cd.gate & Cd.pass_f_light).sum())} ({time.time() - t0:.0f} s)')
    G = Cd[Cd.gate & Cd.pass_f_light]
    rows = [(i, dict(variant=r.variant, tsec=r.tsec, node=r.node, need_true_m0=bool(r.need_true_m0))) for i, r in G.iterrows()]
    res = []
    if rows:
        nch = max(1, min(len(rows), a.workers * 16)); ch = [rows[k::nch] for k in range(nch)]
        with Pool(a.workers, initializer=winit, initargs=(a.run,)) as p:
            for k, r in enumerate(p.imap_unordered(chunk_full, ch)):
                res += r
                if k % 8 == 0:
                    log(f'{a.run}: full chunk {k + 1}/{nch} ({time.time() - t0:.0f} s)')
    if res:
        F = pd.DataFrame(res).set_index('idx')
        F = F.drop(columns=[c for c in F.columns if c in Cd.columns])
        Cd = Cd.join(F)
    Cd['run'] = a.run; Cd['year'] = Y
    Cd.to_parquet(f'{NV}/data/bat_{Y}_{a.run}' + ('_test' if a.limit else '') + '.parquet')
    # live time per period: navail >= 3 and outside the e' time masks
    nav = np.load(f'{rdir}/navail.npy')
    mk = M2.mask('bw5') | M2.mask('reg') | M2.mask('reg3') | head2_m55_mask()
    live = (nav >= 3) & ~mk[:len(nav)]
    nd = len(nav) // 86400
    per = np.repeat(period_vec(np.arange(nd) * 86400 + 43200), 86400)[:len(nav)]   # episodes are whole days
    nvar = 4 if a.run.startswith('polB') else 1
    out = {p: float(live[per == p].sum()) / (365.25 * 86400) for p in ('ETS', 'between', 'other')}
    json.dump(dict(year=Y, run=a.run, n_catalogues=nvar, live_yr_per_catalogue=out), open(f'{NV}/data/live_{Y}_{a.run}.json', 'w'), indent=1)
    log(f'{a.run}: done {len(Cd)} rows, live {out}, {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
