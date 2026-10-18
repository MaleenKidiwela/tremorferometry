#!/usr/bin/env python
"""Amendment Q1: 2011-12 leave-one-out re-test on station sets hg9 / era2007 / era2016 / max2011, plain vs weighted sum.
Uses the E3 machinery (cat_e3) with an extended station list; E3 alignment shifts reused for the stacks.
  bank      -> data/q1/bank.npz (+ weights), tables/q1_bank.csv
  mf        -> data/q1/mf/NAME__SET.npz (plain 'real'/'A1_r'/'A2_r' and weighted 'w_real'/'w_A1_r'/'w_A2_r' peaks)
  detect    -> data/q1/rows.parquet, tables/q1_thresholds.csv, q1_refind.csv
  battery Y -> data/q1/bat_Y.parquet (SCAN_YEAR=Y)
  summary   -> tables/q1_summary.csv, q1_verdict.json, figures/q1_loo.png
Sets in 2011-12 (station-era table; SQM/HDW/OTR are not broadband in 2011-12 -> absent):
  era2007 -> {PGC, VGZ, SNB, NLLB, GNW}; era2016 -> {PGC, VGZ, LRIV, SNB, PFB, NLLB, DOSE, GNW, SP2};
  max2011 -> broadband with data on >= 7 of the 8 vetted-event days (17 stations)."""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
import cat_e3 as E3
import hg_common as H

Q1 = f'{K.DATA}/q1'
MFD = f'{Q1}/mf'
os.makedirs(MFD, exist_ok=True)
HG9 = E3.HG9
DENSE = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.YOUB', 'CN.LZB', 'UW.LRIV', 'UW.DOSE', 'UW.GNW', 'UW.SP2',
         'TA.A04D', 'TA.B05D', 'TA.D03D', 'TA.D04D', 'UW.WISH', 'US.NLWA', 'CC.OBSR', 'UW.ELW', 'UW.LEBA', 'UW.RATT']
ALL = DENSE + ['CN.PFB']
SETS = {'dense2011': DENSE,                                                 # Amendment Q1' sets
        'era2007': ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'UW.GNW'],
        'era2010': ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'UW.GNW', 'UW.LRIV', 'CN.PFB'],
        'era2016': ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'UW.LRIV', 'UW.DOSE', 'UW.GNW', 'UW.SP2', 'UW.WISH', 'US.NLWA', 'CN.PFB'],
        'hg9': HG9, 'cn4': ['CN.PGC', 'CN.VGZ', 'CN.NLLB', 'CN.SNB']}
CH = [f'{s}.{c}' for s in ALL for c in 'NEZ']
NEWCONT = ['CN.LZB', 'UW.WISH', 'US.NLWA', 'CC.OBSR', 'TA.D04D']
# point the E3 machinery at the extended station list
E3.ALL, E3.CH, E3.SETS = ALL, CH, SETS
NDRAW, NNOISE, PRE = E3.NDRAW, E3.NNOISE, E3.PRE


def rms_ratio(start):
    """per channel: RMS of the 90-s window / RMS of the 10 min before it (raw velocity), NaN if not computable."""
    out = {}
    for ch in CH:
        s, c = ch.rsplit('.', 1); ci = 'NEZ'.index(c)
        X, D = E3.cont(s)
        w = np.array(X[ci, start:start + H.L], np.float64); pre = np.array(X[ci, start - 600:start], np.float64)
        ok = np.isfinite(w).all() and np.isfinite(pre).sum() > 300
        out[ch] = float(np.sqrt(np.mean(w ** 2)) / np.sqrt(np.nanmean(pre ** 2))) if ok and np.nanmean(pre ** 2) > 0 else np.nan
    return out


def run_bank():
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    TT = H.templates().set_index('template')
    NT = pd.read_csv(f'{H.MF}/noise_times.csv')
    bank, wts, meta = {}, {}, []

    def add(name, starts, kind, members, extra=None):
        Wn = [E3.cut(int(s)) for s in starts]
        st = E3.stack(Wn) if len(Wn) > 1 else Wn[0]
        rr = [rms_ratio(int(s)) for s in starts]
        w = {ch: float(np.nanmean([r[ch] for r in rr])) if np.isfinite([r[ch] for r in rr]).any() else np.nan for ch in CH}
        bank[name] = st; wts[name] = w
        meta.append(dict(name=name, kind=kind, members=members, n_members=len(starts),
                         n_ch=sum(v is not None for v in st.values()), primary=kind in ('cluster', 'hg_recut'), **(extra or {})))
    for c, g in V.groupby('cluster'):
        starts = (g.start0 + g.shift_s).astype(int).values
        add(f'C{c}', starts, 'cluster', ','.join(map(str, g.cand_id)))
        for k, r in enumerate(g.itertuples()):
            if r.role == 'vetted' and len(g) >= 2:
                idx = [j for j in range(len(g)) if j != k]
                add(f'C{c}_loo{r.cand_id}', starts[idx], 'loo', ','.join(str(g.cand_id.values[j]) for j in idx),
                    dict(left_out=r.cand_id))
    for k in (2, 4):
        add(f'HG_T{k}', [H.tidx(TT.loc[k].time) - PRE], 'hg_recut', '')
    for j in range(NNOISE):
        add(f'N{j:02d}', [int(NT.i[j])], 'noise', str(int(NT.i[j])), dict(i0=int(NT.i[j])))
    arr = {nm: np.array([st[ch] if st[ch] is not None else np.full(H.L, np.nan) for ch in CH]) for nm, st in bank.items()}
    warr = {f'w__{nm}': np.array([wts[nm][ch] for ch in CH]) for nm in bank}
    np.savez(f'{Q1}/bank.npz', chans=np.array(CH), **arr, **warr)
    M = pd.DataFrame(meta); M.to_csv(f'{K.TB}/q1_bank.csv', index=False)
    print(M.to_string())


def mf_job(job):
    name, setname, controls = job
    f = f'{MFD}/{name}__{setname}.npz'
    if os.path.isfile(f):
        return f'{name} {setname} exists'
    B = np.load(f'{Q1}/bank.npz')
    W = B[name]; wv = B[f'w__{name}']
    chs = [ch for k, ch in enumerate(CH) if ch.rsplit('.', 1)[0] in SETS[setname] and np.isfinite(W[k]).all()
           and np.isfinite(wv[k])]
    wk = np.array([wv[CH.index(ch)] for ch in chs]); wk = wk * len(wk) / wk.sum() if len(wk) else wk
    names = ['real'] + ([f'A1_{r}' for r in range(NDRAW)] + [f'A2_{r}' for r in range(NDRAW)] if controls else [])
    S = {n: np.zeros(E3.MFN, np.float32) for n in names}
    SW = {n: np.zeros(E3.MFN, np.float32) for n in names}
    NL = {n: np.zeros(E3.MFN, np.int16) for n in names}
    offs, signs = draws18()
    for j, ch in enumerate(chs):
        k = CH.index(ch); si = ALL.index(ch.rsplit('.', 1)[0])
        cc, live = E3.ccser(ch, W[k])
        ww = np.float32(wk[j])
        S['real'] += cc; SW['real'] += ww * cc; NL['real'] += live
        if controls:
            for r in range(NDRAW):
                sh = int(offs[r, si]) * 86400
                rc = np.roll(cc, sh)
                S[f'A1_{r}'] += rc; SW[f'A1_{r}'] += ww * rc; NL[f'A1_{r}'] += np.roll(live, sh)
                sg = np.float32(signs[r, k])
                S[f'A2_{r}'] += sg * cc; SW[f'A2_{r}'] += (sg * ww) * cc; NL[f'A2_{r}'] += live
        del cc, live
    out = {'name': name, 'set': setname, 'chans': np.array(chs), 'weights': wk}
    for n in names:
        for tag, A in (('', S), ('w_', SW)):
            pi, pv, pn, med, mad = E3.peaks(A[n], NL[n])
            out.update({f'{tag}{n}_pi': pi, f'{tag}{n}_pv': pv, f'{tag}{n}_med': med, f'{tag}{n}_mad': mad})
            A[n] = None
        NL[n] = None
    np.savez_compressed(f + '.part.npz', **out); os.replace(f + '.part.npz', f)
    return f'{name} {setname} done ({len(chs)} ch, controls {controls})'


def draws18():
    rng = np.random.default_rng(E3.SEED + 18)
    offs = []
    for _ in range(NDRAW):
        while True:
            o = rng.choice(np.arange(30, E3.NDF - 30), size=len(ALL), replace=False)
            if np.min(np.diff(np.sort(o))) >= 2:
                break
        offs.append(o.astype(int))
    signs = []
    for _ in range(NDRAW):
        while True:
            sg = rng.choice([-1.0, 1.0], size=len(CH))
            if 0 < (sg > 0).sum() < len(sg):
                break
        signs.append(sg)
    return np.array(offs), np.array(signs)


def run_mf():
    M = pd.read_csv(f'{K.TB}/q1_bank.csv')
    jobs = [(r.name, s, bool(r.primary)) for r in M.itertuples() for s in SETS]
    jobs.sort(key=lambda j: (not j[2], j[1] != 'dense2011'))
    from multiprocessing import Pool
    with Pool(6, maxtasksperchild=3) as p:
        for k, r in enumerate(p.imap_unordered(mf_job, jobs)):
            K.log(f'q1 mf {k + 1}/{len(jobs)}: {r}', 'q1')


def run_detect():
    M = pd.read_csv(f'{K.TB}/q1_bank.csv')
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    ndays = {2011: (pd.Timestamp('2012-01-01') - H.T0).days, 2012: 366}
    noise = M[M.kind == 'noise']
    TH, rows, refind = [], [], []
    for s in SETS:
        for stat, tag in (('plain', ''), ('weighted', 'w_')):
            zs = {2011: [], 2012: []}
            for r in noise.itertuples():
                z = np.load(f'{MFD}/{r.name}__{s}.npz')
                med, mad = float(z[f'{tag}real_med']), float(z[f'{tag}real_mad'])
                pi, pv = z[f'{tag}real_pi'], z[f'{tag}real_pv']
                m = np.abs(pi - int(r.i0)) > 90; yy = E3.year_of(pi)
                for y in (2011, 2012):
                    zs[y].append((pv[m & (yy == y)] - med) / mad)
            zstar = {}
            for y in (2011, 2012):
                allz = np.sort(np.concatenate(zs[y]))[::-1]
                for rate in (0.25, 0.5, 1.0):
                    n = int(round(rate * len(noise) * ndays[y]))
                    zstar[(y, rate)] = float(allz[n - 1]) if 0 < n <= len(allz) else float(allz[-1])
                    TH.append(dict(set=s, stat=stat, year=y, rate=rate, z_star=zstar[(y, rate)]))
            for r in M[M.primary].itertuples():
                z = np.load(f'{MFD}/{r.name}__{s}.npz')
                med, mad = float(z[f'{tag}real_med']), float(z[f'{tag}real_mad'])
                for src in ['real'] + [f'A1_{q}' for q in range(NDRAW)] + [f'A2_{q}' for q in range(NDRAW)]:
                    pi, pv = z[f'{tag}{src}_pi'], z[f'{tag}{src}_pv']
                    zz = (pv - med) / mad; yy = E3.year_of(pi)
                    m = zz >= np.where(yy == 2011, zstar[(2011, 0.5)], zstar[(2012, 0.5)])
                    if m.any():
                        rows.append(pd.DataFrame(dict(set=s, stat=stat, template=r.name, source=src if src == 'real' else src[:2],
                                                      draw=-1 if src == 'real' else int(src[3:]), p=pi[m], zmad=zz[m])))
            for e in V[V.role == 'vetted'].itertuples():
                ye = int(E3.year_of([e.i])[0])
                best = (-np.inf, '')
                for t in M[M.kind.isin(['cluster', 'loo', 'hg_recut'])].itertuples():
                    if str(e.cand_id) in str(t.members).split(','):
                        continue
                    if t.kind == 'loo' and int(float(t.left_out)) != e.cand_id:
                        continue
                    z = np.load(f'{MFD}/{t.name}__{s}.npz')
                    med, mad = float(z[f'{tag}real_med']), float(z[f'{tag}real_mad'])
                    pi, pv = z[f'{tag}real_pi'], z[f'{tag}real_pv']
                    m = np.abs(pi + PRE - e.i) <= 30
                    if m.any():
                        zz = float(((pv[m] - med) / mad).max())
                        if zz > best[0]:
                            best = (zz, t.name)
                refind.append(dict(set=s, stat=stat, cand_id=e.cand_id, time_pt=e.time_pt, best_template=best[1], best_z=best[0],
                                   z_star=zstar[(ye, 0.5)], refound=bool(best[0] >= zstar[(ye, 0.5)])))
    R = pd.concat(rows, ignore_index=True)
    R['i'] = R.p + PRE
    R['epoch'] = (H.T0 - pd.Timestamp('1970-01-01')).total_seconds() + R.i
    R['year'] = pd.to_datetime(R.epoch, unit='s').dt.year
    R['period'] = K.period_of(R.epoch.values)
    R.to_parquet(f'{Q1}/rows.parquet')
    pd.DataFrame(TH).to_csv(f'{K.TB}/q1_thresholds.csv', index=False)
    F = pd.DataFrame(refind); F.to_csv(f'{K.TB}/q1_refind.csv', index=False)
    print(F.groupby(['set', 'stat']).refound.sum().unstack())
    K.log(f'q1 detect: rows {len(R)}, unique times {R.i.nunique()}', 'q1')


def run_battery(year):
    import cat_e2
    cat_e2.run_battery(year, rows_file=f'{Q1}/rows.parquet', out_file=f'{Q1}/bat_{year}.parquet', tag='q1')


def run_summary():
    R = pd.read_parquet(f'{Q1}/rows.parquet')
    B = pd.concat([pd.read_parquet(f'{Q1}/bat_{y}.parquet') for y in (2011, 2012)], ignore_index=True).drop(columns=['epoch'])
    B['grid_edge'] = B.grid_edge.where(B.grid_edge.notna(), B.grid_edge_node) if 'grid_edge' in B else B.grid_edge_node
    D = R.merge(B, on='i', how='left')
    for c in ('n_far', 'n6'):
        D[c] = D[c].fillna(0)
    D = K.apply_rules(D)
    D.to_parquet(f'{Q1}/q1_rows_scored.parquet')
    F = pd.read_csv(f'{K.TB}/q1_refind.csv')
    yb = K.BETWEEN_YR
    out = []
    for (s, stat), g in D[D.battery_e2].groupby(['set', 'stat']):
        o = dict(set=s, stat=stat)
        for per in ('ETS', 'between'):
            y = g[g.period == per]
            o[f'real_{per}'] = E3.dedupe_count(y[y.source == 'real'].i)
            cc = [E3.dedupe_count(y[(y.source == a) & (y.draw == q)].i) for a in ('A1', 'A2') for q in range(NDRAW)]
            o[f'ctrl_mean_{per}'] = float(np.mean(cc)); o[f'ctrl_max_{per}'] = int(np.max(cc))
        o['ctrl_rate_between_per_yr'] = o['ctrl_mean_between'] / yb
        f = F[(F.set == s) & (F.stat == stat)]
        o['loo_refound'] = int(f.refound.sum()); o['loo_n'] = len(f)
        out.append(o)
    T = pd.DataFrame(out)
    for s in SETS:                                   # sets with no battery passes at all
        for stat in ('plain', 'weighted'):
            if not ((T.set == s) & (T.stat == stat)).any():
                f = F[(F.set == s) & (F.stat == stat)]
                T.loc[len(T)] = dict(set=s, stat=stat, loo_refound=int(f.refound.sum()), loo_n=len(f), ctrl_rate_between_per_yr=0.0)
    T = T.sort_values(['set', 'stat'])
    T.to_csv(f'{K.TB}/q1_summary.csv', index=False)
    get = lambda s, st, c: float(T[(T.set == s) & (T.stat == st)][c].iloc[0])
    V = {}
    V['weighted_ge_plain'] = {s: get(s, 'weighted', 'loo_refound') >= get(s, 'plain', 'loo_refound') for s in SETS}
    stat_used = 'weighted' if (V['weighted_ge_plain']['hg9'] and V['weighted_ge_plain']['dense2011']) else 'plain'
    V['stat_for_later'] = stat_used
    V['hg9_ge7'] = get('hg9', stat_used, 'loo_refound') >= 7
    V['dense2011_pass'] = get('dense2011', stat_used, 'loo_refound') >= 8 and get('dense2011', stat_used, 'ctrl_rate_between_per_yr') <= 1
    V['era2016_pass'] = get('era2016', stat_used, 'loo_refound') >= 7 and get('era2016', stat_used, 'ctrl_rate_between_per_yr') <= 1
    V['era2007_pass'] = get('era2007', stat_used, 'loo_refound') >= 5 and get('era2007', stat_used, 'ctrl_rate_between_per_yr') <= 1
    json.dump(V, open(f'{K.TB}/q1_verdict.json', 'w'), indent=1)
    pd.set_option('display.width', 250)
    print(T.round(2).to_string(index=False)); print(V)
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4))
    xs = {s: k for k, s in enumerate(SETS)}
    for stat, mk, dx in (('plain', 'o', -0.12), ('weighted', 's', 0.12)):
        f = F[F.stat == stat]
        ax.scatter([xs[s] + dx for s in f.set], f.best_z - f.z_star, marker=mk, label=stat, alpha=0.7)
    ax.axhline(0, color='k', lw=0.6); ax.set_xticks(range(len(SETS))); ax.set_xticklabels(list(SETS))
    ax.set_ylabel('best leave-one-out MF z minus threshold (MAD)'); ax.legend(); fig.tight_layout()
    fig.savefig(f'{K.FIG}/q1_loo.png', dpi=130)


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'battery':
        run_battery(int(sys.argv[2]))
    else:
        {'bank': run_bank, 'mf': run_mf, 'detect': run_detect, 'summary': run_summary}[m]()
