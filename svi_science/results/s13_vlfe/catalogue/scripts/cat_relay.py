#!/usr/bin/env python
"""Amendment T2: relay (era-matched) template banks for a production year. Run with PROD_YEAR=Y PILOT_SFX=_v2 and the scan-v2
environment of that year (SCAN_YEAR, SCAN_RUN, SCAN_BP, SCAN_GRID, SCAN_GAIN_FALLBACK).
  bank     -> data/prodY/relay/bank.npz, tables/prodY_relay_bank.csv, prodY_relay_events.csv
  mf       -> data/prodY/relay/mf/NAME.npz (weighted; primary relay templates with 8 + 8 scrambled controls)
  loo      -> tables/prodY_relay_loo.csv; data/prodY/relay/rows.parquet (tier candidates of the relay bank)
  battery  -> data/prodY/relay/bat.parquet
  verdict  -> tables/prodY_relay_verdict.json"""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
import cat_pilot as P
import cat_e3 as E3

Y = P.PILOT
RD = f'{P.PD}/relay'
MFD = f'{RD}/mf'
os.makedirs(MFD, exist_ok=True)


def stack_frame(starts):
    W = [P.cut_frame(int(s)) for s in starts]
    out = {}
    for ch in P.CH:
        ws = [w[ch] for w in W if w[ch] is not None]
        if not ws or len(ws) < 0.5 * len(W):
            out[ch] = None; continue
        s = np.mean(ws, 0); s = s - s.mean(); n = np.linalg.norm(s); out[ch] = s / n if n > 0 else None
    return out, W


def sumcc(a, b):
    return sum(float(a[c] @ b[c]) for c in P.CH if a[c] is not None and b[c] is not None)


def run_bank():
    E = pd.read_csv(f'{K.TB}/{P.NAME}{Y}_events{P.SFX}.csv')
    V = E[(E.period == 'ETS') & (E.lfe_status == 'with LFEs') & (E.vr_thrust >= 0.12)].copy().reset_index(drop=True)
    V['rid'] = [f'R{k}' for k in range(len(V))]
    V['start0'] = V.i - P.PRE
    if len(V) < 3:
        V.to_csv(f'{K.TB}/{P.NAME}{Y}_relay_events.csv', index=False)
        json.dump(dict(year=Y, n_vetted=len(V), accepted=False, reason='< 3 vetted relay events'),
                  open(f'{K.TB}/{P.NAME}{Y}_relay_verdict.json', 'w'), indent=1)
        print('fewer than 3 vetted relay events'); return
    V['cluster'] = E3.complete_linkage(V.lat_best.values, V.lon_best.values, 15.0)
    bank, wts, meta = {}, {}, []
    V['shift_s'] = 0
    for c, g in V.groupby('cluster'):
        starts = list(g.start0.astype(int).values)
        st, W = stack_frame(starts)
        if len(starts) >= 2:                                      # +-5 s alignment to the stack of the others (2 passes)
            sh = np.zeros(len(starts), int)
            for _ in range(2):
                for k, s0 in enumerate(starts):
                    ref, _w = stack_frame([starts[j] + sh[j] for j in range(len(starts)) if j != k])
                    sh[k] = max(range(-5, 6), key=lambda d: sumcc(P.cut_frame(s0 + d), ref))
            starts = [s0 + d for s0, d in zip(starts, sh)]
            V.loc[g.index, 'shift_s'] = sh
            st, W = stack_frame(starts)
        rr = [P.rms_ratio_frame(s) for s in starts]
        name = f'RC{c}'
        bank[name] = st
        wts[name] = {ch: float(np.nanmean([r[ch] for r in rr])) if np.isfinite([r[ch] for r in rr]).any() else np.nan for ch in P.CH}
        meta.append(dict(name=name, kind='relay', members=','.join(g.rid), n_members=len(g), primary=True, off=P.PRE))
        for k, r in enumerate(g.itertuples()):
            if len(g) < 2:
                continue
            idx = [j for j in range(len(starts)) if j != k]
            st2, _ = stack_frame([starts[j] for j in idx])
            nm = f'RC{c}_loo{r.rid}'
            bank[nm] = st2
            wts[nm] = {ch: float(np.nanmean([rr[j][ch] for j in idx])) if np.isfinite([rr[j][ch] for j in idx]).any() else np.nan for ch in P.CH}
            meta.append(dict(name=nm, kind='relay_loo', members=','.join(g.rid.values[j] for j in idx), n_members=len(idx),
                             primary=False, off=P.PRE, left_out=r.rid))
    arr = {nm: np.array([st[ch] if st[ch] is not None else np.full(P.L, np.nan) for ch in P.CH]) for nm, st in bank.items()}
    warr = {f'w__{nm}': np.array([wts[nm][ch] for ch in P.CH]) for nm in bank}
    np.savez(f'{RD}/bank.npz', chans=np.array(P.CH), **arr, **warr)
    pd.DataFrame(meta).to_csv(f'{K.TB}/{P.NAME}{Y}_relay_bank.csv', index=False)
    V.to_csv(f'{K.TB}/{P.NAME}{Y}_relay_events.csv', index=False)
    print(V[['rid', 'time_pt', 'cluster', 'lat_best', 'lon_best', 'zmad', 'shift_s']].to_string()); print(pd.DataFrame(meta).to_string())


def mf_job(job):
    name, controls = job
    f = f'{MFD}/{name}.npz'
    if os.path.isfile(f):
        return f'{name} exists'
    B = np.load(f'{RD}/bank.npz'); W = B[name]; wv = B[f'w__{name}']
    chs = [ch for k, ch in enumerate(P.CH) if np.isfinite(W[k]).all() and np.isfinite(wv[k])]
    wk = np.array([wv[P.CH.index(ch)] for ch in chs]); wk = wk * len(wk) / wk.sum()
    names = ['real'] + ([f'A1_{r}' for r in range(P.NDRAW)] + [f'A2_{r}' for r in range(P.NDRAW)] if controls else [])
    S = {n: np.zeros(P.M, np.float32) for n in names}; NL = {n: np.zeros(P.M, np.int16) for n in names}
    offs, signs = P.draws()
    for j, ch in enumerate(chs):
        k = P.CH.index(ch); si = P.STAS.index(ch.rsplit('.', 1)[0]); ww = np.float32(wk[j])
        cc, live = P.ccser(ch, W[k])
        S['real'] += ww * cc; NL['real'] += live
        if controls:
            for r in range(P.NDRAW):
                sh = int(offs[r, si]) * 86400
                S[f'A1_{r}'] += ww * np.roll(cc, sh); NL[f'A1_{r}'] += np.roll(live, sh)
                S[f'A2_{r}'] += (ww * np.float32(signs[r, k])) * cc; NL[f'A2_{r}'] += live
    out = {'name': name, 'chans': np.array(chs), 'weights': wk}
    for n in names:
        pi, pv, pn, med, mad = E3.peaks(S[n], NL[n])
        out.update({f'{n}_pi': pi, f'{n}_pv': pv, f'{n}_med': med, f'{n}_mad': mad}); S[n] = None; NL[n] = None
    np.savez_compressed(f + '.part.npz', **out); os.replace(f + '.part.npz', f)
    return f'{name} done ({len(chs)} ch)'


def run_mf():
    Mt = pd.read_csv(f'{K.TB}/{P.NAME}{Y}_relay_bank.csv')
    jobs = sorted([(r.name, bool(r.primary)) for r in Mt.itertuples()], key=lambda j: not j[1])
    from multiprocessing import Pool
    with Pool(P.NW, maxtasksperchild=3) as p:
        for k, r in enumerate(p.imap_unordered(mf_job, jobs)):
            K.log(f'relay {Y} mf {k + 1}/{len(jobs)}: {r}', 'relay')


def run_loo():
    Mt = pd.read_csv(f'{K.TB}/{P.NAME}{Y}_relay_bank.csv')
    V = pd.read_csv(f'{K.TB}/{P.NAME}{Y}_relay_events.csv')
    TH = pd.read_csv(f'{K.TB}/{P.NAME}{Y}_noise_threshold.csv').set_index('rate').z_star
    zs = TH[0.5]
    out = []
    for e in V.itertuples():
        best = (-np.inf, '')
        for t in Mt.itertuples():
            mem = str(t.members).split(',')
            if e.rid in mem:
                continue
            if t.kind == 'relay_loo' and t.left_out != e.rid:
                continue
            z = np.load(f'{MFD}/{t.name}.npz')
            pi, pv = z['real_pi'], z['real_pv']
            m = np.abs(pi + P.PRE - e.i) <= 30
            if m.any():
                v = float(((pv[m] - float(z['real_med'])) / float(z['real_mad'])).max())
                if v > best[0]:
                    best = (v, t.name)
        out.append(dict(rid=e.rid, time_pt=e.time_pt, cluster=e.cluster, best_template=best[1], best_z=best[0], z_star=zs,
                        refound=bool(best[0] >= zs)))
    F = pd.DataFrame(out); F.to_csv(f'{K.TB}/{P.NAME}{Y}_relay_loo.csv', index=False)
    # tier candidates of the relay bank (primary relay templates; >= 5.5 MAD; agreement among relay + 2011 bank rows)
    rows = []
    for t in Mt[Mt.primary].itertuples():
        z = np.load(f'{MFD}/{t.name}.npz'); med, mad = float(z['real_med']), float(z['real_mad'])
        for src in ['real'] + [f'A1_{q}' for q in range(P.NDRAW)] + [f'A2_{q}' for q in range(P.NDRAW)]:
            zz = (z[f'{src}_pv'] - med) / mad; m = zz >= P.LEVEL
            if m.any():
                rows.append(pd.DataFrame(dict(template=t.name, kind='relay', source=src if src == 'real' else src[:2],
                                              draw=-1 if src == 'real' else int(src[3:]), i=z[f'{src}_pi'][m] + P.PRE, zmad=zz[m])))
    R = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=['template', 'kind', 'source', 'draw', 'i', 'zmad'])
    B11 = pd.read_parquet(f'{P.PD}/rows.parquet')
    B11 = B11[B11.template.isin(P.VALID) & (B11.zmad >= P.LEVEL)][['template', 'kind', 'source', 'draw', 'i', 'zmad']]
    A = pd.concat([R, B11], ignore_index=True)
    key = {k: np.sort(g.i.values) for k, g in A.groupby(['source', 'draw', 'template'])}
    ag = []
    for r in R.itertuples():
        ag.append(','.join(k[2] for k, a in key.items() if k[0] == r.source and k[1] == r.draw and k[2] != r.template
                           and np.any(np.abs(a[np.searchsorted(a, r.i - 30):np.searchsorted(a, r.i + 31)] - r.i) <= 30)))
    R['agree'] = ag
    R['epoch'] = P.EP0 + R.i; R['year'] = Y; R['period'] = P.periods(R.epoch.values)
    R.to_parquet(f'{RD}/rows.parquet')
    print(F.to_string()); K.log(f'relay {Y} loo: {int(F.refound.sum())}/{len(F)}; tier candidate rows {len(R)}', 'relay')


def run_battery():
    import cat_e2
    cat_e2.run_battery(Y, rows_file=f'{RD}/rows.parquet', out_file=f'{RD}/bat.parquet', tag=f'relay{Y}')


def run_verdict():
    R = pd.read_parquet(f'{RD}/rows.parquet')
    B = pd.read_parquet(f'{RD}/bat.parquet').drop(columns=['epoch'])
    B['grid_edge'] = B.grid_edge.where(B.grid_edge.notna(), B.grid_edge_node) if 'grid_edge' in B else B.grid_edge_node
    D = R.merge(B, on='i', how='left')
    for c in ('n_far', 'n6'):
        D[c] = D[c].fillna(0)
    D = K.apply_rules(D)
    D['tier'] = (D.agree != '') & D.battery_e2
    PY = P.period_years()

    def dd(x):
        t = np.sort(x.i.values); return int(1 + np.sum(np.diff(t) > 30)) if len(t) else 0
    x = D[D.tier & (D.period == 'between')]
    cc = [dd(x[(x.source == a) & (x.draw == q)]) for a in ('A1', 'A2') for q in range(P.NDRAW)]
    rate = float(np.mean(cc)) / PY['between'] if PY['between'] > 0 else np.nan
    F = pd.read_csv(f'{K.TB}/{P.NAME}{Y}_relay_loo.csv')
    frac = float(F.refound.mean()) if len(F) else 0.0
    V = dict(year=Y, n_vetted=len(F), loo_refound=int(F.refound.sum()), loo_frac=frac, scrambled_rate_between_per_yr=rate,
             real_ETS=dd(D[D.tier & (D.period == 'ETS') & (D.source == 'real')]),
             real_between=dd(D[D.tier & (D.period == 'between') & (D.source == 'real')]))
    try:                                                     # Amendment Y7: noise-template false rate reported, scrambled not gated
        Vy = json.load(open(f'{K.TB}/{P.NAME}{Y}_verdict{P.SFX}.json'))
        V['noise_template_tier_rate_between_per_tpl_yr'] = Vy['counts'].get('noise_pair_rate_between_per_yr')
    except Exception:
        V['noise_template_tier_rate_between_per_tpl_yr'] = None
    V['accepted'] = bool(len(F) >= 3 and frac >= 0.85)
    V['rule'] = 'Amendment Y7: LOO >= 0.85 (>= 3 vetted); scrambled rate reported only'
    json.dump(V, open(f'{K.TB}/{P.NAME}{Y}_relay_verdict.json', 'w'), indent=1)
    print(V)


if __name__ == '__main__':
    {'bank': run_bank, 'mf': run_mf, 'loo': run_loo, 'battery': run_battery, 'verdict': run_verdict}[sys.argv[1]]()
