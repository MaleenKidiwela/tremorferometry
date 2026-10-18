#!/usr/bin/env python
"""Amendment Y9 step 4: similarity-matched noise banks (2011-2013). Run with the year's production env (PROD_YEAR, PILOT_SFX=_v2,
PROD_QC for 2013) and, for the battery, the year's scan-v2 env.
  build   -> data/y9nb/<run><Y>/banks.npz + tables/y9nb_<Y>_banks.csv (achieved mean CC per bank, real-bank target)
  mf      -> data/y9nb/<run><Y>/mf/<bank>_s<j>.npz (siblings; the base reuses the production noise MF)
  detect  -> data/y9nb/<run><Y>/rows.parquet (agreement events: >= 2 bank members within 30 s, MF >= 5.5)
  battery -> data/y9nb/<run><Y>/bat.parquet (scan-v2 battery, every agreement event)
  summary (no env needed) -> tables/y9nb_summary.csv, tables/y9nb_verdict.json"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

NB, NSIB, LAG = 8, 10, 30
RUNS = {2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq'}   # Y9b: 2014-17


def paths(y):
    d = f'{K.DATA}/y9nb/{RUNS[y]}{y}'
    os.makedirs(f'{d}/mf', exist_ok=True)
    return d


def pair_cc(A, B):
    """mean over channels of the max normalised CC within +-LAG s; A, B (nch, L) with NaN rows for missing channels."""
    v = []
    for a, b in zip(A, B):
        if not (np.isfinite(a).all() and np.isfinite(b).all()):
            continue
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na == 0 or nb == 0:
            continue
        c = np.correlate(a, b, mode='full') / (na * nb); m = len(a) - 1
        v.append(c[m - LAG:m + LAG + 1].max())
    return float(np.mean(v)) if v else np.nan


def bank_cc(T):
    n = len(T)
    return float(np.nanmean([pair_cc(T[i], T[j]) for i in range(n) for j in range(i + 1, n)]))


def build():
    import cat_pilot as P
    y = P.PILOT; d = paths(y)
    B = np.load(f'{P.PD}/bank.npz'); Mt = pd.read_csv(f'{K.TB}/{P.NAME}{y}_bank.csv')
    real = [t for t in [f'T{k}' for k in range(1, 6)] + [f'C{k}' for k in range(6)]]
    target = bank_cc([B[t] for t in real])
    noise = Mt[Mt.kind == 'noise'].sort_values('name').head(NB)
    rng = np.random.default_rng(P.SEED + 99)
    wins = [(max(P.tidx(a_), 700), min(P.tidx(b_), P.N - 700)) for a_, b_ in [w[:2] for w in P.CFG['wins']]]
    wins = [w for w in wins if w[1] - w[0] > 3600] or [(700, 60 * 86400)]
    need = [ch for ch in P.CH if sum(np.isfinite(B[t][P.CH.index(ch)]).all() for t in real) >= 6]
    used = set(int(x) for x in Mt[Mt.kind == 'noise'].i0)
    out, meta = {}, []
    for r in noise.itertuples():
        base = B[r.name]; Q, qi = [], []
        while len(Q) < NSIB:
            a, b = wins[int(rng.integers(0, len(wins)))]
            i = int(rng.integers(a + 600, b - P.L))
            if any(abs(i - x) < 180 for x in used | set(qi)):
                continue
            w = P.cut_frame(i)
            qm = {ch: P.qc_mask(ch, P.N) for ch in need}
            if all(w[ch] is not None or (qm[ch] is not None and qm[ch][i]) for ch in need):
                Q.append(np.array([w[ch] if w[ch] is not None else np.full(P.L, np.nan) for ch in P.CH])); qi.append(i)

        def mk(a):
            S = []
            for q in Q:
                s = a * base + np.sqrt(1 - a * a) * q
                nr = np.linalg.norm(np.nan_to_num(s), axis=1, keepdims=True)
                S.append(np.where(nr > 0, s / np.where(nr > 0, nr, 1), np.nan))
            return [base] + S
        lo, hi = 0.0, 1.0
        for _ in range(14):                                  # bisection on the bank mean pairwise CC
            a = 0.5 * (lo + hi); c = bank_cc(mk(a))
            lo, hi = (a, hi) if c < target else (lo, a)
        a = 0.5 * (lo + hi); T = mk(a); c = bank_cc(T)
        for j, s in enumerate(T[1:]):
            out[f'{r.name}_s{j}'] = s
        out[f'w__{r.name}'] = B[f'w__{r.name}']
        meta.append(dict(year=y, bank=r.name, base_i0=int(r.i0), a=a, mean_cc=c, target_cc=target,
                         base_sib_cc=float(np.mean([pair_cc(T[0], s) for s in T[1:]])), q_i0=','.join(map(str, qi))))
        K.log(f'Y9nb build {y} {r.name}: a {a:.3f} cc {c:.3f} target {target:.3f}', 'y9nb')
    np.savez(f'{d}/banks.npz', chans=np.array(P.CH), **out)
    pd.DataFrame(meta).to_csv(f'{K.TB}/y9nb_{y}_banks.csv', index=False)
    print(pd.DataFrame(meta).round(3).to_string())


def mf_job(job):
    import cat_pilot as P
    import cat_e3 as E3
    y = P.PILOT; d = paths(y); name, base = job
    f = f'{d}/mf/{name}.npz'
    if os.path.isfile(f):
        return f'{name} exists'
    Bk = np.load(f'{d}/banks.npz'); W = Bk[name]; wv = Bk[f'w__{base}']
    chs = [ch for k, ch in enumerate(P.CH) if np.isfinite(W[k]).all() and (P.stat() == 'plain' or np.isfinite(wv[k]))]
    wk = np.array([wv[P.CH.index(ch)] for ch in chs])
    wk = np.ones(len(chs)) if (P.stat() == 'plain' or not np.isfinite(wk).all() or wk.sum() <= 0) else wk * len(wk) / wk.sum()
    S = np.zeros(P.M, np.float32); NL = np.zeros(P.M, np.int16)
    for j, ch in enumerate(chs):
        cc, live = P.ccser(ch, W[P.CH.index(ch)])
        S += np.float32(wk[j]) * cc; NL += live
    pi, pv, pn, med, mad = E3.peaks(S, NL)
    np.savez_compressed(f + '.part.npz', real_pi=pi, real_pv=pv, real_med=med, real_mad=mad, chans=np.array(chs))
    os.replace(f + '.part.npz', f)
    return f'{name} done ({len(chs)} ch)'


def mf():
    import cat_pilot as P
    y = P.PILOT
    meta = pd.read_csv(f'{K.TB}/y9nb_{y}_banks.csv')
    jobs = [(f'{b}_s{j}', b) for b in meta.bank for j in range(NSIB)]
    from multiprocessing import Pool
    with Pool(P.NW, maxtasksperchild=4) as p:
        for k, r in enumerate(p.imap_unordered(mf_job, jobs)):
            K.log(f'Y9nb mf {y} {k + 1}/{len(jobs)}: {r}', 'y9nb')


def detect():
    import cat_pilot as P
    y = P.PILOT; d = paths(y)
    meta = pd.read_csv(f'{K.TB}/y9nb_{y}_banks.csv')
    rows = []
    for m in meta.itertuples():
        excl = np.array([m.base_i0] + [int(x) for x in str(m.q_i0).split(',')])
        members = [(m.bank, f'{P.MFD}/{m.bank}.npz')] + [(f'{m.bank}_s{j}', f'{d}/mf/{m.bank}_s{j}.npz') for j in range(NSIB)]
        for nm, f in members:
            z = np.load(f); zz = (z['real_pv'] - float(z['real_med'])) / float(z['real_mad']); pi = z['real_pi']
            ok = (zz >= P.LEVEL) & (np.min(np.abs(pi[:, None] - excl[None, :]), axis=1) > 90)
            rows.append(pd.DataFrame(dict(bank=m.bank, member=nm, i=pi[ok].astype(np.int64), zmad=zz[ok])))
    R = pd.concat(rows, ignore_index=True).sort_values(['bank', 'i'])
    ev = []
    for b, g in R.groupby('bank'):
        g = g.assign(grp=(g.i.diff().fillna(1e9) > 30).cumsum())
        for _, h in g.groupby('grp'):
            r = h.loc[h.zmad.idxmax()]
            ev.append(dict(bank=b, member=r.member, i_frame=int(r.i), zmad=float(r.zmad), n_members=int(h.member.nunique())))
    E = pd.DataFrame(ev)
    E['epoch'] = P.EP0 + E.i_frame; E['year'] = y; E['period'] = P.periods(E.epoch.values)
    E.to_parquet(f'{d}/events_all.parquet')
    A = E[E.n_members >= 2].reset_index(drop=True); A['i'] = np.arange(len(A)); A['template'] = A.member
    A.to_parquet(f'{d}/rows.parquet')
    K.log(f'Y9nb detect {y}: {len(E)} events, {len(A)} with agreement', 'y9nb')


def battery():
    import cat_pilot as P
    import cat_e2
    y = P.PILOT; d = paths(y)
    cat_e2.run_battery(y, rows_file=f'{d}/rows.parquet', out_file=f'{d}/bat.parquet', tag='y9nb', gate_all=True)


def summary():
    import cat_ychecks as YC
    real_ag = pd.read_csv(f'{K.TB}/r20_footing_raw.csv').set_index('year')
    lad = pd.read_csv(f'{K.TB}/y9_ladder_table.csv')
    out = []
    for y, run in RUNS.items():
        d = f'{K.DATA}/y9nb/{run}{y}'
        if not os.path.isfile(f'{d}/bat.parquet'):
            continue
        TB = YC.between_years(y, run)
        A = pd.read_parquet(f'{d}/rows.parquet'); b = K.apply_rules(pd.read_parquet(f'{d}/bat.parquet'))
        A = A.merge(b[['i', 'pass_a', 'pass_b', 'pass_c', 'pass_e2', 'pass_f', 'a_testable', 'travelling', 'z_best', 'z_loso']].rename(columns={'z_best': 'z_best_b'}), on='i', how='left')
        meas = dict(a=A.a_testable.fillna(False).astype(bool), b=A.travelling.notna(), c=A.z_best_b.notna(), e=pd.Series(True, index=A.index), f=A.z_loso.notna())
        pas = dict(a=A.pass_a, b=A.pass_b, c=A.pass_c, e=A.pass_e2, f=A.pass_f)
        A['level'] = sum((meas[k] & pas[k].fillna(False).astype(bool)).astype(int) for k in 'abcef')
        E = pd.read_parquet(f'{d}/events_all.parquet')
        for bk in sorted(A.bank.unique()):
            a = A[A.bank == bk]; e = E[E.bank == bk]
            r = dict(year=y, bank=bk, between_yr=TB, agree_between_per_yr=float(((a.period == 'between')).sum() / TB),
                     real_bank_agree_between_per_yr=float(real_ag.loc[y, 'real_bank_agree_per_yr']),
                     events_between_per_yr=float((e.period == 'between').sum() / TB))
            for L in range(6):
                r[f'between_level_ge{L}'] = int(((a.period == 'between') & (a.level >= L)).sum())
                r[f'ETS_level_ge{L}'] = int(((a.period == 'ETS') & (a.level >= L)).sum())
            out.append(r)
    S = pd.DataFrame(out); S['agree_ratio_to_real'] = S.agree_between_per_yr / S.real_bank_agree_between_per_yr
    S.to_csv(f'{K.TB}/y9nb_summary.csv', index=False)
    real5 = {y: int(lad[(lad.year == str(y)) & (lad.period == 'between') & (lad.level_ge == 5)].real.iloc[0]) for y in sorted(S.year.unique())}
    # 2026-10-16 fix: the ladder's 'pooled' row grew to 2010-2021 with the Y11 extension; use the 2010-2018 years explicitly
    l18 = lad[(lad.year.isin([str(y) for y in range(2010, 2019)])) & (lad.period == 'between') & (lad.level_ge == 5)]
    pooled_real5_all = int(l18.real.sum())
    live_all = float(l18.live_yr.sum())
    from scipy.stats import chi2
    def gw(n):
        return (0.0 if n == 0 else float(chi2.ppf(.025, 2 * n) / 2), float(chi2.ppf(.975, 2 * n + 2) / 2))
    acc = S.agree_ratio_to_real.between(1 / 1.5, 1.5); S['accepted'] = acc
    S.to_csv(f'{K.TB}/y9nb_summary.csv', index=False)
    live13 = float(sum(S.groupby('year').between_yr.first()))
    per_year = []
    for y, g in S.groupby('year'):
        ga = g[g.accepted]; n = int(ga.between_level_ge5.sum()); T = float(ga.between_yr.sum()); lo, hi = gw(n); ly = float(g.between_yr.iloc[0])
        rl, rh = gw(real5[y])
        per_year.append(dict(year=int(y), n_accepted=len(ga), ctrl_rate=n / T if T else np.nan, ctrl_lo=lo / T if T else np.nan, ctrl_hi=hi / T if T else np.nan,
                             real=real5[y], live_yr=ly, real_rate=real5[y] / ly, real_lo=rl / ly, real_hi=rh / ly,
                             expected_ctrl=(n / T * ly) if T else np.nan))
    PY = pd.DataFrame(per_year); PY.to_csv(f'{K.TB}/y9nb_per_year.csv', index=False); print(PY.round(2).to_string(index=False))
    def pooled(x):
        n = int(x.between_level_ge5.sum()); T = float(x.between_yr.sum()); lo, hi = gw(n)
        return dict(n_banks=len(x), level5=n, bank_years=T, rate=n / T, rate_lo=lo / T, rate_hi=hi / T,
                    scaled_2010_18=n / T * live_all, scaled_lo=lo / T * live_all, scaled_hi=hi / T * live_all)
    A_, All_ = pooled(S[acc]), pooled(S)
    r13 = int(sum(real5.values())); lo13, hi13 = gw(r13); loa, hia = gw(pooled_real5_all)
    sc = A_['scaled_2010_18']
    V = dict(rule='Y9: accepted banks only (between-ETS agreement within x1.5 of the real bank); A if scaled control >= 20, B if <= 5',
             accepted=A_, all_banks_column=All_, agree_ratio_range=[float(S.agree_ratio_to_real.min()), float(S.agree_ratio_to_real.max())],
             years=sorted(int(x) for x in S.year.unique()), real_years=dict(level5=r13, live_yr=live13, rate=r13 / live13, rate_lo=lo13 / live13, rate_hi=hi13 / live13),
             real_2010_18=dict(level5=pooled_real5_all, live_yr=live_all, lo=loa, hi=hia),
             outcome=('A' if sc >= 20 else 'B' if sc <= 5 else 'intermediate (neither A nor B)'))
    json.dump(V, open(f'{K.TB}/y9nb_verdict.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250)
    print(S.round(2).to_string(index=False)); print(json.dumps(V, indent=1, default=float))


if __name__ == '__main__':
    {'build': build, 'mf': mf, 'detect': detect, 'battery': battery, 'summary': summary}[sys.argv[1]]()
