#!/usr/bin/env python
"""S13 v2 stacker (PLAN.md Amendment 5, 2026-10-14 PT). Same physics, window QC, template fit and Amendment-4 channel
guard v2 (s13_common.guard2, Amendment 4 v2) as s13_stack2.py, which is imported read-only and not edited. Changes:
  * station-outer loop (one station's continuous data in memory at a time) so that +-k-day controls fit in memory;
  * +-k-day shifted stacks for k = 1..KMAX (both signs; Merlin step 3), same QC and weights, thrust fit;
  * neighbour counts fixed (Merlin step 6): pooled rows only (unit 'pooled' of the row file, or of --nb-pool, else the
    unit's own rows), self excluded, +-2/10/20/60 s, plus a kernel-weighted sum (autocorrelation of the 20-50 s band-pass
    impulse response, normalised to 1 at lag 0, summed over neighbours within +-150 s);
  * extra thrust lag grids as columns, each with its own scrambled-control z: 'c6' = -14..+2 s, 'wide' = -40..+40 s
    (the headline grid stays -8..+8 s);
  * instrument vetoes (svi_science/instrument_vetoes.csv) blanked per component-day before the band-pass;
  * control draws are generated in the same order from the same seed as s13_stack2.py, so with --kmax 1 every m0/z equals
    the s13_stack2.py value (checked on svi_2025-09, see NOTES.md).

usage: s13_stack2_v2.py ROWS.parquet TAG STATIONS [--ndraw 500] [--ndraw-unit 200] [--kmax 5] [--nb-pool ROWS.parquet]
Writes tables/<TAG>_units.csv, tables/<TAG>_stacomp.csv, data/stacks/<TAG>.pkl, marker STACK_<TAG>_DONE."""
import os, sys, argparse, json, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G
from s13_stack import (NWIN, TREL, FITM, MECHS, LAGS, stack_windows, templates, elem_templates, fit, fit_free)
from s13_stack2 import to_idx, gather, SEP
from s13_stack import DIL
from scipy.ndimage import binary_dilation

# Instrument vetoes (addendum to Amendment 5, 2026-10-14 ~7:35 pm PT): rows of svi_science/instrument_vetoes.csv with
# decision == 'veto' blank the named component(s) of that station on every day in [start, end] (inclusive) BEFORE the
# band-pass, so data, scrambled controls and +-k-day stacks lose the same windows (window QC drops windows > 5 % masked;
# a component left with < 5 windows is then dropped by guard2). Component from the channels field: ??Z/??N/??E -> Z/N/E;
# ??1/??2 (unrotated) -> N and E; nothing recognisable -> all three.
VETO_FILE = f'{C.ROOT}/svi_science/instrument_vetoes.csv'


def vetoes(netsta):
    import re
    if not os.path.isfile(VETO_FILE):
        return []
    V = pd.read_csv(VETO_FILE)
    V = V[(V.decision.astype(str).str.strip() == 'veto') & (V.netsta == netsta)]
    out = []
    for r in V.itertuples():
        ch = set(re.findall(r'\b[A-Z]{2}([ZNE12])\b', str(r.channels)))
        comps = set()
        for c in ch:
            comps |= {'Z': {2}, 'N': {0}, 'E': {1}, '1': {0, 1}, '2': {0, 1}}[c]
        out.append((int(pd.Timestamp(r.start).value // 86400e9), int(pd.Timestamp(r.end).value // 86400e9),
                    sorted(comps or {0, 1, 2})))
    return out


def load_runs(netsta, dords, logf=None, tag=''):
    """s13_stack2.load_runs with instrument-veto blanking (identical when the station has no veto)."""
    VV = vetoes(netsta)
    runs, cur = [], [dords[0]]
    for d in dords[1:]:
        if d == cur[-1] + 1:
            cur.append(d)
        else:
            runs.append(cur); cur = [d]
    runs.append(cur)
    L = sum(len(r) * 86400 for r in runs) + SEP * (len(runs) + 1)
    B = np.full((3, L), np.nan, np.float32)
    pos, i, nveto = {}, SEP, 0
    for r in runs:
        X = np.full((3, len(r) * 86400), np.nan)
        got = False
        for k, d in enumerate(r):
            pos[d] = i + k * 86400
            day = pd.Timestamp(d * 86400, unit='s')
            f = f'{C.DISP}/{netsta}/{day.year}/{day.dayofyear:03d}.npz'
            if os.path.isfile(f):
                X[:, k * 86400:(k + 1) * 86400] = np.load(f)['d']; got = True
                for a0, a1, comps in VV:
                    if a0 <= d <= a1:
                        X[comps, k * 86400:(k + 1) * 86400] = np.nan; nveto += len(comps)
        if got:
            bad = ~np.isfinite(X)
            Y = G.bandpass(np.where(bad, 0.0, X), 1.0)
            for c in range(3):
                if bad[c].any():
                    Y[c, binary_dilation(bad[c], structure=np.ones(2 * DIL + 1, bool))] = np.nan
            B[:, i:i + X.shape[1]] = Y
        i += X.shape[1] + SEP
    if nveto:
        C.log(f'{tag}: instrument veto {netsta}: {nveto} component-days blanked', logf)
    return B, pos

GRIDS = {'': LAGS, 'c6_': np.arange(-14, 3), 'wide_': np.arange(-40, 41)}
# Amendment 6 (2026-10-14 PT): optional --grids selects which of GRIDS_ALL are fitted (default = the three above, so every
# earlier call is unchanged); 'fix2_' = the single lag +2 s (no scan). The headline '' grid is always fitted.
GRIDS_ALL = dict(GRIDS, fix2_=np.array([2]))
NB_W = (2, 10, 20, 60)


def band_ac():
    x = np.zeros(4001); x[2000] = 1.0
    h = G.bandpass(x, 1.0)
    ac = np.correlate(h, h, 'full')
    ac = ac / ac.max()
    mid = len(ac) // 2
    return ac[mid - 150: mid + 151]           # lags -150..+150 s


def nb_counts(to, pool):
    pool = np.sort(pool)
    self_in = np.isin(to, pool).astype(int)
    out = {}
    for w in NB_W:
        out[f'n_nb{w}'] = float((np.searchsorted(pool, to + w, 'right') - np.searchsorted(pool, to - w, 'left') - self_in).mean())
    ac = band_ac()
    tot = 0.0
    for t in to:
        lo, hi = np.searchsorted(pool, t - 150, 'left'), np.searchsorted(pool, t + 150, 'right')
        d = np.round(pool[lo:hi] - t).astype(int)
        tot += ac[d + 150].sum() - (1.0 if t in pool[lo:hi] else 0.0)
    out['nb_kern'] = tot / len(to)
    return out


def tfit_prep(T, sts, lags):
    """rolled+sliced templates per lag: dict L -> dict sta -> (3, nfit)."""
    return {L: {s: np.roll(T[s], L, axis=1)[:, FITM] for s in sts} for L in lags}


def fit_many(D, TF, wts, sts, lags):
    """D: dict sta -> (n, 3, nfit) stacks (FITM slice). Same rule as s13_stack.fit (first maximum of explained energy).
    Returns m0 (n,), lag (n,), vr (n,)."""
    n = D[sts[0]].shape[0]
    num = sum((wts[s][None, :, None] * D[s] ** 2).sum(axis=(1, 2)) for s in sts)
    sd = np.zeros((len(lags), n)); ss = np.zeros(len(lags))
    for i, L in enumerate(lags):
        for s in sts:
            wt = wts[s][:, None] * TF[L][s]
            sd[i] += np.einsum('ct,nct->n', wt, D[s])
            ss[i] += (wt * TF[L][s]).sum()
    expl = sd ** 2 / ss[:, None]
    k = np.argmax(expl, axis=0)
    j = np.arange(n)
    m0 = sd[k, j] / ss[k]
    with np.errstate(all='ignore'):
        vr = np.where(num > 0, expl[k, j] / num, np.nan)
    return m0, np.asarray(lags)[k], vr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rows'); ap.add_argument('tag'); ap.add_argument('stations')
    ap.add_argument('--ndraw', type=int, default=500); ap.add_argument('--ndraw-unit', type=int, default=200)
    ap.add_argument('--kmax', type=int, default=5); ap.add_argument('--nb-pool', default='')
    ap.add_argument('--grids', default=',c6_,wide_')
    a = ap.parse_args()
    global GRIDS
    GRIDS = {g: GRIDS_ALL[g] for g in [''] + [x for x in a.grids.split(',') if x]}
    lf = f'{C.LOGD}/stack2_v2.log'
    rng = np.random.default_rng(C.SEED)
    R = pd.read_parquet(a.rows)
    coords = C.station_coords()
    K = list(range(1, a.kmax + 1))
    # units and control draws (same order / seed as s13_stack2.py)
    UN = []
    for uname, U in R.groupby('unit', sort=False):
        full = bool(U.full.iloc[0]); nd = a.ndraw if (full and len(U) <= 3000) else a.ndraw_unit
        to = U.to.values; P = U.P.values
        dayd = np.floor(to / 86400) * 86400
        draws = [dayd + rng.uniform(0, 86400, len(to)) for _ in range(nd)]
        UN.append(dict(name=uname, U=U, full=full, nd=nd, to=to, P=P, draws=draws,
                       src=(U.lat.mean(), U.lon.mean(), U.dep.mean())))
    # neighbour pool
    if a.nb_pool:
        Q = pd.read_parquet(a.nb_pool); poolall = Q[Q.unit == 'pooled'].to.values
    elif (R.unit == 'pooled').any():
        poolall = R[R.unit == 'pooled'].to.values
    else:
        poolall = None
    for u in UN:
        u['nb'] = nb_counts(u['to'], poolall if poolall is not None else u['to'])
    days0 = sorted({int(d) for d in np.floor(R.to / 86400)})
    dords = sorted({d + k for d in days0 for k in range(-(a.kmax + 1), a.kmax + 2)})
    stations = []
    S = {}                                   # (unit index, sta) -> dict
    for s in a.stations.split(','):
        if s not in coords.index:
            continue
        B, pos = load_runs(s, dords, lf, a.tag)
        fr = np.mean([np.isfinite(B[2, pos[d]:pos[d] + 86400]).mean() for d in days0])
        if not fr > 0.2:
            del B; continue
        stations.append(s)
        for j, u in enumerate(UN):
            to = u['to']
            Wd = gather(B, to_idx(to, pos))
            Wc = np.concatenate([gather(B, to_idx(u['draws'][i], pos)) for i in range(min(20, u['nd']))])
            allw = np.concatenate([Wd, Wc])
            with np.errstate(all='ignore'):
                rms = np.sqrt(np.nanmean(allw ** 2, axis=2))
                thr = 3 * np.nanmedian(rms, axis=0)
            stk, nused = stack_windows(Wd, thr)
            stkP, _ = stack_windows(Wd, thr, u['P'])
            ctrl = np.stack([stack_windows(gather(B, to_idx(dr, pos)), thr)[0] for dr in u['draws']])
            sh = {}
            for k in K:
                sh[k] = stack_windows(gather(B, to_idx(to + 86400 * k, pos)), thr)[0]
                sh[-k] = stack_windows(gather(B, to_idx(to - 86400 * k, pos)), thr)[0]
            S[(j, s)] = dict(stk=stk, stkP=stkP, nused=nused, ctrlF=ctrl[:, :, FITM].copy(),
                             sig=ctrl[:, :, FITM].std(axis=0).mean(axis=1),
                             c_lo=np.quantile(ctrl, .025, axis=0), c_hi=np.quantile(ctrl, .975, axis=0), sh=sh)
            del ctrl
        del B
        C.log(f'{a.tag}: station {s} done', lf)
    C.log(f'{a.tag}: {len(R)} rows, {len(UN)} units, stations {stations}', lf)
    rows, srows, store = [], [], {}
    for j, u in enumerate(UN):
        uname, full, nd = u['name'], u['full'], u['nd']
        sts = [s for s in stations if S[(j, s)]['nused'][2] >= 5]
        if not sts:
            continue
        sig = {s: S[(j, s)]['sig'] for s in sts}
        wts, wdrop, wreason = C.guard2(sig, {s: S[(j, s)]['nused'] for s in sts})   # Amendment 4 v2 channel guard
        live = [s for s in sts if wts[s].sum() > 0]
        if not live:
            continue
        gmsg = ' '.join(f'{s}.{cn}:{wreason[s][c]}' for s in sts for c, cn in enumerate('NEZ') if wdrop[s][c])
        if gmsg:
            C.log(f'{a.tag} {uname} guard {C.GUARD_RULE}: dropped [{gmsg}]', lf)
        F1 = lambda key: {s: S[(j, s)][key][None, :, FITM] for s in sts}
        Dd, DP = F1('stk'), F1('stkP')
        Dc = {s: S[(j, s)]['ctrlF'] for s in sts}
        res = {}
        for mname, mech in (MECHS.items() if full else [('thrust', MECHS['thrust'])]):
            T = templates(u['src'], sts, coords, mech)
            grids = GRIDS.items() if mname == 'thrust' else [('', LAGS)]
            for gp, lags in grids:
                TF = tfit_prep(T, sts, lags)
                m0, L, vr = (x[0] for x in fit_many(Dd, TF, wts, sts, lags))
                cm = fit_many(Dc, TF, wts, sts, lags)[0]
                r = dict(m0=m0, lag=int(L), vr=vr, ctrl_mean=cm.mean(), ctrl_sd=cm.std(), ctrl_p99=np.quantile(cm, .99),
                         z=(m0 - cm.mean()) / cm.std(), p=(1 + (cm >= m0).sum()) / (1 + len(cm)))
                if gp == '':
                    shs = {k: fit_many({s: S[(j, s)]['sh'][k][None, :, FITM] for s in sts}, TF, wts, sts, lags)[0][0]
                           for k in sorted(S[(j, sts[0])]['sh'])}
                    r['m0_shift_p1'] = shs.get(1, np.nan); r['m0_shift_m1'] = shs.get(-1, np.nan)
                    r['m0_Pw'] = fit_many(DP, TF, wts, sts, lags)[0][0]
                    if mname == 'thrust':
                        for k, v in shs.items():
                            r[f'm0_k{k:+d}'] = v
                        zs = np.array([(v - cm.mean()) / cm.std() for v in shs.values()])
                        r['zk_rms'] = float(np.sqrt(np.mean(zs ** 2)))
                        Tthr, Lthr = T, int(L)
                for kk, v in r.items():
                    res[f'{mname}_{gp}{kk}'] = v
        stk_ = {s: S[(j, s)]['stk'] for s in sts}
        en = {s: wts[s] * (np.roll(Tthr[s], Lthr, axis=1)[:, FITM] ** 2).sum(axis=1) for s in sts}
        etot = sum(v.sum() for v in en.values())
        tshare = {s: en[s] / etot if etot > 0 else en[s] * np.nan for s in sts}
        msc = max(((s, cn, tshare[s][c]) for s in sts for c, cn in enumerate('NEZ')), key=lambda x: x[2])
        loso = {}
        if len(live) > 1:
            for s0 in live:
                loso[s0] = fit(stk_, Tthr, {s: (np.zeros(3) if s == s0 else wts[s]) for s in sts})[0]
        lmin = min(loso, key=loso.get) if loso else ''
        lmax = max(loso, key=loso.get) if loso else ''
        row = dict(unit=uname, n=len(u['to']), **u['nb'], src_lat=u['src'][0], src_lon=u['src'][1], src_dep=u['src'][2],
                   n_sta=len(live), stations=','.join(live), guard_rule=C.GUARD_RULE, guard_dropped=gmsg,
                   max_share=msc[2], max_share_ch=f'{msc[0]}.{msc[1]}', loso_m0_min=loso.get(lmin, np.nan),
                   loso_m0_max=loso.get(lmax, np.nan), loso_min_without=lmin, loso_max_without=lmax, full=full, ndraw=nd,
                   kmax=a.kmax, date0=pd.Timestamp(u['to'].min(), unit='s').strftime('%Y-%m-%d'),
                   date1=pd.Timestamp(u['to'].max(), unit='s').strftime('%Y-%m-%d'))
        if full:
            mfree, vrfree = fit_free(stk_, elem_templates(u['src'], sts, coords), wts, Lthr)
            row.update(vr_free=vrfree, mt_free=json.dumps(list(np.round(mfree, 1))))
        row.update(res)
        rows.append(row)
        for s in sts:
            for c, cn in enumerate('NEZ'):
                srows.append(dict(unit=uname, sta=s, comp=cn, n_used=int(S[(j, s)]['nused'][c]),
                                  sigma_ctrl_nm=sig[s][c] * 1e9, tpl_amp_per1e12_nm=np.abs(Tthr[s][c]).max() * 1e21,
                                  tpl_share=float(tshare[s][c]), guard_dropped=bool(wdrop[s][c]),
                                  drop_reason=wreason[s][c]))
        store[uname] = dict(stk={s: S[(j, s)]['stk'] for s in sts}, c_lo={s: S[(j, s)]['c_lo'] for s in sts},
                            c_hi={s: S[(j, s)]['c_hi'] for s in sts}, tpl=Tthr, lag=Lthr,
                            shp={s: S[(j, s)]['sh'][1] for s in sts}, shm={s: S[(j, s)]['sh'][-1] for s in sts},
                            sig={s: sig[s] for s in sts}, wts={s: wts[s] for s in sts})
        C.log(f'{a.tag} {uname} n={len(u["to"])} M0={res["thrust_m0"]:.3g} z={res["thrust_z"]:.2f} lag={res["thrust_lag"]} '
              f'VR={res["thrust_vr"]:.2f} sta={len(sts)}', lf)
    pd.DataFrame(rows).to_csv(f'{C.TB}/{a.tag}_units.csv', index=False)
    pd.DataFrame(srows).to_csv(f'{C.TB}/{a.tag}_stacomp.csv', index=False)
    os.makedirs(f'{C.DATA}/stacks', exist_ok=True)
    with open(f'{C.DATA}/stacks/{a.tag}.pkl', 'wb') as fh:
        pickle.dump(dict(trel=TREL, units=store), fh)
    open(f'{C.OUT}/STACK_{a.tag}_DONE', 'w').write('ok\n')


if __name__ == '__main__':
    main()
