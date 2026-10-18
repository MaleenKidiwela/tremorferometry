#!/usr/bin/env python
"""S13 scale-up stacker (same physics/QC/fit as the pilot s13_stack.py; PLAN.md §3, Amendment 1; channel guard
Amendment 4 v2, 2026-10-14: flat/dead/missing station-components dropped, no weight cap, see s13_common.guard2).

usage: s13_stack2.py ROWS.parquet TAG STATIONS [--ndraw 500] [--ndraw-unit 200]
ROWS columns: unit, to (origin, epoch s), P, lat, lon, dep (km), full (bool: alternatives + ndraw draws; else thrust
only + ndraw-unit draws). Days are loaded as contiguous runs (each needed day +-1 d), separated by NaN blocks, so sparse
sets (between-ETS bursts) work too. Writes tables/<TAG>_units.csv, tables/<TAG>_stacomp.csv, data/stacks/<TAG>.pkl,
marker STACK_<TAG>_DONE."""
import os, sys, argparse, json, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G
from s13_stack import (W0, NWIN, TREL, FITM, MECHS, DIL, stack_windows, templates, elem_templates, fit, fit_free)

SEP = 2 * DIL + NWIN + 10


def load_runs(netsta, dords):
    """dords: sorted day ordinals (epoch day numbers). Returns band-passed array (3, L) and dict dord -> start index."""
    runs, cur = [], [dords[0]]
    for d in dords[1:]:
        if d == cur[-1] + 1:
            cur.append(d)
        else:
            runs.append(cur); cur = [d]
    runs.append(cur)
    L = sum(len(r) * 86400 for r in runs) + SEP * (len(runs) + 1)
    B = np.full((3, L), np.nan, np.float32)
    pos, i = {}, SEP
    for r in runs:
        X = np.full((3, len(r) * 86400), np.nan)
        got = False
        for k, d in enumerate(r):
            pos[d] = i + k * 86400
            day = pd.Timestamp(d * 86400, unit='s')
            f = f'{C.DISP}/{netsta}/{day.year}/{day.dayofyear:03d}.npz'
            if os.path.isfile(f):
                X[:, k * 86400:(k + 1) * 86400] = np.load(f)['d']; got = True
        if got:
            bad = ~np.isfinite(X)
            Y = G.bandpass(np.where(bad, 0.0, X), 1.0)
            for c in range(3):
                if bad[c].any():
                    Y[c, binary_dilation(bad[c], structure=np.ones(2 * DIL + 1, bool))] = np.nan
            B[:, i:i + X.shape[1]] = Y
        i += X.shape[1] + SEP
    return B, pos


def to_idx(to, pos):
    d = np.floor(to / 86400).astype(np.int64)
    base = np.array([pos.get(int(x), -10 ** 12) for x in d], np.int64)
    return base + np.round(to - d * 86400).astype(np.int64) + W0


def gather(B, idx):
    n = len(idx)
    out = np.full((n, 3, NWIN), np.nan)
    ok = (idx >= 0) & (idx + NWIN <= B.shape[1])
    if ok.any():
        ii = idx[ok][:, None] + np.arange(NWIN)[None, :]
        out[ok] = B[:, ii].transpose(1, 0, 2).astype(np.float64)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rows'); ap.add_argument('tag'); ap.add_argument('stations')
    ap.add_argument('--ndraw', type=int, default=500); ap.add_argument('--ndraw-unit', type=int, default=200)
    a = ap.parse_args()
    lf = f'{C.LOGD}/stack2.log'
    rng = np.random.default_rng(C.SEED)
    R = pd.read_parquet(a.rows)
    coords = C.station_coords()
    dords = sorted({int(d) + k for d in np.floor(R.to / 86400) for k in (-1, 0, 1)})
    B, pos = {}, {}
    for s in a.stations.split(','):
        if s not in coords.index:
            continue
        B[s], pos[s] = load_runs(s, dords)
    # keep stations with Z data on >= 20 % of the needed days' samples
    need = sorted({int(d) for d in np.floor(R.to / 86400)})
    stations = []
    for s in B:
        fr = np.mean([np.isfinite(B[s][2, pos[s][d]:pos[s][d] + 86400]).mean() for d in need])
        if fr > 0.2:
            stations.append(s)
    C.log(f'{a.tag}: {len(R)} rows, {R.unit.nunique()} units, stations {stations}', lf)
    allt = np.sort(R.to.values)
    rows, srows, store = [], [], {}
    for uname, U in R.groupby('unit', sort=False):
        full = bool(U.full.iloc[0]); nd = a.ndraw if (full and len(U) <= 3000) else a.ndraw_unit
        to = U.to.values; P = U.P.values
        nb = (np.searchsorted(allt, to + 10) - np.searchsorted(allt, to - 10) - 1).mean()
        src = (U.lat.mean(), U.lon.mean(), U.dep.mean())
        dayd = np.floor(to / 86400) * 86400
        draws = [dayd + rng.uniform(0, 86400, len(to)) for _ in range(nd)]
        stk, stkP, ctrl, shp, shm, nused, thr = {}, {}, {}, {}, {}, {}, {}
        for s in stations:
            Wd = gather(B[s], to_idx(to, pos[s]))
            Wc = np.concatenate([gather(B[s], to_idx(draws[i], pos[s])) for i in range(min(20, nd))])
            allw = np.concatenate([Wd, Wc])
            with np.errstate(all='ignore'):
                rms = np.sqrt(np.nanmean(allw ** 2, axis=2))
                thr[s] = 3 * np.nanmedian(rms, axis=0)
            stk[s], nused[s] = stack_windows(Wd, thr[s])
            stkP[s], _ = stack_windows(Wd, thr[s], P)
            ctrl[s] = np.stack([stack_windows(gather(B[s], to_idx(dr, pos[s])), thr[s])[0] for dr in draws])
            shp[s], _ = stack_windows(gather(B[s], to_idx(to + 86400, pos[s])), thr[s])
            shm[s], _ = stack_windows(gather(B[s], to_idx(to - 86400, pos[s])), thr[s])
        sts = [s for s in stations if nused[s][2] >= 5]
        if not sts:
            continue
        sig = {s: ctrl[s][:, :, FITM].std(axis=0).mean(axis=1) for s in sts}
        wts, wdrop, wreason = C.guard2(sig, {s: nused[s] for s in sts})   # PLAN.md Amendment 4 v2 channel guard
        live = [s for s in sts if wts[s].sum() > 0]
        if not live:
            continue
        gmsg = ' '.join(f'{s}.{cn}:{wreason[s][c]}' for s in sts for c, cn in enumerate('NEZ') if wdrop[s][c])
        if gmsg:
            C.log(f'{a.tag} {uname} guard {C.GUARD_RULE}: dropped [{gmsg}]', lf)
        stk_, stkP_, shp_, shm_ = ({s: d[s] for s in sts} for d in (stk, stkP, shp, shm))
        res = {}
        for mname, mech in (MECHS.items() if full else [('thrust', MECHS['thrust'])]):
            T = templates(src, sts, coords, mech)
            m0, L, vr = fit(stk_, T, wts)
            cm = np.array([fit({s: ctrl[s][i] for s in sts}, T, wts)[0] for i in range(nd)])
            res[mname] = dict(m0=m0, lag=L, vr=vr, ctrl_mean=cm.mean(), ctrl_sd=cm.std(), ctrl_p99=np.quantile(cm, .99),
                              z=(m0 - cm.mean()) / cm.std(), p=(1 + (cm >= m0).sum()) / (1 + len(cm)),
                              m0_shift_p1=fit(shp_, T, wts)[0], m0_shift_m1=fit(shm_, T, wts)[0], m0_Pw=fit(stkP_, T, wts)[0])
            if mname == 'thrust':
                Tthr, Lthr = T, L
        # template-energy share per station-component, (tpl/sigma)^2 over the fit window, and leave-one-station-out M0
        en = {s: wts[s] * (np.roll(Tthr[s], Lthr, axis=1)[:, FITM] ** 2).sum(axis=1) for s in sts}
        etot = sum(v.sum() for v in en.values())
        tshare = {s: en[s] / etot if etot > 0 else en[s] * np.nan for s in sts}
        msc = max(((s, cn, tshare[s][c]) for s in sts for c, cn in enumerate('NEZ')), key=lambda x: x[2])
        loso = {}
        if len(live) > 1:
            for s0 in live:
                w0 = {s: (np.zeros(3) if s == s0 else wts[s]) for s in sts}
                loso[s0] = fit(stk_, Tthr, w0)[0]
        lmin = min(loso, key=loso.get) if loso else ''
        lmax = max(loso, key=loso.get) if loso else ''
        row = dict(unit=uname, n=len(U), n_nb10=nb, src_lat=src[0], src_lon=src[1], src_dep=src[2], n_sta=len(live),
                   stations=','.join(live), guard_rule=C.GUARD_RULE, guard_dropped=gmsg,
                   max_share=msc[2], max_share_ch=f'{msc[0]}.{msc[1]}', loso_m0_min=loso.get(lmin, np.nan),
                   loso_m0_max=loso.get(lmax, np.nan), loso_min_without=lmin, loso_max_without=lmax, full=full, ndraw=nd, date0=pd.Timestamp(to.min(), unit='s').strftime('%Y-%m-%d'),
                   date1=pd.Timestamp(to.max(), unit='s').strftime('%Y-%m-%d'))
        if full:
            mfree, vrfree = fit_free(stk_, elem_templates(src, sts, coords), wts, Lthr)
            row.update(vr_free=vrfree, mt_free=json.dumps(list(np.round(mfree, 1))))
        for mname, r in res.items():
            for k, v in r.items():
                row[f'{mname}_{k}'] = v
        rows.append(row)
        for s in sts:
            for c, cn in enumerate('NEZ'):
                srows.append(dict(unit=uname, sta=s, comp=cn, n_used=int(nused[s][c]), sigma_ctrl_nm=sig[s][c] * 1e9,
                                  tpl_amp_per1e12_nm=np.abs(Tthr[s][c]).max() * 1e21,
                                  tpl_share=float(tshare[s][c]), guard_dropped=bool(wdrop[s][c]),
                                  drop_reason=wreason[s][c]))
        store[uname] = dict(stk=stk_, c_lo={s: np.quantile(ctrl[s], .025, axis=0) for s in sts},
                            c_hi={s: np.quantile(ctrl[s], .975, axis=0) for s in sts}, tpl=Tthr, lag=Lthr,
                            shp=shp_, shm=shm_, wts=wts)
        C.log(f'{a.tag} {uname} n={len(U)} M0={res["thrust"]["m0"]:.3g} z={res["thrust"]["z"]:.2f} '
              f'VR={res["thrust"]["vr"]:.2f} sta={len(sts)}', lf)
    pd.DataFrame(rows).to_csv(f'{C.TB}/{a.tag}_units.csv', index=False)
    pd.DataFrame(srows).to_csv(f'{C.TB}/{a.tag}_stacomp.csv', index=False)
    os.makedirs(f'{C.DATA}/stacks', exist_ok=True)
    with open(f'{C.DATA}/stacks/{a.tag}.pkl', 'wb') as fh:
        pickle.dump(dict(trel=TREL, units=store), fh)
    open(f'{C.OUT}/STACK_{a.tag}_DONE', 'w').write('ok\n')


if __name__ == '__main__':
    main()
