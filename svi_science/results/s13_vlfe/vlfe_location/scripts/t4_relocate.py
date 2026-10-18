#!/usr/bin/env python
"""Test 4 (PLAN.md + NOTES amendment A1): re-locate tier VLFEs with alternative Green's functions on the production data
windows (data/win_YEAR.npz from t2_zsurface.py), same grid / Slab2 depths / weights / on-flags; statistic = weighted
VR of a thrust template, free M0 >= 0, one common lag (+-10 s). Variants:
  ws      production whole-space dynamic templates x production gains g                 (0) baseline
  ws_st   whole-space dynamic x stack terms (calibrated on S13 pooled ETS stacks)      (c0)
  ok      half-space Okada static x band-passed step delayed r/beta x production g     (b)
  ok_g1   same, g = 1                                                                   (b1, column)
  ok_st   Okada x stack terms                                                           (c)
  wsst    whole-space STATIC x2 x step (diagnostic column: dynamic vs free-surface)
-> tables/t4_station_terms.csv, tables/t4_reloc.csv, tables/t4_summary.csv"""
import os, sys, glob, pickle
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S13 = os.path.dirname(OUT)
sys.path.insert(0, f'{S13}/scripts')
import s13_green as G
import scan_common as S
from s13_stack import TREL, FITM
from pyrocko.modelling import okada_ext

VS, VP, RHO = G.VS, G.VP, G.RHO
MU = VS ** 2 * RHO; LAM = VP ** 2 * RHO - 2 * MU
STRIKE, DIP = 315.0, 12.0
LAGS = np.arange(-10, 11)
TS = S.TS


def step_bp(delay, out_t, tau=1.0):
    from scipy.special import erf
    h = 0.5 * (1 + erf((TS - delay) / tau))
    hb = G.bandpass(h[None, :], 4.0)[0]
    return np.interp(out_t, TS, hb)


def static_okada(src, sll):
    """unit-moment static displacement N,E,Up (m per N m) at stations sll (k,2) for a thrust point-like patch."""
    n, e = G.ned_offset(src[0], src[1], sll[:, 0], sll[:, 1])
    L = 500.0
    patch = np.array([[0.0, 0.0, src[2] * 1e3, STRIKE, DIP, -L / 2, L / 2, -L / 2, L / 2]])
    dis = np.array([[0.0, 1.0, 0.0]])                       # pure up-dip slip = rake 90
    rec = np.c_[n, e, np.zeros_like(n)]
    u = okada_ext.okada(patch, dis, rec, LAM, MU, 1)[:, :3]
    m0 = MU * L * L * 1.0
    return np.c_[u[:, 0], u[:, 1], -u[:, 2]] / m0          # N, E, Up


def static_ws(src, sll):
    M = G.mt_ned(STRIKE, DIP, 90, 1.0)
    t = np.array([-1.0, 2000.0])
    out = []
    for la, lo in sll:
        n, e = G.ned_offset(src[0], src[1], la, lo)
        u = S.synth_fast(M, (0, 0, src[2] * 1e3), (n, e, 0.0), np.arange(-20, 200, 0.5))
        out.append(u[:, -1])
    return np.array(out)


def static_templates(src, sll, out_t, kind):
    U = static_okada(src, sll) if kind == 'ok' else static_ws(src, sll)
    n, e = G.ned_offset(src[0], src[1], sll[:, 0], sll[:, 1])
    r = np.sqrt(n ** 2 + e ** 2 + (src[2] * 1e3) ** 2)
    T = np.zeros((len(sll), 3, len(out_t)))
    for k in range(len(sll)):
        h = step_bp(r[k] / VS, out_t)
        T[k] = U[k][:, None] * h[None, :]
    return T * 1e21                                         # same scaling as the production grid (nm per 1e12 N m)


def dyn_templates(src, sll, out_t):
    M = G.mt_ned(STRIKE, DIP, 90, 1.0)
    return np.array([S.template(M, src, (la, lo), out_t=out_t) for la, lo in sll]) * 1e21


# ------------------------------------------------------------------------------- stack terms
def stack_list():
    out = []
    for f in sorted(glob.glob(f'{S13}/tables/*_units.csv')):
        b = os.path.basename(f)[:-10]
        parts = b.split('_', 2)
        if len(parts) != 3:
            continue
        pre, reg, ep = parts
        if pre not in ('scalef', 'refit', 'scale') or reg not in ('svi', 'wa') or ep == 'between':
            continue
        out.append((pre, reg, ep))
    pref = {'scalef': 0, 'refit': 1, 'scale': 2}
    best = {}
    for pre, reg, ep in out:
        k = (reg, ep)
        if k not in best or pref[pre] < pref[best[k]]:
            best[k] = pre
    return [(p, r, e) for (r, e), p in sorted(best.items())]


def stack_terms(stas, coords, kinds=('ws', 'ok'), zmin=4.0):
    rows = []
    acc = {k: {} for k in kinds}
    for pre, reg, ep in stack_list():
        tag = f'{pre}_{reg}_{ep}'
        U = pd.read_csv(f'{S13}/tables/{tag}_units.csv'); U = U[U.unit == 'pooled']
        f = f'{S13}/data/stacks/{tag}.pkl'
        if U.empty or not os.path.isfile(f) or not np.isfinite(U.thrust_z.iloc[0]) or U.thrust_z.iloc[0] < zmin:
            continue
        u = U.iloc[0]
        st = pickle.load(open(f, 'rb'))['units']['pooled']
        Sc = pd.read_csv(f'{S13}/tables/{tag}_stacomp.csv'); Sc = Sc[Sc.unit == 'pooled']
        src = (u.src_lat, u.src_lon, u.src_dep)
        chans = []
        for s in st['stk']:
            if s not in coords:
                continue
            for c, cn in enumerate('NEZ'):
                r = Sc[(Sc.sta == s) & (Sc.comp == cn)]
                if r.empty or not np.isfinite(r.sigma_ctrl_nm.iloc[0]) or r.sigma_ctrl_nm.iloc[0] <= 0 or r.n_used.iloc[0] < 5:
                    continue
                if 'guard_dropped' in r and bool(r.guard_dropped.iloc[0]):
                    continue
                d = st['stk'][s][c]
                if not np.isfinite(d[FITM]).all():
                    continue
                chans.append((s, c, d * 1e9, r.sigma_ctrl_nm.iloc[0]))
        if len(chans) < 4:
            continue
        sll = np.array([coords[s] for s, _, _, _ in chans])
        for kind in kinds:
            Tall = dyn_templates(src, sll, TREL) if kind == 'ws' else static_templates(src, sll, TREL, 'ok')
            T = np.array([Tall[k, c] for k, (s, c, _, _) in enumerate(chans)])          # (nch, 401)
            D = np.array([x for _, _, x, _ in chans]); w = np.array([1 / sg ** 2 for *_, sg in chans])
            best = None
            for L in range(-8, 9):
                Tl = np.roll(T, L, axis=1)[:, FITM]
                sd = (w * (Tl * D[:, FITM]).sum(1)).sum(); ss = (w * (Tl * Tl).sum(1)).sum()
                if best is None or sd * sd / ss > best[0]:
                    best = (sd * sd / ss, sd / ss, L)
            _, m0, L = best
            if m0 <= 0:
                continue
            Tl = np.roll(T, L, axis=1)[:, FITM]
            for k, (s, c, x, sg) in enumerate(chans):
                tt = (Tl[k] ** 2).sum()
                if tt <= 0:
                    continue
                ratio = (Tl[k] * x[FITM]).sum() / (m0 * tt)
                var = sg ** 2 / (m0 ** 2 * tt)
                a = acc[kind].setdefault((s, c), [0.0, 0.0, []])
                a[0] += ratio / var; a[1] += 1 / var; a[2].append(ep)
                rows.append(dict(kind=kind, stack=tag, sta=s, comp='NEZ'[c], m0=m0 / 1e12 * 1e12, lag=L, ratio=ratio, sd=np.sqrt(var)))
    terms = []
    for kind in kinds:
        for (s, c), (num, den, eps) in acc[kind].items():
            t = num / den; sd = 1 / np.sqrt(den)
            snr = t / sd
            use = t if (t > 0 and snr >= 2) else (0.0 if (t < 0 and snr <= -2) else 1.0)
            terms.append(dict(kind=kind, sta=s, comp='NEZ'[c], term=t, sd=sd, snr=snr, n_stacks=len(eps),
                              used=float(np.clip(use, 0.25, 4.0)) if use > 0 else 0.0))
    return pd.DataFrame(terms), pd.DataFrame(rows)


# ------------------------------------------------------------------------------- relocation
def relocate(x, sig, won, sv, Tg):
    """x (nsta,3,211) nm, Tg (nodes,nsta,3,181) effective templates. -> best node, lag, VR, VR per node (best lag)."""
    ok = won & np.isfinite(sig) & (sig > 0) & np.isfinite(x).all(-1)
    ok &= (np.abs(Tg).sum(axis=(0, 3)) > 0)                       # channels with a zero template (term 0) drop out
    idx = np.argwhere(ok)
    if len(idx) < 3:
        return None
    w = np.array([1 / sig[i, c] ** 2 for i, c in idx])
    X = np.array([x[i, c] * sv[i, c] for i, c in idx])            # (nch, 211)
    T = np.stack([Tg[:, i, c] for i, c in idx], axis=1)           # (nodes, nch, 181)
    Wn = np.lib.stride_tricks.sliding_window_view(X, 181, axis=1)  # (nch, 31, 181)
    Wn = Wn[:, 15 + LAGS, :]                                       # (nch, nlag, 181)
    num = np.einsum('c,nct,clt->nl', w, T, Wn)
    den = np.einsum('c,nct->n', w, T * T)
    dd = np.einsum('c,clt->l', w, Wn * Wn)
    vr = np.where(num > 0, num ** 2 / (den[:, None] * dd[None, :]), 0.0)
    vb = vr.max(1)
    k = int(np.argmax(vb)); l = int(np.argmax(vr[k]))
    return k, int(LAGS[l]), float(vb[k]), vb, len(idx)


def main():
    years = [int(a) for a in sys.argv[1:]] or list(range(2010, 2019))
    E = pd.read_csv(f'{OUT}/tables/t2_zsurface_all.csv')
    Wd = {y: np.load(f'{OUT}/data/win_{y}.npz') for y in years}
    Z = np.load(f'{OUT}/data/zsurf_2011.npz')
    nla, nlo, ndep = Z['lat'], Z['lon'], Z['dep']
    nodes = list(zip(nla, nlo, ndep))
    allst = sorted({s for y in years for s in Wd[y]['stas']})
    coords = {}
    for y in years:
        for s, ll in zip(Wd[y]['stas'], Wd[y]['sll']):
            coords[s] = tuple(ll)
    sll = np.array([coords[s] for s in allst])
    # station terms
    ft = f'{OUT}/tables/t4_station_terms.csv'
    if os.path.isfile(ft):
        TM = pd.read_csv(ft)
    else:
        TM, TR = stack_terms(allst, coords)
        TM.to_csv(ft, index=False); TR.to_csv(f'{OUT}/tables/t4_station_terms_per_stack.csv', index=False)
    print(TM.round(3).to_string())
    # template banks over all stations (nodes, nsta_all, 3, 181)
    fb = f'{OUT}/data/t4_templates.npz'
    if os.path.isfile(fb):
        B = dict(np.load(fb))
    else:
        B = dict(ok=np.array([static_templates(n, sll, S.TAU, 'ok') for n in nodes]),
                 wsst=np.array([static_templates(n, sll, S.TAU, 'ws') for n in nodes]))
        np.savez_compressed(fb, **B)
    rows = []
    nx, ny = S.kmxy(nla, nlo)
    for y in years:
        W = Wd[y]; stas = list(W['stas']); si = [allst.index(s) for s in stas]
        Tws = W['T_thrust']                                         # production (nodes, nsta, 3, 181), x1e21
        term = {}
        for kind in ('ws', 'ok'):
            t = np.ones((len(stas), 3))
            for r in TM[TM.kind == kind].itertuples():
                if r.sta in stas:
                    t[stas.index(r.sta), 'NEZ'.index(r.comp)] = r.used
            term[kind] = t
        Ey = E[E.year == y].set_index('eid')
        for j, eid in enumerate(W['eid']):
            r = Ey.loc[eid]
            g = W['g'][j]
            V = {'ws': Tws * g[None, :, :, None], 'ws_st': Tws * term['ws'][None, :, :, None],
                 'ok': B['ok'][:, si] * g[None, :, :, None], 'ok_g1': B['ok'][:, si],
                 'ok_st': B['ok'][:, si] * term['ok'][None, :, :, None], 'wsst': B['wsst'][:, si] * g[None, :, :, None]}
            o = dict(eid=eid, year=y, period=r.period, episode=r.episode, z_stored=r.z_stored, node_prod=int(r.node_stored),
                     c10_lat=r.c10_lat, c10_lon=r.c10_lon, lat_prod=r.lat_best, lon_prod=r.lon_best)
            for nm, Tg in V.items():
                res = relocate(W['x'][j], W['sig'][j], W['won'][j], W['sv'][j], Tg)
                if res is None:
                    continue
                k, L, vr, vb, nch = res
                o[f'{nm}_node'] = k; o[f'{nm}_lat'] = nla[k]; o[f'{nm}_lon'] = nlo[k]; o[f'{nm}_lag'] = L
                o[f'{nm}_vr'] = vr; o[f'{nm}_nch'] = nch
                o[f'{nm}_n_vr002'] = int((vb >= vr - 0.02).sum())
                if np.isfinite(r.c10_lat):
                    cx, cy = S.kmxy(r.c10_lat, r.c10_lon)
                    o[f'{nm}_dE'] = nx[k] - cx; o[f'{nm}_dN'] = ny[k] - cy; o[f'{nm}_dist'] = float(np.hypot(nx[k] - cx, ny[k] - cy))
            rows.append(o)
    R = pd.DataFrame(rows)
    R.to_csv(f'{OUT}/tables/t4_reloc.csv', index=False)
    px, py = S.kmxy(R.lat_prod, R.lon_prod); cx, cy = S.kmxy(R.c10_lat, R.c10_lon)
    R['prod_dist'] = np.hypot(px - cx, py - cy); R['prod_dE'] = px - cx; R['prod_dN'] = py - cy
    R.to_csv(f'{OUT}/tables/t4_reloc.csv', index=False)
    summ = []
    vet = set(pd.read_csv(f'{S13}/catalogue/tables/e3_events.csv').epoch.round().astype(int))
    Ev = E.set_index('eid').epoch.round().astype(int)
    R['vetted'] = [any(abs(Ev[e] - v) <= 30 for v in vet) for e in R.eid]
    for subset, g0 in (('2011-2013 ETS with LFE', R[(R.year.between(2011, 2013)) & (R.period == 'ETS') & np.isfinite(R.c10_lat)]),
                       ('vetted+near-miss 2011-12', R[R.vetted & np.isfinite(R.c10_lat)]),
                       ('2010-2018 ETS with LFE', R[(R.period == 'ETS') & np.isfinite(R.c10_lat)])):
        for nm in ('prod', 'ws', 'ws_st', 'ok', 'ok_g1', 'ok_st', 'wsst'):
            if f'{nm}_dist' not in g0:
                continue
            g = g0[np.isfinite(g0[f'{nm}_dist'])]
            d = g[f'{nm}_dist']
            o = dict(subset=subset, variant=nm, n=len(g), median_dist=d.median(), mean_dist=d.mean(),
                     mean_dE=g[f'{nm}_dE'].mean(), mean_dN=g[f'{nm}_dN'].mean(), sd_dE=g[f'{nm}_dE'].std(), sd_dN=g[f'{nm}_dN'].std(),
                     frac_closer_than_ws=float((d < g.ws_dist).mean()) if nm != 'ws' else np.nan,
                     frac_farther_than_ws=float((d > g.ws_dist).mean()) if nm != 'ws' else np.nan,
                     same_node_as_prod=float((g[f'{nm}_node'] == g.node_prod).mean()) if nm != 'prod' else 1.0,
                     median_vr=g[f'{nm}_vr'].median() if nm != 'prod' else np.nan,
                     median_n_vr002=g[f'{nm}_n_vr002'].median() if nm != 'prod' else np.nan)
            summ.append(o)
    Sm = pd.DataFrame(summ); Sm.to_csv(f'{OUT}/tables/t4_summary.csv', index=False)
    pd.set_option('display.width', 250)
    print(Sm.round(3).to_string())


if __name__ == '__main__':
    main()
