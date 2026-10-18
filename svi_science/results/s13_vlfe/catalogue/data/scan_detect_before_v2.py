#!/usr/bin/env python
"""VLFE scan detector (PLAN_SCAN.md §4-5), GPU (cupy) FFT matched filter on the source grid.

usage: scan_detect.py RUN [--bp DIR] [--ndraw 50] [--shift S]
  RUN        output name (data/scan/runs/RUN/)
  --bp       band-passed year directory (scan_bp.py output), default data/scan/bp_2011
  --ndraw    polarity-scrambled template draws computed in the same pass (control b); 0 for injection / control a
  --shift    control (a): station k (PILOT order) uses its data from day (d + shift*k) mod 365

Outputs (per run): cand.parquet (variant 0 = real signs, 1..ndraw = scrambled; positive and negative peaks, |z| >= 4,
local extrema within +-150 s over time of the max/min over grid nodes), madscale.npy (variants, nodes, days),
signs.npy, and for variant 0: zmax.npy, zmin.npy, argmax.npy, m0_centroid.npy (node 289 = S13 2011-08 centroid),
navail.npy (stations available per second), daycfg.parquet (weights/gains per station-component-day). Marker DONE.
"""
import os, sys, argparse, time, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S

ZTHR = 4.0
DECL = 150
NF = 131072
PAD = 300                                    # segment = day +- 300 s
I0 = PAD + S.FIT0                            # correlation index of origin u=0 (240)


def gain_table():
    G = pd.read_csv(f'{C.TB}/sc_gain.csv')
    Nm = pd.read_csv(S.tab('scan_noise_month.csv'))          # S.tab: year-specific (Amendment S8)
    excl = {(r.sta, int(r.month)) for r in Nm[Nm.exclude].itertuples()}
    return G, excl


def day_cfg(G, excl, stas, days):
    """w-independent config: gain g (0 if excluded) per (station, comp, day of source data)."""
    g = np.zeros((len(stas), 3, S.NDAY))
    for i, s in enumerate(stas):
        for d in range(S.NDAY):
            day = S.Y0 + pd.Timedelta(days=d)
            if (s, day.month) in excl:
                continue
            ds = str(day)
            for c, cn in enumerate('NEZ'):
                r = G[(G.sta == s) & (G.comp == cn) & (G.epoch_start <= ds) & (G.epoch_end >= ds)]
                if len(r) and bool(r.keep.iloc[0]):
                    g[i, c, d] = r.g.iloc[0]
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('--bp', default=f'{S.SD}/bp_2011')
    ap.add_argument('--ndraw', type=int, default=50); ap.add_argument('--shift', type=int, default=0)
    ap.add_argument('--days', default='')
    a = ap.parse_args()
    import cupy as cp
    out = f'{S.SD}/runs/{a.run}'
    os.makedirs(out, exist_ok=True)
    lf = f'{C.LOGD}/scan_detect.log'
    TP = np.load(f'{S.SD}/templates_grid.npz')
    stas = list(TP['stations'])
    nst = len(stas)
    T = TP['thrust'] * 1e21                                  # (nodes, sta, 3, 181) nm per 1e12 N m (float32-safe units)
    nn = T.shape[0]
    ndet = nn - 1                                            # last node = centroid, excluded from max
    E = (T ** 2).sum(axis=3)                                 # (nodes, sta, 3)
    kern = np.zeros((nn, nst, 3, NF), np.float32); kern[..., :S.NT] = T
    Kf = cp.conj(cp.fft.rfft(cp.asarray(kern), axis=-1)).astype(cp.complex64)     # (nodes, sta, 3, NF/2+1)
    del kern
    Ec = cp.asarray(E.astype(np.float32))
    Gt, excl = gain_table()
    gcfg = day_cfg(Gt, excl, stas, None)
    offs = np.array([a.shift * k if s in S.PILOT[:7] else 0 for k, s in enumerate(stas)])
    BP = {s: np.load(f'{a.bp}/{s}.npy', mmap_mode='r') for s in stas}
    DS = {s: np.load(f'{a.bp}/{s}_daystats.npz') for s in stas}
    mask = S.mask_vector('head')
    rng = np.random.default_rng(20261014)
    nv = 1 + a.ndraw
    signs = np.ones((nv, nst, 3), np.float32)
    rhos = [1.0]
    if a.ndraw:                                   # Amendment S5.2: draws with |rho| <= 0.3 vs the true template
        wbar = np.zeros((nst, 3))
        for i, s in enumerate(stas):
            sg_ = DS[s]['sig'] * 1e9
            for c in range(3):
                ok = (gcfg[i, c] > 0) & np.isfinite(sg_[:, c]) & (sg_[:, c] > 0)
                if ok.any():
                    wbar[i, c] = np.median(1.0 / sg_[ok, c] ** 2) * np.median(gcfg[i, c][ok]) ** 2
        Ew = wbar * E[-1]
        k = 1
        while k < nv:
            s_ = rng.choice([-1.0, 1.0], size=(nst, 3))
            rho = (s_ * Ew).sum() / Ew.sum()
            if abs(rho) <= 0.3:
                signs[k] = s_; rhos.append(float(rho)); k += 1
    np.save(f'{out}/rho.npy', np.array(rhos))
    np.save(f'{out}/signs.npy', signs)
    days = range(S.NDAY) if not a.days else [int(x) for x in a.days.split(',')]
    zmax_y = np.full(S.NDAY * 86400, np.nan, np.float32); zmin_y = zmax_y.copy()
    arg_y = np.full(S.NDAY * 86400, -1, np.int16); m0c_y = zmax_y.copy(); navail_y = np.zeros(S.NDAY * 86400, np.int8)
    madscale = np.full((nv, nn, S.NDAY), np.nan, np.float32)
    cands = []
    cfgrows = []
    t_start = time.time()
    for d in days:
        seg = np.full((nst, 3, 86400 + 2 * PAD), np.nan, np.float32)
        w = np.zeros((nst, 3)); g = np.zeros((nst, 3))
        for i, s in enumerate(stas):
            ds = (d + offs[i]) % S.NDAY
            lo, hi = ds * 86400 - PAD, (ds + 1) * 86400 + PAD
            a0, b0 = max(lo, 0), min(hi, S.NDAY * 86400)
            seg[i, :, a0 - lo:b0 - lo] = BP[s][:, a0:b0] * 1e9                    # nm
            sig = DS[s]['sig'][ds]
            for c in range(3):
                gg = gcfg[i, c, ds]
                if gg > 0 and np.isfinite(sig[c]) and sig[c] > 0:
                    w[i, c] = 1.0 / (sig[c] * 1e9) ** 2; g[i, c] = gg
            cfgrows += [dict(day=d, src_day=int(ds), sta=s, comp='NEZ'[c], g=g[i, c], sigma_nm=sig[c] * 1e9,
                             w_on=bool(w[i, c] > 0)) for c in range(3)]
        act = w > 0
        if act.sum() == 0:
            continue
        fin = np.isfinite(seg)
        # availability: window [u+240, u+420] fully finite (per station-component)
        bad = np.concatenate([np.zeros((nst, 3, 1)), np.cumsum(~fin, axis=2)], axis=2)
        u = np.arange(86400)
        avail = (bad[:, :, u + I0 + S.NT] - bad[:, :, u + I0]) == 0              # (sta, 3, 86400)
        avail &= act[:, :, None]
        sta_av = avail.any(axis=1)                                                 # (sta, 86400)
        nav = sta_av.sum(axis=0).astype(np.int8)
        if a.ndraw >= 0:
            navail_y[d * 86400:(d + 1) * 86400] = nav
        D = cp.asarray(np.where(fin, seg, 0.0).astype(np.float32))
        Df = cp.fft.rfft(D, n=NF, axis=-1)                                         # (sta, 3, F)
        av = cp.asarray(avail.astype(np.float32))
        wg = cp.asarray((w * g).astype(np.float32)); wg2 = cp.asarray((w * g * g).astype(np.float32))
        coef = av * wg[:, :, None]                                                 # (sta, 3, U)
        cw2 = av * wg2[:, :, None]
        sg = cp.asarray(signs)                                                     # (v, sta, 3)
        zsel = cp.zeros((nst, 3), cp.float32); zsel[:, 2] = 1
        msk = cp.asarray(mask[d * 86400:(d + 1) * 86400])
        sub = cp.arange(0, 86400, 4)
        zmx = cp.full((nv, 86400), -cp.inf, cp.float32); zmn = cp.full((nv, 86400), cp.inf, cp.float32)
        amx = cp.zeros((nv, 86400), cp.int32); amn = cp.zeros((nv, 86400), cp.int32)
        m0mx = cp.zeros((nv, 86400), cp.float32); m0mn = cp.zeros((nv, 86400), cp.float32)
        zZmx = cp.zeros(86400, cp.float32)
        CH = 29
        for n0 in range(0, nn, CH):
            n1 = min(n0 + CH, nn)
            corr = cp.fft.irfft(Kf[n0:n1] * Df[None], n=NF, axis=-1)[..., I0:I0 + 86400]    # c[i] = sum_j k[j] d[i+j], i = u + I0
            X = corr * coef[None]                                                  # (n, sta, 3, U)
            den = cp.einsum('nsc,scu->nu', Ec[n0:n1], cw2)
            num = cp.einsum('vsc,nscu->vnu', sg, X)
            dd = cp.where(den > 0, den, cp.nan)[None]
            m0 = num / dd                                                          # moment estimate (N m)
            st = num / cp.sqrt(dd)                                                 # Amendment S3: availability-normalised
            ok = (~msk[sub]) & cp.isfinite(st[0, 0, sub])
            if float(ok.mean()) < 0.1:
                ok = cp.isfinite(st[0, 0, sub])
            sts = cp.where(ok[None, None], st[:, :, sub], cp.nan)
            med = cp.nanmedian(sts, axis=2, keepdims=True)
            mad = 1.4826 * cp.nanmedian(cp.abs(sts - med), axis=2)                 # (v, n)
            madscale[:, n0:n1, d] = cp.asnumpy(mad)
            z = st / mad[:, :, None]
            numZ = cp.einsum('sc,nscu->nu', zsel, X)
            denZ = cp.einsum('nsc,scu->nu', Ec[n0:n1] * zsel[None], cw2)
            stZ = numZ / cp.sqrt(cp.where(denZ > 0, denZ, cp.nan))
            okZ = (~msk[sub]) & cp.isfinite(stZ[0, sub])
            if float(okZ.mean()) < 0.1:
                okZ = cp.isfinite(stZ[0, sub])
            sZ = cp.where(okZ[None], stZ[:, sub], cp.nan)
            madZ = 1.4826 * cp.nanmedian(cp.abs(sZ - cp.nanmedian(sZ, axis=1, keepdims=True)), axis=1)
            zZall = stZ / madZ[:, None]
            if n1 == nn:                                                           # centroid node: store, exclude
                m0c_y[d * 86400:(d + 1) * 86400] = cp.asnumpy(m0[0, -1]) * 1e12
                z = z[:, :-1]; m0 = m0[:, :-1]; zZall = zZall[:-1]
                if z.shape[1] == 0:
                    continue
            zf = cp.nan_to_num(z, nan=-cp.inf)
            k = cp.argmax(zf, axis=1); zm = cp.take_along_axis(zf, k[:, None], 1)[:, 0]
            upd = zm > zmx
            zmx = cp.where(upd, zm, zmx); amx = cp.where(upd, k + n0, amx)
            m0mx = cp.where(upd, cp.take_along_axis(m0, k[:, None], 1)[:, 0], m0mx)
            zZmx = cp.where(upd[0], zZall[k[0], cp.arange(86400)], zZmx)
            zn = cp.nan_to_num(z, nan=cp.inf)
            k2 = cp.argmin(zn, axis=1); zmi = cp.take_along_axis(zn, k2[:, None], 1)[:, 0]
            upd2 = zmi < zmn
            zmn = cp.where(upd2, zmi, zmn); amn = cp.where(upd2, k2 + n0, amn)
            m0mn = cp.where(upd2, cp.take_along_axis(m0, k2[:, None], 1)[:, 0], m0mn)
            del corr, X, num, den, m0, z, zf, zn, st, sts, stZ, zZall
        zmx_h = cp.asnumpy(zmx); zmn_h = cp.asnumpy(zmn); amx_h = cp.asnumpy(amx); amn_h = cp.asnumpy(amn)
        m0mx_h = cp.asnumpy(m0mx) * 1e12; m0mn_h = cp.asnumpy(m0mn) * 1e12; zZ_h = cp.asnumpy(zZmx)   # back to N m
        zmax_y[d * 86400:(d + 1) * 86400] = np.where(np.isfinite(zmx_h[0]), zmx_h[0], np.nan)
        zmin_y[d * 86400:(d + 1) * 86400] = np.where(np.isfinite(zmn_h[0]), zmn_h[0], np.nan)
        arg_y[d * 86400:(d + 1) * 86400] = amx_h[0]
        from scipy.ndimage import maximum_filter1d, minimum_filter1d
        for v in range(nv):
            for sgn, zz, aa, mm in ((1, zmx_h[v], amx_h[v], m0mx_h[v]), (-1, zmn_h[v], amn_h[v], m0mn_h[v])):
                x = np.where(np.isfinite(zz), zz, -np.inf if sgn > 0 else np.inf)
                if sgn > 0:
                    pk = (x >= ZTHR) & (x == maximum_filter1d(x, 2 * DECL + 1, mode='nearest'))
                else:
                    pk = (x <= -ZTHR) & (x == minimum_filter1d(x, 2 * DECL + 1, mode='nearest'))
                for i in np.where(pk)[0]:
                    cands.append((v, sgn, d * 86400 + int(i), int(aa[i]), float(zz[i]), float(mm[i]),
                                  float(zZ_h[i]) if (v == 0 and sgn > 0) else np.nan, int(nav[i])))
        if d % 20 == 0:
            C.log(f'{a.run}: day {d} done, {len(cands)} raw peaks, {time.time() - t_start:.0f} s', lf)
    Cd = pd.DataFrame(cands, columns=['variant', 'sign', 'tsec', 'node', 'z', 'm0', 'zZ', 'navail'])
    # decluster across day boundaries (+-150 s, keep max |z|) per variant and sign
    keep = []
    for (v, sgn), g_ in Cd.groupby(['variant', 'sign']):
        g_ = g_.sort_values('tsec')
        t = g_.tsec.values; zz = np.abs(g_.z.values); idx = g_.index.values
        kk = np.ones(len(g_), bool)
        for j in range(1, len(g_)):
            if t[j] - t[j - 1] <= DECL:
                if zz[j] > zz[j - 1]:
                    kk[j - 1] = False
                else:
                    kk[j] = False
        keep += list(idx[kk])
    Cd = Cd.loc[sorted(keep)].reset_index(drop=True)
    Cd['shift'] = a.shift
    Cd.to_parquet(f'{out}/cand.parquet')
    np.save(f'{out}/madscale.npy', madscale)
    np.save(f'{out}/zmax.npy', zmax_y); np.save(f'{out}/zmin.npy', zmin_y); np.save(f'{out}/argmax.npy', arg_y)
    np.save(f'{out}/m0_centroid.npy', m0c_y); np.save(f'{out}/navail.npy', navail_y)
    pd.DataFrame(cfgrows).to_parquet(f'{out}/daycfg.parquet')
    json.dump(dict(run=a.run, bp=a.bp, ndraw=a.ndraw, shift=a.shift, offsets=offs.tolist(), stations=stas,
                   n_cand=len(Cd), seconds=time.time() - t_start), open(f'{out}/meta.json', 'w'), indent=1)
    open(f'{out}/DONE', 'w').write('ok\n')
    C.log(f'{a.run}: DONE {len(Cd)} candidates, {time.time() - t_start:.0f} s', lf)


if __name__ == '__main__':
    main()
