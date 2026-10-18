#!/usr/bin/env python
"""VLFE scan pilot 2 detector (PLAN_SCAN.md S7.5): the pilot-1 detector (scan_detect.py, §4 + S3 + S4; unchanged maths)
preceded by the single-station step cleaner, plus z maximised over footprint-interior nodes (column).

Pass 1 (cleaner), per UTC day of the run timeline: z_s(t, n) = sum_c w g (T * d) / sqrt(sum_c w g^2 |T|^2) for one station
s (S3 statistic restricted to one station), normalised by its own daily MAD (unmasked, every 4th sample). For each
station and second, n* = argmax_n |z_s|; flag if |z_s(t, n*)| >= 8 and every other available station has
|z_s'(t, n*)| < 2.5, with >= 2 other stations available. Flags dilated +-600 s -> the station's three components NaN.
Pass 2: detection on the cleaned data (all variants of the run). --noclean skips pass 1 (code check only).

usage: scan_detect2.py RUN [--bp DIR] [--ndraw 50] [--shift S] [--days d1,d2,..] [--noclean]
Outputs data/scan/runs/RUN/: as scan_detect.py (cand.parquet, madscale.npy, signs.npy, rho.npy, zmax.npy, zmin.npy,
argmax.npy, m0_centroid.npy, navail.npy, daycfg.parquet, meta.json) + zmax_fp.npy, argmax_fp.npy, clean_iv.parquet
(sta, t0, t1 in run-timeline seconds, dilated), clean_flags.csv (flagged seconds per station). Marker DONE."""
import os, sys, argparse, time, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S
from scan_detect import gain_table, day_cfg, ZTHR, DECL, NF, PAD, I0

ZS_HI, ZS_LO, CL_DIL = 8.0, 2.5, 600


def footprint_nodes(lat, lon):
    from matplotlib.path import Path
    P = pd.read_csv(f'{C.TB}/scan_footprint_polygon.csv')
    return Path(np.c_[P.lon, P.lat]).contains_points(np.c_[lon, lat])


def build_seg(stas, BP, DS, gcfg, offs, d, extra_nan=None):
    nst = len(stas)
    seg = np.full((nst, 3, 86400 + 2 * PAD), np.nan, np.float32)
    w = np.zeros((nst, 3)); g = np.zeros((nst, 3)); rows = []
    for i, s in enumerate(stas):
        ds = (d + offs[i]) % S.NDAY
        lo, hi = ds * 86400 - PAD, (ds + 1) * 86400 + PAD
        a0, b0 = max(lo, 0), min(hi, S.NDAY * 86400)
        seg[i, :, a0 - lo:b0 - lo] = BP[s][:, a0:b0] * 1e9
        sig = DS[s]['sig'][ds]
        for c in range(3):
            gg = gcfg[i, c, ds]
            if gg > 0 and np.isfinite(sig[c]) and sig[c] > 0:
                w[i, c] = 1.0 / (sig[c] * 1e9) ** 2; g[i, c] = gg
            rows.append(dict(day=d, src_day=int(ds), sta=s, comp='NEZ'[c], g=g[i, c], sigma_nm=sig[c] * 1e9,
                             w_on=bool(w[i, c] > 0)))
    if extra_nan is not None:                                  # cleaner mask on the run timeline (sta, 86400+2PAD)
        seg[np.broadcast_to(extra_nan[:, None, :], seg.shape)] = np.nan
    return seg, w, g, rows


def availability(seg, act):
    fin = np.isfinite(seg)
    nst = seg.shape[0]
    bad = np.concatenate([np.zeros((nst, 3, 1)), np.cumsum(~fin, axis=2)], axis=2)
    u = np.arange(86400)
    avail = (bad[:, :, u + I0 + S.NT] - bad[:, :, u + I0]) == 0
    avail &= act[:, :, None]
    return fin, avail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('--bp', default=f'{S.SD}/bp_2011g')
    ap.add_argument('--ndraw', type=int, default=50); ap.add_argument('--shift', type=int, default=0)
    ap.add_argument('--days', default=''); ap.add_argument('--noclean', action='store_true')
    a = ap.parse_args()
    import cupy as cp
    from scipy.ndimage import maximum_filter1d, minimum_filter1d, binary_dilation
    out = f'{S.SD}/runs/{a.run}'
    os.makedirs(out, exist_ok=True)
    lf = f'{C.LOGD}/scan2_detect.log'
    TP = np.load(f'{S.SD}/templates_grid.npz')
    stas = list(TP['stations']); nst = len(stas)
    T = TP['thrust'] * 1e21
    nn = T.shape[0]
    fpn = footprint_nodes(TP['lat'], TP['lon']); fpn[-1] = False             # centroid node excluded
    fp_c = cp.asarray(fpn)
    E = (T ** 2).sum(axis=3)
    kern = np.zeros((nn, nst, 3, NF), np.float32); kern[..., :S.NT] = T
    Kf = cp.conj(cp.fft.rfft(cp.asarray(kern), axis=-1)).astype(cp.complex64)
    del kern
    Ec = cp.asarray(E.astype(np.float32))
    Gt, excl = gain_table()
    gcfg = day_cfg(Gt, excl, stas, None)
    offs = np.array([a.shift * k if s in S.PILOT[:7] else 0 for k, s in enumerate(stas)])
    BP = {s: np.load(f'{a.bp}/{s}.npy', mmap_mode='r') for s in stas}
    DS = {s: np.load(f'{a.bp}/{s}_daystats.npz') for s in stas}
    mask = S.mask_vector('head')
    days = list(range(S.NDAY)) if not a.days else [int(x) for x in a.days.split(',')]
    t_start = time.time()
    # ---------------- pass 1: single-station step cleaner ----------------
    flag_y = np.zeros((nst, S.NDAY * 86400), bool)
    if not a.noclean:
        for d in days:
            seg, w, g, _ = build_seg(stas, BP, DS, gcfg, offs, d)
            act = w > 0
            if act.sum() == 0:
                continue
            fin, avail = availability(seg, act)
            sta_av = avail.any(axis=1)                                        # (sta, U)
            nav = sta_av.sum(axis=0)
            D = cp.asarray(np.where(fin, seg, 0.0).astype(np.float32))
            Df = cp.fft.rfft(D, n=NF, axis=-1)
            av = cp.asarray(avail.astype(np.float32))
            coef = av * cp.asarray((w * g).astype(np.float32))[:, :, None]
            cw2 = av * cp.asarray((w * g * g).astype(np.float32))[:, :, None]
            msk = cp.asarray(mask[d * 86400:(d + 1) * 86400]); sub = cp.arange(0, 86400, 4)
            best = cp.zeros((nst, 86400), cp.float32); oth = cp.zeros((nst, 86400), cp.float32)
            sidx = cp.arange(nst)[None, :, None]
            for n0 in range(0, nn - 1, 29):
                n1 = min(n0 + 29, nn - 1)
                corr = cp.fft.irfft(Kf[n0:n1] * Df[None], n=NF, axis=-1)[..., I0:I0 + 86400]
                num = (corr * coef[None]).sum(axis=2)                                 # (n, sta, U)
                den = cp.einsum('nsc,scu->nsu', Ec[n0:n1], cw2)
                st = num / cp.sqrt(cp.where(den > 0, den, cp.nan))
                ok = (~msk[sub])[None, None] & cp.isfinite(st[:, :, sub])
                sts = cp.where(ok, st[:, :, sub], cp.nan)
                madv = 1.4826 * cp.nanmedian(cp.abs(sts - cp.nanmedian(sts, axis=2, keepdims=True)), axis=2)
                z = cp.abs(st / madv[:, :, None])
                z = cp.nan_to_num(z, nan=0.0)                                         # unavailable -> 0
                i1 = cp.argmax(z, axis=1)                                             # (n, U)
                m1 = cp.take_along_axis(z, i1[:, None], 1)[:, 0]
                z2 = cp.where(sidx == i1[:, None], -1.0, z)
                m2 = z2.max(axis=1)
                others = cp.where(sidx == i1[:, None], m2[:, None], m1[:, None])     # (n, sta, U)
                kb = cp.argmax(z, axis=0)                                             # best node per (sta, U)
                zb = cp.take_along_axis(z, kb[None], 0)[0]
                ob = cp.take_along_axis(others, kb[None], 0)[0]
                upd = zb > best
                best = cp.where(upd, zb, best); oth = cp.where(upd, ob, oth)
                del corr, num, den, st, sts, z, z2, others
            fl = (cp.asnumpy(best) >= ZS_HI) & (cp.asnumpy(oth) < ZS_LO) & (nav[None, :] - 1 >= 2) & sta_av
            flag_y[:, d * 86400:(d + 1) * 86400] = fl
            if d % 40 == 0:
                C.log(f'{a.run}: clean pass day {d}, flagged s so far {int(flag_y.sum())}, {time.time() - t_start:.0f} s', lf)
    cl_y = np.zeros_like(flag_y)
    ivrows = []
    for i in range(nst):
        if flag_y[i].any():
            cs = np.concatenate([[0], np.cumsum(flag_y[i], dtype=np.int64)])
            n_ = flag_y.shape[1]; ii = np.arange(n_)
            cl_y[i] = (cs[np.minimum(ii + CL_DIL + 1, n_)] - cs[np.maximum(ii - CL_DIL, 0)]) > 0     # +-600 s dilation
            r0, r1 = S.runs(cl_y[i])
            ivrows += [dict(sta=stas[i], t0=int(x), t1=int(y)) for x, y in zip(r0, r1)]
    pd.DataFrame(ivrows, columns=['sta', 't0', 't1']).to_parquet(f'{out}/clean_iv.parquet')
    pd.DataFrame(dict(sta=stas, flagged_s=flag_y.sum(axis=1), cleaned_s=cl_y.sum(axis=1),
                      n_intervals=[sum(1 for r in ivrows if r['sta'] == s) for s in stas])).to_csv(f'{out}/clean_flags.csv', index=False)
    del flag_y
    C.log(f'{a.run}: cleaner done, {len(ivrows)} intervals, {time.time() - t_start:.0f} s', lf)
    # ---------------- pass 2: detection on cleaned data (pilot-1 maths) ----------------
    rng = np.random.default_rng(20261014)
    nv = 1 + a.ndraw
    signs = np.ones((nv, nst, 3), np.float32); rhos = [1.0]
    if a.ndraw:
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
    np.save(f'{out}/rho.npy', np.array(rhos)); np.save(f'{out}/signs.npy', signs)
    zmax_y = np.full(S.NDAY * 86400, np.nan, np.float32); zmin_y = zmax_y.copy(); zfp_y = zmax_y.copy()
    arg_y = np.full(S.NDAY * 86400, -1, np.int16); argfp_y = arg_y.copy(); m0c_y = zmax_y.copy()
    navail_y = np.zeros(S.NDAY * 86400, np.int8)
    madscale = np.full((nv, nn, S.NDAY), np.nan, np.float32)
    cands = []; cfgrows = []
    for d in days:
        lo = d * 86400 - PAD
        idx = np.arange(lo, lo + 86400 + 2 * PAD) % (S.NDAY * 86400)
        seg, w, g, rows = build_seg(stas, BP, DS, gcfg, offs, d, extra_nan=cl_y[:, idx])
        cfgrows += rows
        act = w > 0
        if act.sum() == 0:
            continue
        fin, avail = availability(seg, act)
        nav = avail.any(axis=1).sum(axis=0).astype(np.int8)
        navail_y[d * 86400:(d + 1) * 86400] = nav
        D = cp.asarray(np.where(fin, seg, 0.0).astype(np.float32))
        Df = cp.fft.rfft(D, n=NF, axis=-1)
        av = cp.asarray(avail.astype(np.float32))
        coef = av * cp.asarray((w * g).astype(np.float32))[:, :, None]
        cw2 = av * cp.asarray((w * g * g).astype(np.float32))[:, :, None]
        sg = cp.asarray(signs)
        zsel = cp.zeros((nst, 3), cp.float32); zsel[:, 2] = 1
        msk = cp.asarray(mask[d * 86400:(d + 1) * 86400]); sub = cp.arange(0, 86400, 4)
        zmx = cp.full((nv, 86400), -cp.inf, cp.float32); zmn = cp.full((nv, 86400), cp.inf, cp.float32)
        amx = cp.zeros((nv, 86400), cp.int32); amn = cp.zeros((nv, 86400), cp.int32)
        m0mx = cp.zeros((nv, 86400), cp.float32); m0mn = cp.zeros((nv, 86400), cp.float32)
        zZmx = cp.zeros(86400, cp.float32)
        zfp = cp.full(86400, -cp.inf, cp.float32); afp = cp.zeros(86400, cp.int32)
        for n0 in range(0, nn, 29):
            n1 = min(n0 + 29, nn)
            corr = cp.fft.irfft(Kf[n0:n1] * Df[None], n=NF, axis=-1)[..., I0:I0 + 86400]
            X = corr * coef[None]
            den = cp.einsum('nsc,scu->nu', Ec[n0:n1], cw2)
            num = cp.einsum('vsc,nscu->vnu', sg, X)
            dd = cp.where(den > 0, den, cp.nan)[None]
            m0 = num / dd
            st = num / cp.sqrt(dd)
            ok = (~msk[sub]) & cp.isfinite(st[0, 0, sub])
            if float(ok.mean()) < 0.1:
                ok = cp.isfinite(st[0, 0, sub])
            sts = cp.where(ok[None, None], st[:, :, sub], cp.nan)
            med = cp.nanmedian(sts, axis=2, keepdims=True)
            mad = 1.4826 * cp.nanmedian(cp.abs(sts - med), axis=2)
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
            fpc = fp_c[n0:n1]
            if n1 == nn:
                m0c_y[d * 86400:(d + 1) * 86400] = cp.asnumpy(m0[0, -1]) * 1e12
                z = z[:, :-1]; m0 = m0[:, :-1]; zZall = zZall[:-1]; fpc = fpc[:-1]
                if z.shape[1] == 0:
                    continue
            zf = cp.nan_to_num(z, nan=-cp.inf)
            k = cp.argmax(zf, axis=1); zm = cp.take_along_axis(zf, k[:, None], 1)[:, 0]
            upd = zm > zmx
            zmx = cp.where(upd, zm, zmx); amx = cp.where(upd, k + n0, amx)
            m0mx = cp.where(upd, cp.take_along_axis(m0, k[:, None], 1)[:, 0], m0mx)
            zZmx = cp.where(upd[0], zZall[k[0], cp.arange(86400)], zZmx)
            if bool(fpc.any()):
                zf0 = cp.where(fpc[:, None], zf[0], -cp.inf)
                kf = cp.argmax(zf0, axis=0); zmf = cp.take_along_axis(zf0, kf[None], 0)[0]
                u2 = zmf > zfp
                zfp = cp.where(u2, zmf, zfp); afp = cp.where(u2, kf + n0, afp)
            zn = cp.nan_to_num(z, nan=cp.inf)
            k2 = cp.argmin(zn, axis=1); zmi = cp.take_along_axis(zn, k2[:, None], 1)[:, 0]
            upd2 = zmi < zmn
            zmn = cp.where(upd2, zmi, zmn); amn = cp.where(upd2, k2 + n0, amn)
            m0mn = cp.where(upd2, cp.take_along_axis(m0, k2[:, None], 1)[:, 0], m0mn)
            del corr, X, num, den, m0, z, zf, zn, st, sts, stZ, zZall
        zmx_h = cp.asnumpy(zmx); zmn_h = cp.asnumpy(zmn); amx_h = cp.asnumpy(amx); amn_h = cp.asnumpy(amn)
        m0mx_h = cp.asnumpy(m0mx) * 1e12; m0mn_h = cp.asnumpy(m0mn) * 1e12; zZ_h = cp.asnumpy(zZmx)
        sl = slice(d * 86400, (d + 1) * 86400)
        zmax_y[sl] = np.where(np.isfinite(zmx_h[0]), zmx_h[0], np.nan)
        zmin_y[sl] = np.where(np.isfinite(zmn_h[0]), zmn_h[0], np.nan)
        arg_y[sl] = amx_h[0]
        zf_h = cp.asnumpy(zfp); zfp_y[sl] = np.where(np.isfinite(zf_h), zf_h, np.nan); argfp_y[sl] = cp.asnumpy(afp)
        for v in range(nv):
            for sgn, zz, aa, mm in ((1, zmx_h[v], amx_h[v], m0mx_h[v]), (-1, zmn_h[v], amn_h[v], m0mn_h[v])):
                x = np.where(np.isfinite(zz), zz, -np.inf if sgn > 0 else np.inf)
                if sgn > 0:
                    pk = (x >= ZTHR) & (x == maximum_filter1d(x, 2 * DECL + 1, mode='nearest'))
                else:
                    pk = (x <= -ZTHR) & (x == minimum_filter1d(x, 2 * DECL + 1, mode='nearest'))
                for i in np.where(pk)[0]:
                    cands.append((v, sgn, d * 86400 + int(i), int(aa[i]), float(zz[i]), float(mm[i]),
                                  float(zZ_h[i]) if (v == 0 and sgn > 0) else np.nan, int(nav[i]),
                                  float(zf_h[i]) if v == 0 else np.nan))
        if d % 20 == 0:
            C.log(f'{a.run}: day {d} done, {len(cands)} raw peaks, {time.time() - t_start:.0f} s', lf)
    Cd = pd.DataFrame(cands, columns=['variant', 'sign', 'tsec', 'node', 'z', 'm0', 'zZ', 'navail', 'z_fp'])
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
    np.save(f'{out}/zmax_fp.npy', zfp_y); np.save(f'{out}/argmax_fp.npy', argfp_y)
    np.save(f'{out}/m0_centroid.npy', m0c_y); np.save(f'{out}/navail.npy', navail_y)
    pd.DataFrame(cfgrows).to_parquet(f'{out}/daycfg.parquet')
    json.dump(dict(run=a.run, bp=a.bp, ndraw=a.ndraw, shift=a.shift, offsets=offs.tolist(), stations=stas,
                   cleaned=not a.noclean, n_clean_iv=len(ivrows), n_cand=len(Cd), seconds=time.time() - t_start),
              open(f'{out}/meta.json', 'w'), indent=1)
    open(f'{out}/DONE', 'w').write('ok\n')
    C.log(f'{a.run}: DONE {len(Cd)} candidates, {time.time() - t_start:.0f} s', lf)


if __name__ == '__main__':
    main()
