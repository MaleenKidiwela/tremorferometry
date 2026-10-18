#!/usr/bin/env python
"""s11_r3_kernel -- S11 round 3 (PLAN §6.1, §6.3). Same physics and random numbers as s11_kernel_mc (imported, not edited):
  * new variants ell10 (10 km everywhere) and two3_20 (3 km for z < 3 km, 20 km below) -> kernel_coef_r3.parquet;
  * ell40 re-run that additionally tallies the 25-45 km time residence in 5-km ALONG-STRIKE bins (exact segment-plane
    intervals, all zones) -> kernel_asbins_r3.parquet (a_d per bin; sum over bins = a_d) + reproduction check vs stored ell40.
Axis = mean Slab2 strike at the region's family epicentres, origin = family centroid (frame of s11_kernel_mc.geometry).
    S11_REGION=... python s11_r3_kernel.py [--workers 12] [--npart 100000] [--ext] [--pairs N]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_kernel_mc as K
import s11_r3_common as R

VAR3 = {"ell40": (40.0, 40.0, 0), "ell10": (10.0, 10.0, 0), "two3_20": (3.0, 20.0, 0)}
AS_EDGES = np.arange(-120.0, 120.0 + 1e-9, 5.0)          # inner edges; outer bins open-ended
NAS = len(AS_EDGES) + 1
DEEP = np.flatnonzero(np.isin(K.REP, K_GROUPS := C.GROUPS["d"]))      # fine layers inside 25-45 km
AX = {"v": None, "c": None}


def strike_axis():
    """mean Slab2 strike at the family epicentres (unit vector in the kernel km frame, northward) + family centroid"""
    from scipy.interpolate import RegularGridInterpolator
    F, S, fl = K.geometry()
    if C.REGION == "SVI":
        F2, _, fl2 = K.geometry(f"{C.TB}/ext/families.csv"); F.update(F2); fl = pd.concat([fl, fl2], ignore_index=True)
    z = pd.read_csv(f"{C.ROOT}/data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz", header=None, names=["lon", "lat", "dep"]).dropna()
    z["lon"] = np.where(z.lon > 180, z.lon - 360, z.lon)
    lo, la = np.sort(z.lon.unique()), np.sort(z.lat.unique())
    g = z.pivot_table(index="lat", columns="lon", values="dep").reindex(index=la, columns=lo).to_numpy()
    f = RegularGridInterpolator((la, lo), g)
    gr = []
    for r in fl.drop_duplicates("fam").itertuples():
        d = 0.05
        gx = (f((r.lat, r.lon + d)) - f((r.lat, r.lon - d))) / (2 * d * 111.19 * np.cos(np.radians(r.lat)))
        gy = (f((r.lat + d, r.lon)) - f((r.lat - d, r.lon))) / (2 * d * 111.19)
        gr.append([float(gx), float(gy)])
    m = np.nanmean(np.array(gr), 0); dip = m / np.linalg.norm(m)
    st = np.array([-dip[1], dip[0]]); st = st if st[1] > 0 else -st
    P = np.array([F[k][:2] for k in sorted(F)]); c = P.mean(0)
    xf = {k: float((F[k][:2] - c) @ st) for k in F}
    xs = {k: float((S[k][:2] - c) @ st) for k in S}
    return st, c, xf, xs


def seg_bins(p0, p1):
    """deep-layer (25-45 km) time residence per along-strike bin for straight segments p0->p1 (n,3) -> (n, NAS) s"""
    d = p1 - p0; L = np.sqrt((d * d).sum(1)); n = len(L)
    z0, dz = p0[:, 2], d[:, 2]
    flat = np.abs(dz) < 1e-12
    a0 = (p0[:, :2] - AX["c"][None, :]) @ AX["v"]; da = d[:, :2] @ AX["v"]
    E = AS_EDGES
    with np.errstate(divide="ignore", invalid="ignore"):
        se = (E[None, :] - a0[:, None]) / da[:, None]                 # (n, NE) s at each inner edge
    pos = da > 1e-12; neg = da < -1e-12; zer = ~(pos | neg)
    lo_b = np.full((n, NAS), -np.inf); hi_b = np.full((n, NAS), np.inf)
    # da > 0: bin j = [se[j-1], se[j]] (se[-1] = -inf, se[NE] = +inf)
    sp = np.where(pos[:, None], se, 0.0)
    lo_b[pos, 1:] = sp[pos]; hi_b[pos, :-1] = sp[pos]
    sn = np.where(neg[:, None], se, 0.0)                               # da < 0: s decreases with a
    hi_b[neg, 1:] = sn[neg]; lo_b[neg, :-1] = sn[neg]
    if zer.any():                                                      # no along-strike motion: whole segment in a0's bin
        j = np.searchsorted(E, a0[zer], side="right")
        lo_b[zer] = np.inf; hi_b[zer] = -np.inf
        lo_b[np.flatnonzero(zer), j] = -np.inf; hi_b[np.flatnonzero(zer), j] = np.inf
    out = np.zeros((n, NAS))
    for k in DEEP:
        with np.errstate(divide="ignore", invalid="ignore"):
            sa = (K.FZ0[k] - z0) / dz
            sb = ((K.FZ1[k] if np.isfinite(K.FZ1[k]) else 1e9) - z0) / dz
        lo = np.clip(np.minimum(sa, sb), 0, 1); hi = np.clip(np.maximum(sa, sb), 0, 1); hi = np.maximum(hi, lo)
        infl = (z0 >= K.FZ0[k]) & (z0 < K.FZ1[k])
        lo = np.where(flat, 0.0, lo); hi = np.where(flat, np.where(infl, 1.0, 0.0), hi)
        ov = np.maximum(0.0, np.minimum(hi[:, None], hi_b) - np.maximum(lo[:, None], lo_b))
        out += ov * (L / K.BETA[k])[:, None]
    return out


def leg_bins(p, rec, image=False):
    n = len(p)
    if not image or rec[2] <= 0:
        return seg_bins(p, np.tile(rec, (n, 1)))
    im = np.array([rec[0], rec[1], -rec[2]])
    d = im[None, :] - p
    s = p[:, 2] / (p[:, 2] - im[2])
    pm = p + s[:, None] * d
    return seg_bins(p, pm) + seg_bins(pm, np.tile(rec, (n, 1)))


def simulate(src, rec, rxy, var, npart, seed):
    """s11_kernel_mc.simulate (identical random-number use) + along-strike deep tally"""
    ell_top, ell_bot, _ = VAR3[var]
    ELL = np.where(K.TOPL, ell_top, ell_bot)
    rng = np.random.default_rng(seed)
    NB, BIN, NCAT, Z_TWO = K.NB, K.BIN, K.NCAT, K.Z_TWO
    Td, Ld, _ = K.seg_cats(src[None, :], rec[None, :], src, rxy)
    bd = seg_bins(src[None, :], rec[None, :])
    if rec[2] > 0:
        Tg, Lg, _ = K.leg_to(src[None, :], rec, src, rxy, image=True)
        rS = 0.5 * (Td[0] + Tg[0]); rSb = 0.5 * (bd[0] + leg_bins(src[None, :], rec, image=True)[0])
    else:
        rS = Td[0]; rSb = bd[0]
    tS = float(rS.sum())
    tmax = tS + K.BIN_EDGES[-1]
    n = npart
    pos = np.tile(src, (n, 1)).astype(float)
    u = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); r = np.sqrt(1 - u * u)
    dirv = np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
    cum = np.zeros((n, NCAT)); cumb = np.zeros((n, NAS)); T = np.zeros(n); alive = np.ones(n, bool)
    sumc = np.zeros(NB); sumc2 = np.zeros(NB); sumR = np.zeros((NB, NCAT)); sumRb = np.zeros((NB, NAS)); sumtau = np.zeros(NB)
    chk = []
    for leg in range(400):
        ia = np.flatnonzero(alive)
        if not len(ia):
            break
        p = pos[ia].copy(); dv = dirv[ia].copy(); tau_rem = rng.exponential(1.0, len(ia))
        act = np.ones(len(ia), bool)
        for it in range(60):
            j = np.flatnonzero(act)
            if not len(j):
                break
            z = p[j, 2]; dz = dv[j, 2]; top = (z < Z_TWO - 1e-9) | ((z <= Z_TWO + 1e-9) & (dz < 0))
            ell_c = np.where(top, ell_top, ell_bot)
            with np.errstate(divide="ignore", invalid="ignore"):
                d_surf = np.where(dz < 0, z / -dz, np.inf)
                d_int = np.where(top, np.where(dz > 0, (Z_TWO - z) / dz, np.inf), np.where(dz < 0, (z - Z_TWO) / -dz, np.inf))
            d_b = np.minimum(d_surf, d_int)
            if ell_top == ell_bot:
                d_b = d_surf
            opt_b = d_b / ell_c
            stop = tau_rem[j] <= opt_b
            step = np.where(stop, tau_rem[j] * ell_c, d_b)
            p1 = p[j] + step[:, None] * dv[j]
            p1[:, 2] = np.maximum(p1[:, 2], 0.0)
            Tc, Lk, Ls = K.seg_cats(p[j], p1, src, rxy)
            Tb = seg_bins(p[j], p1)
            cum[ia[j]] += Tc; cumb[ia[j]] += Tb; T[ia[j]] += Tc.sum(1)
            tau_rem[j] = np.where(stop, 0.0, tau_rem[j] - opt_b)
            hit_s = ~stop & (d_surf <= d_b + 1e-12)
            dv[j[hit_s], 2] *= -1.0
            p[j] = p1
            act[j[stop]] = False
        pos[ia] = p; dirv[ia] = dv
        for image in (False, True):
            Tr, Lr, Rr = K.leg_to(pos[ia], rec, src, rxy, image=image)
            opt = (Lr / ELL[None, :]).sum(1)
            tr = Tr.sum(1)
            lap = T[ia] + tr - tS
            b = np.floor(lap / BIN).astype(int)
            ok = (lap >= 0) & (b < NB)
            if ok.any():
                c = np.exp(-opt[ok]) / (4 * np.pi * np.maximum(Rr[ok], 0.5) ** 2)
                bb = b[ok]
                sumc += np.bincount(bb, c, NB); sumc2 += np.bincount(bb, c * c, NB)
                sumtau += np.bincount(bb, c * lap[ok], NB)
                Rtot = cum[ia[ok]] + Tr[ok]
                for k in range(NCAT):
                    sumR[:, k] += np.bincount(bb, c * Rtot[:, k], NB)
                Rb = cumb[ia[ok]] + leg_bins(pos[ia[ok]], rec, image=image)
                for k in range(NAS):
                    sumRb[:, k] += np.bincount(bb, c * Rb[:, k], NB)
                dsum = Rtot[:, [i for i in range(NCAT) if i // 3 in C.GROUPS["d"]]].sum(1)
                chk.append(float(np.max(np.abs(Rb.sum(1) - dsum))))
        u = rng.uniform(-1, 1, len(ia)); ph = rng.uniform(0, 2 * np.pi, len(ia)); r = np.sqrt(1 - u * u)
        dirv[ia] = np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
        alive[ia[T[ia] > tmax]] = False
    with np.errstate(invalid="ignore", divide="ignore"):
        Rb_ = sumR / sumc[:, None]; Rbb = sumRb / sumc[:, None]; taub = sumtau / sumc; neff = sumc ** 2 / sumc2
    return dict(Rb=Rb_, Rbb=Rbb, taub=taub, sumc=sumc, neff=neff, rS=rS, rSb=rSb, tS=tS, legs=leg + 1,
                binchk=float(max(chk) if chk else 0.0))


def _pair(args):
    i, fam, ns, src, rec, npart, variants = args
    out = {}
    for var in variants:
        seed = C.SEED + 1000 * i + (500000 if VAR3[var][2] else 0)       # = s11_kernel_mc seeds (common random numbers)
        out[var] = simulate(src, rec, rec[:2], var, npart, seed)
    return i, fam, ns, out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--npart", type=int, default=100000); ap.add_argument("--pairs", type=int, default=0)
    ap.add_argument("--ext", action="store_true"); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    lg = R.log(f"kernel{'_ext' if a.ext else ''}{a.tag}"); t0 = time.time()
    st, c0, xf, xs = strike_axis(); AX["v"], AX["c"] = st, c0
    R.dump(dict(strike=st.tolist(), az_deg=float(np.degrees(np.arctan2(st[0], st[1]))), centroid_km=c0.tolist(),
                edges=AS_EDGES.tolist(), fam_x=xf, sta_x=xs, written=R.pt_now()), f"{R.TB3}/along_strike_axis.json")
    lg(f"axis: strike {st.round(3).tolist()} az {np.degrees(np.arctan2(st[0], st[1])):.1f} deg, centroid {c0.round(2).tolist()}; "
       f"family x {min(xf.values()):.1f}..{max(xf.values()):.1f} km; deep fine layers {DEEP.tolist()}")
    F, S, fl = K.geometry(f"{C.TB}/ext/families.csv" if a.ext else None)
    jobs = []; i = 600 if a.ext else 0                       # identical job order and seeds to s11_kernel_mc
    if C.REGION != "SVI":
        ser = pd.read_csv(f"{C.TB}/series.csv"); ser = ser[(ser.n_data > 0) & (ser.causal25 > 1.5)]
        for fam, ns in sorted(set(zip(ser.fam.astype(int), ser.netsta))):
            jobs.append((i, fam, ns, F[fam], S[ns], a.npart, list(VAR3))); i += 1
    for fam in (sorted(F) if C.REGION == "SVI" else []):
        for ns in C.STATIONS:
            jobs.append((i, fam, ns, F[fam], S[ns], a.npart, list(VAR3))); i += 1
    if a.pairs:
        jobs = [jobs[j] for j in np.random.default_rng(0).choice(len(jobs), a.pairs, replace=False)]
    lg(f"{len(jobs)} pairs x {list(VAR3)}, {a.npart} particles, {a.workers} workers")
    import multiprocessing as mp
    res = {}
    with mp.get_context("fork").Pool(a.workers) as p:
        for k, (i, fam, ns, out) in enumerate(p.imap_unordered(_pair, jobs)):
            res[(fam, ns)] = out
            if (k + 1) % 50 == 0 or k + 1 == len(jobs):
                lg(f"  {k + 1}/{len(jobs)} pairs {time.time() - t0:.0f}s")
    W = K.lapse_weights(f"{C.TB}/store_ext" if a.ext else (f"{C.TB}/store" if C.REGION != "SVI" else None))
    rows, brows = [], []
    for key, Eb in W.items():
        fam, ns, ep = key.split("|"); fam = int(fam)
        if (fam, ns) not in res:
            continue
        for v in VAR3:
            k = res[(fam, ns)][v]
            for wn, win in C.WINS.items():
                wb = K.window_weights(Eb, win, "3c")
                if wb.sum() <= 0:
                    continue
                av = K.coef(k, wb)
                r = dict(series=key, fam=fam, netsta=ns, epoch=ep, site=C.SITE[ns], variant=v, win=wn, comp="3c",
                         unit_sum=float(np.nansum(av)), a_s=float(C.group_coef(av, "s")), a_m=float(C.group_coef(av, "m")),
                         a_d=float(C.group_coef(av, "d")), a_b=float(C.group_coef(av, "b")),
                         neff_min_win=float(np.nanmin(np.where(wb > 0, k["neff"], np.inf))), binchk=k["binchk"])
                for ci, cn in enumerate(C.CATS):
                    r[f"c_{cn}"] = float(av[ci])
                rows.append(r)
                if v == "ell40":
                    kb = dict(Rb=k["Rbb"], rS=k["rSb"], taub=k["taub"], sumc=k["sumc"])
                    ab = K.coef(kb, wb)
                    b = dict(series=key, fam=fam, netsta=ns, epoch=ep, win=wn, a_d=r["a_d"], a_d_bins_sum=float(np.nansum(ab)))
                    for j in range(NAS):
                        b[f"b{j}"] = float(ab[j])
                    brows.append(b)
    KC = pd.DataFrame(rows); KB = pd.DataFrame(brows)
    sfx = ("_ext" if a.ext else "") + a.tag
    KC.to_parquet(f"{R.TB3}/kernel_coef_r3{sfx}.parquet", index=False); KB.to_parquet(f"{R.TB3}/kernel_asbins_r3{sfx}.parquet", index=False)
    # ---- checks (PLAN §6.1/§6.3 stop rules)
    old = pd.read_parquet(f"{C.TB}/kernel_coef{'_ext' if a.ext else ''}.parquet")
    old = old[(old.variant == "ell40") & (old.comp == "3c")][["series", "win", "a_s", "a_d"]]
    m = KC[KC.variant == "ell40"].merge(old, on=["series", "win"], suffixes=("", "_old"))
    rep = float(np.nanmax(np.abs(np.r_[m.a_s - m.a_s_old, m.a_d - m.a_d_old]))) if len(m) else np.nan
    ut = float((KC.unit_sum - 1).abs().max())
    prim = KC[(KC.win == C.HEAD)].groupby(["variant", "fam", "netsta"]).a_d.first().reset_index()
    neg = prim.groupby("variant").a_d.apply(lambda x: float((x < 0).mean())).to_dict()
    bsum = float(np.nanmax(np.abs(KB.a_d_bins_sum - KB.a_d))) if len(KB) else np.nan
    summ = KC[KC.win == C.HEAD].groupby("variant")[["a_s", "a_m", "a_d", "a_b"]].median().round(4).to_dict("index")
    lg(f"unit test max |sum-1| {ut:.2e}; a_d<0 fraction {neg}; ell40 reproduction max |diff| {rep:.2e} over {len(m)} rows; "
       f"along-strike bins: max |sum - a_d| {bsum:.2e}; bin path check {KC.binchk.max():.1e} s")
    lg(f"medians 2-5 s: {summ}")
    stop = (ut > 0.02) or any(v > 0.10 for v in neg.values()) or not (rep <= 1e-6) or not (bsum <= 1e-6)
    mk = f"{R.R3}/KERNEL{sfx}_{'STOP' if stop else 'DONE'}"
    open(mk, "w").write(f"{R.pt_now()}\nunit {ut:.3e} neg {neg} repro {rep:.3e} bins {bsum:.3e}\n")
    lg(f"wrote {mk} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
