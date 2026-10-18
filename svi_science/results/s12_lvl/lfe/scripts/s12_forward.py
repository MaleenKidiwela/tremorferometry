#!/usr/bin/env python
"""s12_forward -- S12-LVL (LFE side), task 1: forward calibration of a Gosselin-sized LVL velocity change.

Target model (Gosselin et al. 2020): Vs inside a ~3.5 km low-velocity layer (LVL) just below the plate interface rises by
+0.05-0.1 km/s after ETS (~+1.5-3 % of an LVL Vs ~3.0-3.4 km/s). Nominal: dVs/Vs = +2 % in a layer of thickness H = 3.5 km
whose top passes through the LFE hypocentre, dipping 10.2 deg toward N52E (plane fit to Slab2 under the SVI families),
LVL Vs = 3.0 km/s (Vp/Vs 2.3 -> Vp 6.9).

Part A (coda): the S11 S-pinned Monte-Carlo kernel (s11_kernel_mc.py, imported, not edited) re-run with an extra tally of
the time each scattered path spends inside the dipping LVL slab (several slab placements tallied at once), and the LVL
slowness put into the travel time. Coefficient a_LVL(window) = sum_b W_b tau_b (R_LVL,b - r_S,LVL) / sum_b W_b tau_b^2 with
the S11 energy weights from each series' grand stack, so  dv/v_apparent = a_LVL * dVs/Vs_LVL  (S-anchored stretch).
Part B (ray theory, flat layers, LVL inserted below the source): lag after direct S (or P) of SxS / PxP reflections and
S-to-P / P-to-S conversions at the LVL base and at an oceanic Moho 7 km below the top; their change for the target.
Part C (source side): P/S amplitude ratio for a source inside the LVL (far-field amplitude ~ (rho c^5)^-1/2 at the source).

    python s12_forward.py mc [--workers 10] [--npart 8000] [--pairs N] [--tag T]
    python s12_forward.py coef [--tag T]
    python s12_forward.py rays
    python s12_forward.py summary
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json, glob
import numpy as np
import pandas as pd
ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, f"{ROOT}/svi_science/scripts")
import s11_common as C
import s11_kernel_mc as KM

OUT = f"{ROOT}/svi_science/results/s12_lvl/lfe"
TB, FIG, LOG = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/logs"
GX, GY = 0.141, 0.110                     # slab-plane depth gradient (km/km, east, north): dip 10.2 deg toward N52E
H_NOM, VS_LVL, VPVS_LVL = 3.5, 3.0, 2.3
DELTA = 0.02                               # nominal dVs/Vs inside the LVL
# slab placements tallied at once: (name, top offset above the source km (+ = top above source), thickness km)
SLABS = [("nom", 0.0, 3.5), ("thin", 0.0, 2.0), ("thick", 0.0, 5.0), ("srcin", 1.0, 3.5), ("below1", -1.0, 3.5),
         ("foot15", 0.0, 3.5), ("foot25", 0.0, 3.5)]
# Gosselin footprint: nominal slab restricted to horizontal distance <= R of CN.PGC / CN.SNB / CN.VGZ (sum of the three
# cylinders; their overlap is small at R 15 km and double-counted at R 25 km -- labelled). Set in cmd_mc from the geometry.
FOOT = {"foot15": 15.0, "foot25": 25.0}
FOOT_XY = []
# MC variants: (ell_top, ell_bot, LVL Vs or None = CVM speed in the LVL)
VAR = {"ell40_cvm": (40.0, 40.0, None), "ell40": (40.0, 40.0, VS_LVL), "two3_20": (3.0, 20.0, VS_LVL),
       "ell20": (20.0, 20.0, VS_LVL)}
WINS = {"1-2": (1.0, 2.0), "2-3": (2.0, 3.0), "3-4": (3.0, 4.0), "4-6": (4.0, 6.0), "2-5": (2.0, 5.0)}
NB, BIN, EDGES = KM.NB, KM.BIN, KM.BIN_EDGES
NCAT = KM.NCAT


def log(m):
    s = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"
    print(s, flush=True)
    open(f"{LOG}/forward.log", "a").write(s + "\n")


# ------------------------------------------------------------------------------ LVL slab tally
def slab_frac(p0, p1, src, top_off, H):
    """fraction (s_lo, s_hi) of straight segments p0->p1 inside the dipping slab whose top is top_off above src"""
    d = p1 - p0
    zt0 = src[2] - top_off + GX * (p0[:, 0] - src[0]) + GY * (p0[:, 1] - src[1])   # slab top depth below p0
    h0 = p0[:, 2] - zt0                                                              # depth below the slab top
    dh = d[:, 2] - (GX * d[:, 0] + GY * d[:, 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        sa = (0.0 - h0) / dh; sb = (H - h0) / dh
    lo = np.clip(np.minimum(sa, sb), 0, 1); hi = np.clip(np.maximum(sa, sb), 0, 1)
    flat = np.abs(dh) < 1e-12
    inside = (h0 >= 0) & (h0 <= H)
    lo = np.where(flat, 0.0, lo); hi = np.where(flat, np.where(inside, 1.0, 0.0), np.maximum(hi, lo))
    return lo, hi


def seg_all(p0, p1, src, rxy, vs_lvl):
    """CVM categories (KM.seg_cats) + LVL length per slab placement + time correction for LVL slowness (nominal slab)"""
    T, Lk, L = KM.seg_cats(p0, p1, src, rxy)
    Ls = np.zeros((len(L), len(SLABS)))
    extra = np.zeros(len(L))
    for j, (nm, off, H) in enumerate(SLABS):
        lo, hi = slab_frac(p0, p1, src, off, H)
        if nm in FOOT:
            acc = np.zeros(len(L))
            for cxy in FOOT_XY:
                cl, ch = KM._quad_iv(p0, p1 - p0, np.r_[cxy, 0.0], FOOT[nm], horiz=True)
                acc += np.maximum(0.0, np.minimum(hi, ch) - np.maximum(lo, cl))
            Ls[:, j] = L * acc
            continue
        Ls[:, j] = L * (hi - lo)
        if nm == "nom" and vs_lvl is not None:
            zm = p0[:, 2] + 0.5 * (lo + hi) * (p1[:, 2] - p0[:, 2])
            extra = Ls[:, j] * (1.0 / vs_lvl - 1.0 / KM.beta_at(zm))
    return T, Lk, L, Ls, extra


def leg_all(p, rec, src, rxy, vs_lvl, image=False):
    n = len(p)
    if not image or rec[2] <= 0:
        return seg_all(p, np.tile(rec, (n, 1)), src, rxy, vs_lvl)
    im = np.array([rec[0], rec[1], -rec[2]])
    d = im[None, :] - p
    s = p[:, 2] / (p[:, 2] - im[2])
    pm = p + s[:, None] * d
    a = seg_all(p, pm, src, rxy, vs_lvl); b = seg_all(pm, np.tile(rec, (n, 1)), src, rxy, vs_lvl)
    return tuple(x + y for x, y in zip(a, b))


def simulate(src, rec, var, npart, seed):
    """copy of KM.simulate with the LVL tallies; returns per lapse bin the mean time in each slab placement"""
    ell_top, ell_bot, vs_lvl = VAR[var]
    ELL = np.where(KM.TOPL, ell_top, ell_bot)
    vs_t = vs_lvl if vs_lvl is not None else float(KM.beta_at(src[2] + 0.5 * H_NOM))   # LVL residence-time speed
    rng = np.random.default_rng(seed)
    rxy = rec[:2]
    Td, Ld, _, Lsd, exd = seg_all(src[None, :], rec[None, :], src, rxy, vs_lvl)
    if rec[2] > 0:
        Tg, Lg, _, Lsg, exg = leg_all(src[None, :], rec, src, rxy, vs_lvl, image=True)
        rS = 0.5 * (Td[0] + Tg[0]); rSl = 0.5 * (Lsd[0] + Lsg[0]); exS = 0.5 * (exd[0] + exg[0])
    else:
        rS = Td[0]; rSl = Lsd[0]; exS = exd[0]
    tS = float(rS.sum() + exS)
    tmax = tS + EDGES[-1]
    n = npart
    pos = np.tile(src, (n, 1)).astype(float)
    u = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); r = np.sqrt(1 - u * u)
    dirv = np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
    cum = np.zeros((n, NCAT)); cumL = np.zeros((n, len(SLABS))); T = np.zeros(n); alive = np.ones(n, bool)
    sumc = np.zeros(NB); sumR = np.zeros((NB, NCAT)); sumL = np.zeros((NB, len(SLABS))); sumtau = np.zeros(NB)
    sumc2 = np.zeros(NB)
    Z_TWO = KM.Z_TWO
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
                d_int = np.where(top, np.where(dz > 0, (Z_TWO - z) / dz, np.inf),
                                 np.where(dz < 0, (z - Z_TWO) / -dz, np.inf))
            d_b = np.minimum(d_surf, d_int)
            if ell_top == ell_bot:
                d_b = d_surf
            opt_b = d_b / ell_c
            stop = tau_rem[j] <= opt_b
            step = np.where(stop, tau_rem[j] * ell_c, d_b)
            p1 = p[j] + step[:, None] * dv[j]
            p1[:, 2] = np.maximum(p1[:, 2], 0.0)
            Tc, Lk, Ls_, Lsl, ex = seg_all(p[j], p1, src, rxy, vs_lvl)
            cum[ia[j]] += Tc; cumL[ia[j]] += Lsl; T[ia[j]] += Tc.sum(1) + ex
            tau_rem[j] = np.where(stop, 0.0, tau_rem[j] - opt_b)
            hit_s = ~stop & (d_surf <= d_b + 1e-12)
            dv[j[hit_s], 2] *= -1.0
            p[j] = p1
            act[j[stop]] = False
        pos[ia] = p; dirv[ia] = dv
        for image in (False, True):
            Tr, Lr, Rr, Lsr, exr = leg_all(pos[ia], rec, src, rxy, vs_lvl, image=image)
            opt = (Lr / ELL[None, :]).sum(1)
            tr = Tr.sum(1) + exr
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
                Ltot = cumL[ia[ok]] + Lsr[ok]
                for k in range(len(SLABS)):
                    sumL[:, k] += np.bincount(bb, c * Ltot[:, k], NB)
        u = rng.uniform(-1, 1, len(ia)); ph = rng.uniform(0, 2 * np.pi, len(ia)); r = np.sqrt(1 - u * u)
        dirv[ia] = np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
        alive[ia[T[ia] > tmax]] = False
    with np.errstate(invalid="ignore", divide="ignore"):
        Rb = sumR / sumc[:, None]; Lb = sumL / sumc[:, None] / vs_t; taub = sumtau / sumc; neff = sumc ** 2 / sumc2
    return dict(Rb=Rb, Rl=Lb, taub=taub, sumc=sumc, neff=neff, rS=rS, rSl=rSl / vs_t, tS=tS)


def _pair(args):
    i, fam, ns, src, rec, npart, variants = args
    out = {}
    for var in variants:
        out[var] = simulate(src, rec, var, npart, C.SEED + 1000 * i + 777)
    return i, fam, ns, out


def cmd_mc(a):
    F, S, fl = KM.geometry()
    try:
        Fx, _, flx = KM.geometry(f"{C.TB}/ext/families.csv")
        F.update(Fx)
    except Exception as e:                                     # pragma: no cover
        log(f"ext families not added: {e}")
    variants = list(VAR)
    FOOT_XY.extend([np.asarray(S[ns][:2], float) for ns in ("CN.PGC", "CN.SNB", "CN.VGZ")])
    jobs = []; i = 0
    for fam in sorted(F):
        for ns in C.STATIONS:
            jobs.append((i, fam, ns, F[fam], S[ns], a.npart, variants)); i += 1
    if a.pairs:
        jobs = [jobs[j] for j in np.random.default_rng(0).choice(len(jobs), a.pairs, replace=False)]
    log(f"MC: {len(jobs)} pairs x {len(variants)} variants, {a.npart} particles, {a.workers} workers; slabs {SLABS}")
    import multiprocessing as mp
    res = {}; t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as p:
        for k, (i, fam, ns, out) in enumerate(p.imap_unordered(_pair, jobs)):
            res[(fam, ns)] = out
            if (k + 1) % 50 == 0 or k + 1 == len(jobs):
                log(f"  {k + 1}/{len(jobs)} pairs {time.time() - t0:.0f}s")
    keys = sorted(res)
    np.savez_compressed(f"{TB}/fwd_mc{a.tag}.npz", fam=np.array([k[0] for k in keys]), netsta=np.array([k[1] for k in keys]),
                        variants=np.array(variants), slabs=np.array([s[0] for s in SLABS]), bin_edges=EDGES,
                        **{f"{f}": np.array([[res[k][v][f] for v in variants] for k in keys])
                           for f in ("Rb", "Rl", "taub", "sumc", "neff", "rS", "rSl", "tS")})
    json.dump({"src": {str(k): v.tolist() for k, v in F.items()}, "sta": {k: v.tolist() for k, v in S.items()}},
              open(f"{TB}/fwd_geometry.json", "w"))
    open(f"{OUT}/FWD_MC{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"MC done {time.time() - t0:.0f}s")


def energy_weights():
    W = {}
    from scipy.signal import hilbert
    for d in (f"{C.TB}/store", f"{C.TB}/store_ext"):
        for p in sorted(glob.glob(f"{d}/*.npz")):
            z = np.load(p); tau = z["tau"]
            for key, G, tS in zip(z["keys"], z["G1"], z["tS"]):
                e = np.abs(hilbert(np.nan_to_num(G.astype(np.float64)), axis=1)) ** 2
                tr = tau - tS
                Eb = np.array([[e[c][(tr >= EDGES[b]) & (tr < EDGES[b + 1])].mean() for b in range(NB)] for c in range(3)])
                Eb[~np.isfinite(G).all(1)] = np.nan
                W[str(key)] = Eb
    return W


def cmd_coef(a):
    z = np.load(f"{TB}/fwd_mc{a.tag}.npz")
    V = list(z["variants"]); SL = list(z["slabs"])
    idx = {(int(f), str(n)): i for i, (f, n) in enumerate(zip(z["fam"], z["netsta"]))}
    W = energy_weights()
    ser = pd.concat([pd.read_csv(f"{C.TB}/series.csv"), pd.read_csv(f"{C.TB}/series_ext.csv")], ignore_index=True)
    ser["key"] = ser.fam.astype(str) + "|" + ser.netsta + "|" + ser.epoch
    good = set(ser[(ser.n_data > 0) & (ser.causal25_3c > 1.5)].key)
    rows = []
    for key, Eb in W.items():
        fam, ns, ep = key.split("|"); fam = int(fam)
        if (fam, ns) not in idx:
            continue
        i = idx[(fam, ns)]
        for vi, v in enumerate(V):
            taub, Rb, Rl, rS, rSl = z["taub"][i, vi], z["Rb"][i, vi], z["Rl"][i, vi], z["rS"][i, vi], z["rSl"][i, vi]
            k = dict(taub=taub, Rb=Rb, rS=rS, sumc=z["sumc"][i, vi])
            for wn, win in WINS.items():
                for comp in ("3c", "Z", "H"):
                    wb = KM.window_weights(Eb, win, comp)
                    if wb.sum() <= 0:
                        continue
                    ok = np.isfinite(taub) & (k["sumc"] > 0) & (wb > 0)
                    den = (wb * taub ** 2)[ok].sum()
                    if den <= 0:
                        continue
                    r = dict(series=key, fam=fam, netsta=ns, epoch=ep, variant=v, win=wn, comp=comp, causal=key in good,
                             tS_mc=float(z["tS"][i, vi]))
                    acvm = ((wb * taub)[ok, None] * (Rb[ok] - rS[None, :])).sum(0) / den
                    r["a_d"] = float(C.group_coef(acvm, "d")); r["a_s"] = float(C.group_coef(acvm, "s"))
                    r["unit_sum_cvm"] = float(acvm.sum())
                    for j, sn in enumerate(SL):
                        r[f"a_lvl_{sn}"] = float(((wb * taub)[ok] * (Rl[ok, j] - rSl[j])).sum() / den)
                    rows.append(r)
    K = pd.DataFrame(rows)
    K.to_parquet(f"{TB}/fwd_coef{a.tag}.parquet", index=False)
    g = K[K.causal & (K.comp == "3c")].groupby(["variant", "win"])
    sm = g.agg(n=("a_lvl_nom", "size"), a_lvl_med=("a_lvl_nom", "median"), a_lvl_q10=("a_lvl_nom", lambda x: x.quantile(.1)),
               a_lvl_q90=("a_lvl_nom", lambda x: x.quantile(.9)), a_thin=("a_lvl_thin", "median"),
               a_thick=("a_lvl_thick", "median"), a_srcin=("a_lvl_srcin", "median"), a_below1=("a_lvl_below1", "median"),
               a_d_med=("a_d", "median"), unit_cvm=("unit_sum_cvm", "median")).reset_index()
    sm["dvv_target_pct_med"] = 100 * DELTA * sm.a_lvl_med
    sm.to_csv(f"{TB}/fwd_coef_summary{a.tag}.csv", index=False)
    log("coef summary (3C, causal series):\n" + sm.round(4).to_string())
    open(f"{OUT}/FWD_COEF{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


# ------------------------------------------------------------------------------ Part B: flat-layer ray theory
def model_layers(zs, H, vs_lvl=VS_LVL, vpvs=VPVS_LVL, moho=7.0):
    """layers as (ztop, zbot, vp, vs): CVM above zs, LVL zs..zs+H, oceanic lower crust to zs+moho (vp 6.9, vs 3.9), mantle"""
    vm = pd.read_csv(C.VEL)
    L = []
    for r in vm.itertuples():
        zt, zb = r.ztop_km, min(r.zbot_km, zs)
        if zt < zs:
            L.append((zt, zb, r.vp, r.vs))
    L.append((zs, zs + H, vs_lvl * vpvs, vs_lvl))
    L.append((zs + H, zs + moho, 6.9, 3.9))
    L.append((zs + moho, zs + moho + 50, 8.1, 4.6))
    return L


def legs_time(p, segs):
    """segs: list of (thickness, v) traversed once each -> (X, T, T_in_each)"""
    X = 0.0; T = 0.0; Ti = []
    for h, v in segs:
        q = np.sqrt(max(1 - (p * v) ** 2, 1e-12))
        X += h * p * v / q; t = h / (v * q); T += t; Ti.append(t)
    return X, T, Ti


def shoot(segs, x):
    vmax = max(v for h, v in segs if h > 0)
    lo, hi = 0.0, (1 - 1e-9) / vmax
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if legs_time(m, segs)[0] < x:
            lo = m
        else:
            hi = m
    return legs_time(0.5 * (lo + hi), segs), 0.5 * (lo + hi)


_PG = np.sin(np.linspace(0, np.pi / 2, 1500))


def shoot_fast(segs, x):
    """vectorised version of shoot(): X(p) on a dense p grid, p interpolated at x, then exact T, Ti at that p"""
    h = np.array([a for a, b in segs]); v = np.array([b for a, b in segs])
    pmax = (1 - 1e-7) / v[h > 0].max()
    p = pmax * _PG
    q = np.sqrt(np.maximum(1 - (p[:, None] * v[None, :]) ** 2, 1e-14))
    X = (h[None, :] * p[:, None] * v[None, :] / q).sum(1)
    p0 = float(np.interp(x, X, p))
    return legs_time(p0, segs), p0


def phase_segs(L, zs, zr, kind):
    """kind: 'S','P' direct up; 'SrS_b','PrP_b','SrP_b','PrS_b' reflect/convert at LVL base; '_m' at Moho"""
    up = lambda w: [(min(zb, zs) - max(zt, zr), (vp if w == "P" else vs)) for zt, zb, vp, vs in L
                    if min(zb, zs) - max(zt, zr) > 0]
    if kind in ("S", "P"):
        return up(kind), []
    base = kind.split("_")[1]
    d_w, u_w = kind[0], kind[2]
    zref = L[[i for i, l in enumerate(L) if abs(l[0] - zs) < 1e-9][0]][1] if base == "b" else \
        L[[i for i, l in enumerate(L) if abs(l[0] - zs) < 1e-9][0] + 1][1]
    below = [(min(zb, zref) - max(zt, zs), vp, vs) for zt, zb, vp, vs in L if min(zb, zref) - max(zt, zs) > 0]
    down = [(h, vp if d_w == "P" else vs) for h, vp, vs in below]
    upb = [(h, vp if u_w == "P" else vs) for h, vp, vs in below]
    lvl_idx = [0, len(down)]                                   # first of down legs and first of up-below legs are LVL
    return down + upb + up(u_w), lvl_idx


def cmd_rays(a):
    geo = json.load(open(f"{TB}/fwd_geometry.json")) if os.path.exists(f"{TB}/fwd_geometry.json") else None
    if geo is None:
        F, S, _ = KM.geometry(); geo = {"src": {str(k): v.tolist() for k, v in F.items()}, "sta": {k: v.tolist() for k, v in S.items()}}
    rows = []
    for fam, src in geo["src"].items():
        src = np.array(src)
        for ns, rec in geo["sta"].items():
            rec = np.array(rec)
            x = float(np.hypot(*(src[:2] - rec[:2]))); zs = float(src[2]); zr = float(rec[2])
            r = dict(fam=int(fam), netsta=ns, epi_km=x, zs=zs)
            for H in (2.0, 3.5, 5.0):
                for vsl in (2.7, 3.0, 3.4):
                    L = model_layers(zs, H, vsl)
                    tt = {}
                    for kind in ("S", "P", "SrS_b", "PrP_b", "SrP_b", "PrS_b", "SrS_m", "PrP_m"):
                        segs, li = phase_segs(L, zs, zr, kind)
                        (X, T, Ti), p = shoot(segs, x)
                        tlvl = sum(Ti[j] for j in li) if li else 0.0
                        tt[kind] = (T, tlvl, p)
                    tag = f"H{H:g}_v{vsl:g}"
                    for kind in ("SrS_b", "SrP_b", "PrS_b", "SrS_m"):
                        r[f"lagS_{kind}_{tag}"] = tt[kind][0] - tt["S"][0]
                        r[f"tlvl_{kind}_{tag}"] = tt[kind][1]
                    for kind in ("PrP_b", "PrP_m"):
                        r[f"lagP_{kind}_{tag}"] = tt[kind][0] - tt["P"][0]
                        r[f"tlvl_{kind}_{tag}"] = tt[kind][1]
                    if tag == "H3.5_v3":
                        r["inc_S_deg"] = float(np.degrees(np.arcsin(min(tt["S"][2] * L[-4][3], 1))))  # at the source
                        r["sp"] = tt["S"][0] - tt["P"][0]
            rows.append(r)
    R = pd.DataFrame(rows); R.to_csv(f"{TB}/fwd_rays.csv", index=False)
    t = "H3.5_v3"
    log(f"rays: SrS(LVL base) lag after S: median {R[f'lagS_SrS_b_{t}'].median():.2f} s (range {R[f'lagS_SrS_b_{t}'].min():.2f}-"
        f"{R[f'lagS_SrS_b_{t}'].max():.2f}); time in LVL {R[f'tlvl_SrS_b_{t}'].median():.2f} s -> change for +2 %: "
        f"{-1000 * DELTA * R[f'tlvl_SrS_b_{t}'].median():.0f} ms")
    open(f"{OUT}/FWD_RAYS_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    os.makedirs(TB, exist_ok=True); os.makedirs(LOG, exist_ok=True)
    ap = argparse.ArgumentParser(); ap.add_argument("cmd")
    ap.add_argument("--workers", type=int, default=10); ap.add_argument("--npart", type=int, default=8000)
    ap.add_argument("--pairs", type=int, default=0); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    {"mc": cmd_mc, "coef": cmd_coef, "rays": cmd_rays}[a.cmd](a)
