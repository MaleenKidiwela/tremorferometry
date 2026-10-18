#!/usr/bin/env python
"""s11_kernel_mc -- S-pinned differential TIME-residence kernels by Monte-Carlo radiative transfer for the S11 pilot
(rules: svi_science/results/s11_depth/PLAN.md D1; adapted from fault_tomography/kernels/residence_kernel_mc.py, which is
not edited).

Physics. The stretch is measured about the picked S after per-event realignment on S, so the observable at lapse tau past
S is the coda's EXCESS delay over the direct ray: delta_tau = -(R - r_S) . delta, delta = dbeta/beta per category, where
R = time-residence vector (s) of the coda path and r_S that of the direct S. A stretch fitted over a window weights lapse
tau by energy x tau^2, so  eps_apparent = -dv/v_apparent = -a . delta  with
        a = sum_b W_b tau_b (R_b - r_S) / sum_b W_b tau_b^2        (sum_k a_k = 1: homogeneous change -> dv/v = delta).
MC: isotropic S-energy random walk from the hypocentre, isotropic scattering, mean free path ell(z), perfectly reflecting
free surface, straight legs; next-event estimation at every scattering to the receiver at its sensor depth AND to its
free-surface image (equal weights); weight = exp(-optical depth)/(4 pi R^2). Every segment is cut EXACTLY at the union of
the CVM v1.7 (south_VI) layer edges and the reporting layer edges; time residence = length/beta(layer). Zones are exact
line-sphere / line-cylinder intervals: near-source sphere r <= 10 km (priority), near-station vertical cylinder r <= 10 km,
rest. 13 lapse bins of 0.5 s past S (0-6.5 s). Direct S = straight ray (for a borehole: mean of direct and surface ghost).

Variants: ell 20 / 40 (primary) / 80 km; two-layer (5 km for z < 3 km, 40 km below); second seed at ell 40.
Outputs (tables/): kernel_mc.npz (per pair x variant: R_b, tau_b, sum_c, n_eff, r_S, t_S), kernel_coef.parquet (per
series x variant x window x component: 30 coefficients, a_s, a_m, a_d, a_b, unit sum, shares), kernel_unit_test.csv,
kernel_summary.csv. Marker KERNEL_DONE, or KERNEL_STOP if the unit test / a_d sign rule fails.

    python s11_kernel_mc.py [--workers 8] [--npart 8000] [--pairs N]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C

R_SRC, R_STA = 10.0, 10.0
BIN = 0.5
NB = 13
BIN_EDGES = np.arange(NB + 1) * BIN
VARIANTS = {"ell20": (20.0, 20.0, 0), "ell40": (40.0, 40.0, 0), "ell80": (80.0, 80.0, 0),
            "two5_40": (5.0, 40.0, 0), "ell40_s2": (40.0, 40.0, 1)}
Z_TWO = 3.0                                     # two-layer ell boundary (km)
NCAT = len(C.CATS)

# ------------------------------------------------------------------------------ velocity model and fine layers
vm = pd.read_csv(C.VEL)
cvm_edges = np.r_[vm.ztop_km.to_numpy(), np.inf]
edges = np.unique(np.r_[C.LAYERS, cvm_edges, Z_TWO])
FZ0, FZ1 = edges[:-1], edges[1:]
K = len(FZ0)
mid = np.where(np.isfinite(FZ1), 0.5 * (FZ0 + FZ1), FZ0 + 1.0)
BETA = np.array([vm.vs.to_numpy()[np.searchsorted(vm.ztop_km.to_numpy(), z, side="right") - 1] for z in mid])
REP = np.searchsorted(C.LAYERS, mid, side="right") - 1          # reporting layer of each fine layer
TOPL = FZ1 <= Z_TWO + 1e-9                                       # fine layers inside the two-layer top


def beta_at(z):
    return BETA[np.clip(np.searchsorted(FZ0, z, side="right") - 1, 0, K - 1)]


# ------------------------------------------------------------------------------ exact segment geometry
def _clip_iv(a, b):
    lo = np.clip(np.minimum(a, b), 0.0, 1.0); hi = np.clip(np.maximum(a, b), 0.0, 1.0)
    return lo, np.maximum(hi, lo)


def _quad_iv(p0, d, c, R, horiz=False):
    """s-interval in [0,1] where |p0 + s d - c| <= R (horizontal components only if horiz)"""
    if horiz:
        p0 = p0[:, :2]; d = d[:, :2]; c = c[:2]
    q = p0 - c[None, :]
    A = (d * d).sum(1); B = 2 * (d * q).sum(1); Cc = (q * q).sum(1) - R * R
    disc = B * B - 4 * A * Cc
    ok = (disc >= 0) & (A > 1e-12)
    sq = np.sqrt(np.where(ok, disc, 0.0))
    s1 = np.where(ok, (-B - sq) / np.where(A > 1e-12, 2 * A, 1.0), 1.0)
    s2 = np.where(ok, (-B + sq) / np.where(A > 1e-12, 2 * A, 1.0), 0.0)
    # degenerate (no horizontal motion / zero length): inside -> whole segment
    deg = A <= 1e-12
    inside = Cc <= 0
    s1 = np.where(deg, np.where(inside, 0.0, 1.0), s1); s2 = np.where(deg, np.where(inside, 1.0, 0.0), s2)
    lo = np.clip(s1, 0, 1); hi = np.clip(s2, 0, 1)
    return lo, np.maximum(hi, lo)


def seg_cats(p0, p1, src, rxy):
    """physical straight segments p0->p1 (n,3), z >= 0 on both ends.
    -> time per category (n, 30) [s], length per fine layer (n, K) [km], total length (n,)"""
    d = p1 - p0
    L = np.sqrt((d * d).sum(1))
    z0, dz = p0[:, 2], d[:, 2]
    flat = np.abs(dz) < 1e-12
    with np.errstate(divide="ignore", invalid="ignore"):
        sa = (FZ0[None, :] - z0[:, None]) / dz[:, None]
        sb = (np.where(np.isfinite(FZ1), FZ1, 1e9)[None, :] - z0[:, None]) / dz[:, None]
    lo, hi = _clip_iv(sa, sb)
    infl = (z0[:, None] >= FZ0[None, :]) & (z0[:, None] < FZ1[None, :])
    lo = np.where(flat[:, None], 0.0, lo); hi = np.where(flat[:, None], np.where(infl, 1.0, 0.0), hi)
    sl, sh = _quad_iv(p0, d, src, R_SRC)
    cl, ch = _quad_iv(p0, d, np.r_[rxy, 0.0], R_STA, horiz=True)
    ov = lambda a0, a1, b0, b1: np.maximum(0.0, np.minimum(a1, b1) - np.maximum(a0, b0))
    f_lay = hi - lo                                                   # (n,K) fraction of segment in fine layer
    f_src = ov(lo, hi, sl[:, None], sh[:, None])
    bl = np.maximum(sl, cl); bh = np.minimum(sh, ch)
    f_both = ov(lo, hi, bl[:, None], np.maximum(bh, bl)[:, None])
    f_sta = ov(lo, hi, cl[:, None], ch[:, None]) - f_both
    f_rest = f_lay - f_src - f_sta
    Lk = L[:, None] * f_lay
    tk = L[:, None] / BETA[None, :]
    T = np.zeros((len(L), NCAT))
    for zi, f in enumerate((f_src, f_sta, f_rest)):
        tz = tk * f                                                   # (n, K)
        for k in range(K):
            T[:, REP[k] * 3 + zi] += tz[:, k]
    return T, Lk, L


def leg_to(p, rec, src, rxy, image=False):
    """receiver leg from p (n,3) to rec (3,) or to its free-surface image (path reflected at z = 0)."""
    n = len(p)
    if not image or rec[2] <= 0:
        T, Lk, L = seg_cats(p, np.tile(rec, (n, 1)), src, rxy)
        return T, Lk, L
    im = np.array([rec[0], rec[1], -rec[2]])
    d = im[None, :] - p
    s = p[:, 2] / (p[:, 2] - im[2])                                   # crossing of z = 0
    pm = p + s[:, None] * d
    T1, L1, l1 = seg_cats(p, pm, src, rxy)
    T2, L2, l2 = seg_cats(pm, np.tile(rec, (n, 1)), src, rxy)
    return T1 + T2, L1 + L2, l1 + l2


# ------------------------------------------------------------------------------ MC for one pair and variant
def simulate(src, rec, rxy, var, npart, seed):
    ell_top, ell_bot, _ = VARIANTS[var]
    ELL = np.where(TOPL, ell_top, ell_bot)                           # per fine layer
    rng = np.random.default_rng(seed)
    # direct S (+ ghost for a borehole)
    Td, Ld, _ = seg_cats(src[None, :], rec[None, :], src, rxy)
    if rec[2] > 0:
        Tg, Lg, _ = leg_to(src[None, :], rec, src, rxy, image=True)
        rS = 0.5 * (Td[0] + Tg[0])
    else:
        rS = Td[0]
    tS = float(rS.sum())
    tmax = tS + BIN_EDGES[-1]
    n = npart
    pos = np.tile(src, (n, 1)).astype(float)
    u = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); r = np.sqrt(1 - u * u)
    dirv = np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
    cum = np.zeros((n, NCAT)); T = np.zeros(n); alive = np.ones(n, bool)
    sumc = np.zeros(NB); sumc2 = np.zeros(NB); sumR = np.zeros((NB, NCAT)); sumtau = np.zeros(NB)
    lenchk = []
    for leg in range(400):
        ia = np.flatnonzero(alive)
        if not len(ia):
            break
        # ---- free path with optical-depth stepping (boundaries: z = Z_TWO and the free surface)
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
                d_b = d_surf                                          # single ell: only the surface matters
            opt_b = d_b / ell_c
            stop = tau_rem[j] <= opt_b
            step = np.where(stop, tau_rem[j] * ell_c, d_b)
            p1 = p[j] + step[:, None] * dv[j]
            p1[:, 2] = np.maximum(p1[:, 2], 0.0)
            Tc, Lk, Ls = seg_cats(p[j], p1, src, rxy)
            if leg == 0 and it == 0:
                lenchk.append(np.max(np.abs(Lk.sum(1) - Ls) / np.maximum(Ls, 1e-9)))
            cum[ia[j]] += Tc; T[ia[j]] += Tc.sum(1)
            tau_rem[j] = np.where(stop, 0.0, tau_rem[j] - opt_b)
            hit_s = ~stop & (d_surf <= d_b + 1e-12)
            dv[j[hit_s], 2] *= -1.0
            p[j] = p1
            act[j[stop]] = False
        pos[ia] = p; dirv[ia] = dv
        # ---- next-event estimation to the receiver and its image
        for image in (False, True):
            Tr, Lr, Rr = leg_to(pos[ia], rec, src, rxy, image=image)
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
        # ---- scatter (isotropic) and kill
        u = rng.uniform(-1, 1, len(ia)); ph = rng.uniform(0, 2 * np.pi, len(ia)); r = np.sqrt(1 - u * u)
        dirv[ia] = np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
        alive[ia[T[ia] > tmax]] = False
    with np.errstate(invalid="ignore", divide="ignore"):
        Rb = sumR / sumc[:, None]; taub = sumtau / sumc; neff = sumc ** 2 / sumc2
    return dict(Rb=Rb, taub=taub, sumc=sumc, neff=neff, rS=rS, tS=tS, legs=leg + 1,
                lenchk=float(max(lenchk) if lenchk else 0.0))


# ------------------------------------------------------------------------------ inputs
def geometry_region():
    """PLAN §5 (WA): families (Slab2-replaced depths) and stations (FDSN sensor depths) from the region inventory"""
    fl = pd.read_csv(f"{C.MTB}/families.csv")
    st = pd.read_csv(f"{C.MTB}/stations.csv")
    lat0, lon0 = float(fl.lat.mean()), float(fl.lon.mean())
    pj = lambda la, lo: np.array([(lo - lon0) * 111.19 * np.cos(np.radians(lat0)), (la - lat0) * 111.19])
    S = {r.netsta: np.r_[pj(r.lat, r.lon), r.sensor_depth_m / 1000.0] for r in st.itertuples()}
    F = {int(r.fam): np.r_[pj(r.lat, r.lon), r.depth_km] for r in fl.itertuples()}
    return F, S, fl


def geometry(fams_csv=None):
    if C.REGION != "SVI":
        return geometry_region()
    fam = pd.read_csv(f"{C.ROOT}/data/catalogue/families_v1_3.csv")
    fam = fam[fam.region == "S_Vancouver_Is"].copy(); fam["fam"] = pd.to_numeric(fam.family_id, errors="coerce")
    fl = pd.read_csv(fams_csv or f"{C.MTB}/families.csv").merge(fam[["fam", "lat", "lon", "depth_km"]], on="fam", how="left")
    inv = pd.read_csv(f"{C.ROOT}/svi_science/results/backward_scope/inventory_stations.csv")
    inv["netsta"] = inv.net + "." + inv.sta; inv = inv.drop_duplicates("netsta").set_index("netsta")
    dep = pd.read_csv(f"{C.TB}/sensor_depths.csv"); dep["netsta"] = dep.Network + "." + dep.Station
    dep = dep.set_index("netsta").sensor_depth_m / 1000.0
    lat0, lon0 = 48.5, -123.6
    pj = lambda la, lo: np.array([(lo - lon0) * 111.19 * np.cos(np.radians(lat0)), (la - lat0) * 111.19])
    S = {ns: np.r_[pj(inv.loc[ns, "lat"], inv.loc[ns, "lon"]), dep[ns]] for ns in C.STATIONS}
    F = {int(r.fam): np.r_[pj(r.lat, r.lon), r.depth_km] for r in fl.itertuples()}
    return F, S, fl


def _pair(args):
    i, fam, ns, src, rec, npart = args
    out = {}
    for vi, var in enumerate(VARIANTS):
        seed = C.SEED + 1000 * i + (500000 if VARIANTS[var][2] else 0)   # common random numbers across ell variants
        out[var] = simulate(src, rec, rec[:2], var, npart, seed)
    return i, fam, ns, out


# ------------------------------------------------------------------------------ lapse weights from grand stacks
def lapse_weights(store=None):
    from scipy.signal import hilbert
    if store:                                      # Revision 2: grand stacks from the S11 store of the added families
        import glob
        ks, Gs, ts = [], [], []
        for p in sorted(glob.glob(f"{store}/*.npz")):
            zz = np.load(p); ks += list(zz["keys"]); Gs += list(zz["G1"]); ts += list(zz["tS"])
        z = dict(keys=np.array(ks), G1=np.array(Gs), tS=np.array(ts), tau=S_TAU())
    else:
        z = np.load(f"{C.MTB}/grandstacks.npz")
    tau = z["tau"]
    W = {}
    for key, G, tS in zip(z["keys"], z["G1"], z["tS"]):
        e = np.abs(hilbert(np.nan_to_num(G), axis=1)) ** 2                    # (3, NW)
        tr = tau - tS
        Eb = np.array([[e[c][(tr >= BIN_EDGES[b]) & (tr < BIN_EDGES[b + 1])].mean() for b in range(NB)] for c in range(3)])
        Eb[~np.isfinite(G).all(1)] = np.nan
        W[str(key)] = Eb
    return W


def S_TAU():
    import s10_common as S10
    return S10.TAU


def coef(k, wb):
    """k = simulate() output, wb (NB,) weights -> a (30,)"""
    ok = np.isfinite(k["taub"]) & (k["sumc"] > 0) & (wb > 0)
    num = ((wb * k["taub"])[ok, None] * (k["Rb"][ok] - k["rS"][None, :])).sum(0)
    den = (wb * k["taub"] ** 2)[ok].sum()
    return num / den if den > 0 else np.full(NCAT, np.nan)


def window_weights(Eb, win, comp):
    """Eb (3, NB) energy per bin per component -> NB weights for a stretch over win with components comp
    (each component z-scored on the window, as the joint 3C stretch does)"""
    inw = np.array([(BIN_EDGES[b] >= win[0] - 1e-9) and (BIN_EDGES[b + 1] <= win[1] + 1e-9) for b in range(NB)])
    cs = {"3c": (0, 1, 2), "Z": (0,), "H": (1, 2)}[comp]
    w = np.zeros(NB)
    for c in cs:
        if not np.isfinite(Eb[c]).all():
            continue
        e = np.where(inw, Eb[c], 0.0); w += e / max(e.sum(), 1e-30)
    return w


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--npart", type=int, default=8000); ap.add_argument("--pairs", type=int, default=0)
    ap.add_argument("--tag", default="")
    ap.add_argument("--ext", action="store_true", help="Revision 2: the added families (tables/ext/families.csv, store_ext)")
    a = ap.parse_args()
    C.ensure(); log = C.Log(f"{C.LOGD}/kernel{a.tag}.log")
    F, S, fl = geometry(f"{C.TB}/ext/families.csv" if a.ext else None)
    log(f"fine layers {K}: edges {np.round(edges[:-1], 2).tolist()} beta {BETA.round(3).tolist()}")
    out_of_25_45 = fl[(fl.depth_km < 25) | (fl.depth_km > 45)]
    log(f"families {len(F)}, depth {fl.depth_km.min():.1f}-{fl.depth_km.max():.1f} km; outside 25-45 km: "
        f"{out_of_25_45.fam.tolist()} (kept at their hypocentre depth; flagged)")
    jobs = []; i = 600 if a.ext else 0          # independent MC seeds for the added families
    if C.REGION != "SVI":                        # PLAN §5: only (family, station) pairs with a causal measured series
        ser = pd.read_csv(f"{C.TB}/series.csv")
        ser = ser[(ser.n_data > 0) & (ser.causal25 > 1.5)]
        pairs = sorted(set(zip(ser.fam.astype(int), ser.netsta)))
        log(f"region {C.REGION}: {len(pairs)} (family, station) pairs with a causal series")
        for fam, ns in pairs:
            jobs.append((i, fam, ns, F[fam], S[ns], a.npart)); i += 1
    for fam in (sorted(F) if C.REGION == "SVI" else []):
        for ns in C.STATIONS:
            jobs.append((i, fam, ns, F[fam], S[ns], a.npart)); i += 1
    if a.pairs:
        jobs = [jobs[j] for j in np.random.default_rng(0).choice(len(jobs), a.pairs, replace=False)]
    log(f"{len(jobs)} pairs x {len(VARIANTS)} variants, {a.npart} particles, {a.workers} workers")
    import multiprocessing as mp
    res = {}; t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as p:
        for k, (i, fam, ns, out) in enumerate(p.imap_unordered(_pair, jobs)):
            res[(fam, ns)] = out
            if (k + 1) % 50 == 0 or k + 1 == len(jobs):
                log(f"  {k + 1}/{len(jobs)} pairs {time.time() - t0:.0f}s")
    keys = sorted(res); V = list(VARIANTS)
    np.savez_compressed(f"{C.TB}/kernel_mc{a.tag}.npz", fam=np.array([k[0] for k in keys]), netsta=np.array([k[1] for k in keys]),
                        variants=np.array(V), cats=np.array(C.CATS), bin_edges=BIN_EDGES,
                        Rb=np.array([[res[k][v]["Rb"] for v in V] for k in keys]),
                        taub=np.array([[res[k][v]["taub"] for v in V] for k in keys]),
                        sumc=np.array([[res[k][v]["sumc"] for v in V] for k in keys]),
                        neff=np.array([[res[k][v]["neff"] for v in V] for k in keys]),
                        rS=np.array([[res[k][v]["rS"] for v in V] for k in keys]),
                        tS=np.array([[res[k][v]["tS"] for v in V] for k in keys]),
                        legs=np.array([[res[k][v]["legs"] for v in V] for k in keys]),
                        lenchk=np.array([[res[k][v]["lenchk"] for v in V] for k in keys]))
    # ---- per-series coefficients
    W = lapse_weights(f"{C.TB}/store_ext" if a.ext else (f"{C.TB}/store" if C.REGION != "SVI" else None))
    rows = []
    for key, Eb in W.items():
        fam, ns, ep = key.split("|"); fam = int(fam)
        if (fam, ns) not in res:
            continue
        for v in V:
            k = res[(fam, ns)][v]
            for wn, win in C.WINS.items():
                for comp in ("3c", "Z", "H"):
                    if comp != "3c" and wn != C.HEAD:
                        continue
                    wb = window_weights(Eb, win, comp)
                    if wb.sum() <= 0:
                        continue
                    av = coef(k, wb)
                    r = dict(series=key, fam=fam, netsta=ns, epoch=ep, site=C.SITE[ns], variant=v, win=wn, comp=comp,
                             tS_mc=k["tS"], unit_sum=float(np.nansum(av)),
                             a_s=float(C.group_coef(av, "s")), a_m=float(C.group_coef(av, "m")),
                             a_d=float(C.group_coef(av, "d")), a_b=float(C.group_coef(av, "b")),
                             neff_min_win=float(np.nanmin(np.where(wb > 0, k["neff"], np.inf))),
                             hypo_km=float(np.linalg.norm(S[ns] - F[fam])), epi_km=float(np.linalg.norm(S[ns][:2] - F[fam][:2])),
                             fam_depth=float(F[fam][2]), sensor_depth_km=float(S[ns][2]))
                    a3 = av.reshape(len(C.LAYER_NAMES), 3)
                    r["a_s_sta_share"] = float(a3[C.GROUPS["s"], 1].sum() / r["a_s"]) if r["a_s"] != 0 else np.nan
                    r["a_d_src_share"] = float(a3[C.GROUPS["d"], 0].sum() / r["a_d"]) if r["a_d"] != 0 else np.nan
                    for ci, cn in enumerate(C.CATS):
                        r[f"c_{cn}"] = float(av[ci])
                    rows.append(r)
    KC = pd.DataFrame(rows)
    KC.to_parquet(f"{C.TB}/kernel_coef{a.tag}.parquet", index=False)
    # ---- unit test + summary
    UT = KC.assign(dev=(KC.unit_sum - 1).abs())
    ut = UT.groupby(["variant", "win", "comp"]).agg(n=("dev", "size"), max_dev=("dev", "max"),
                                                    fail=("dev", lambda x: int((x > 0.02).sum()))).reset_index()
    lc = max(res[k][v]["lenchk"] for k in res for v in V)
    ut["max_len_conservation_err"] = lc
    ut.to_csv(f"{C.TB}/kernel_unit_test{a.tag}.csv", index=False)
    prim = KC[(KC.variant == "ell40") & (KC.win == C.HEAD) & (KC.comp == "3c")]
    pp = prim.groupby(["fam", "netsta"]).first().reset_index()
    neg = float((pp.a_d < 0).mean())
    sm = []
    for (v, wn, comp), g in KC.groupby(["variant", "win", "comp"]):
        sm.append(dict(variant=v, win=wn, comp=comp, n=len(g), a_s_med=g.a_s.median(), a_m_med=g.a_m.median(),
                       a_d_med=g.a_d.median(), a_d_q10=g.a_d.quantile(.1), a_d_q90=g.a_d.quantile(.9), a_b_med=g.a_b.median(),
                       a_d_neg_frac=float((g.a_d < 0).mean()), a_s_sta_share_med=g.a_s_sta_share.median(),
                       a_d_src_share_med=g.a_d_src_share.median(), unit_maxdev=float((g.unit_sum - 1).abs().max())))
    SM = pd.DataFrame(sm); SM.to_csv(f"{C.TB}/kernel_summary{a.tag}.csv", index=False)
    # ell spread of a_s per pair (2-5 s 3C) and family-depth spread of a_d per site
    h = KC[(KC.win == C.HEAD) & (KC.comp == "3c")].groupby(["fam", "netsta", "variant"]).a_s.first().unstack()
    ellv = [c for c in ("ell20", "ell40", "ell80", "two5_40") if c in h]
    spread = (h[ellv].max(1) / h[ellv].min(1).where(h[ellv].min(1) > 0))
    sd = pp.groupby("site").a_d.agg(lambda x: x.max() / x.min() if x.min() > 0 else np.nan)
    seedd = KC[(KC.win == C.HEAD) & (KC.comp == "3c") & KC.variant.isin(["ell40", "ell40_s2"])].groupby(
        ["fam", "netsta", "variant"]).a_d.first().unstack()
    log(f"unit test: max |sum-1| {UT.dev.max():.2e}, failures {int((UT.dev > 0.02).sum())}; segment length conservation max rel err {lc:.1e}")
    log(f"primary (ell40, 2-5 s, 3C): a_s median {pp.a_s.median():.3f}, a_m {pp.a_m.median():.3f}, a_d {pp.a_d.median():.3f} "
        f"(10-90% {pp.a_d.quantile(.1):.3f}..{pp.a_d.quantile(.9):.3f}), a_b {pp.a_b.median():.3f}; a_d < 0 for {100 * neg:.1f}% of pairs")
    log(f"ell spread of a_s (max/min over {ellv}): median x{spread.median():.2f}, 90% x{spread.quantile(.9):.2f}")
    log(f"family-depth spread of a_d per site (max/min): {sd.round(2).to_dict()}")
    if "ell40_s2" in seedd:
        log(f"second seed: |a_d(s2) - a_d(s1)| median {np.nanmedian(np.abs(seedd.ell40_s2 - seedd.ell40)):.4f}")
    stop = (UT.dev > 0.02).any() or neg > 0.10
    mk = f"{C.OUT}/KERNEL{a.tag}_{'STOP' if stop else 'DONE'}"
    open(mk, "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n") + f"unit_maxdev {UT.dev.max():.3e} a_d_neg {neg:.3f}\n")
    log(f"wrote {mk}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
