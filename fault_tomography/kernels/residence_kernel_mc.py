#!/usr/bin/env python
"""STEP 0 of Merlin's 2026-09-02 3-5 s ruling: S-PINNED DIFFERENTIAL RESIDENCE KERNELS by Monte-Carlo radiative
transfer, for every production (station, cell) pair. Diagnostic + candidate G_diff operators. Nothing in
res_catalog is touched (writes res_catalog_kernels/).

Physics (Merlin F1-F3): the stacks are aligned on the matched-filter detection (direct S pinned at stack +1.0 s) and
the stretch is computed about the S, so the measurement is the coda's EXCESS delay over the direct ray:
    eps_apparent = -[ t * <K_coda, delta> - t_S * <K_ray, delta> ] / (t - t_S),   delta = dbeta/beta
i.e. the effective kernel is K_eff = [t K_coda - t_S K_ray]/(t - t_S), where K are PATH-RESIDENCE fractions (length
of path inside a voxel / total path length; sum = 1), t the (energy-weighted) total lapse from origin of the coda
window and t_S the direct-S lapse. For a homogeneous delta, K_eff sums to 1 (eps = delta). For a change confined to
the +-1.5 km interface slab the net sensitivity is [t f_coda - t_S f_ray]/(t - t_S) with f = in-slab path fraction —
the small difference f_coda - f_ray gets amplified by t/(t - t_S) ~ 8-18 and can be negative.

Monte Carlo: isotropic S-energy random walk from the source (mean free path ell, isotropic scattering, no intrinsic
absorption — it cancels in the kernel shape), perfectly reflecting free surface (z -> |z|), next-event estimation to
the surface receiver at every scattering event. Every contributing path's full residence (all legs + final leg) is
accumulated per lapse bin (0.5 s past S, 0.5..4.5 s). Single-scatter (order-1) contributions kept separately.
Window kernels combine bins with weights = OBSERVED fleet-median coda energy profile (data/window_snr_profile.csv,
rms_rel^2 by 0.5-s bin) [primary] and with the MC's own energy [secondary]; the direct ray is rasterized exactly.
Windows: W24 = stack 2-4 s = 1-3 s past S; W35 = stack 3-5 s = 2-4 s past S (Flaw 3 bookkeeping fix).

Outputs (OUT): residence_diag.csv — per pair x ell: t_S, t_bar, f_ray_iface, f_coda_iface (all-order + single),
  f_*_top5, net_iface (K_eff on the slab), net_top5, net_srccell (slab within 15 km of the source), naive (= f_coda)
  G_diff_w24.npz / G_diff_w35.npz — K_eff projected to the 872 cells (nearest cell within 15 km of the interface
  voxel column), stored as G = -K_eff rows (eps convention as G.npz), ell = 40 km primary, plus per-ell nets.
Caveats: isotropic scattering + isotropic source (no radiation pattern), uniform beta (straight rays), interface =
  Slab2 input constraints (+-1.5 km), 5 x 5 x 1 km voxels, 1-km path sampling.
Usage: [NPART=8000] [WORKERS=8] [PAIRS=all|N] python fault_tomography/kernels/residence_kernel_mc.py
"""
import os, sys, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from concurrent.futures import ProcessPoolExecutor

RD = "fault_tomography/inversion/res_catalog"
OUT = "fault_tomography/inversion/res_catalog_kernels"
BETA = 3.5
ELLS = (20.0, 40.0, 80.0)
NPART = int(os.environ.get("NPART", "8000"))
DS = float(os.environ.get("DS", "1.0"))                       # path sampling (km)
DX, DZ, ZMAX = float(os.environ.get("DX", "5.0")), float(os.environ.get("DZ", "1.0")), 70.0
HSLAB, TOP, RCELL = float(os.environ.get("HSLAB", "1.5")), 5.0, 15.0
SEED_OFFSET = int(os.environ.get("SEED_OFFSET", "0"))         # Merlin seed-reproducibility test
HG_G = float(os.environ.get("HG_G", "0"))                     # Henyey-Greenstein anisotropy; ell = ell* (1-g), ell* = ELLS
OUTTAG = os.environ.get("OUTTAG", "")                         # suffix for variant outputs (never overwrite the primary)
BIN_EDGES = np.arange(0.5, 4.5 + 1e-9, 0.5)          # lapse past S: 8 bins
NB = len(BIN_EDGES) - 1
BIN_C = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])
WIN = {"w24": (1.0, 3.0), "w35": (2.0, 4.0)}
TMAX_PAST_S = BIN_EDGES[-1]

# ---------------- geometry ----------------
gz = np.load(f"{RD}/G.npz", allow_pickle=True)
cells = gz["cell"]; CXY = gz["cxy"]; CZ = gz["depth_km"].astype(float); M = len(cells)
lat0, lon0 = float(gz["lat0"]), float(gz["lon0"])
ptag, pcell = gz["pair_tag"], gz["pair_cell"]
cellpos = {c: i for i, c in enumerate(cells)}
pcell_idx = np.array([cellpos[c] for c in pcell])
pairs = pd.read_csv(f"{RD}/pairs.csv")
sta_ll = pairs.groupby("tag")[["sta_lat", "sta_lon"]].first()


def proj(lat, lon):
    return np.column_stack([(np.asarray(lon) - lon0) * 111 * np.cos(np.radians(lat0)),
                            (np.asarray(lat) - lat0) * 111])


def unproj(x, y):
    return lat0 + y / 111.0, lon0 + x / (111.0 * np.cos(np.radians(lat0)))


SXY = proj(sta_ll.sta_lat.values, sta_ll.sta_lon.values)
sxy_of = {t: SXY[i] for i, t in enumerate(sta_ll.index)}

# voxel grid (projected km): covers cells + stations + 150 km margin
allx = np.concatenate([CXY[:, 0], SXY[:, 0]]); ally = np.concatenate([CXY[:, 1], SXY[:, 1]])
X0, Y0 = allx.min() - 150, ally.min() - 150
NX = int(np.ceil((allx.max() + 150 - X0) / DX)); NY = int(np.ceil((ally.max() + 150 - Y0) / DX))
NZ = int(ZMAX / DZ)
NV = NX * NY * NZ
# Slab2 depth per column (same source + filters as build_G_captured.py)
sl = pd.read_csv("data/cas_slab2_input_04-18.csv", usecols=["lat", "lon", "depth", "etype"])
sl = sl[~sl.etype.isin(["BA", "TO"])]; sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna()
cx = X0 + (np.arange(NX) + 0.5) * DX; cy = Y0 + (np.arange(NY) + 0.5) * DX
GX, GY = np.meshgrid(cx, cy, indexing="ij")
glat, glon = unproj(GX.ravel(), GY.ravel())
vs = griddata(sl[["lon", "lat"]].values, sl.depth.values, np.column_stack([glon, glat]), method="linear")
vsn = griddata(sl[["lon", "lat"]].values, sl.depth.values, np.column_stack([glon, glat]), method="nearest")
SLAB = np.where(np.isfinite(vs), vs, vsn).reshape(NX, NY)              # slab depth per column (km)
# nearest cell per column (within RCELL)
dcell, icell = cKDTree(CXY).query(np.column_stack([GX.ravel(), GY.ravel()]))
COLCELL = np.where(dcell <= RCELL, icell, -1).reshape(NX, NY)
zc = (np.arange(NZ) + 0.5) * DZ
IFACE = (np.abs(zc[None, None, :] - SLAB[:, :, None]) <= HSLAB)          # (NX,NY,NZ)
TOP5 = np.zeros((NX, NY, NZ), bool); TOP5[:, :, zc <= TOP] = True
IFACE_F = IFACE.ravel(); TOP5_F = TOP5.ravel()
CELL_F = np.repeat(COLCELL.ravel(), NZ)                                    # per voxel

# observed coda energy profile (fleet median rms_rel per 0.5-s stack-time bin) -> weights per lapse-past-S bin
try:
    pr = pd.read_csv("data/window_snr_profile.csv")
    prof = pr.groupby("bin_t0").rms_rel.median()
    W_OBS = np.array([float(prof.get(round(b + 0.75, 2), np.nan)) ** 2 for b in BIN_C])  # stack t0 = (bin start past S) + 1.0
    W_OBS = np.where(np.isfinite(W_OBS), W_OBS, np.nanmedian(W_OBS))
except Exception:
    W_OBS = np.ones(NB)


def hg_phase(mu, g):
    """Henyey-Greenstein phase function value (per steradian) at cos(theta)=mu."""
    if g == 0:
        return np.full_like(mu, 1.0 / (4 * np.pi))
    return (1 - g ** 2) / (4 * np.pi * (1 + g ** 2 - 2 * g * mu) ** 1.5)


def scatter_dirs(din, g, rng):
    """new directions after scattering from incoming unit vectors din (n,3): HG with anisotropy g (0 = isotropic)."""
    n = len(din)
    if g == 0:
        u = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); r = np.sqrt(1 - u ** 2)
        return np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
    xi = rng.uniform(0, 1, n)
    mu = (1 + g ** 2 - ((1 - g ** 2) / (1 - g + 2 * g * xi)) ** 2) / (2 * g)
    ph = rng.uniform(0, 2 * np.pi, n); sn = np.sqrt(np.clip(1 - mu ** 2, 0, None))
    # orthonormal frame around din
    a = np.where(np.abs(din[:, 2:3]) < 0.9, np.array([[0, 0, 1.0]]), np.array([[1.0, 0, 0]]))
    e1 = np.cross(din, a); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(din, e1)
    return mu[:, None] * din + sn[:, None] * (np.cos(ph)[:, None] * e1 + np.sin(ph)[:, None] * e2)


def vox_of(p):
    """flat voxel id (-1 outside domain) for points p (n,3) with z >= 0."""
    ix = np.floor((p[:, 0] - X0) / DX).astype(np.int64); iy = np.floor((p[:, 1] - Y0) / DX).astype(np.int64)
    iz = np.floor(p[:, 2] / DZ).astype(np.int64)
    ok = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY) & (iz >= 0) & (iz < NZ)
    v = np.where(ok, (ix * NY + iy) * NZ + iz, -1)
    return v


def rasterize(p0, p1, pid):
    """Sample straight segments p0->p1 (n,3) every DS km; z folded (|z|) = free-surface reflection.
    Returns (pid_s, vox_s, len_s) — out-of-domain samples get vox -1 (kept for total-length bookkeeping)."""
    d = p1 - p0; L = np.sqrt((d ** 2).sum(1))
    n = np.maximum(1, np.ceil(L / DS)).astype(np.int64)
    tot = int(n.sum())
    rep = np.repeat(np.arange(len(L)), n)
    k = np.arange(tot) - np.repeat(np.cumsum(n) - n, n)
    frac = (k + 0.5) / n[rep]
    p = p0[rep] + frac[:, None] * d[rep]
    p[:, 2] = np.abs(p[:, 2])
    return pid[rep], vox_of(p), (L / n)[rep]


class Acc:
    """Per-bin accumulators of contribution-weighted residence."""
    def __init__(self):
        self.tot = np.zeros(NB); self.iface = np.zeros(NB); self.top5 = np.zeros(NB); self.src = np.zeros(NB)
        self.cells = np.zeros((NB, M)); self.energy = np.zeros(NB)

    def add(self, b_of_p, c_of_p, pid_s, vox_s, len_s, srcmask_vox):
        """samples (pid, vox, len) weighted by c[pid], binned by b[pid] (b<0 = no contribution)."""
        b = b_of_p[pid_s]; ok = b >= 0
        if not ok.any():
            return
        b = b[ok]; w = c_of_p[pid_s[ok]] * len_s[ok]; v = vox_s[ok]
        self.tot += np.bincount(b, w, NB)
        ind = v >= 0
        if ind.any():
            bi, wi, vi = b[ind], w[ind], v[ind]
            fi = IFACE_F[vi]
            self.iface += np.bincount(bi, wi * fi, NB)
            self.top5 += np.bincount(bi, wi * TOP5_F[vi], NB)
            self.src += np.bincount(bi, wi * (fi & srcmask_vox[vi]), NB)
            cf = CELL_F[vi]; kk = fi & (cf >= 0)
            if kk.any():
                self.cells += np.bincount(bi[kk] * M + cf[kk], wi[kk], NB * M).reshape(NB, M)


def simulate(src, rec, ell, rng):
    """returns dict of per-bin residence fractions (all-order, single-scatter) + energies."""
    t_S = np.linalg.norm(rec - src) / BETA
    Lmax = BETA * (t_S + TMAX_PAST_S)
    ell_t = ell; ell = ell_t * (1 - HG_G)                                # scattering mfp at fixed transport mfp
    # source-cell mask: interface voxels within RCELL of the source horizontally
    dxs = (GX - src[0]) ** 2 + (GY - src[1]) ** 2
    srcmask = np.repeat((dxs <= RCELL ** 2).ravel(), NZ)
    acc_all, acc_ss = Acc(), Acc()
    n = NPART
    pos = np.tile(src, (n, 1)).astype(float)
    u = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n); r = np.sqrt(1 - u ** 2)
    dirv = np.column_stack([r * np.cos(ph), r * np.sin(ph), u])
    L = np.zeros(n); alive = np.ones(n, bool); order = np.zeros(n, int)
    hist = []                                                          # list of (pid, vox, len) arrays
    for leg in range(200):
        ia = np.where(alive)[0]
        if not len(ia):
            break
        s = rng.exponential(ell, len(ia))
        p1 = pos[ia] + s[:, None] * dirv[ia]
        ps, vs_, ls = rasterize(pos[ia], p1, ia)
        hist.append((ps, vs_, ls))
        crossed = p1[:, 2] < 0
        p1[:, 2] = np.abs(p1[:, 2]); dirv[ia[crossed], 2] *= -1
        pos[ia] = p1; L[ia] += s; order[ia] += 1
        # next-event estimation to the receiver
        dr = rec[None, :] - pos[ia]; Rr = np.sqrt((dr ** 2).sum(1))
        tp = (L[ia] + Rr) / BETA - t_S                                  # lapse past S at receiver
        b = np.floor((tp - BIN_EDGES[0]) / 0.5).astype(int)
        b[(tp < BIN_EDGES[0]) | (tp >= BIN_EDGES[-1])] = -1
        mu_r = (dirv[ia] * (dr / np.maximum(Rr, 1e-9)[:, None])).sum(1)  # cos(angle incoming -> receiver)
        c = hg_phase(mu_r, HG_G) * np.exp(-Rr / ell) / np.maximum(Rr, 0.5) ** 2 / n
        b_of_p = np.full(n, -1); c_of_p = np.zeros(n)
        b_of_p[ia] = b; c_of_p[ia] = np.where(b >= 0, c, 0.0)
        contrib = b >= 0
        if contrib.any():
            for acc, sel in ((acc_all, np.ones(len(ia), bool)), (acc_ss, order[ia] == 1)):
                m = contrib & sel
                if not m.any():
                    continue
                bb = b_of_p.copy(); bb[ia[~m]] = -1
                acc.energy += np.bincount(b[m], c[m], NB)
                # final leg
                pf, vf, lf = rasterize(pos[ia[m]], np.tile(rec, (m.sum(), 1)), ia[m])
                acc.add(bb, c_of_p, pf, vf, lf, srcmask)
                for (ps_, vs2, ls2) in hist:                            # full history of contributing particles
                    acc.add(bb, c_of_p, ps_, vs2, ls2, srcmask)
        dirv[ia] = scatter_dirs(dirv[ia], HG_G, rng)                     # re-scatter (HG; g=0 isotropic)
        alive[ia[L[ia] > Lmax]] = False
        # drop dead particles' history to bound memory
        if leg % 5 == 4:
            hist = [(p_, v_, l_) for (p_, v_, l_) in hist]
    # direct ray
    pr_, vr, lr = rasterize(src[None, :], rec[None, :], np.array([0]))
    okv = vr >= 0
    ray = dict(tot=lr.sum(), iface=(lr[okv] * IFACE_F[vr[okv]]).sum(), top5=(lr[okv] * TOP5_F[vr[okv]]).sum(),
               src=(lr[okv] * (IFACE_F[vr[okv]] & srcmask[vr[okv]])).sum())
    cf = CELL_F[vr[okv]]; kk = IFACE_F[vr[okv]] & (cf >= 0)
    ray["cells"] = np.bincount(cf[kk], lr[okv][kk], M) if kk.any() else np.zeros(M)
    return t_S, acc_all, acc_ss, ray


def window_kernel(t_S, acc, ray, win, wobs=True):
    """combine bins -> (t_bar, f_coda_*, net_*, K_eff cells)"""
    lo, hi = win
    inwin = (BIN_C > lo) & (BIN_C < hi)
    wb = (W_OBS if wobs else acc.energy) * inwin
    good = (acc.tot > 0) & (wb > 0)
    if not good.any():
        return None
    wb = wb * good; wb = wb / wb.sum()
    fb = lambda a: np.nansum(wb * np.where(good, a / np.maximum(acc.tot, 1e-30), 0.0))
    f_if, f_t5, f_sc = fb(acc.iface), fb(acc.top5), fb(acc.src)
    kc = (wb[:, None] * acc.cells / np.maximum(acc.tot, 1e-30)[:, None]).sum(0)
    t_bar = t_S + (wb * BIN_C).sum()
    g = lambda fc, fr: (t_bar * fc - t_S * fr) / (t_bar - t_S)
    rt = max(ray["tot"], 1e-9)
    net_if = g(f_if, ray["iface"] / rt); net_t5 = g(f_t5, ray["top5"] / rt); net_sc = g(f_sc, ray["src"] / rt)
    kcells = g(kc, ray["cells"] / rt)
    return dict(t_bar=t_bar, f_coda_iface=f_if, f_coda_top5=f_t5, f_coda_src=f_sc, f_ray_iface=ray["iface"] / rt,
                f_ray_top5=ray["top5"] / rt, net_iface=net_if, net_top5=net_t5, net_src=net_sc), kcells


def one_pair(i):
    rng = np.random.RandomState(1234 + i + SEED_OFFSET)
    src = np.array([CXY[pcell_idx[i], 0], CXY[pcell_idx[i], 1], CZ[pcell_idx[i]]])
    rec = np.array([sxy_of[ptag[i]][0], sxy_of[ptag[i]][1], 0.0])
    rows, kcells = [], {}
    for ell in ELLS:
        t_S, acc_all, acc_ss, ray = simulate(src, rec, ell, rng)
        for wname, win in WIN.items():
            for label, acc in (("all", acc_all), ("single", acc_ss)):
                for wobs in (True, False):
                    res = window_kernel(t_S, acc, ray, win, wobs)
                    if res is None:
                        continue
                    r, kc = res
                    r.update(pair=i, tag=ptag[i], cell=pcell[i], ell=ell, window=wname, orders=label,
                             weights="obs" if wobs else "mc", t_S=t_S, offset_km=float(np.linalg.norm(rec[:2] - src[:2])),
                             depth=src[2], e_win=float(acc.energy[(BIN_C > win[0]) & (BIN_C < win[1])].sum()))
                    rows.append(r)
                    if ell == 40.0 and label == "all" and wobs:
                        kcells[wname] = kc.astype(np.float32)
    return rows, kcells


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    npairs = len(ptag)
    sel = os.environ.get("PAIRS", "all")
    idx = np.arange(npairs) if sel == "all" else np.random.RandomState(0).choice(npairs, int(sel), replace=False)
    workers = int(os.environ.get("WORKERS", "8"))
    print(f"grid {NX}x{NY}x{NZ} voxels ({DX} km x {DZ} km) | interface voxels {IFACE_F.sum():,} | pairs {len(idx)} "
          f"| NPART {NPART} | ells* {ELLS} | HG g {HG_G} | HSLAB {HSLAB} | seed+{SEED_OFFSET} | tag '{OUTTAG}' "
          f"| W_OBS {np.round(W_OBS / W_OBS.max(), 3)}", flush=True)
    t0 = time.time(); rows = []; K24 = np.zeros((npairs, M), np.float32); K35 = np.zeros((npairs, M), np.float32)
    done = 0
    with ProcessPoolExecutor(workers) as ex:
        for r, kc in ex.map(one_pair, idx, chunksize=4):
            rows += r
            i = r[0]["pair"]
            if "w24" in kc: K24[i] = kc["w24"]
            if "w35" in kc: K35[i] = kc["w35"]
            done += 1
            if done % 200 == 0:
                print(f"  {done}/{len(idx)} pairs, {time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/residence_diag{OUTTAG}.csv", index=False)
    for wname, K in (("w24", K24), ("w35", K35)):
        np.savez_compressed(f"{OUT}/G_diff_{wname}{OUTTAG}.npz", G=-K, cell=cells, cxy=CXY, depth_km=CZ,
                            pair_tag=ptag, pair_cell=pcell, pairs_done=idx, ell=40.0, window=wname,
                            hg_g=HG_G, hslab=HSLAB, dx=DX, dz=DZ, npart=NPART, seed_offset=SEED_OFFSET,
                            note="G = -K_eff in ABSOLUTE residence units (row sums != -1 BY DESIGN: interface rows sum "
                                 "to ~-0.2; never renormalize rows). K_eff = S-pinned differential residence kernel, "
                                 "MC-RTE all orders, obs-profile lapse weights, interface +-HSLAB km, cells<=15 km")
    # summary
    d = df[(df.orders == "all") & (df.weights == "obs")]
    print("\n=== SUMMARY (all orders, observed-profile weights) ===")
    for (w, ell), g in d.groupby(["window", "ell"]):
        print(f"{w} ell {ell:3.0f}: t_bar-t_S {np.median(g.t_bar - g.t_S):.2f} s | f_ray_iface {g.f_ray_iface.median():.3f} "
              f"| f_coda_iface {g.f_coda_iface.median():.3f} | NET iface median {g.net_iface.median():+.3f} "
              f"IQR {g.net_iface.quantile(.25):+.3f}..{g.net_iface.quantile(.75):+.3f} neg {100 * (g.net_iface < 0).mean():.0f}% "
              f"| net top5 {g.net_top5.median():+.3f} (neg {100 * (g.net_top5 < 0).mean():.0f}%) | net srccell {g.net_src.median():+.3f}")
    d = df[(df.orders == "single") & (df.weights == "obs") & (df.ell == 40.0)]
    for w, g in d.groupby("window"):
        print(f"SINGLE-SCATTER {w} ell 40: net iface median {g.net_iface.median():+.3f} IQR {g.net_iface.quantile(.25):+.3f}.."
              f"{g.net_iface.quantile(.75):+.3f} neg {100 * (g.net_iface < 0).mean():.0f}% | naive f_coda {g.f_coda_iface.median():.3f} "
              f"| net top5 {g.net_top5.median():+.3f}")
    print(f"wrote {OUT}/residence_diag.csv, G_diff_w24.npz, G_diff_w35.npz  ({time.time() - t0:.0f}s)")
