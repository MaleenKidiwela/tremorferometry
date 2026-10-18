"""s11_r7_common -- shared helpers for S11 round 7 (PLAN.md §11). Region via S11_REGION (SVI default, WA_NW).
Read-only on every earlier S11 table; outputs under s11_depth/round7/<svi|wa_nw>/. eps = -dv/v; cell d = dv/v in %.
Contents: paths; a generic weighted linear model with per-series demeaning and exact null-space constraints (LinJ);
the E4 design (U1 units, D ⊥ PROJ, M, B) generalised to several deep unknowns per episode (one per zone, or S/M/D);
neighbour (ΔD) functionals for any group; a 1-D layered ray tracer (P and S travel time per depth group)."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_v2 as V
import s11_r3_fit as LF
import s11_r4_common as R4
import s11_r5_common as R5

BASE = R4.BASE
SUB = R4.SUB
R7 = f"{BASE}/round7"
OUT7 = f"{R7}/{SUB}"
TB7, LOG7, FIG7 = f"{OUT7}/tables", f"{OUT7}/logs", f"{R7}/figures"
for _d in (TB7, LOG7, FIG7):
    os.makedirs(_d, exist_ok=True)
U = C.episodes()
SEED = C.SEED
LSEL, NOM = "two3_20", "ell40"
WIN3 = LF.WIN3
GROUP_EDGES = {"s": (0.0, 5.0), "m": (5.0, 25.0), "d": (25.0, 45.0), "b": (45.0, np.inf)}


def log(name):
    return C.Log(f"{LOG7}/{name}.log")


def dump(obj, path):
    json.dump(obj, open(path, "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))


def pt_now():
    return C.pt_now()


def tbdir(test=False):
    d = f"{TB7}/test" if test else TB7
    os.makedirs(d, exist_ok=True)
    return d


def majors():
    return R4.majors()


# ------------------------------------------------------------------------------ generic linear model
class LinJ:
    """y = o_series + G th; weights w; th constrained by A th = 0 (exact null space); estimable functionals from a
    pseudo-inverse (column-scaled). names: tuples (group, zone/label, key)."""
    def __init__(self, G, names, series, w, A=None, rcond=1e-10):
        self.names = list(names); self.idx = {nm: i for i, nm in enumerate(self.names)}
        self.groups = np.array([nm[0] for nm in self.names])
        codes, uniq = pd.factorize(pd.Series(series).astype(str))
        self.codes = codes; self.nser = len(uniq); self.w = np.asarray(w, float)
        Wg = np.bincount(codes, self.w); self.Wg = Wg
        Gm = np.zeros((len(uniq), G.shape[1])); np.add.at(Gm, codes, self.w[:, None] * G)
        Gt = G - (Gm / Wg[:, None])[codes]
        self.Gt = Gt
        N = Gt.T @ (Gt * self.w[:, None])
        self.N = N
        P = N.shape[0]
        if A is not None and len(A):
            from scipy.linalg import null_space
            Z = null_space(np.asarray(A, float))
        else:
            Z = np.eye(P)
        Mz = Z.T @ N @ Z
        sc = np.sqrt(np.clip(np.diag(Mz), 1e-300, None))
        Mzs = Mz / sc[:, None] / sc[None, :]
        Mi = np.linalg.pinv(Mzs, rcond=rcond, hermitian=True) / sc[:, None] / sc[None, :]
        self.Ci = Z @ Mi @ Z.T
        self.Op = self.Ci @ (Gt.T * self.w[None, :])
        self.peff = float(np.trace(self.Ci @ N))

    def demean(self, y):
        m = np.bincount(self.codes, self.w * y) / self.Wg
        return y - m[self.codes]

    def solve(self, y, scale_chi2=True):
        th = self.Op @ y
        r = self.demean(y) - self.Gt @ th
        dof = max(len(y) - self.nser - self.peff, 1)
        chi2 = float(np.sum(self.w * r ** 2) / dof)
        Cm = self.Ci * (max(chi2, 1.0) if scale_chi2 else 1.0)
        return dict(th=th, C=Cm, chi2=chi2, idx=self.idx, names=self.names, groups=self.groups, resid=r, n=len(y))


def episodes_of(cells):
    return sorted(cells.episode.unique(), key=lambda e: U.mid[e])


def proj_rows(names_idx, P, grp, zone, eps_, proj=R5.PROJ):
    """constraint rows on the episode terms of (grp, zone): sum-to-zero + ⊥ proj at the episode mids"""
    iD = [names_idx[(grp, zone, e)] for e in eps_]
    A = []
    a = np.zeros(P); a[iD] = 1.0; A.append(a)
    if proj:
        em = U.mid.reindex(eps_); te, de = V.tdec(em.to_numpy()); bb = V.basis(te, de, proj)
        for tm in proj:
            if np.any(bb[tm] != 0):
                a = np.zeros(P); a[iD] = bb[tm]; A.append(a)
    return A


def e4_blocks(cells, ep_cols, units=True, MB=True, uterms=None, usmall=None, rt=None):
    """E4 design for coda cells. ep_cols: list of (group, zone, coefficient array) -> one episode term per (group, zone, e).
    Returns (G, names, eps). Units: U1 (coefficient 1) with uterms (>= 40 cells) / usmall (10-39); M, B: a_m, a_b x RT."""
    c = cells.reset_index(drop=True)
    mid = c.episode.map(U.mid); t, doy = V.tdec(mid)
    n = len(c); r = np.arange(n)
    eps_ = episodes_of(c); ie = {e: i for i, e in enumerate(eps_)}; col_e = np.array([ie[e] for e in c.episode])
    cols, names = [], []
    for grp, z, coef in ep_cols:
        Gd = np.zeros((n, len(eps_))); Gd[r, col_e] = coef
        cols.append(Gd); names += [(grp, z, e) for e in eps_]
    if units:
        old = (R5.STERMS, R5.SMALL)
        if uterms is not None:
            R5.STERMS, R5.SMALL = list(uterms), list(usmall if usmall is not None else ["c1", "s1"])
        try:
            lab, terms = R5.units_of(c)
            ST = list(R5.STERMS)
        finally:
            R5.STERMS, R5.SMALL = old
        allt = sorted(set(sum([list(v) for v in terms.values()], [])))
        bS = V.basis(t, doy, allt) if allt else {}
        for u_, tl in sorted(terms.items()):
            if not tl:
                continue
            m = lab == u_
            G = np.zeros((n, len(tl)))
            for j, tm in enumerate(tl):
                G[m, j] = bS[tm][m]
            cols.append(G); names += [("U", u_, tm) for tm in tl]
    if MB:
        RT = rt if rt is not None else R5.RT
        bR = V.basis(t, doy, RT)
        cols.append(np.column_stack([c.a_m.to_numpy(float) * bR[tm] for tm in RT])); names += [("M", "", tm) for tm in RT]
        cols.append(np.column_stack([c.a_b.to_numpy(float) * bR[tm] for tm in RT])); names += [("B", "", tm) for tm in RT]
    return np.column_stack(cols), names, eps_


def e4_model(cells, ep_cols, proj=R5.PROJ, **kw):
    G, names, eps_ = e4_blocks(cells, ep_cols, **kw)
    idx = {nm: i for i, nm in enumerate(names)}; P = len(names)
    A = []
    for grp, z, _ in ep_cols:
        A += proj_rows(idx, P, grp, z, eps_, proj)
    c = cells.reset_index(drop=True)
    L = LinJ(G, names, c.series.to_numpy(), 1.0 / c.sig.to_numpy() ** 2, A)
    L.eps = eps_
    return L


def dvec_g(F, eps_, e, grp, z=""):
    """ΔX(e) functional: X(e) minus the mean of the model episodes within ±3 yr (same as s11_r3_fit.dvec)"""
    if (grp, z, e) not in F["idx"]:
        return None
    oth = [x for x in eps_ if x != e and abs((U.mid[x] - U.mid[e]).days) <= WIN3 and (grp, z, x) in F["idx"]]
    if not oth:
        return None
    g = np.zeros(len(F["th"]))
    g[F["idx"][(grp, z, e)]] = 1.0
    for x in oth:
        g[F["idx"][(grp, z, x)]] -= 1.0 / len(oth)
    return g


def val_sig(F, g):
    return float(g @ F["th"]), float(np.sqrt(max(g @ F["C"] @ g, 0.0)))


def e4_reference(H, inj, setup):
    """round-5 E4 ΔD (L-SEL) on the same cells, as the round-6 size stage builds it -> (cells, dict e -> (dD, sig))"""
    import s11_r6_size as Z
    src = Z.fit_source(H, inj, "E4", setup)
    return src


# ------------------------------------------------------------------------------ ray tracing (1-D layers, flat earth)
class Rays:
    """P and S rays in the region's 1-D profile; travel time per depth group for a source depth and epicentral distance"""
    def __init__(self, vel_csv=None):
        vm = pd.read_csv(vel_csv or C.VEL)
        zt = vm.ztop_km.to_numpy(float)
        ed = sorted(set(zt.tolist()) | {5.0, 25.0, 45.0})
        self.z0 = np.array(ed)
        j = np.searchsorted(zt, self.z0, side="right") - 1
        self.vp = vm.vp.to_numpy(float)[j]; self.vs = vm.vs.to_numpy(float)[j]
        self.mid_group = None

    def _layers(self, zs):
        z0 = self.z0[self.z0 < zs]
        z1 = np.r_[z0[1:], zs]
        h = z1 - z0
        return z0, z1, h

    def trace(self, zs, x, phase="S"):
        """-> (total time s, {group: time s}, ray parameter); straight vertical if x ~ 0"""
        z0, z1, h = self._layers(zs)
        v = (self.vs if phase == "S" else self.vp)[: len(z0)]
        vmax = v.max()
        f = lambda p: np.sum(h * p * v / np.sqrt(np.clip(1 - (p * v) ** 2, 1e-15, None)))
        lo, hi = 0.0, (1.0 - 1e-12) / vmax
        if x <= 1e-6:
            p = 0.0
        else:
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if f(mid) < x:
                    lo = mid
                else:
                    hi = mid
            p = 0.5 * (lo + hi)
        tl = h / (v * np.sqrt(np.clip(1 - (p * v) ** 2, 1e-15, None)))
        out = {}
        zm = 0.5 * (z0 + z1)
        for g, (a, b) in GROUP_EDGES.items():
            m = (zm >= a) & (zm < b)
            out[g] = float(tl[m].sum())
        return float(tl.sum()), out, p


def ray_table(pairs, lg=None):
    """pairs: DataFrame fam, netsta -> columns tS, tP, tS_<g>, tP_<g> (s) from s11_kernel_mc.geometry()"""
    import s11_kernel_mc as KM
    F, S, _ = KM.geometry()
    R = Rays()
    rows = []
    for f, ns in pairs[["fam", "netsta"]].drop_duplicates().itertuples(index=False):
        f = int(f)
        if f not in F or ns not in S or not np.all(np.isfinite(F[f])) or not np.all(np.isfinite(S[ns])):
            continue
        src, rec = F[f], S[ns]
        x = float(np.hypot(src[0] - rec[0], src[1] - rec[1])); zs = float(src[2]) - float(max(rec[2], 0.0))
        tS, gS, _ = R.trace(zs, x, "S"); tP, gP, _ = R.trace(zs, x, "P")
        r = dict(fam=f, netsta=ns, x_km=x, zs_km=zs, tS=tS, tP=tP)
        for g in GROUP_EDGES:
            r[f"tS_{g}"] = gS[g]; r[f"tP_{g}"] = gP[g]
        rows.append(r)
    T = pd.DataFrame(rows)
    if lg is not None and len(T):
        lg(f"rays: {len(T)} family-station pairs; tS {T.tS.median():.2f} s, tS-tP {np.median(T.tS - T.tP):.2f} s (median); "
           f"S time share s/m/d/b " + "/".join(f"{np.median(T[f'tS_{g}'] / T.tS):.2f}" for g in GROUP_EDGES))
    return T
