"""s11_dist_model -- cells and models for the S11 distant-station lever (distant/PLAN.md §6-§8).
Near cells = the region's headline n >= 20 cells exactly as s11_q3.cells_for (SVI tag _all = 115 families; WA tag "").
Far cells = the same rule on distant/<region>/tables (s11_invert.load_cells / prepare run with C.TB pointed there).
Model: d = o_series + a_s S_site(t) + a_d D + a_m M(t) + a_b B(t); S per site = 2 harmonics + slow (s11_v2 'smooth');
D = slow + 1 harmonic ('smooth') or free per episode ('temporal'); M, B = slow + 1 harmonic per region (optional).
Read-only on every near table."""
from __future__ import annotations
import os, sys
from contextlib import contextmanager
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_dist_common as D
import s11_invert as I
import s11_v2 as V
import s11_q3 as Q3

NEAR_TAG = "_all" if D.REG == "SVI" else ""
NMIN, CC2 = 20, 0.35
VARS4 = ("ell20", "ell40", "ell80", "two5_40")
RT = V.SLOW + ["c1", "s1"]
MAJ = list(Q3.MAJ15)
U = C.episodes()
FAR_TAG = os.environ.get("S11_DIST_TAG", "")          # test runs only (tables with a _test suffix)


@contextmanager
def tb(path):
    old = C.TB; C.TB = path
    try:
        yield
    finally:
        C.TB = old


def _kern(var, far):
    if far:
        K = pd.read_parquet(f"{D.DTB}/kernel_coef{FAR_TAG}.parquet")
        return K[(K.variant == var) & (K.comp == "3c")]
    return I.kernels(var)


def with_variant(cells, var, far, win=C.HEAD):
    """copy of cells (same index/order) with a_s/a_m/a_d/a_b of a kernel variant ('AVG4' = mean of the 4 variants)"""
    out = cells.copy()
    if var == "AVG4":
        parts = [with_variant(cells, v, far, win)[["a_s", "a_m", "a_d", "a_b"]].to_numpy() for v in VARS4]
        out[["a_s", "a_m", "a_d", "a_b"]] = np.mean(parts, 0); return out
    base = cells[["fam", "netsta", "epoch", "episode"]].reset_index(drop=True)
    k = I.attach_kernels(base, _kern(var, far), win)       # merge keeps the row order of base (left join)
    out[["a_s", "a_m", "a_d", "a_b"]] = k[["a_s", "a_m", "a_d", "a_b"]].to_numpy()
    return out


def near_cells(win=C.HEAD, control=False):
    if control:
        c = I.load_cells("raw", C.HEAD, NEAR_TAG, control=True, nmin=NMIN, ccmin=CC2)
        c = I.prepare("raw", C.HEAD, NEAR_TAG, cells=c)
    else:
        c = I.prepare("raw", win, NEAR_TAG, nmin=NMIN, ccmin=CC2)
    c["far"] = False
    return c


def site_map(far_stas):
    XY = D.station_coords()
    near = dict(C.SITE)
    out = {}
    for ns in far_stas:
        if ns in near:
            out[ns] = near[ns]; continue
        best = None
        for m, s in near.items():
            if m in XY.index and D.km(XY.lat[ns], XY.lon[ns], XY.lat[m], XY.lon[m]) < D.SITE_KM:
                best = s; break
        out[ns] = best
    rest = [ns for ns in far_stas if out[ns] is None]
    lab = {ns: i for i, ns in enumerate(rest)}
    for i, a in enumerate(rest):
        for b in rest[i + 1:]:
            if D.km(XY.lat[a], XY.lon[a], XY.lat[b], XY.lon[b]) < D.SITE_KM:
                la, lb = lab[a], lab[b]
                for k in lab:
                    if lab[k] == lb:
                        lab[k] = la
    grp = {}
    for ns, g in lab.items():
        grp.setdefault(g, []).append(ns)
    for g, mem in grp.items():
        name = "+".join(sorted(s.split(".")[1] for s in mem))
        for ns in mem:
            out[ns] = name
    return out


def far_cells(win=C.HEAD, control=False):
    tg = FAR_TAG
    if not os.path.exists(f"{D.DTB}/episode_raw{tg}.parquet") or not os.path.exists(f"{D.DTB}/kernel_coef{tg}.parquet"):
        return pd.DataFrame()
    with tb(D.DTB):                         # = s11_invert.prepare, with the far tag also on the kernel table
        if control:
            c = I.load_cells("raw", C.HEAD, tg, control=True, nmin=NMIN, ccmin=CC2)
        else:
            c = I.load_cells("raw", win, tg, nmin=NMIN, ccmin=CC2)
        c = I.attach_kernels(c, _kern("ell40", True), win if not control else C.HEAD)
        c = c[np.isfinite(c.a_d) & np.isfinite(c.a_s)]
        inf, _ = I.control_inflation("raw", win if not control else C.HEAD, tg)
        c["infl"] = [inf.get((f, n, e), 1.0) for f, n, e in zip(c.fam, c.netsta, c.epoch)]
        c["sig_raw"] = c.sig; c["sig"] = c.sig * c.infl
        key = c.fam.astype(str) + "|" + c.netsta + "|" + c.epoch
        c = c[key.map(key.value_counts()) >= 3].reset_index(drop=True)
    if not len(c):
        return c
    sm = site_map(sorted(c.netsta.unique()))
    c["site"] = c.netsta.map(sm); c["far"] = True
    return c


def all_cells(win=C.HEAD, control=False, variant="ell40"):
    n = near_cells(win, control); f = far_cells(win, control)
    if variant != "ell40":
        n = with_variant(n, variant, False, win)
        if len(f):
            f = with_variant(f, variant, True, win)
    return pd.concat([n, f], ignore_index=True) if len(f) else n


# ------------------------------------------------------------------------------ model
def fit(cells, kind="smooth", MB=True, prior=1.0, d=None, mid=None, Sc=False):
    """Sc (temporal only, PLAN Amendment 1): an episode-locked REGIONAL shallow term Sc(e) x a_s, free per episode"""
    m = cells.episode.map(U.mid) if mid is None else mid
    t, doy = V.tdec(m)
    Sterms = V.MODELS["smooth"][0]
    bS = V.basis(t, doy, Sterms); bR = V.basis(t, doy, RT)
    names = []
    if kind == "smooth":
        names += [("D", "", tm) for tm in RT]; eps_ = None
    else:
        eps_ = sorted(cells.episode.unique(), key=lambda e: U.mid[e]); names += [("D", "", e) for e in eps_]
    sites = sorted(cells.site.unique())
    names += [("S", s, tm) for s in sites for tm in Sterms]
    if MB:
        names += [("M", "", tm) for tm in RT] + [("B", "", tm) for tm in RT]
    if Sc and kind != "smooth":
        names += [("Sc", "", e) for e in eps_]
    idx = {nm: i for i, nm in enumerate(names)}
    G = np.zeros((len(cells), len(names)))
    ad, as_, am, ab = (cells[c].to_numpy() for c in ("a_d", "a_s", "a_m", "a_b"))
    if kind == "smooth":
        for tm in RT:
            G[:, idx[("D", "", tm)]] = ad * bR[tm]
    else:
        G[np.arange(len(cells)), [idx[("D", "", e)] for e in cells.episode]] = ad
    site = cells.site.to_numpy()
    for s in sites:
        ms = site == s
        for tm in Sterms:
            G[ms, idx[("S", s, tm)]] = as_[ms] * bS[tm][ms]
    if MB:
        for tm in RT:
            G[:, idx[("M", "", tm)]] = am * bR[tm]; G[:, idx[("B", "", tm)]] = ab * bR[tm]
    if Sc and kind != "smooth":
        G[np.arange(len(cells)), [idx[("Sc", "", e)] for e in cells.episode]] = as_
    dd = cells.d.to_numpy() if d is None else d
    w = 1 / cells.sig.to_numpy() ** 2
    Gt, dt, ns = I.demean(cells, G, dd, w)
    N = Gt.T @ (Gt * w[:, None]); Ci = np.linalg.inv(N + np.eye(len(names)) / prior ** 2)
    th = Ci @ (Gt.T @ (w * dt)); r = dt - Gt @ th
    chi = float(np.sum(w * r ** 2) / max(len(dd) - ns - np.trace(Ci @ N), 1))
    return dict(th=th, C=Ci * max(chi, 1.0), Ci=Ci, idx=idx, names=names, eps=eps_, chi2=chi, n=len(dd), nser=ns,
                ncell=cells.groupby("episode").size(), U=U)


def Qs(F):
    """smooth-model quantities"""
    q = V.Q(F)
    for k in ("M", "B"):
        if (k, "", "c1") in F["idx"]:
            A, sA, pk, _ = V.amp(F, k)
            q[f"{k}_annual"] = A; q[f"{k}_annual_sig"] = sA; q[f"{k}_peak"] = pk
    return q


def q3tab(F):
    T, comp = Q3.q3(F)
    return T, comp or {}


def comp_vec(F):
    gs = [Q3.dvec(F, e) for e in MAJ if ("D", "", e) in F["idx"]]
    gs = [g for g in gs if g is not None]
    return np.mean(gs, 0) if gs else None


def diff_fit(F, dth):
    Fd = dict(F); Fd["th"] = dth
    return Fd
