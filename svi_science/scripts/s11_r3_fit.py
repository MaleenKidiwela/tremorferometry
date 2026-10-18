"""s11_r3_fit -- linear-operator versions of the S11 smooth (s11_v2.fit) and temporal (s11_q3.fit) models, plus the
along-strike zone model (PLAN §6.3). Same design, weights, demeaning (series offsets), prior and chi2 scaling as the
originals; the operator Op = Ci Gt^T W maps any cell-data vector (natural or injected difference) to parameters, so
injection responses cost one matrix product. Checked against the original fits (check_against_original)."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I
import s11_v2 as V

U = C.episodes()


def _demean_G(cells, G, w):
    key = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    codes, uniq = pd.factorize(key)
    W = np.bincount(codes, w)
    Gm = np.zeros((len(uniq), G.shape[1])); np.add.at(Gm, codes, w[:, None] * G)
    return G - (Gm / W[:, None])[codes], codes, W, len(uniq)


class Lin:
    """generic: cells (with a_s, a_d, sig, site, episode, d), design G and parameter names"""
    def __init__(self, cells, G, names, prior=1.0, extra=None):
        self.cells = cells.reset_index(drop=True); self.names = names; self.idx = {n: i for i, n in enumerate(names)}
        self.w = 1.0 / self.cells.sig.to_numpy() ** 2
        self.Gt, self.codes, self.Wg, self.nser = _demean_G(self.cells, G, self.w)
        self.N = self.Gt.T @ (self.Gt * self.w[:, None])
        self.prior = prior
        self.Ci = np.linalg.inv(self.N + np.eye(len(names)) / prior ** 2)
        self.Op = self.Ci @ (self.Gt.T * self.w[None, :])
        self.peff = float(np.trace(self.Ci @ self.N))
        self.extra = extra or {}

    def demean_d(self, d):
        dm = np.bincount(self.codes, self.w * d) / self.Wg
        return d - dm[self.codes]

    def solve(self, d=None):
        d = self.cells.d.to_numpy() if d is None else d
        th = self.Op @ d                         # Gt^T W M d = 0 (Gt weighted-demeaned) -> no need to demean d
        r = self.demean_d(d) - self.Gt @ th
        chi2 = float(np.sum(self.w * r ** 2) / max(len(d) - self.nser - self.peff, 1))
        F = dict(th=th, C=self.Ci * max(chi2, 1.0), idx=self.idx, names=self.names, chi2=chi2, n=len(d), nser=self.nser,
                 resid=r, wr2=self.w * r ** 2)
        F.update(self.extra)
        return F

    def delta(self, d):
        """parameter response to a data difference (no chi2; for injections)"""
        F = dict(th=self.Op @ d, C=self.Ci, idx=self.idx, names=self.names); F.update(self.extra); return F


def smooth(cells, model="smooth", prior=1.0, mid=None):
    Sterms, Dterms = V.MODELS[model]
    m = cells.episode.map(U.mid) if mid is None else mid
    t, doy = V.tdec(m)
    bS, bD = V.basis(t, doy, Sterms), V.basis(t, doy, Dterms)
    names = [("D", "", tm) for tm in Dterms]
    sites = sorted(cells.site.unique())
    names += [("S", s, tm) for s in sites for tm in Sterms]
    idx = {nm: i for i, nm in enumerate(names)}
    G = np.zeros((len(cells), len(names)))
    ad, as_ = cells.a_d.to_numpy(), cells.a_s.to_numpy(); site = cells.site.to_numpy()
    for tm in Dterms:
        G[:, idx[("D", "", tm)]] = ad * bD[tm]
    for s in sites:
        ms = site == s
        for tm in Sterms:
            G[ms, idx[("S", s, tm)]] = as_[ms] * bS[tm][ms]
    return Lin(cells, G, names, prior)


def temporal(cells, prior=1.0, zcols=None):
    """temporal model (s11_q3.fit); zcols = list of zone coefficient columns -> D_z(e) per zone (PLAN §6.3)"""
    m = cells.episode.map(U.mid); t, doy = V.tdec(m)
    Sterms = V.MODELS["smooth"][0]; bS = V.basis(t, doy, Sterms)
    eps_ = sorted(cells.episode.unique(), key=lambda e: U.mid[e])
    zs = [""] if zcols is None else list(zcols)
    names = [("D", z, e) for z in zs for e in eps_]
    sites = sorted(cells.site.unique())
    names += [("S", s, tm) for s in sites for tm in Sterms]
    idx = {n: i for i, n in enumerate(names)}
    G = np.zeros((len(cells), len(names)))
    r = np.arange(len(cells))
    for z in zs:
        col = cells.a_d.to_numpy() if z == "" else cells[z].to_numpy()
        G[r, [idx[("D", z, e)] for e in cells.episode]] = col
    for s in sites:
        ms = (cells.site == s).to_numpy()
        for tm in Sterms:
            G[ms, idx[("S", s, tm)]] = cells.a_s.to_numpy()[ms] * bS[tm][ms]
    return Lin(cells, G, names, prior, extra=dict(eps=eps_, U=U, ncell=cells.groupby("episode").size(), zs=zs))


# ------------------------------------------------------------------------------ Q3 functionals
WIN3 = 3 * 365.25


def dvec(F, e, z=""):
    g = np.zeros(len(F["th"]))
    oth = [x for x in F["eps"] if x != e and abs((U.mid[x] - U.mid[e]).days) <= WIN3]
    if not oth or ("D", z, e) not in F["idx"]:
        return None
    g[F["idx"][("D", z, e)]] = 1.0
    for x in oth:
        g[F["idx"][("D", z, x)]] -= 1.0 / len(oth)
    return g


def q3_table(F, majors, z=""):
    rows, gs = [], []
    for e in F["eps"]:
        g = dvec(F, e, z)
        if g is None:
            continue
        rows.append(dict(episode=e, zone=z, mid=U.mid[e], major=e in majors, n_cells=int(F["ncell"].get(e, 0)),
                         dD=float(g @ F["th"]), sig=float(np.sqrt(g @ F["C"] @ g))))
        if e in majors:
            gs.append(g)
    T = pd.DataFrame(rows)
    comp = None
    if gs:
        gc = np.mean(gs, 0)
        comp = dict(composite=float(gc @ F["th"]), composite_sig=float(np.sqrt(gc @ F["C"] @ gc)), n_majors=len(gs), g=gc)
    return T, comp


def comp_vec(F, majors, z=""):
    gs = [dvec(F, e, z) for e in majors]
    gs = [g for g in gs if g is not None]
    return np.mean(gs, 0) if gs else None


def Qsmooth(F):
    A, sA, pk, cs = V.amp(F)
    g = V.q2vec(F)
    return dict(Q1=A, Q1_sig=sA, Q1_peak=pk, Q2=float(g @ F["th"]), Q2_sig=float(np.sqrt(g @ F["C"] @ g)), chi2=F.get("chi2", np.nan),
                n_cells=F.get("n", np.nan))


def check_against_original(cells, lg):
    """sanity: operator fits reproduce s11_v2.fit and s11_q3.fit on the natural data"""
    import s11_q3 as Q3
    a = smooth(cells).solve(); b = V.fit(cells, "smooth")
    t1 = float(np.max(np.abs(a["th"] - b["th"])))
    c = temporal(cells).solve(); d = Q3.fit(cells)
    t2 = float(np.max(np.abs(c["th"] - d["th"])))
    lg(f"operator check: smooth max |dth| {t1:.1e}, temporal max |dth| {t2:.1e}")
    return max(t1, t2)
