"""s11_redpy_common -- shared code for the S11 REDPy-inspired pilot (rules: svi_science/results/s11_depth/redpy/PLAN.md).

Loading of S10 pilot windows, realignment, 3C S-anchored stretch, node building, all-pairs graph, reference-stack method.
eps = -dv/v (positive = slower)."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_measure as M

OUT = f"{S.ROOT}/svi_science/results/s11_depth/redpy"
TB, FIG, LOGD = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/logs"
PIL = S.OUT
EPS = np.linspace(-0.01, 0.01, 401)
DE = EPS[1] - EPS[0]
RAIL = 0.0099
WIN = (2.0, 5.0)                       # lapse after picked S
CLW = (-1.0, 5.0)                      # clustering window about picked S
# Weaver et al. (2011) a-priori stretch error scale for 2-8 Hz, lapse 2-5 s (relative weights; rescaled by misclosure)
_T, _WC = 1 / 6.0, 2 * np.pi * 5.0
WEAVER = np.sqrt(6 * np.sqrt(np.pi / 2) * _T / (_WC ** 2 * (WIN[1] ** 3 - WIN[0] ** 3)))
EPI = pd.read_csv(f"{S.ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")


# ------------------------------------------------------------------------------------------ data
def load_series(fam, ns, pmin=0.5):
    """-> dict per epoch: g (rows DataFrame, time-sorted), W (n,3,NW) realigned, G (3,NW), tS, Wc (control windows realigned,
    family's own c1 rows), gc (control rows), Wca (all-family control windows realigned, for the clustering threshold)"""
    rows = pd.read_parquet(f"{PIL}/tables/rows.parquet")
    rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet"); rs = rs[rs.netsta == ns][["row_id", "epoch"]]
    ser = pd.read_csv(f"{PIL}/tables/series.csv"); ser = ser[(ser.fam == fam) & (ser.netsta == ns)].set_index("epoch")
    z = np.load(f"{PIL}/win/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}
    zc = np.load(f"{PIL}/win/{ns}__c1f.npz"); posc = {r: i for i, r in enumerate(zc["row_id"])}
    d = rows[(rows.fam == fam) & (rows.series == "data") & (rows.P >= pmin)].merge(rs, on="row_id")
    d = d[d.row_id.isin(pos)]
    c = rows[(rows.series == "c1") & rows.row_id.isin(posc)].copy()
    c["epoch"] = S.epoch_of(ns, c.date)
    out = {}
    for ep, g in d.groupby("epoch"):
        if ep not in ser.index or len(g) < 30:
            continue
        g = g.sort_values("t").reset_index(drop=True)
        W = z["W"][[pos[r] for r in g.row_id]].astype(np.float64)
        ok = np.isfinite(W[:, 0, 0]) & (np.abs(np.nan_to_num(W[:, 0])).sum(1) > 0)
        g, W = g[ok].reset_index(drop=True), W[ok]
        tS = float(ser.loc[ep, "tS_pick"]); w = g.P.to_numpy()
        G = wstack(W, w)
        Wa = M.realign(W, G, tS)[0].astype(np.float64)
        G1 = wstack(Wa, w)
        cc_ = c[c.epoch == ep]
        Wca = zc["W"][[posc[r] for r in cc_.row_id]].astype(np.float64)
        okc = np.isfinite(Wca[:, 0, 0]) & (np.abs(np.nan_to_num(Wca[:, 0])).sum(1) > 0)
        cc_, Wca = cc_[okc].reset_index(drop=True), Wca[okc]
        if len(Wca):
            Wca = M.realign(Wca, G, tS)[0].astype(np.float64)
        own = (cc_.fam == fam).to_numpy()
        g["date"] = pd.to_datetime(g.date)
        g["cluster"] = M.clusters(g.date.dt.strftime("%Y-%m-%d"), EPI)
        gc = cc_[own].reset_index(drop=True); gc["date"] = pd.to_datetime(gc.date)
        gc["cluster"] = M.clusters(gc.date.dt.strftime("%Y-%m-%d"), EPI) if len(gc) else []
        out[ep] = dict(g=g, W=Wa, G=G1, tS=tS, Wca=Wca, Wc=Wca[own], gc=gc)
    return out


def wstack(W, w):
    """nan-aware weighted mean over events per component -> (3, NW)"""
    out = np.full(W.shape[1:], np.nan)
    for c in range(W.shape[1]):
        ok = np.isfinite(W[:, c, 0])
        if ok.any() and w[ok].sum() > 0:
            out[c] = (w[ok, None] * W[ok, c]).sum(0) / w[ok].sum()
    return out


# ------------------------------------------------------------------------------------------ single-event cc
def feat(W, tS, win=CLW):
    """normalised 3C feature vectors on tS+win: per component demeaned, unit norm, / sqrt(3); missing comps = 0"""
    m = (S.TAU >= tS + win[0]) & (S.TAU <= tS + win[1])
    X = W[:, :, m].copy()
    X = np.nan_to_num(X - np.nanmean(X, 2, keepdims=True))
    n = np.linalg.norm(X, axis=2, keepdims=True); n[n == 0] = 1
    return (X / n / np.sqrt(3)).reshape(len(W), -1)


# ------------------------------------------------------------------------------------------ 3C stretch
class Ref3:
    """3C S-anchored stretch bank: ref r (3,NW) stretched r(tau'/(1+eps)) about S; cc = mean of component cc's"""
    def __init__(self, r, tS, win=WIN):
        from scipy.interpolate import CubicSpline
        tr = S.TAU - tS
        self.idx = np.flatnonzero((tr >= win[0]) & (tr <= win[1])); x = tr[self.idx]
        self.comps, self.B = [], []
        for c in range(3):
            if not np.isfinite(r[c]).all():
                continue
            B = CubicSpline(tr, r[c])(x[None, :] / (1.0 + EPS[:, None]))
            B = B - B.mean(1, keepdims=True); B /= np.linalg.norm(B, axis=1, keepdims=True)
            self.comps.append(c); self.B.append(B)

    def measure(self, C):
        C = np.asarray(C)
        if C.ndim == 2:
            C = C[None]
        cc = np.zeros((C.shape[0], len(EPS))); nc = np.zeros((C.shape[0], 1))
        for c, B in zip(self.comps, self.B):
            X = C[:, c][:, self.idx]
            okr = np.isfinite(X).all(1)
            X = np.nan_to_num(X - np.nanmean(X, 1, keepdims=True))
            nx = np.linalg.norm(X, axis=1, keepdims=True); nx[nx == 0] = 1
            cc[okr] += (X[okr] / nx[okr]) @ B.T; nc[okr, 0] += 1
        cc = cc / np.maximum(nc, 1)
        i = np.argmax(cc, 1); ii = np.clip(i, 1, len(EPS) - 2); r = np.arange(len(i))
        y0, y1, y2 = cc[r, ii - 1], cc[r, ii], cc[r, ii + 1]; d2 = y0 - 2 * y1 + y2
        delta = np.where((d2 < 0) & (i == ii), 0.5 * (y0 - y2) / np.where(d2 == 0, -1, d2), 0.0)
        return EPS[i] + np.clip(delta, -1, 1) * DE, cc[r, i]


def inject(W, tS, a):
    """stretch windows about S: cur(tau') = x(tau'/(1+a)) -> reads eps = +a"""
    from scipy.interpolate import CubicSpline
    tr = S.TAU - tS; out = W.copy()
    for c in range(3):
        ok = np.isfinite(W[:, c, 0])
        if ok.any():
            out[ok, c] = CubicSpline(tr, W[ok, c], axis=1)(tr / (1.0 + a))
    return out


# ------------------------------------------------------------------------------------------ nodes
def make_nodes(g, kind, N=None, min_day=3):
    """-> list of (label, event index array) ; g time-sorted with 'date', 'cluster'"""
    nodes = []
    if kind == "day":
        for d, q in g.groupby("date"):
            if len(q) >= min_day:
                nodes.append((str(d.date()), q.index.to_numpy()))
    elif kind == "episode":
        for cl, q in g.groupby("cluster", sort=False):
            if len(q) >= 10:
                nodes.append((cl, q.index.to_numpy()))
    elif kind == "chunk":
        for cl, q in g.groupby("cluster", sort=False):
            idx = q.index.to_numpy()
            if len(idx) < N:
                if len(idx) >= N / 2:
                    nodes.append((f"{cl}#0", idx))
                continue
            k = len(idx) // N; parts = [idx[i * N:(i + 1) * N] for i in range(k)]; rem = idx[k * N:]
            if len(rem) >= N / 2:
                parts.append(rem)
            elif len(rem):
                parts[-1] = np.r_[parts[-1], rem]
            nodes += [(f"{cl}#{i}", p) for i, p in enumerate(parts)]
    nodes.sort(key=lambda x: g.t.iloc[x[1]].median())
    return nodes


def node_stacks(W, w, nodes):
    return np.array([wstack(W[ix], w[ix]) for _, ix in nodes])


# ------------------------------------------------------------------------------------------ graph
def pairs(St, tS):
    n = len(St); pi, pj, dd, cc = [], [], [], []
    for i in range(n - 1):
        e, c = Ref3(St[i], tS).measure(St[i + 1:])
        pi += [i] * (n - 1 - i); pj += list(range(i + 1, n)); dd += list(e); cc += list(c)
    return np.array(pi, int), np.array(pj, int), np.array(dd), np.array(cc)


def solve(n, pi, pj, d, cc, huber=4):
    """v_j - v_i = d; a-priori sigma Weaver(cc); scale by chi2 of the plain LS, then Huber IRLS.
    -> v (zero mean over connected nodes), sigma_v (misclosure-scaled), chi2_red (vs rescaled sigma), scale, n_pairs"""
    keep = np.abs(d) < RAIL
    pi, pj, d, cc = pi[keep], pj[keep], d[keep], np.clip(cc[keep], 0.02, 0.999)
    if len(pi) < 2:
        return np.full(n, np.nan), np.full(n, np.nan), np.nan, np.nan, 0
    sig = np.sqrt(1 - cc ** 2) / (2 * cc) * WEAVER
    A = np.zeros((len(pi), n)); A[np.arange(len(pi)), pi] = -1; A[np.arange(len(pi)), pj] = 1
    dof = max(len(pi) - n + 1, 1)

    def ls(w):
        N = A.T @ (A * w[:, None]); Ni = np.linalg.pinv(N, rcond=1e-10); return Ni @ (A.T @ (w * d)), Ni
    w0 = 1 / sig ** 2; v, Ni = ls(w0)
    scale = np.sqrt(max(np.sum(w0 * (A @ v - d) ** 2) / dof, 1e-12)); sig = sig * scale; w0 = 1 / sig ** 2; w = w0.copy()
    for _ in range(huber):
        v, Ni = ls(w); r = (A @ v - d) / sig
        w = w0 * np.where(np.abs(r) <= 1.345, 1.0, 1.345 / np.maximum(np.abs(r), 1e-12))
    v, Ni = ls(w)
    chi2 = float(np.sum(w * (A @ v - d) ** 2) / dof)
    deg = np.bincount(np.r_[pi, pj], minlength=n)
    sv = np.sqrt(np.clip(np.diag(Ni), 0, None)) * np.sqrt(max(chi2, 1.0))
    v = v - v[deg > 0].mean(); v[deg == 0] = np.nan; sv[deg == 0] = np.nan
    return v, sv, chi2, scale, int(len(pi))


def triplet_misclosure(n, pi, pj, d, rng, nmax=20000):
    keep = np.abs(d) < RAIL
    D = np.full((n, n), np.nan); D[pi[keep], pj[keep]] = d[keep]
    if n < 3:
        return np.nan
    T = np.sort(rng.integers(0, n, (nmax * 3, 3)), 1); T = T[(T[:, 0] < T[:, 1]) & (T[:, 1] < T[:, 2])][:nmax]
    m = D[T[:, 0], T[:, 1]] + D[T[:, 1], T[:, 2]] - D[T[:, 0], T[:, 2]]
    m = m[np.isfinite(m)]
    return float(np.sqrt(np.mean(m ** 2))) if len(m) > 10 else np.nan


def graph_run(W, w, nodes, tS):
    St = node_stacks(W, w, nodes)
    pi, pj, d, cc = pairs(St, tS)
    v, sv, chi2, scale, npair = solve(len(nodes), pi, pj, d, cc)
    return dict(v=v, sv=sv, chi2=chi2, scale=scale, npair=npair, pi=pi, pj=pj, d=d, cc=cc, St=St)


def ref_run(W, w, nodes, tS, rng=None, nboot=0):
    """reference-stack method: node vs leave-node-out stack of the node set's events"""
    allix = np.unique(np.concatenate([ix for _, ix in nodes]))
    Wsum = np.zeros(W.shape[1:]); wsum = np.zeros(3)
    for c in range(3):
        ok = np.isfinite(W[allix, c, 0])
        Wsum[c] = (w[allix][ok, None] * W[allix][ok, c]).sum(0); wsum[c] = w[allix][ok].sum()
    eps, cc, sb = [], [], []
    for _, ix in nodes:
        st = wstack(W[ix], w[ix]); sub = np.zeros_like(Wsum); ws = np.zeros(3)
        for c in range(3):
            ok = np.isfinite(W[ix, c, 0]); sub[c] = (w[ix][ok, None] * W[ix][ok, c]).sum(0); ws[c] = w[ix][ok].sum()
        ref = (Wsum - sub) / np.maximum(wsum - ws, 1e-9)[:, None]
        ref[wsum - ws <= 0] = np.nan
        R = Ref3(ref, tS); e, c_ = R.measure(st); eps.append(e[0]); cc.append(c_[0])
        if nboot:
            bs = []
            for _ in range(nboot):
                b = rng.choice(ix, len(ix)); bs.append(wstack(W[b], w[b]))
            sb.append(float(np.std(R.measure(np.array(bs))[0], ddof=1)))
        else:
            sb.append(np.nan)
    eps = np.array(eps); bad = np.abs(eps) >= RAIL
    return dict(eps=np.where(bad, np.nan, eps), eps_raw=eps, cc=np.array(cc), sb=np.array(sb))


def ensure():
    for d in (TB, FIG, LOGD):
        os.makedirs(d, exist_ok=True)
