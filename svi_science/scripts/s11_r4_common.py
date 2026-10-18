"""s11_r4_common -- shared helpers for S11 round 4 (PLAN.md §8). Region via S11_REGION (SVI default, WA_NW).
Read-only on every earlier S11 table; outputs under s11_depth/round4/<svi|wa_nw>/. eps = -dv/v; cell d = dv/v in %.
Contents: paths; kernel coefficients (near + far, all variants, a_s/a_m/a_d/a_b); B2 cell builders (near, far); the combined
model (D, Sc per episode with exact sum-to-zero, S per site smooth, R per station, M, B) as a linear operator."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I
import s11_v2 as V
import s11_r3_common as R

BASE = f"{C.ROOT}/svi_science/results/s11_depth"
SUB = "svi" if C.REGION == "SVI" else C.REGION.lower()
R4 = f"{BASE}/round4"
OUT4 = f"{R4}/{SUB}"
TB4, LOG4, FIG4 = f"{OUT4}/tables", f"{OUT4}/logs", f"{R4}/figures"
for _d in (TB4, LOG4, FIG4):
    os.makedirs(_d, exist_ok=True)
DIST = f"{BASE}/distant/{SUB}/tables"
V6 = R.V6
NODE_MIN, CTRL_NODE_MIN = 50, 10
RMIN = 20                     # receiver term only for stations with >= RMIN cells (PLAN §8.1)
U = C.episodes()
SLOW = V.SLOW
STERMS = V.MODELS["smooth"][0]
RT = SLOW + ["c1", "s1"]
RTERMS = ["c1", "s1", "lin"]


def log(name):
    return C.Log(f"{LOG4}/{name}.log")


def dump(obj, path):
    json.dump(obj, open(path, "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))


def majors():
    return R.majors()


def lapse_sel():
    v = json.load(open(f"{R.TB3}/kerr_verdict.json"))
    return v["lapse"]["sel"]


# ------------------------------------------------------------------------------ kernels (near + far, every variant)
_K = {}


def coef_table():
    if "K" in _K:
        return _K["K"]
    parts = [R.kernel_table().assign(far=False)]
    p = f"{DIST}/kernel_coef.parquet"
    if os.path.exists(p):
        k = pd.read_parquet(p); k = k[k.comp == "3c"]
        parts.append(k[["series", "fam", "netsta", "variant", "win", "a_s", "a_m", "a_d", "a_b"]].assign(far=True))
    p = f"{TB4}/kernel_coef_far_r4.parquet"
    if os.path.exists(p):
        k = pd.read_parquet(p); k = k[k.variant.isin(["ell10", "two3_20"])]
        parts.append(k[["series", "fam", "netsta", "variant", "win", "a_s", "a_m", "a_d", "a_b"]].assign(far=True))
    K = pd.concat(parts, ignore_index=True).drop_duplicates(["series", "variant", "win"])
    _K["K"] = K
    return K


def _one(cells, variant, win):
    K = coef_table(); k = K[(K.variant == variant) & (K.win == win)]
    ks = k.set_index("series")[["a_s", "a_m", "a_d", "a_b"]]
    key = (cells.fam.astype(int).astype(str) + "|" + cells.netsta + "|" + cells.epoch).to_numpy()
    out = ks.reindex(key).to_numpy(float)
    miss = ~np.isfinite(out[:, 2])
    if miss.any():
        alt = k.groupby(["fam", "netsta"])[["a_s", "a_m", "a_d", "a_b"]].mean()
        for i in np.flatnonzero(miss):
            kk = (int(cells.fam.iloc[i]), cells.netsta.iloc[i])
            if kk in alt.index:
                out[i] = alt.loc[kk].to_numpy(float)
    return out


def coefs(cells, setup, win=C.HEAD):
    """setup: variant name | 'AVG6' | ('avg', [variants]) -> (n, 4) a_s, a_m, a_d, a_b"""
    if setup == "AVG6":
        setup = ("avg", V6)
    if isinstance(setup, str):
        return _one(cells, setup, win)
    kind, arg = setup
    assert kind == "avg"
    return np.nanmean([_one(cells, v, win) for v in arg], 0)


def with_setup(cells, setup, win=C.HEAD):
    c = cells.copy()
    c[["a_s", "a_m", "a_d", "a_b"]] = coefs(c, setup, win)
    return c


# ------------------------------------------------------------------------------ sites
def far_site_map(stas):
    import s11_dist_model as Mo
    return Mo.site_map(sorted(stas))


# ------------------------------------------------------------------------------ cells
def r4_episode(far=False):
    if far:
        p = f"{TB4}/r4_episode_far.parquet"
        return pd.read_parquet(p) if os.path.exists(p) else pd.DataFrame()
    E = [pd.read_parquet(f"{TB4}/r4_episode.parquet")]
    if C.REGION == "SVI" and os.path.exists(f"{TB4}/r4_episode_ext.parquet"):
        E.append(pd.read_parquet(f"{TB4}/r4_episode_ext.parquet"))
    return pd.concat(E, ignore_index=True)


def inflate(cells, ctrl):
    """= s11_r3_q3.inflate: per series sqrt(max(1, mean (d/sig)^2 of its control values)), >= 3 values; else station median; else 1"""
    import s11_r3_q3 as RQ
    return RQ.inflate(cells, ctrl)


def b2_cells(far=False, ccmin=R.CCMIN, ctrl="g"):
    """B2 graph cells (nodes >= 50 events, b1_cc >= ccmin) of the region, near or far; columns d, sig (inflated), cc, site,
    far, dc/sc (control value and its sigma of the chosen control kind: g = B2 control nodes, b1 = reference-path control)"""
    E = r4_episode(far)
    if not len(E):
        return pd.DataFrame()
    E = E.copy()
    ce, cs_ = {"g": ("g_c2_eps", "g_c2_sig"), "b1": ("b1_c2_eps", "b1_c2_sig")}[ctrl]
    E["dc"] = -100 * E[ce] if ce in E else np.nan
    E["sc"] = 100 * E[cs_] if cs_ in E else np.nan
    if far:
        ser = pd.read_csv(f"{DIST}/series.csv"); ok = set(zip(ser.fam[ser.causal25 > 1.5], ser.netsta[ser.causal25 > 1.5], ser.epoch[ser.causal25 > 1.5]))
        E = E[[(f, n, e) in ok for f, n, e in zip(E.fam, E.netsta, E.epoch)]]
    q = E[(E.n >= NODE_MIN) & (E.b1_cc >= ccmin)].copy()
    q["d"] = -100 * q.g_eps; q["sig"] = 100 * q.g_sig; q["cc"] = q.b1_cc
    if far:
        sm = far_site_map(set(E.netsta)); q["site"] = q.netsta.map(sm)
    else:
        q["site"] = q.netsta.map(C.SITE)
    q = q[np.isfinite(q.d) & np.isfinite(q.sig) & (q.sig > 0) & q.site.notna()].reset_index(drop=True)
    E2 = E.copy()
    q["sig_raw"] = q.sig; q["infl"] = inflate(q, E2); q["sig"] = q.sig_raw * q.infl
    key = q.fam.astype(str) + "|" + q.netsta + "|" + q.epoch
    q = q[key.map(key.value_counts()) >= 3].reset_index(drop=True)
    q["series"] = q.fam.astype(str) + "|" + q.netsta + "|" + q.epoch
    q[["a_s", "a_m", "a_d", "a_b"]] = coefs(q, "ell40")
    q = q[np.isfinite(q.a_d) & np.isfinite(q.a_s)].reset_index(drop=True)
    q["far"] = far
    q["mid"] = q.episode.map(U.mid)
    return q


def b0_far_cells():
    import s11_dist_model as Mo
    f = Mo.far_cells()
    if not len(f):
        return f
    f["series"] = f.fam.astype(str) + "|" + f.netsta + "|" + f.epoch
    f[["a_s", "a_m", "a_d", "a_b"]] = coefs(f, "ell40")
    f["mid"] = f.episode.map(U.mid)
    return f


def control_cells(cells):
    """control version of a cell set (same cells, value dc, sigma sc x the same inflation)"""
    c = cells[np.isfinite(cells.dc) & np.isfinite(cells.sc) & (cells.sc > 0)].copy()
    c["d"] = c.dc; c["sig"] = c.sc * c.infl
    key = c.fam.astype(str) + "|" + c.netsta + "|" + c.epoch
    return c[key.map(key.value_counts()) >= 3].reset_index(drop=True)


def identifiable(cells):
    """keep cells of episodes with cells at >= 2 sites (PLAN §8.1)"""
    ns = cells.groupby("episode").site.nunique()
    ok = cells.episode.map(ns) >= 2
    return cells[ok.to_numpy()].reset_index(drop=True), int((~ok).sum())


# ------------------------------------------------------------------------------ combined model
def _demean_G(cells, G, w):
    key = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    codes, uniq = pd.factorize(key)
    W = np.bincount(codes, w)
    Gm = np.zeros((len(uniq), G.shape[1])); np.add.at(Gm, codes, w[:, None] * G)
    return G - (Gm / W[:, None])[codes], codes, W, len(uniq)


class Comb:
    """d = o_series + a_d D(e) [+ a_s Sc(e)] + a_s S_site(t) [+ R_sta(t)] [+ a_m M(t) + a_b B(t)];
    constraint=True: no prior on D, Sc and exact sum-to-zero; else prior prior_ep on them. apparent=True: every kernel
    coefficient = 1 (kernel-free version). Q3 functionals use names ("D", "", e) as s11_r3_fit."""
    def __init__(self, cells, Sc=True, Rr=True, MB=True, prior=1.0, constraint=True, prior_ep=1.0, apparent=False, rmin=RMIN, Ssite=True):
        c = cells.reset_index(drop=True); self.cells = c
        mid = c.episode.map(U.mid); t, doy = V.tdec(mid)
        bS = V.basis(t, doy, STERMS); bR = V.basis(t, doy, RT); bRR = V.basis(t, doy, RTERMS)
        eps_ = sorted(c.episode.unique(), key=lambda e: U.mid[e])
        sites = sorted(c.site.unique())
        cnt = c.netsta.value_counts()
        stas = sorted([s for s in cnt.index if cnt[s] >= rmin]) if Rr else []
        names = [("D", "", e) for e in eps_]
        if Sc:
            names += [("Sc", "", e) for e in eps_]
        if not Ssite:
            sites = []
        names += [("S", s, tm) for s in sites for tm in STERMS]
        names += [("R", s, tm) for s in stas for tm in RTERMS]
        if MB:
            names += [("M", "", tm) for tm in RT] + [("B", "", tm) for tm in RT]
        idx = {n: i for i, n in enumerate(names)}
        n = len(c); G = np.zeros((n, len(names)))
        one = np.ones(n)
        ad, as_, am, ab = (one, one, one, one) if apparent else (c.a_d.to_numpy(), c.a_s.to_numpy(), c.a_m.to_numpy(), c.a_b.to_numpy())
        r = np.arange(n)
        G[r, [idx[("D", "", e)] for e in c.episode]] = ad
        if Sc:
            G[r, [idx[("Sc", "", e)] for e in c.episode]] = as_
        site = c.site.to_numpy(); ns = c.netsta.to_numpy()
        for s in sites:
            m = site == s
            for tm in STERMS:
                G[m, idx[("S", s, tm)]] = as_[m] * bS[tm][m]
        for s in stas:
            m = ns == s
            for tm in RTERMS:
                G[m, idx[("R", s, tm)]] = bRR[tm][m]
        if MB:
            for tm in RT:
                G[:, idx[("M", "", tm)]] = am * bR[tm]; G[:, idx[("B", "", tm)]] = ab * bR[tm]
        isep = np.array([k[0] in ("D", "Sc") for k in names])
        prec = np.where(isep, 0.0 if constraint else 1.0 / prior_ep ** 2, 1.0 / prior ** 2)
        self.w = 1.0 / c.sig.to_numpy() ** 2
        self.Gt, self.codes, self.Wg, self.nser = _demean_G(c, G, self.w)
        self.N = self.Gt.T @ (self.Gt * self.w[:, None])
        M = self.N + np.diag(prec)
        self.singular = False
        if constraint:
            A = np.zeros((1 + int(Sc), len(names)))
            A[0, [idx[("D", "", e)] for e in eps_]] = 1.0
            if Sc:
                A[1, [idx[("Sc", "", e)] for e in eps_]] = 1.0
            from scipy.linalg import null_space
            Z = null_space(A)
            Mz = Z.T @ M @ Z
            cond = np.linalg.cond(Mz)
            self.cond = float(cond)
            if not np.isfinite(cond) or cond > 1e12:
                self.singular = True; Mzi = np.linalg.pinv(Mz, hermitian=True)
            else:
                Mzi = np.linalg.inv(Mz)
            self.Ci = Z @ Mzi @ Z.T
        else:
            self.cond = float(np.linalg.cond(M)); self.Ci = np.linalg.inv(M)
        self.Op = self.Ci @ (self.Gt.T * self.w[None, :])
        self.peff = float(np.trace(self.Ci @ self.N))
        self.names, self.idx, self.eps = names, idx, eps_
        self.extra = dict(eps=eps_, U=U, ncell=c.groupby("episode").size(), zs=[""])
        self.stas_R = stas

    def demean_d(self, d):
        dm = np.bincount(self.codes, self.w * d) / self.Wg
        return d - dm[self.codes]

    def solve(self, d=None):
        d = self.cells.d.to_numpy() if d is None else d
        th = self.Op @ d
        r = self.demean_d(d) - self.Gt @ th
        chi2 = float(np.sum(self.w * r ** 2) / max(len(d) - self.nser - self.peff, 1))
        F = dict(th=th, C=self.Ci * max(chi2, 1.0), idx=self.idx, names=self.names, chi2=chi2, n=len(d), nser=self.nser, resid=r)
        F.update(self.extra)
        return F

    def delta(self, x):
        F = dict(th=self.Op @ x, C=self.Ci, idx=self.idx, names=self.names); F.update(self.extra); return F


def q3(F, majors_):
    import s11_r3_fit as LF
    return LF.q3_table(F, majors_)


def summary(F, T, comp, label=""):
    if not len(T):
        T = pd.DataFrame(columns=["episode", "major", "dD", "sig"]).astype({"major": bool, "dD": float, "sig": float})
    mj = T[T.major]
    return dict(label=label, n_cells=F["n"], n_series=F["nser"], n_episodes=len(F["eps"]), chi2=F["chi2"],
                n_majors=int(len(mj)), sig_major_median=float(mj.sig.median()) if len(mj) else np.nan,
                n_major_sig_le_0p1=int((mj.sig <= 0.10).sum()) if len(mj) else 0,
                composite=comp["composite"] if comp else np.nan, composite_sig=comp["composite_sig"] if comp else np.nan,
                chi2_majors=float(np.mean((mj.dD / mj.sig) ** 2)) if len(mj) else np.nan,
                n_beyond_2sig=int((np.abs(mj.dD / mj.sig) > 2).sum()) if len(mj) else 0,
                n_beyond_3sig=int((np.abs(mj.dD / mj.sig) >= 3).sum()) if len(mj) else 0)


def fmt(r):
    return ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items())


def pt_now():
    return C.pt_now()
