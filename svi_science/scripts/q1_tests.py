"""Q1 tests: slip budget vs front kinematics for per-patch residual shares. Rules: results/q1_patch_budget/PLAN.md.

Core object: a 'design' = one catalogue/variant/patch definition, an ordered list of history slots (episodes) and a set
of target slots (big episodes). Count vectors (with their event-time offsets, tremor extents and directions) can be
permuted among slots of the same (era, target) stratum: the order-shuffle control.
"""
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q1_patch_budget"
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
NCAL = int(sys.argv[2]) if len(sys.argv) > 2 else 200
GAP = 90.0
rng = np.random.default_rng(20261007)
EPS = 0.5


# ----------------------------------------------------------------------------------------------- data
def load():
    CL = pd.read_parquet(f"{OUT}/counts_long.parquet")
    PA = pd.read_csv(f"{OUT}/patches.csv", dtype={"patch": str})
    EP = pd.read_csv(f"{OUT}/episodes.csv", parse_dates=["w0", "w1"])
    TR = pd.read_csv(f"{OUT}/tremor_ep.csv")
    DT = pd.read_csv(f"{OUT}/det_patch_ep.csv", dtype={"patch": str})
    S3 = pd.read_csv(f"{ROOT}/svi_science/results/s3/episodes_all.csv")
    KIN = pd.read_csv(f"{ROOT}/svi_science/results/s3/episodes_kinematics.csv")
    GEO = pd.read_csv(f"{ROOT}/svi_science/results/s1/episode_geodetic.csv")
    return CL, PA, EP, TR, DT, S3, KIN, GEO


class Design:
    def __init__(self, D, cat, variant, patchdef, hist_eps, target_eps, extra_cat=None):
        """hist_eps: DataFrame rows of episodes (catalogue, episode, w0, w1, era, major) in time order.
        target_eps: list of (catalogue, episode) keys that are targets."""
        CL, PA = D["CL"], D["PA"]
        H = hist_eps.sort_values("w0").reset_index(drop=True)
        self.H = H
        pa = PA[(PA.catalogue == (cat if cat != "JOINT" else "SVI")) & (PA.patchdef == patchdef)].copy()
        self.patches = pa.patch.tolist(); self.pos = pa.as_km.to_numpy(); P = len(self.patches)
        self.pa = pa
        A = len(H)
        C = np.zeros((P, A)); O = np.full((P, A), np.nan)
        for a, r in H.iterrows():
            v = variant[r.catalogue] if isinstance(variant, dict) else variant
            q = CL[(CL.catalogue == r.catalogue) & (CL.variant == v) & (CL.patchdef == patchdef) & (CL.episode == r.episode)]
            q = q.set_index("patch").reindex(self.patches)
            C[:, a] = q.n.fillna(0).to_numpy()
            O[:, a] = ((pd.to_datetime(q.t_med) - r.w0) / pd.Timedelta(days=1)).to_numpy()
        self.C, self.O = C, O
        self.N = C.sum(0)
        self.w0 = ((H.w0 - pd.Timestamp("2000-01-01")) / pd.Timedelta(days=1)).to_numpy()
        self.w1 = ((H.w1 - pd.Timestamp("2000-01-01")) / pd.Timedelta(days=1)).to_numpy()
        self.mid = (self.w0 + self.w1) / 2
        self.era = H.era.to_numpy()
        keys = list(zip(H.catalogue, H.episode))
        tset = set(target_eps)
        self.is_t = np.array([k in tset for k in keys])
        self.tslots = np.where(self.is_t)[0]
        self.keys = keys
        # era sums (all episodes of the era, as long-term reference for the activity threshold)
        self.eras = sorted(set(self.era))
        self.Sera = {e: C[:, self.era == e].sum(1) for e in self.eras}
        # residuals: per target vector, LOO within era among targets; pair-out matrices
        T = self.tslots; nt = len(T)
        self.rL = np.zeros((P, nt)); self.pL = np.zeros((P, nt))
        self.rP = np.zeros((P, nt, nt)); self.pP = np.zeros((P, nt, nt))
        for i in range(nt):
            same = [k for k in range(nt) if self.era[T[k]] == self.era[T[i]]]
            for j in range(nt):
                excl = {i, j}
                oth = [k for k in same if k not in excl]
                s = C[:, T[oth]].sum(1); s = np.clip(s / s.sum(), 1e-4, None); s /= s.sum()
                E = self.N[T[i]] * s
                rr = np.log((C[:, T[i]] + EPS) / (E + EPS)); pr = (C[:, T[i]] - E) / np.sqrt(E)
                self.rP[:, i, j] = rr; self.pP[:, i, j] = pr
                if j == i:
                    self.rL[:, i] = rr; self.pL[:, i] = pr
            # if other era: pair-out == LOO (other episode not in this era's pool)
            for j in range(nt):
                if self.era[T[j]] != self.era[T[i]]:
                    self.rP[:, i, j] = self.rL[:, i]; self.pP[:, i, j] = self.pL[:, i]
        self.D = D; self.cat = cat

    # ------------------------------------------------------------------ permutations
    def identity(self):
        return np.arange(len(self.H))

    def order_shuffle(self):
        perm = np.arange(len(self.H))
        for e in self.eras:
            for t in (True, False):
                idx = np.where((self.era == e) & (self.is_t == t))[0]
                perm[idx] = rng.permutation(idx)
        return perm

    # ------------------------------------------------------------------ predictors for slot assignment perm (perm[slot] = vector)
    def predictors(self, perm, thr=(3, 0.25), patch_perm=None):
        """Returns DataFrame rows (patch, target-slot) with r, pr, T1, censored, T2 (+ its responses), T3."""
        C = self.C[:, perm] if patch_perm is None else np.take_along_axis(self.C, patch_perm, 1)
        O = self.O[:, perm] if patch_perm is None else np.take_along_axis(self.O, patch_perm, 1)
        N = C.sum(0)
        P = C.shape[0]
        tpos = {s: i for i, s in enumerate(self.tslots)}       # vector index in target arrays
        vt = perm[self.tslots] if patch_perm is None else None
        rows = []
        tslots = self.tslots
        for ii, j in enumerate(tslots):
            v = tpos[perm[j]] if patch_perm is None else ii
            ej = self.era[j]
            # long-term shares for threshold, excluding the target vector
            act = np.zeros((P, j), bool)
            for e in self.eras:
                cols = np.where(self.era[:j] == e)[0]
                if not len(cols):
                    continue
                S = self.Sera[e] - (C[:, j] if e == ej else 0)
                s = S / max(S.sum(), 1)
                act[:, cols] = (C[:, cols] >= thr[0]) & (C[:, cols] >= thr[1] * N[cols][None, :] * s[:, None])
            T1 = np.full(P, np.nan); cens = np.ones(P, bool)
            for p in range(P):
                a = np.where(act[p])[0]
                if len(a):
                    a = a[-1]; tl = self.w0[a] + (O[p, a] if np.isfinite(O[p, a]) else (self.w1[a] - self.w0[a]) / 2)
                    T1[p] = max(self.w0[j] - tl, 1.0); cens[p] = False
                else:
                    T1[p] = max(self.w0[j] - self.w0[0], 1.0)
            # earlier targets >= 90 d
            prev = [k for k in tslots if k < j and (self.mid[j] - self.mid[k]) >= GAP]
            if patch_perm is None:
                rL = self.rL[:, v]; pL = self.pL[:, v]
                if prev:
                    kk = prev[-1]; w = tpos[perm[kk]]
                    T2 = self.rP[:, w, v]; r2 = self.rP[:, v, w]
                    T3 = sum(self.rP[:, tpos[perm[k]], v] for k in prev)
                else:
                    T2 = r2 = T3 = np.full(P, np.nan)
            else:
                rL, pL, T2, r2, T3 = self._pp_resid(C, ii, prev)
            for p in range(P):
                rows.append((p, ii, j, v, rL[p], pL[p], T1[p], cens[p], T2[p], r2[p], T3[p]))
        return pd.DataFrame(rows, columns=["p", "ti", "slot", "v", "r", "pr", "T1", "cens", "T2", "r2", "T3"])

    def _pp_resid(self, C, ii, prev):
        """residuals recomputed from a per-patch-shuffled matrix (sum-to-one broken; caller's control)."""
        T = self.tslots; j = T[ii]
        same = [k for k in T if self.era[k] == self.era[j]]
        N = C.sum(0)

        def res(i, excl):
            oth = [k for k in same if k not in excl]
            s = C[:, oth].sum(1); s = np.clip(s / s.sum(), 1e-4, None); s /= s.sum()
            E = N[i] * s
            return np.log((C[:, i] + EPS) / (E + EPS)), (C[:, i] - E) / np.sqrt(E)
        rL, pL = res(j, {j})
        if prev:
            kk = prev[-1]
            T2 = res(kk, {kk, j})[0]; r2 = res(j, {kk, j})[0]
            T3 = sum(res(k, {k, j})[0] for k in prev)
        else:
            nan = np.full(C.shape[0], np.nan); T2 = r2 = T3 = nan
        return rL, pL, T2, r2, T3


# ----------------------------------------------------------------------------------------------- regression helpers
def fe_slope(r, x, g):
    m = np.isfinite(r) & np.isfinite(x)
    r, x, g = r[m], x[m], g[m]
    if len(r) < 5 or np.std(x) == 0:
        return np.nan, np.nan, 0
    x = (x - x.mean()) / x.std()
    df = pd.DataFrame(dict(r=r, x=x, g=g))
    rt = r - df.groupby("g").r.transform("mean").to_numpy(); xt = x - df.groupby("g").x.transform("mean").to_numpy()
    sxx = (xt ** 2).sum()
    if sxx <= 1e-12:
        return np.nan, np.nan, len(r)
    b = (xt * rt).sum() / sxx
    r2 = b ** 2 * sxx / (rt ** 2).sum()
    return float(b), float(r2), int(len(r))


def fe_multi_r2(r, X, g):
    m = np.isfinite(r) & np.all(np.isfinite(X), 1)
    r, X, g = r[m], X[m], g[m]
    df = pd.DataFrame(X); df["g"] = g; df["r"] = r
    Xt = (df.drop(columns=["g", "r"]) - df.groupby("g")[list(range(X.shape[1]))].transform("mean")).to_numpy()
    rt = r - df.groupby("g").r.transform("mean").to_numpy()
    keep = Xt.std(0) > 1e-9
    Xt = Xt[:, keep]
    if Xt.shape[1] == 0:
        return 0.0
    beta, *_ = np.linalg.lstsq(Xt, rt, rcond=None)
    return float(1 - ((rt - Xt @ beta) ** 2).sum() / (rt ** 2).sum())


def budget_stats(df, resp="r"):
    out = {}
    g = df.ti.to_numpy()
    nc = df[~df.cens]
    out["T1"] = fe_slope(nc[resp].to_numpy(), np.log10(nc.T1.to_numpy()), nc.ti.to_numpy())
    out["T1_cens_kept"] = fe_slope(df[resp].to_numpy(), np.log10(df.T1.to_numpy()), g)
    out["T2"] = fe_slope(df.r2.to_numpy() if resp == "r" else df[resp].to_numpy(), df.T2.to_numpy(), g)
    out["T3"] = fe_slope(df[resp].to_numpy(), df.T3.to_numpy(), g)
    return out


def ptest(obs, null, side):
    null = np.asarray([x for x in null if np.isfinite(x)])
    if not np.isfinite(obs) or not len(null):
        return np.nan, np.nan
    one = ((null >= obs).sum() + 1) / (len(null) + 1) if side > 0 else ((null <= obs).sum() + 1) / (len(null) + 1)
    two = ((np.abs(null - np.median(null)) >= abs(obs - np.median(null))).sum() + 1) / (len(null) + 1)
    return float(one), float(two)


SIDES = {"T1": +1, "T1_cens_kept": +1, "T2": -1, "T3": -1}


def run_budget(des, thr=(3, 0.25), resp="r", nperm=NPERM, per_patch=False, drop_close=False):
    def stat(perm=None, pp=None):
        df = des.predictors(des.identity() if perm is None else perm, thr=thr, patch_perm=pp)
        if drop_close:
            close = [ii for ii, j in enumerate(des.tslots) if any((des.mid[j] - des.mid[k]) < GAP for k in des.tslots if k < j)]
            df = df[~df.ti.isin(close)]
        return budget_stats(df, resp), df
    obs, df0 = stat()
    nulls = {k: [] for k in obs}
    for _ in range(nperm):
        s, _ = stat(perm=des.order_shuffle())
        for k in s:
            nulls[k].append(s[k][0])
    res = {}
    for k, (b, r2, n) in obs.items():
        p1, p2 = ptest(b, nulls[k], SIDES[k])
        res[k] = dict(beta=b, r2=r2, n=n, p_order=p1, p_order_2s=p2, null_mean=float(np.nanmean(nulls[k])),
                      null_q01=float(np.nanquantile(nulls[k], 0.01)), null_q99=float(np.nanquantile(nulls[k], 0.99)))
    if per_patch:
        pn = {k: [] for k in obs}
        for _ in range(min(nperm, 1000)):
            pp = np.tile(np.arange(len(des.H)), (len(des.patches), 1))
            for e in des.eras:
                for t in (True, False):
                    idx = np.where((des.era == e) & (des.is_t == t))[0]
                    for p in range(len(des.patches)):
                        pp[p, idx] = rng.permutation(idx)
            s, _ = stat(pp=pp)
            for k in s:
                pn[k].append(s[k][0])
        for k in res:
            res[k]["p_perpatch"] = ptest(res[k]["beta"], pn[k], SIDES[k])[0]
    return res, df0


# ----------------------------------------------------------------------------------------------- designs
def hist_table(EP, cat, eras=None):
    h = EP[EP.catalogue == cat].copy()
    if eras is not None:
        h = h[h.era.isin(eras)]
    return h


def make(D, cat, variant, patchdef, targets="big", hist_eras=None, add_bostock=False):
    EP, S3 = D["EP"], D["S3"]
    if cat in ("SVI", "NW"):
        reg = cat
        if targets == "big":
            tl = S3[(S3.region == reg) & S3.big].episode.tolist()
        else:
            tl = S3[(S3.region == reg) & (S3.n_headline >= 100)].episode.tolist()
        H = hist_table(EP, cat, hist_eras)
        tk = [(cat, e) for e in tl]
        if add_bostock:
            Hb = hist_table(EP, "BOSTOCK")
            H = pd.concat([Hb, H], ignore_index=True)
            variant = {cat: variant, "BOSTOCK": "all"}
        return Design(D, cat, variant, patchdef, H, tk)
    if cat == "BOSTOCK":
        H = hist_table(EP, "BOSTOCK")
        return Design(D, cat, variant, "group", H, list(zip(H.catalogue, H.episode)))
    if cat == "JOINT":       # Bostock 10 + our SVI big 10, big-only history
        Hb = hist_table(EP, "BOSTOCK")
        big = S3[(S3.region == "SVI") & S3.big].episode.tolist()
        Hs = EP[(EP.catalogue == "SVI") & EP.episode.isin(big)]
        H = pd.concat([Hb, Hs], ignore_index=True)
        return Design(D, "JOINT", {"BOSTOCK": "all", "SVI": variant}, "group", H, list(zip(H.catalogue, H.episode)))


# ----------------------------------------------------------------------------------------------- kinematics, detectability, episode level
def covariates(des, D, region, dir_source="tremor"):
    TR, KIN, DT = D["TR"], D["KIN"], D["DT"]
    T = des.tslots; nt = len(T)
    eps = [des.keys[j][1] for j in T]
    tr = TR[TR.region == region].set_index("episode").reindex(eps)
    p05, p95 = tr.trem_p05.to_numpy(), tr.trem_p95.to_numpy()
    if dir_source == "tremor":
        dr = tr.trem_dir.fillna(0).to_numpy()
    else:
        k = KIN[KIN.region == region].set_index("episode").reindex(eps)
        dr = np.where(k.direction == "NW-ward", 1, np.where(k.direction == "SE-ward", -1, 0))
    x = des.pos; xmin, xmax = x.min(), x.max()
    det = DT[(DT.catalogue == region) & (DT.patchdef == des.pa.patchdef.iloc[0])].copy()
    det["d"] = np.log(det.mean_sta * det.frac_obs.clip(lower=1e-3))
    dm = det.pivot_table(index="patch", columns="episode", values="d").reindex(index=des.patches, columns=eps).to_numpy()
    dm = dm - np.nanmean(dm, 1, keepdims=True)
    return dict(p05=p05, p95=p95, dir=dr, det=dm, eps=eps, xmin=xmin, xmax=xmax)


def kin_rows(des, cov, perm_ext=None, perm_dir=None):
    """per (patch, target) kinematic predictors; perm_* permute episode attributes among targets."""
    P = len(des.patches); nt = len(des.tslots); x = des.pos
    pe = np.arange(nt) if perm_ext is None else perm_ext
    pd_ = np.arange(nt) if perm_dir is None else perm_dir
    lo = cov["p05"][pe] - 5; hi = cov["p95"][pe] + 5
    inr = ((x[:, None] >= lo[None]) & (x[:, None] <= hi[None])).astype(float)
    kout = np.maximum(0, np.maximum(lo[None] - x[:, None], x[:, None] - hi[None]))
    dr = cov["dir"][pd_]
    u = np.where(dr[None] > 0, (x[:, None] - cov["xmin"]) / (cov["xmax"] - cov["xmin"]),
                 np.where(dr[None] < 0, (cov["xmax"] - x[:, None]) / (cov["xmax"] - cov["xmin"]), np.nan))
    # direction mode: mean pair-out residual over other targets with same direction, >= 90 d apart
    dm = np.full((P, nt), np.nan)
    for i in range(nt):
        if dr[i] == 0:
            continue
        ks = [k for k in range(nt) if k != i and dr[k] == dr[i] and abs(des.mid[des.tslots[k]] - des.mid[des.tslots[i]]) >= GAP]
        if ks:
            dm[:, i] = np.mean([des.rP[:, k, i] for k in ks], 0)
    return inr, kout, u, dm


def flat(M):
    return M.reshape(-1)        # patch-major, matches rows ordering p*nt + ti? (we build r the same way)


def kin_tests(des, cov, nperm=NPERM):
    P = len(des.patches); nt = len(des.tslots)
    r = des.rL.reshape(-1); g = np.tile(np.arange(nt), P)
    inr, kout, u, dm = kin_rows(des, cov)
    out = {}
    obs = dict(K1_inreach=fe_slope(r, inr.reshape(-1), g), K1_kmout=fe_slope(r, kout.reshape(-1), g),
               K2_dirmode=fe_slope(r, dm.reshape(-1), g), K3_u=fe_slope(r, u.reshape(-1), g),
               D1_det=fe_slope(r, cov["det"].reshape(-1), g))
    nul = {k: [] for k in obs}
    for _ in range(nperm):
        pe = rng.permutation(nt); pdir = rng.permutation(nt)
        a, bkm, _, _ = kin_rows(des, cov, perm_ext=pe)
        _, _, uu, dd = kin_rows(des, cov, perm_dir=pdir)
        nul["K1_inreach"].append(fe_slope(r, a.reshape(-1), g)[0]); nul["K1_kmout"].append(fe_slope(r, bkm.reshape(-1), g)[0])
        nul["K2_dirmode"].append(fe_slope(r, dd.reshape(-1), g)[0]); nul["K3_u"].append(fe_slope(r, uu.reshape(-1), g)[0])
        dp = cov["det"].copy()
        for p in range(P):
            ok = np.isfinite(dp[p]); dp[p, ok] = rng.permutation(dp[p, ok])
        nul["D1_det"].append(fe_slope(r, dp.reshape(-1), g)[0])
    side = dict(K1_inreach=1, K1_kmout=-1, K2_dirmode=1, K3_u=0, D1_det=1)
    for k, (b, r2, n) in obs.items():
        p1, p2 = ptest(b, nul[k], side[k] if side[k] else 1)
        out[k] = dict(beta=b, r2=r2, n=n, p=(p2 if side[k] == 0 else p1), p_2s=p2)
    return out, dict(inr=inr, kout=kout, u=u, dm=dm)


def h4_and_l2(des, cov, K, nperm=NPERM):
    """H4: R2 gain of budget {T1,T2,T3} over kinematics {inreach, kmout, dirmode, u}; L2: inreach ~ log T1."""
    P = len(des.patches); nt = len(des.tslots)

    def build(df):
        df = df.sort_values(["p", "ti"])
        r = df.r.to_numpy(); g = df.ti.to_numpy()
        pi, vi = df.p.to_numpy(), df.v.to_numpy()      # kinematics travel with the count vector
        Kx = np.column_stack([K["inr"][pi, vi], K["kout"][pi, vi], np.nan_to_num(K["dm"][pi, vi]), np.nan_to_num(K["u"][pi, vi], nan=0.5)])
        T1 = np.where(df.cens, np.nan, np.log10(df.T1.to_numpy()))
        Bx = np.column_stack([T1, np.nan_to_num(df.T2.to_numpy()), np.nan_to_num(df.T3.to_numpy())])
        full = np.column_stack([Kx, Bx])
        r2k = fe_multi_r2(np.where(np.isfinite(T1), r, np.nan), Kx, g)
        r2f = fe_multi_r2(r, full, g)
        r2b = fe_multi_r2(r, Bx, g)
        l2 = fe_slope(K["inr"][pi, vi], T1, g)
        return r2f - r2k, r2k, r2f, r2b, l2
    obs = build(des.predictors(des.identity()))
    nd, nl = [], []
    for _ in range(nperm):
        o = build(des.predictors(des.order_shuffle()))
        nd.append(o[0]); nl.append(o[4][0])
    return dict(H4_dR2=obs[0], R2_kin=obs[1], R2_full=obs[2], R2_budget=obs[3], p_H4=ptest(obs[0], nd, 1)[0],
                L2_beta=obs[4][0], L2_n=obs[4][2], p_L2=ptest(obs[4][0], nl, 1)[0], p_L2_2s=ptest(obs[4][0], nl, 1)[1])


def episode_level(des, cov, D, region):
    T = des.tslots; eps = cov["eps"]
    P = len(des.patches)
    U = (des.pL ** 2).sum(0) / (P - 1)
    tr = D["TR"][D["TR"].region == region].set_index("episode").reindex(eps)
    det = D["DT"][(D["DT"].catalogue == region) & (D["DT"].patchdef == des.pa.patchdef.iloc[0])]
    msta = det.groupby("episode").mean_sta.mean().reindex(eps).to_numpy()
    X = dict(logN=np.log10(des.N[T]), abs_trem_speed=np.abs(tr.trem_speed.to_numpy()),
             trem_width=(tr.trem_p95 - tr.trem_p05).to_numpy(), mean_sta=msta)
    if region == "SVI":
        geo = D["GEO"].set_index("episode").reindex(eps)
        X["gnss_D_mm"] = geo.D.to_numpy()
    out = {}
    for k, x in X.items():
        m = np.isfinite(x) & np.isfinite(U)
        if m.sum() < 5:
            continue
        rho = stats.spearmanr(U[m], x[m])[0]
        nul = [stats.spearmanr(U[m], rng.permutation(x[m]))[0] for _ in range(5000)]
        out[k] = dict(rho=float(rho), n=int(m.sum()), p=float((np.sum(np.abs(nul) >= abs(rho)) + 1) / 5001))
    return out, dict(eps=eps, U=U.tolist(), **{k: v.tolist() for k, v in X.items()})


# ----------------------------------------------------------------------------------------------- calibration
def calibrate(des, nsim=NCAL, nperm=200):
    """no-memory multinomial episodes (era pooled shares, real N and dates) -> false-positive rate of order-shuffle test."""
    C0 = des.C.copy(); hits = {k: 0 for k in SIDES}; betas = {k: [] for k in SIDES}
    for _ in range(nsim):
        Cs = np.zeros_like(C0)
        for e in des.eras:
            cols = np.where(des.era == e)[0]; s = C0[:, cols].sum(1); s = s / s.sum()
            for a in cols:
                Cs[:, a] = rng.multinomial(int(C0[:, a].sum()), s)
        d2 = object.__new__(Design); d2.__dict__.update(des.__dict__)
        d2.C = Cs; d2.O = np.where(Cs > 0, np.nan_to_num(des.O, nan=5.0), np.nan)
        d2.Sera = {e: Cs[:, d2.era == e].sum(1) for e in d2.eras}
        # recompute residual arrays
        tmp = Design.__new__(Design); tmp.__dict__.update(d2.__dict__)
        _recompute_resid(tmp)
        res, _ = run_budget(tmp, nperm=nperm)
        for k in SIDES:
            betas[k].append(res[k]["beta"])
            hits[k] += int(np.isfinite(res[k]["p_order"]) and res[k]["p_order"] < 0.01)
    return {k: dict(fpr_p01=hits[k] / nsim, beta_mean=float(np.nanmean(betas[k])), beta_sd=float(np.nanstd(betas[k]))) for k in SIDES}


def _recompute_resid(self):
    C = self.C; T = self.tslots; nt = len(T); P = C.shape[0]
    self.N = C.sum(0)
    self.rL = np.zeros((P, nt)); self.pL = np.zeros((P, nt)); self.rP = np.zeros((P, nt, nt)); self.pP = np.zeros((P, nt, nt))
    for i in range(nt):
        same = [k for k in range(nt) if self.era[T[k]] == self.era[T[i]]]
        for j in range(nt):
            oth = [k for k in same if k not in {i, j}]
            s = C[:, T[oth]].sum(1); s = np.clip(s / s.sum(), 1e-4, None); s /= s.sum()
            E = self.N[T[i]] * s
            self.rP[:, i, j] = np.log((C[:, T[i]] + EPS) / (E + EPS)); self.pP[:, i, j] = (C[:, T[i]] - E) / np.sqrt(E)
            if j == i:
                self.rL[:, i] = self.rP[:, i, j]; self.pL[:, i] = self.pP[:, i, j]
        for j in range(nt):
            if self.era[T[j]] != self.era[T[i]]:
                self.rP[:, i, j] = self.rL[:, i]; self.pP[:, i, j] = self.pL[:, i]


# ----------------------------------------------------------------------------------------------- main
def main():
    CL, PA, EP, TR, DT, S3, KIN, GEO = load()
    D = dict(CL=CL, PA=PA, EP=EP, TR=TR, DT=DT, S3=S3, KIN=KIN, GEO=GEO)
    results, table = {}, []

    def rec(name, region, res, extra=None):
        results[name] = res if extra is None else {**res, **extra}
        for k, v in res.items():
            if isinstance(v, dict) and "beta" in v:
                table.append(dict(design=name, region=region, predictor=k, **{kk: vv for kk, vv in v.items()}))

    # scrambled-control counts per target (SVI)
    sc = {}
    for v in ["c1_once075", "c1_once050"]:
        big = S3[(S3.region == "SVI") & S3.big].episode.tolist()
        q = CL[(CL.catalogue == "SVI") & (CL.variant == v) & (CL.patchdef == "group") & CL.episode.isin(big)].groupby("episode").n.sum()
        sc[v] = q.reindex(big).fillna(0).astype(int).to_dict()
    results["scrambled_counts_big"] = sc
    print("scrambled", sc, flush=True)

    specs = [
        # name, cat, variant, patchdef, targets, hist_eras, add_bostock, kwargs
        ("SVI_primary", "SVI", "once075", "group", "big", None, False, dict(per_patch=True)),
        ("NW_primary", "NW", "once075", "bin10", "big", ["E3"], False, dict(per_patch=True)),
        ("SVI_S13fixed", "SVI", "S13_once075", "group", "big", None, False, {}),
        ("NW_F9fixed_raw", "NW", "F9_raw075", "bin10", "big", ["E3"], False, {}),
        ("SVI_raw075", "SVI", "raw075", "group", "big", None, False, {}),
        ("NW_raw075", "NW", "raw075", "bin10", "big", ["E3"], False, {}),
        ("SVI_once050", "SVI", "once050", "group", "big", None, False, {}),
        ("NW_once050", "NW", "once050", "bin10", "big", ["E3"], False, {}),
        ("SVI_bin10", "SVI", "once075", "bin10", "big", None, False, {}),
        ("NW_bin10s5", "NW", "once075", "bin10s5", "big", ["E3"], False, {}),
        ("SVI_secondary", "SVI", "once075", "group", "ge100", None, False, {}),
        ("NW_secondary", "NW", "once075", "bin10", "ge100", ["E3"], False, {}),
        ("SVI_hist_bostock", "SVI", "once075", "group", "big", None, True, {}),
        ("NW_hist_E2E3", "NW", "once075", "bin10", "big", ["E2", "E3"], False, {}),
        ("NW_hist_E1E2E3", "NW", "once075", "bin10", "big", ["E1", "E2", "E3"], False, {}),
        ("NW_drop_close", "NW", "once075", "bin10", "big", ["E3"], False, dict(drop_close=True)),
        ("SVI_pearson", "SVI", "once075", "group", "big", None, False, dict(resp="pr")),
        ("NW_pearson", "NW", "once075", "bin10", "big", ["E3"], False, dict(resp="pr")),
        ("SVI_thr1", "SVI", "once075", "group", "big", None, False, dict(thr=(1, 0.0))),
        ("NW_thr1", "NW", "once075", "bin10", "big", ["E3"], False, dict(thr=(1, 0.0))),
        ("SVI_thr5", "SVI", "once075", "group", "big", None, False, dict(thr=(5, 0.5))),
        ("NW_thr5", "NW", "once075", "bin10", "big", ["E3"], False, dict(thr=(5, 0.5))),
    ]
    designs = {}
    for name, cat, var, pdf, tg, he, ab, kw in specs:
        des = make(D, cat, var, pdf, tg, he, ab)
        designs[name] = des
        res, df0 = run_budget(des, **kw)
        region = cat
        rec(name, region, res, dict(n_targets=int(len(des.tslots)), n_patches=len(des.patches), n_history=int(len(des.H)),
                                    n_censored=int(df0.cens.sum())))
        df0["patch"] = [des.patches[p] for p in df0.p]; df0["episode"] = [des.keys[j][1] for j in df0.slot]
        df0.to_csv(f"{OUT}/rows_{name}.csv", index=False)
        print(name, {k: (round(v["beta"], 3) if np.isfinite(v["beta"]) else None, v["p_order"], v.get("p_perpatch")) for k, v in res.items()}, flush=True)

    # Bostock alone and joint 2003-2026
    for name, cat, var in [("BOSTOCK_alone", "BOSTOCK", "all"), ("BOSTOCK_once2s", "BOSTOCK", "once2s"), ("JOINT_2003_2026", "JOINT", "once075")]:
        des = make(D, cat, var, "group")
        designs[name] = des
        res, df0 = run_budget(des)
        rec(name, "SVI", res, dict(n_targets=int(len(des.tslots)), n_patches=len(des.patches), n_censored=int(df0.cens.sum())))
        df0["patch"] = [des.patches[p] for p in df0.p]; df0["episode"] = [des.keys[j][1] for j in df0.slot]
        df0.to_csv(f"{OUT}/rows_{name}.csv", index=False)
        print(name, {k: (round(v["beta"], 3) if np.isfinite(v["beta"]) else None, v["p_order"]) for k, v in res.items()}, flush=True)

    # kinematics / detectability / H4 / L2 / episode level
    epl = {}
    for name, region in [("SVI_primary", "SVI"), ("NW_primary", "NW"), ("SVI_S13fixed", "SVI"), ("NW_F9fixed_raw", "NW")]:
        des = designs[name]
        for dsrc in (["tremor", "lfe"] if name.endswith("primary") else ["tremor"]):
            cov = covariates(des, D, region, dsrc)
            kt, K = kin_tests(des, cov)
            h = h4_and_l2(des, cov, K, nperm=NPERM if name.endswith("primary") else 500)
            key = f"{name}|kin_{dsrc}"
            rec(key, region, kt, dict(H4_L2=h))
            print(key, {k: (round(v["beta"], 3), v["p"]) for k, v in kt.items()}, h, flush=True)
            if dsrc == "tremor":
                np.savez(f"{OUT}/kin_{name}.npz", **K, det=cov["det"], r=des.rL, patches=np.array(des.patches),
                         eps=np.array(cov["eps"]), pos=des.pos)
        if name.endswith("primary"):
            e, ev = episode_level(des, covariates(des, D, region), D, region)
            results[f"{name}|episode_level"] = e; epl[name] = ev
            print(name, "episode level", e, flush=True)
    json.dump(epl, open(f"{OUT}/episode_level_values.json", "w"), indent=1, default=float)

    # partial R2 table (within-episode), primary designs
    pr2 = []
    for name in ["SVI_primary", "NW_primary"]:
        des = designs[name]; K = dict(np.load(f"{OUT}/kin_{name}.npz", allow_pickle=True))
        df = pd.read_csv(f"{OUT}/rows_{name}.csv").sort_values(["p", "ti"])
        r = df.r.to_numpy(); g = df.ti.to_numpy()
        T1 = np.where(df.cens, np.nan, np.log10(df.T1))
        cols = dict(T1=T1, T2=df.T2.to_numpy(), T3=df.T3.to_numpy(), in_reach=K["inr"].reshape(-1), km_outside=K["kout"].reshape(-1),
                    dir_mode=K["dm"].reshape(-1), u_along_path=K["u"].reshape(-1), detectability=K["det"].reshape(-1))
        for k, x in cols.items():
            b_, r2, n = fe_slope(r, x, g)
            pr2.append(dict(design=name, predictor=k, beta=b_, r2_within=r2, n=n))
        kin = np.column_stack([cols["in_reach"], cols["km_outside"], np.nan_to_num(cols["dir_mode"]), np.nan_to_num(cols["u_along_path"], nan=0.5)])
        bud = np.column_stack([np.nan_to_num(T1, nan=np.nanmedian(T1)), np.nan_to_num(cols["T2"]), np.nan_to_num(cols["T3"])])
        pr2.append(dict(design=name, predictor="kinematics_set", r2_within=fe_multi_r2(r, kin, g), n=len(r)))
        pr2.append(dict(design=name, predictor="budget_set", r2_within=fe_multi_r2(r, bud, g), n=len(r)))
        pr2.append(dict(design=name, predictor="all", r2_within=fe_multi_r2(r, np.column_stack([kin, bud, np.nan_to_num(cols["detectability"])]), g), n=len(r)))
    pd.DataFrame(pr2).to_csv(f"{OUT}/partial_r2.csv", index=False)

    # calibration (no-memory multinomial)
    cal = {}
    for name in ["SVI_primary", "NW_primary"]:
        cal[name] = calibrate(designs[name], nsim=NCAL)
        print("calibration", name, cal[name], flush=True)
    results["calibration"] = cal

    json.dump(results, open(f"{OUT}/tests.json", "w"), indent=1, default=float)
    pd.DataFrame(table).to_csv(f"{OUT}/tests_table.csv", index=False)


if __name__ == "__main__":
    main()
