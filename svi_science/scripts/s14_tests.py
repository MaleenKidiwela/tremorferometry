#!/usr/bin/env python
"""s14_tests -- S14 steps C/D: the pre-registered comparisons C1-C7 on the stretching CC (PLAN §5) with the pass rules (§6).
    [S11_REGION=WA_NW] python s14_tests.py [--pilot] [--nperm 1000]
Inputs: episode_raw (S11), s14 store tables (store_episode_*, store_day_*), round-7 front times.
"""
import os, sys, argparse, json, time
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import pandas as pd
import s14_common as C

REGION = os.environ.get("S11_REGION", "SVI").upper()
R7 = {"SVI": f"{C.S11}/round7/svi/tables", "WA_NW": f"{C.S11}/round7/wa_nw/tables"}[REGION]
EP_WINS = ["25", "2035", "3560", "36", "Z25", "H25", "S"]
DAY_WINS = ["25", "2035", "3560", "S"]
COVS = ["lnn", "lnr", "lnd"]


# ------------------------------------------------------------------------------ fast FE machinery
def gdemean(x, codes, ng):
    s = np.bincount(codes, weights=x, minlength=ng); c = np.bincount(codes, minlength=ng)
    return x - (s / np.maximum(c, 1))[codes]


class Design:
    """fixed-effect design on a subset: y, covariates; gamma(term) by Frisch-Waugh"""
    def __init__(self, D, ycol, covs, fe):
        ok = np.isfinite(D[ycol].to_numpy(float))
        for c in covs:
            ok &= np.isfinite(D[c].to_numpy(float))
        self.D = D[ok]
        self.codes, u = pd.factorize(self.D[fe]); self.ng = len(u)
        self.y = gdemean(self.D[ycol].to_numpy(float), self.codes, self.ng)
        X = np.column_stack([gdemean(self.D[c].to_numpy(float), self.codes, self.ng) for c in covs]) if covs else np.zeros((len(self.D), 0))
        self.X = X[:, X.std(0) > 1e-12] if X.shape[1] else X
        if self.X.shape[1]:
            self.P = np.linalg.pinv(self.X.T @ self.X) @ self.X.T
            self.yr = self.y - self.X @ (self.P @ self.y)
        else:
            self.yr = self.y
        self.covs = covs

    def gamma(self, t):
        t = gdemean(np.asarray(t, float), self.codes, self.ng)
        tr = t - self.X @ (self.P @ t) if self.X.shape[1] else t
        den = tr @ tr
        if den <= 0:
            return np.nan, np.nan
        g = float(tr @ self.yr / den)
        r = self.yr - g * tr
        dof = max(len(r) - self.X.shape[1] - 1 - self.ng, 1)
        return g, float(np.sqrt((r @ r / dof) / den))

    def beta(self):
        if not self.X.shape[1]:
            return {}
        return dict(zip(self.covs, self.P @ self.y))


def perm_labels(D, lab, rng, half=45):
    """shuffle a label among cells of the same series within doy blocks of 2*half days (random phase); vectorised"""
    off = rng.integers(0, 2 * half)
    blk = ((D.doy.to_numpy() + off) // (2 * half)).astype(np.int64)
    key = pd.factorize(D.series.to_numpy())[0].astype(np.int64) * 1000 + blk
    lab = np.asarray(lab)
    i_st = np.argsort(key, kind="stable")
    i_rn = np.lexsort((rng.random(len(key)), key))
    out = lab.copy(); out[i_st] = lab[i_rn]
    return out


def matched(D, ycol, lab):
    """within series: pair each label-1 cell with the nearest label-0 cell (|dlnn|<=0.3, |dlnr|<=0.3, no reuse)"""
    d = []
    for s, q in D.groupby("series"):
        a, b = q[q["_lab"] == 1], q[q["_lab"] == 0]
        if not len(a) or not len(b):
            continue
        used = set()
        for i, r in a.iterrows():
            dist = np.hypot(b.lnn - r.lnn, b.lnr - r.lnr)
            ok = (np.abs(b.lnn - r.lnn) <= 0.3) & (np.abs(b.lnr - r.lnr) <= 0.3) & ~b.index.isin(list(used))
            if ok.any():
                j = dist[ok].idxmin(); used.add(j); d.append(r[ycol] - b.loc[j, ycol])
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    return (float(d.mean()), float(d.std(ddof=1) / np.sqrt(len(d))), len(d)) if len(d) > 2 else (np.nan, np.nan, len(d))


# ------------------------------------------------------------------------------ data
def load(pilot):
    E = C.load_episode(REGION, pilot=pilot)
    tag = f"{REGION}{'_pilot' if pilot else ''}"
    SE = pd.read_parquet(f"{C.TB}/store_episode_{tag}.parquet")
    E = E.merge(SE[["fam", "netsta", "epoch", "episode", "rec_cc_S", "rec_eps_S", "next_cc_25", "next_eps_25", "prev_wsum", "wsum",
                    "gap_days", "prev_episode"]], on=["fam", "netsta", "epoch", "episode"], how="left")
    E["z_S"] = C.fz(E.rec_cc_S)
    E["z_next"] = C.fz(E.next_cc_25)
    E["absdeps"] = (E.eps3c_25 - E.groupby("series").eps3c_25.transform("median")).abs() * 100
    E["lnn_c2"] = np.log(E.n_ctrl2.clip(lower=1))
    Dd = pd.read_parquet(f"{C.TB}/store_day_{tag}.parquet")
    return E, Dd, tag


def ycol_of(w):
    return {"Z25": "z_Z25", "H25": "z_H25", "S": "z_S"}.get(w, f"z_{w}")


# ------------------------------------------------------------------------------ class comparisons (C1, C2)
def compare(E, lab, name, rng, nperm, wins=EP_WINS, extra=True):
    rows = []
    D0 = E[np.isfinite(lab)].copy(); D0["_lab"] = lab[np.isfinite(lab)].astype(int)
    for w in wins:
        y = ycol_of(w)
        des = Design(D0, y, COVS, "series")
        if len(des.D) < 20 or des.D._lab.nunique() < 2:
            continue
        L = des.D._lab.to_numpy()
        g, se = des.gamma(L)
        perm = np.array([des.gamma(perm_labels(des.D, L, rng))[0] for _ in range(nperm)])
        p = float((np.sum(np.abs(perm - np.nanmedian(perm)) >= abs(g - np.nanmedian(perm))) + 1) / (np.isfinite(perm).sum() + 1))
        raw = Design(des.D, y, [], "series").gamma(L)
        r = dict(comp=name, win=w, n_cells=len(des.D), n1=int(L.sum()), gamma=g, se=se, perm_p=p, perm_sd=float(np.nanstd(perm)),
                 perm_med=float(np.nanmedian(perm)), raw_gamma=raw[0], raw_se=raw[1], **{f"b_{k}": v for k, v in des.beta().items()})
        if extra:
            mg, ms, mn = matched(des.D, y, None)
            r.update(matched=mg, matched_se=ms, matched_n=mn)
            # co-located / instrument split
            for grp in C.COLOC[REGION]:
                for ns in grp:
                    q = des.D[des.D.netsta == ns]
                    if len(q) > 15 and q._lab.nunique() == 2:
                        dq = Design(q, y, COVS, "series"); r[f"g_{ns}"] = dq.gamma(dq.D._lab.to_numpy())[0]
            for par, nm in ((0, "even"), (1, "odd")):
                q = des.D[des.D.fam % 2 == par]
                if len(q) > 15 and q._lab.nunique() == 2:
                    dq = Design(q, y, COVS, "series"); r[f"g_{nm}"] = dq.gamma(dq.D._lab.to_numpy())[0]
            q = des.D[~des.D.glitch & ~des.D.shift_flag]
            dq = Design(q, y, COVS, "series"); Lq = dq.D._lab.to_numpy()
            gq = dq.gamma(Lq)[0]
            pq = np.array([dq.gamma(perm_labels(dq.D, Lq, rng))[0] for _ in range(max(nperm // 2, 100))])
            r.update(g_noglitch=gq, p_noglitch=float((np.sum(np.abs(pq - np.nanmedian(pq)) >= abs(gq - np.nanmedian(pq))) + 1) / (len(pq) + 1)))
            dq = Design(D0, y, COVS + ["absdeps"], "series"); r["g_reg_eps"] = dq.gamma(dq.D._lab.to_numpy())[0]
            if w in ("25", "2035", "3560", "36"):
                dc = Design(D0, f"cz_{w}", ["lnn_c2", "lnd"], "series")
                if len(dc.D) > 20 and dc.D._lab.nunique() == 2:
                    r["g_c2"], r["se_c2"] = dc.gamma(dc.D._lab.to_numpy())
        rows.append(r)
    return rows


# ------------------------------------------------------------------------------ C3 within episode (day level)
def c3(Dd, E, rng, nperm, pilot):
    U = C.units(REGION)
    F = pd.read_csv(f"{R7}/r7_front_times.csv")
    F = F[F.set == "headline"].copy()
    F["start"] = pd.to_datetime(F.start)
    D = Dd[Dd.episode.notna()].copy()
    D["cls"] = D.episode.map(U.cls)
    D = D[D.cls == "major"].copy()
    D["date"] = pd.to_datetime(D.date)
    D = D.merge(F[["fam", "major", "start", "rel"]], left_on=["fam", "episode"], right_on=["fam", "major"], how="inner")
    D["lnws"] = np.log(D.ws.clip(lower=1e-3)); D["lnr"] = np.log(D.r_day.clip(lower=1e-3))
    D["sm"] = D.series + "#" + D.episode
    D["doy"] = D.date.dt.dayofyear
    for w in DAY_WINS:
        D[f"z_{w}"] = C.fz(D[f"lo_cc_{w}"])           # Amendment 1: leave-day-out reference
        D[f"zg_{w}"] = C.fz(D[f"cc_{w}"])             # grand-stack reference (as round 7), column
    # front-offset derangement among a family's majors (scrambled control)
    fm = F.groupby("fam")
    rels = {f: (q.major.to_numpy(), q.rel.to_numpy()) for f, q in fm}

    def lags(Dx, relmap=None):
        rel = Dx.rel.to_numpy() if relmap is None else np.array([relmap[(f, m)] for f, m in zip(Dx.fam, Dx.episode)])
        return ((Dx.date - Dx.start).dt.days.to_numpy() - rel)

    def contrasts(Dx, lag, y):
        out = {}
        for nm, a, b in (("late_minus_early", lag >= 3, (lag >= 0) & (lag <= 1)), ("pre_minus_early", lag <= -1, (lag >= 0) & (lag <= 1))):
            m = a | b
            q = Dx[m].copy()
            if len(q) < 30:
                out[nm] = (np.nan, np.nan); continue
            des = Design(q, y, ["lnws", "lnr"], "sm")
            lab = a[m][np.isfinite(q[y].to_numpy()) & np.isfinite(q.lnws) & np.isfinite(q.lnr)].astype(float)
            out[nm] = des.gamma(lab)
        return out

    rows = []
    for kind in ("data", "ctrl"):
        Dk = D[D.kind == kind]
        lag0 = lags(Dk)
        for w in DAY_WINS:
            y = f"z_{w}"
            obs = contrasts(Dk, lag0, y)
            scr = {k: [] for k in obs}
            for _ in range(nperm if kind == "data" else max(nperm // 5, 50)):
                relmap = {}
                for f, (mj, rl) in rels.items():
                    p = rng.permutation(len(rl))
                    if len(rl) > 1:
                        while np.any(p == np.arange(len(rl))):
                            p = rng.permutation(len(rl))
                    relmap.update({(f, m): r for m, r in zip(mj, rl[p])})
                o = contrasts(Dk, lags(Dk, relmap), y)
                for k in o:
                    scr[k].append(o[k][0])
            for k, (g, se) in obs.items():
                s = np.array(scr[k], float)
                p = float((np.sum(np.abs(s - np.nanmedian(s)) >= abs(g - np.nanmedian(s))) + 1) / (np.isfinite(s).sum() + 1))
                r = dict(comp=f"C3_{k}", kind=kind, win=w, n_days=len(Dk), gamma=g, se=se, perm_p=p, perm_sd=float(np.nanstd(s)),
                         perm_med=float(np.nanmedian(s)))
                if kind == "data" and w == "25":
                    for grp in C.COLOC[REGION]:
                        for ns in grp:
                            q = Dk[Dk.netsta == ns]
                            if len(q) > 30:
                                r[f"g_{ns}"] = contrasts(q, lags(q), y)[k][0]
                    for par, nm in ((0, "even"), (1, "odd")):
                        q = Dk[Dk.fam % 2 == par]
                        r[f"g_{nm}"] = contrasts(q, lags(q), y)[k][0] if len(q) > 30 else np.nan
                rows.append(r)
    # profile by lag bin (for the figure): REG residual means per bin
    Dk = D[D.kind == "data"].copy(); Dk["lag"] = lags(Dk)
    des = Design(Dk, "z_25", ["lnws", "lnr"], "sm")
    Q = des.D.assign(res=des.yr, lag=lags(des.D))
    prof = Q.assign(lb=pd.cut(Q.lag, [-99, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 6.5, 99], labels=["<=-2", "-1", "0", "1", "2", "3", "4-6", ">=7"])) \
        .groupby("lb", observed=True).agg(res=("res", "mean"), sd=("res", "std"), n=("res", "size"), cc=("lo_cc_25", "median"), lnws=("lnws", "mean"))
    return rows, prof.reset_index(), D


# ------------------------------------------------------------------------------ C4 / C5 / C6 / C7
def resid_table(E):
    out = {}
    for w in EP_WINS:
        des = Design(E, ycol_of(w), COVS, "series")
        out[w] = pd.Series(des.yr, index=des.D.index)
    return out


def coherence(E, res, rng, nperm):
    """C4a cross-family M per station; C4b co-located correlation"""
    Q = E.loc[res.index].assign(res=res.values)
    def M_of(Qx):
        Ms = []
        for ns, q in Qx.groupby("netsta"):
            P = q.pivot_table(index="episode", columns="fam", values="res", aggfunc="mean").to_numpy()
            if P.shape[1] < 3:
                continue
            fin = np.isfinite(P); X = np.where(fin, P, 0.0)
            S_, N_ = X.sum(1, keepdims=True), fin.sum(1, keepdims=True)
            O = (S_ - X) / np.maximum(N_ - fin, 1)                      # mean of the other families
            ok = fin & ((N_ - fin) >= 2)
            k = ok.sum(0)
            for j in np.flatnonzero(k >= 8):
                m = ok[:, j]; x, o = P[m, j], O[m, j]
                if x.std() > 0 and o.std() > 0:
                    Ms.append(np.corrcoef(x, o)[0, 1])
        return float(np.nanmedian(Ms)) if Ms else np.nan, len(Ms)
    M, nM = M_of(Q)
    scr = []
    for _ in range(nperm):
        Qs = Q.copy(); Qs["res"] = perm_labels(Qs, Qs.res.to_numpy(), rng)
        scr.append(M_of(Qs)[0])
    scr = np.array(scr)
    # co-located
    col = {}
    for grp in C.COLOC[REGION]:
        P = Q[Q.netsta.isin(grp)].pivot_table(index=["fam", "episode"], columns="netsta", values="res", aggfunc="mean")
        for a in grp:
            for b in grp:
                if a < b and a in P and b in P:
                    ok = P[a].notna() & P[b].notna()
                    if ok.sum() > 10:
                        col[f"{a}~{b}"] = (float(np.corrcoef(P[a][ok], P[b][ok])[0, 1]), int(ok.sum()))
    return dict(M=M, n_fam_series=nM, M_scr_p99=float(np.nanpercentile(scr, 99)), M_scr_med=float(np.nanmedian(scr)),
                M_p=float((np.sum(scr >= M) + 1) / (len(scr) + 1)), coloc=col)


def seasonal(E, res, rng, nperm):
    Q = E.loc[res.index].assign(res=res.values)
    def amp(doy):
        ph = 2 * np.pi * doy / 365.25
        X = np.column_stack([np.cos(ph), np.sin(ph)])
        b, *_ = np.linalg.lstsq(X, Q.res.to_numpy(), rcond=None)
        return float(np.hypot(*b)), float((np.degrees(np.arctan2(b[1], b[0])) % 360) * 365.25 / 360)
    A, pk = amp(Q.doy.to_numpy())
    codes, u = pd.factorize(Q.series)
    scr = []
    for _ in range(nperm):
        sh = rng.integers(0, 365, len(u))[codes]
        scr.append(amp((Q.doy.to_numpy() + sh) % 365)[0])
    scr = np.array(scr)
    return dict(amp_z=A, peak_doy=pk, scr_p95=float(np.percentile(scr, 95)), p=float((np.sum(scr >= A) + 1) / (len(scr) + 1)))


def longterm(E, res):
    from scipy.stats import theilslopes
    Q = E.loc[res.index].assign(res=res.values)
    s = [theilslopes(q.res, q.year)[0] for _, q in Q.groupby("series") if len(q) >= 8 and np.ptp(q.year) > 3]
    return dict(median_slope_z_per_yr=float(np.median(s)), n_series=len(s), frac_neg=float(np.mean(np.array(s) < 0)))


def cc_only_list(E, res, tag):
    Q = E.loc[res.index].assign(res=res.values)
    g = Q.groupby("series")
    Q["rs"] = g.res.transform(lambda x: 1.4826 * np.median(np.abs(x - np.median(x))))
    Q["dv_a"] = (Q.eps3c_25 - g.eps3c_25.transform("median")).abs()
    Q["flat"] = Q.dv_a < Q.sig3c_25
    Q["low"] = Q.res < -3 * Q.rs; Q["high"] = Q.res > 3 * Q.rs
    L = Q[Q.low & Q.flat]
    ep = L.groupby(["netsta", "episode"]).agg(nfam=("fam", "nunique"), date=("date", "median"), cls=("cls", "first"),
                                             glitch=("glitch", "any"), dcc_shift=("dcc_shift", "median")).reset_index()
    ep = ep[ep.nfam >= 3].sort_values("date")
    # site-common: co-located instruments' mean residual (in robust sigma) for that episode
    sc = []
    for _, r in ep.iterrows():
        grp = next((g_ for g_ in C.COLOC[REGION] if r.netsta in g_), None)
        if grp is None:
            sc.append(np.nan); continue
        q = Q[(Q.episode == r.episode) & Q.netsta.isin([x for x in grp if x != r.netsta])]
        sc.append(float((q.res / q.rs).mean()) if len(q) else np.nan)
    ep["coloc_mean_rsig"] = sc
    ep.to_csv(f"{C.TB}/cconly_episodes_{tag}.csv", index=False)
    Q[Q.low | Q.high][["fam", "netsta", "epoch", "episode", "cls", "date", "n", "cc3c_25", "res", "rs", "low", "high", "flat",
                       "dvv25", "sig25", "glitch", "dcc_shift"]].to_csv(f"{C.TB}/cconly_cells_{tag}.csv", index=False)
    return dict(cells=len(Q), low=int(Q.low.sum()), high=int(Q.high.sum()), low_flat=int((Q.low & Q.flat).sum()),
                high_flat=int((Q.high & Q.flat).sum()), episodes_ge3=len(ep),
                episodes_ge3_sitecommon=int((ep.coloc_mean_rsig < -1).sum()))


def successive(E):
    """C4c: CC between successive episode stacks vs gap and stack sizes (column)"""
    Q = E[np.isfinite(E.z_next) & np.isfinite(E.prev_wsum)].copy()
    Q["lnw1"], Q["lnw0"], Q["lngap"] = np.log(Q.wsum.clip(lower=1e-3)), np.log(Q.prev_wsum.clip(lower=1e-3)), np.log(Q.gap_days.clip(lower=1))
    Q["major"] = (Q.cls == "major").astype(float)
    des = Design(Q, "z_next", ["lnw1", "lnw0", "lngap"], "series")
    b = des.beta()
    g = des.gamma(des.D.major.to_numpy())
    return dict(n=len(des.D), b_lnw1=b.get("lnw1"), b_lnw0=b.get("lnw0"), b_lngap=b.get("lngap"),
                g_major_after=g[0], se_major_after=g[1], median_cc=float(Q.next_cc_25.median()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pilot", action="store_true"); ap.add_argument("--nperm", type=int, default=1000)
    a = ap.parse_args()
    rng = np.random.default_rng(14)
    t0 = time.time()
    E, Dd, tag = load(a.pilot)
    rows = []
    # C1 majors vs bursts
    lab = (E.cls == "major").astype(float).to_numpy()
    rows += compare(E, lab, "C1_major_minus_burst", rng, a.nperm)
    print(f"C1 done {time.time() - t0:.0f} s", flush=True)
    # C2 RF-like post vs pre (bursts only)
    b = E.cls == "burst"
    post = b & E.d_after_major.between(1, 45)
    pre = b & E.d_before_major.between(1, 60) & ~post
    lab = np.where(post, 1.0, np.where(pre, 0.0, np.nan))
    rows += compare(E, lab, "C2_post_minus_pre", rng, a.nperm)
    lab = np.where(post, 1.0, np.where(b, 0.0, np.nan))
    rows += compare(E, lab, "C2b_post_minus_otherbursts", rng, a.nperm, extra=False)
    print(f"C2 done {time.time() - t0:.0f} s", flush=True)
    T = pd.DataFrame(rows)
    # C3
    r3, prof, _ = c3(Dd, E, rng, max(a.nperm // 2, 200), a.pilot)
    T = pd.concat([T, pd.DataFrame(r3)], ignore_index=True)
    prof.to_csv(f"{C.TB}/c3_profile_{tag}.csv", index=False)
    print(f"C3 done {time.time() - t0:.0f} s", flush=True)
    T.to_csv(f"{C.TB}/tests_{tag}.csv", index=False)
    # C4-C7 on covariate-only residuals
    res = resid_table(E)
    extra = {}
    for w in ("25", "2035", "3560", "S"):
        extra[f"C4_coherence_{w}"] = coherence(E, res[w], rng, min(a.nperm, 300))
        extra[f"C5_seasonal_{w}"] = seasonal(E, res[w], rng, a.nperm)
        extra[f"C6_longterm_{w}"] = longterm(E, res[w])
    extra["C4c_successive"] = successive(E)
    extra["C7_cconly"] = cc_only_list(E, res["25"], tag)
    des = Design(E, "z_25", COVS, "series")
    yt = des.y; extra["REG_variance_explained_25"] = float(1 - des.yr.var() / yt.var()); extra["REG_beta_25"] = des.beta()
    json.dump(extra, open(f"{C.TB}/tests_extra_{tag}.json", "w"), indent=1, default=float)
    print(json.dumps(extra, indent=1, default=float)[:3000])
    open(f"{C.OUT}/TESTS_{tag}_DONE", "w").write(C.pt_now() + "\n")
    print(f"done {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
