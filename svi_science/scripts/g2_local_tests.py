#!/usr/bin/env python
"""g2_local_tests (2026-10-11) -- Study G2 tests (PLAN.md §1.8, §2-§4): completeness, Q1 (local earthquakes -> LFEs,
hours and days), Q2 (ETS episodes / between-episode LFE bursts -> local seismicity rate), Q3 (location bins).
Reuses study G's LFE tables (results/g_trigger/tables) and its delayed-test code.  Single process, < 20 GB.

    python svi_science/scripts/g2_local_tests.py  -> results/g2_local/tables/*
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy import stats
import g_trigger_tests as GT

ROOT = "/home/jovyan/tremorferometry"
G2 = f"{ROOT}/svi_science/results/g2_local"
T = f"{G2}/tables"
GTT = f"{ROOT}/svi_science/results/g_trigger/tables"
NS = 10 ** 9; H = 3600 * NS; DAY = 86400 * NS
NR, NSTK = 2000, 10000
RNG = np.random.default_rng(0)
GT.RNG = np.random.default_rng(1)
SPAN = {"SVI": ("2003-01-01", "2026-05-27"), "WA": ("2009-08-06", "2026-05-27"), "WA_pre": ("2007-01-01", "2009-08-05")}
CLASSES = ["crustal", "intraslab", "near_interface"]


def pt(ns):
    return pd.Timestamp(int(ns)).tz_localize("UTC").tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M %Z")


# ------------------------------------------------------------------------------------------- completeness
def maxc(m, bw=0.1):
    m = np.asarray(m); m = m[np.isfinite(m)]
    if len(m) < 50:
        return np.nan
    b = np.round(np.floor(m / bw + 1e-6) * bw, 1)
    u, n = np.unique(b, return_counts=True)
    return float(u[n.argmax()] + 0.2)


def completeness(A):
    rows, Mc = [], {}
    for reg in ("SVI", "WA"):
        X = A[A.in_analysis & (A[f"{reg}_dist_km"] <= 150)]
        y = X.t.dt.year
        vals = []
        for y0 in range(2003, 2027, 3):
            m = X.mag[(y >= y0) & (y < y0 + 3)]
            v = maxc(m); vals.append(v)
            rows.append(dict(region=reg, block=f"{y0}-{min(y0 + 2, 2026)}", n=len(m), Mc_maxc02=v))
        Mc[reg] = max(1.5, float(np.nanmax(vals)))
        rows.append(dict(region=reg, block="Mc_R (max of blocks, >= 1.5)", n=len(X), Mc_maxc02=Mc[reg]))
    pd.DataFrame(rows).to_csv(f"{T}/completeness.csv", index=False)
    print("Mc", Mc, flush=True)
    return Mc


# ------------------------------------------------------------------------------------------- LFE data
def load_lfe(reg):
    r = "WA" if reg.startswith("WA") else "SVI"
    K, days, ep = GT.load_region(r)
    if reg == "WA_pre":
        ep = ep[ep.kind == "wecc"].copy()
        inep = np.zeros(len(days), bool); inpad = np.zeros(len(days), bool)
        for s, e in zip(ep.s, ep.e):
            inep |= (days.dns >= s) & (days.dns < e); inpad |= (days.dns >= s - 2 * DAY) & (days.dns < e + 2 * DAY)
        days["in_ep"] = inep; days["in_pad"] = inpad
    if r == "SVI":
        Ev = pd.read_parquet(f"{GTT}/svi_events.parquet"); pc = "P"
    else:
        Ev = pd.read_parquet(f"{GTT}/wa_events.parquet"); pc = "P_own_WAVE_PNET"
    Ev = Ev.assign(tns=Ev.time_utc.astype("datetime64[ns]").astype("int64"), w=Ev[pc].astype(float))
    W = {}
    for k, sel in (("data", Ev.series == "data"), ("c1", Ev.series == "c1"),
                   ("mz", (Ev.series == "data") & (Ev.n_z2 >= 2)), ("c1mz", (Ev.series == "c1") & (Ev.n_z2 >= 2))):
        x = Ev[sel].sort_values("tns")
        W[k] = (x.tns.to_numpy(), np.r_[0.0, np.cumsum(x.w.to_numpy())])
    # P-weighted daily sums (for the delayed test)
    dP = days.copy()
    for k in ("data", "c1", "mz", "c1mz"):
        t, c = W[k]
        w = np.diff(c)
        dd = (t // DAY) * DAY
        s = pd.Series(w).groupby(dd).sum()
        dP[f"n_{k}"] = s.reindex(dP.dns).fillna(0).to_numpy()
    return K, W, days, dP, ep


def wcount(W, s, e):
    t, c = W
    return c[np.searchsorted(t, e, "left")] - c[np.searchsorted(t, s, "left")]


# ------------------------------------------------------------------------------------------- Q1 hours
WINS = {"H1": (0, H), "H6": (0, 6 * H), "H24": (0, 24 * H)}


def q1_trials(A, reg, Mc):
    r = "WA" if reg.startswith("WA") else "SVI"
    s0, s1 = GT.day_ns([SPAN[reg][0]])[0], GT.day_ns([SPAN[reg][1]])[0] + DAY
    X = A[A.in_analysis & A.gk_main & (A[f"{r}_dist_km"] <= 150) & (A.mag >= Mc[r])].copy()
    X["t0"] = X.t.astype("datetime64[ns]").astype("int64")
    X = X[(X.t0 >= s0) & (X.t0 < s1)].sort_values("t0").reset_index(drop=True)
    X["trial"] = np.arange(len(X))
    for c in ("dist_km", "epi_km", "static_kpa", "dyn_kpa", "zone", "in_strike", "along_km", "downdip_km", "nearest_fid"):
        X[c] = X[f"{r}_{c}"]
    X["beneath"] = (X.cls == "intraslab") & (X.epi_km <= 15)
    X["near"] = X.dist_km <= 30
    X["q3"] = np.where(~X.in_strike, "outside_strike", X.zone)
    return X


def q1_hours(reg, X, K, W, days, big_ns):
    keys_w = ("data", "c1", "mz")
    keys_c = ("data", "c1", "mz", "cand_data", "cand_c1", "trem100")
    D = days.set_index("dns")
    dns = days.dns.to_numpy(); grp = days.grp.to_numpy()
    ok = (~days.in_pad).to_numpy()
    near_big = np.zeros(len(days), bool)
    for t in big_ns:
        near_big |= np.abs(dns - (t // DAY) * DAY) <= 1 * DAY
    okset = ok & ~near_big
    oks = set(dns[okset])
    startok = np.array([okset[i] and (dns[i] - DAY in oks) and (dns[i] + DAY in oks) for i in range(len(dns))])
    rows, rnd, loc = [], {}, {}
    for _, e in X.iterrows():
        t0 = int(e.t0); d0 = (t0 // DAY) * DAY
        touch = np.unique([((t0 - 6 * H) // DAY) * DAY, d0, ((t0 + 6 * H - 1) // DAY) * DAY])
        if not all(d in D.index for d in touch):
            status = "not_observable"
        elif any(D.loc[d, "in_ep"] for d in touch):
            status = "ets"
        elif any(D.loc[d, "in_pad"] for d in touch):
            status = "pad"
        else:
            status = "between"
        g = D.loc[d0, "grp"] if d0 in D.index else None
        row = dict(trial=e.trial, status=status, grp=g)
        for wn, (a, b) in WINS.items():
            for k in keys_w:
                row[f"{wn}_P{k}"] = float(wcount(W[k], t0 + a, t0 + b))
                row[f"{wn}_Pbef_{k}"] = float(wcount(W[k], t0 - b, t0 - a))
            for k in keys_c:
                row[f"{wn}_{k}"] = int(GT.count(K[k], t0 + a, t0 + b))
        rows.append(row)
        tod = t0 - d0
        if status == "between":
            cand = dns[startok & (grp == g) & (np.abs(dns - d0) > 3 * DAY)]
            pick = RNG.choice(cand, NR, replace=len(cand) < NR) + tod
            R = {"n_pool": len(cand)}
            for wn, (a, b) in WINS.items():
                for k in keys_w:
                    R[f"{wn}_P{k}"] = wcount(W[k], pick + a, pick + b).astype(np.float32)
                    R[f"{wn}_Pbef_{k}"] = wcount(W[k], pick - b, pick - a).astype(np.float32)
                for k in keys_c:
                    R[f"{wn}_{k}"] = GT.count(K[k], pick + a, pick + b).astype(np.int32)
            rnd[e.trial] = R
        elif status == "ets":
            L = {}
            for wn in ("H1", "H6"):
                Lw = WINS[wn][1]; nt = int(24 * H // Lw)
                ts = np.concatenate([t0 - Lw * np.arange(1, nt + 1), t0 + Lw * np.arange(1, nt + 1)])
                for k in keys_w:
                    L[f"{wn}_P{k}"] = wcount(W[k], ts, ts + Lw)
                for k in keys_c:
                    L[f"{wn}_{k}"] = GT.count(K[k], ts, ts + Lw)
            loc[e.trial] = L
    TT = X.merge(pd.DataFrame(rows), on="trial")
    # per-trial random-window stats (H6) and the named rule
    ex = []
    for _, r in TT.iterrows():
        src = rnd.get(r.trial) if r.status == "between" else (loc.get(r.trial) if r.status == "ets" else None)
        if src is None:
            ex.append({}); continue
        x = {}
        for k in ("data", "c1", "mz", "cand_data", "cand_c1", "trem100"):
            v = np.asarray(src[f"H6_{k}"])
            x[f"mu_{k}"] = v.mean(); x[f"q99_{k}"] = np.quantile(v, 0.99); x[f"max_{k}"] = v.max()
        v = np.asarray(src["H6_Pdata"]); x["mu_Pdata"] = v.mean()
        x["pois_p"] = stats.poisson.sf(r.H6_data - 1, max(x["mu_data"], 1.0 / NR))
        ex.append(x)
    TT = pd.concat([TT.reset_index(drop=True), pd.DataFrame(ex)], axis=1)
    nT = int(TT.status.isin(["between", "ets"]).sum())
    TT["named_abcd"] = ((TT.H6_data > TT.max_data) & (TT.pois_p < 0.01 / max(nT, 1)) & (TT.H6_mz >= 2)
                        & (TT.H6_c1 <= TT.q99_c1) & (TT.H6_cand_c1 <= TT.q99_cand_c1)).fillna(False)
    TT["H6_excess_P"] = TT.H6_Pdata - TT.mu_Pdata
    return TT, rnd, loc


def stack(TT, src, sel, wn, ets=False):
    S = TT[sel]
    if not len(S):
        return dict(n=0)
    tr = S.trial.to_numpy()
    picks = {t: RNG.integers(0, len(src[t][f"{wn}_Pdata"]), NSTK) for t in tr}
    res = dict(n=len(S))

    def M(fn):
        acc = np.zeros(NSTK)
        for t in tr:
            acc += fn(src[t])[picks[t]]
        return acc
    defs = {
        "S": (lambda s: np.asarray(s[f"{wn}_Pdata"], float), S[f"{wn}_Pdata"].sum()),
        "X": (lambda s: np.asarray(s[f"{wn}_Pdata"], float) - np.asarray(s[f"{wn}_Pc1"], float),
              (S[f"{wn}_Pdata"] - S[f"{wn}_Pc1"]).sum()),
        "Smz": (lambda s: np.asarray(s[f"{wn}_Pmz"], float), S[f"{wn}_Pmz"].sum()),
        "Sc1": (lambda s: np.asarray(s[f"{wn}_Pc1"], float), S[f"{wn}_Pc1"].sum()),
        "N": (lambda s: np.asarray(s[f"{wn}_data"], float), S[f"{wn}_data"].sum()),
        "cand_data": (lambda s: np.asarray(s[f"{wn}_cand_data"], float), S[f"{wn}_cand_data"].sum()),
        "cand_c1": (lambda s: np.asarray(s[f"{wn}_cand_c1"], float), S[f"{wn}_cand_c1"].sum()),
        "trem100": (lambda s: np.asarray(s[f"{wn}_trem100"], float), S[f"{wn}_trem100"].sum()),
    }
    if not ets:
        defs["BA"] = (lambda s: (np.asarray(s[f"{wn}_Pdata"], float) - np.asarray(s[f"{wn}_Pc1"], float))
                      - (np.asarray(s[f"{wn}_Pbef_data"], float) - np.asarray(s[f"{wn}_Pbef_c1"], float)),
                      ((S[f"{wn}_Pdata"] - S[f"{wn}_Pc1"]) - (S[f"{wn}_Pbef_data"] - S[f"{wn}_Pbef_c1"])).sum())
    for k, (fn, o) in defs.items():
        m = M(fn)
        res[f"obs_{k}"] = float(o); res[f"rnd_mean_{k}"] = float(m.mean())
        res[f"p_{k}"] = float((1 + (m >= o - 1e-9).sum()) / (1 + NSTK))
        if k == "S":
            q01 = np.quantile(m, 0.01)
            res["ub99_extra_P_per_eq"] = float((o - q01) / len(S))
            res["ub99_rel"] = float((o - q01) / m.mean()) if m.mean() > 0 else np.nan
    res["PASS"] = bool(res["p_S"] < 0.01 and res["p_X"] < 0.01 and res["p_Smz"] < 0.01)
    return res


def shift_check(TT, reg, W, days, sel, wn="H6", ndraw=2000):
    """Segment-wise circular shift of all selected trial times by one common whole-day offset (keeps clustering)."""
    S = TT[sel]
    if not len(S):
        return {}
    D = days.set_index("dns")
    btw = np.sort(days.dns[~days.in_pad].to_numpy())
    a, b = WINS[wn]

    def val(t0s):
        d = (t0s // DAY) * DAY
        okm = np.isin(d, btw)
        if okm.sum() == 0:
            return np.nan
        tt = t0s[okm]
        return float(np.mean(wcount(W["data"], tt + a, tt + b) - wcount(W["c1"], tt + a, tt + b)))
    obs = val(S.t0.to_numpy())
    out = []
    grp_days = {g: np.sort(days.dns[days.grp == g].to_numpy()) for g in S.grp.dropna().unique()}
    for _ in range(ndraw):
        allt = []
        for g, gd in grp_days.items():
            t = S.t0[S.grp == g].to_numpy()
            span = gd[-1] - gd[0] + DAY
            off = RNG.integers(1, span // DAY) * DAY
            allt.append(gd[0] + (t - gd[0] + off) % span)
        out.append(val(np.concatenate(allt)))
    out = np.array(out)
    return dict(shift_obs_meanX=obs, shift_mean=float(np.nanmean(out)),
                p_shift=float((1 + (out >= obs).sum()) / (1 + np.isfinite(out).sum())))


# ------------------------------------------------------------------------------------------- dose
def dose(S, col):
    S = S.dropna(subset=["mu_data"])
    if len(S) < 5:
        return dict(n=len(S))
    N = S.H6_data.to_numpy(); mu = np.maximum(S.mu_data.to_numpy(), 1.0 / NR); sg = S[col].to_numpy()

    def nll(b):
        return -(stats.poisson.logpmf(N, mu * (1 + b * sg))).sum()
    bs = np.concatenate([[0], np.logspace(-6, 4, 3000)])
    L = np.array([nll(b) for b in bs])
    inn = L <= L.min() + 5.412 / 2
    resp = (N + 0.5) / (mu + 0.5)
    rho = stats.spearmanr(sg, resp).statistic
    perm = np.array([stats.spearmanr(RNG.permutation(sg), resp).statistic for _ in range(2000)])
    p = (1 + (perm >= rho).sum()) / 2001
    return dict(n=len(S), stress=col, beta_hat=float(bs[L.argmin()]), beta_ub99=float(bs[inn].max()),
                beta_lb99=float(bs[inn].min()), rho=float(rho), p_rho=float(p),
                scales=bool(rho > 0 and p < 0.01 and bs[inn].min() > 0))


# ------------------------------------------------------------------------------------------- Q2
def place_random(durs, n_days, pad, rng):
    """Random non-overlapping placement (with pads) of episodes of given lengths in [0, n_days)."""
    occ = np.zeros(n_days, bool); starts = []
    for L in durs:
        for _ in range(200):
            s = int(rng.integers(0, n_days - L + 1))
            if not occ[max(0, s - pad):min(n_days, s + L + pad)].any():
                occ[s:s + L] = True; starts.append((s, L)); break
    return starts


def q2(A, reg, Mc, dP, ep, K):
    r = reg
    s0, s1 = GT.day_ns([SPAN[r][0]])[0], GT.day_ns([SPAN[r][1]])[0] + DAY
    dgrid = np.arange(s0, s1, DAY); nd = len(dgrid)
    X = A[A.in_analysis & A.gk_main & (A[f"{r}_dist_km"] <= 150) & (A.mag >= Mc[r])].copy()
    X["t0"] = X.t.astype("datetime64[ns]").astype("int64")
    X = X[(X.t0 >= s0) & (X.t0 < s1)]
    X["di"] = ((X.t0 - s0) // DAY).astype(int)
    X["q3"] = np.where(~X[f"{r}_in_strike"], "outside_strike", X[f"{r}_zone"])
    X["beneath"] = (X.cls == "intraslab") & (X[f"{r}_epi_km"] <= 15)
    X["near"] = X[f"{r}_dist_km"] <= 30
    sets = {"all": X.index == X.index, **{c: X.cls.eq(c) for c in CLASSES + ["overriding_deep"]}}
    for z in ("updip", "lfe_zone", "downdip", "outside_strike"):
        sets[f"zone_{z}"] = X.q3.eq(z)
        for c in ("crustal", "intraslab"):
            sets[f"zone_{z}_{c}"] = X.q3.eq(z) & X.cls.eq(c)
    sets["beneath_intraslab"] = X.beneath; sets["near30"] = X.near
    cnt = {k: np.bincount(X.di[m], minlength=nd).astype(float) for k, m in sets.items()}
    # episode masks
    E = ep.copy(); E["si"] = ((E.s - s0) // DAY).astype(int); E["ei"] = ((E.e - s0) // DAY).astype(int)
    E = E[(E.ei > 0) & (E.si < nd)]

    def masks(si, ei, pad=2):
        inE = np.zeros(nd, bool); inP = np.zeros(nd, bool)
        for a, b in zip(si, ei):
            inE[max(a, 0):min(b, nd)] = True; inP[max(a - pad, 0):min(b + pad, nd)] = True
        return inE, ~inP
    out = []

    setn = list(cnt); Cm = np.vstack([cnt[k] for k in setn])          # sets x days

    def run(label, mEO_obs, ctrl_masks):
        mE, mO = mEO_obs
        nE_o = Cm[:, mE].sum(1); nO_o = Cm[:, mO].sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):
            RRo = (nE_o / max(mE.sum(), 1)) / (nO_o / max(mO.sum(), 1))
        ctrl = {}
        for cn, (ME, MO) in ctrl_masks.items():
            ME = np.asarray(ME, float); MO = np.asarray(MO, float)
            NE = Cm @ ME.T; NO = Cm @ MO.T                                    # sets x draws
            with np.errstate(divide="ignore", invalid="ignore"):
                ctrl[cn] = (NE / np.maximum(ME.sum(1), 1)) / (NO / np.maximum(MO.sum(1), 1))
        for i, k in enumerate(setn):
            o, nE, nO = RRo[i], nE_o[i], nO_o[i]
            res = dict(region=r, test=label, set=k, n_eq=int(Cm[i].sum()), N_E=int(nE), days_E=int(mE.sum()),
                       N_O=int(nO), days_O=int(mO.sum()), RR=float(o))
            if nE + nO > 0:
                pE = mE.sum() / (mE.sum() + mO.sum())
                lo = stats.beta.ppf(0.01, nE, nO + 1) if nE > 0 else 0.0
                hi = stats.beta.ppf(0.99, nE + 1, nO) if nO > 0 else 1.0
                res["RR_ub99_exact"] = float((hi / (1 - hi)) / (pE / (1 - pE))) if hi < 1 else np.inf
                res["RR_lb99_exact"] = float((lo / (1 - lo)) / (pE / (1 - pE)))
            for cn, V in ctrl.items():
                v = V[i]; v = v[np.isfinite(v)]
                if len(v) and np.isfinite(o):
                    res[f"p_{cn}"] = float((1 + (v >= o - 1e-12).sum()) / (1 + len(v)))
                    res[f"p_low_{cn}"] = float((1 + (v <= o + 1e-12).sum()) / (1 + len(v)))
                    res[f"ctrl_mean_{cn}"] = float(v.mean())
                    q01 = np.quantile(v, 0.01)
                    res[f"RR_ub99_{cn}"] = float(o / q01) if q01 > 0 else np.nan
            ps = [res.get(f"p_{cn}", 1) for cn in ctrl_masks]
            res["pass_both"] = bool(len(ps) >= 2 and max(ps) < 0.01)
            out.append(res)

    def rolled(mE, mO, offs):
        return (np.array([np.roll(mE, o2) for o2 in offs]), np.array([np.roll(mO, o2) for o2 in offs]))
    NDRAW = 10000
    rng = np.random.default_rng(2)
    # (a) all episodes and majors
    for lab, EE in (("ets_all", E), ("ets_major", E[E.major.astype(bool)])):
        durs = (EE.ei - EE.si).to_numpy(); si = EE.si.to_numpy(); ei = EE.ei.to_numpy()
        obs = masks(si, ei)
        shift = rolled(obs[0], obs[1], rng.integers(1, nd, NDRAW))
        pE, pO = [], []
        for _ in range(NDRAW):
            pl = place_random(rng.permutation(durs), nd, 2, rng)
            a, b = masks([x[0] for x in pl], [x[0] + x[1] for x in pl])
            pE.append(a); pO.append(b)
        run(lab, obs, {"shift": shift, "place": (np.array(pE), np.array(pO))})
        # after / before columns (same controls, shifted the same way)
        aft = np.zeros(nd, bool); bef = np.zeros(nd, bool)
        for a, b in zip(si, ei):
            aft[max(b, 0):min(b + 10, nd)] = True; bef[max(a - 10, 0):min(a, nd)] = True
        aft &= ~obs[0]; bef &= ~obs[0]
        run(lab + "_after10", (aft, obs[1] & ~aft), {"shift": rolled(aft, obs[1] & ~aft, rng.integers(1, nd, 2000))})
        run(lab + "_before10", (bef, obs[1] & ~bef), {"shift": rolled(bef, obs[1] & ~bef, rng.integers(1, nd, 2000))})
        if r == "SVI" and lab == "ets_all":
            for per, (a, b) in (("pre2010", ("2003-01-01", "2009-12-31")), ("post2010", ("2010-01-01", "2026-05-27"))):
                pm = (dgrid >= GT.day_ns([a])[0]) & (dgrid <= GT.day_ns([b])[0])
                run(f"{lab}_{per}", (obs[0] & pm, obs[1] & pm),
                    {"shift": (shift[0][:2000] & pm, shift[1][:2000] & pm)})
    # tremor days
    tr = K["trem100"]
    td = np.bincount(((tr[(tr >= s0) & (tr < s1)] - s0) // DAY).astype(int), minlength=nd)
    t0i = int((GT.day_ns(["2010-01-06"])[0] - s0) // DAY)
    valid = np.arange(nd) >= max(t0i, 0)
    mT = valid & (td >= 10); mQ = valid & (td == 0)
    a, b = rolled(mT, mQ, rng.integers(1, nd, 2000))
    run("tremor_days", (mT, mQ), {"shift": (a & valid, b & valid)})
    # LFE-active days (multi-site leg)
    dd = dP.set_index("dns")
    obsd = np.isin(dgrid, dP.dns.to_numpy())
    nmz = dd.n_mz_cnt.reindex(dgrid).fillna(0).to_numpy() if "n_mz_cnt" in dd else None
    act = obsd & (nmz >= 5); sil = obsd & (dd.n_data_cnt.reindex(dgrid).fillna(0).to_numpy() == 0)
    oi = np.where(obsd)[0]
    PA = np.zeros((NDRAW, nd), bool); PB = np.zeros((NDRAW, nd), bool)
    na_, ns_ = int(act.sum()), int(sil.sum())
    for j in range(NDRAW):
        p = rng.permutation(oi)
        PA[j, p[:na_]] = True; PB[j, p[na_:na_ + ns_]] = True
    run("lfe_active_days", (act, sil), {"shift": rolled(act, sil, rng.integers(1, nd, NDRAW)), "place": (PA, PB)})
    del PA, PB
    # (b) between-episode bursts
    obsE = masks(E.si.to_numpy(), E.ei.to_numpy())
    btw_days = obsd & obsE[1]
    bi = np.where(btw_days)[0]
    bursts = []
    for ser in ("data", "c1"):
        t = np.sort(K[ser]); mzt = np.sort(K["mz" if ser == "data" else "c1mz"])
        t = t[(t >= s0) & (t < s1)]
        di = ((t - s0) // DAY).astype(int)
        t = t[btw_days[di]]
        if not len(t):
            continue
        brk = np.where(np.diff(t) >= 24 * H)[0]
        for a, b in zip(np.r_[0, brk + 1], np.r_[brk, len(t) - 1]):
            n = b - a + 1
            nm = GT.count(mzt, t[a], t[b] + 1)
            if n >= 3 and nm >= 2:
                bursts.append(dict(region=r, series=ser, start_PT=pt(t[a]), end_PT=pt(t[b]), n=int(n), n_mz=int(nm),
                                   si=int((t[a] - s0) // DAY), ei=int((t[b] - s0) // DAY) + 2))
    B = pd.DataFrame(bursts)
    if len(B):
        Bd = B[B.series == "data"]
        for lab, ext in (("bursts", 0), ("bursts_plus3", 2)):
            mB = np.zeros(nd, bool)
            for a, b in zip(Bd.si, Bd.ei + ext):
                mB[a:min(b, nd)] = True
            mB &= btw_days; mO = btw_days & ~mB
            pos = np.searchsorted(bi, np.where(mB)[0])
            SH = np.zeros((NDRAW, nd), bool); PL = np.zeros((NDRAW, nd), bool)
            for j in range(NDRAW):
                o2 = int(rng.integers(1, len(bi)))
                SH[j, bi[(pos + o2) % len(bi)]] = True
                for a, b in zip(Bd.si, Bd.ei + ext):
                    L = b - a; k = int(rng.integers(0, len(bi) - L))
                    PL[j, bi[k:k + L]] = True
            run(lab, (mB, mO), {"shift": (SH, btw_days[None] & ~SH), "place": (PL, btw_days[None] & ~PL)})
            del SH, PL
        B.to_csv(f"{T}/q2_bursts_{r}.csv", index=False)
    # detection check
    det = []
    mE, mO = obsE[0], obsE[1]
    Xall = A[A.in_analysis & (A[f"{r}_dist_km"] <= 150)].copy()
    Xall["t0"] = Xall.t.astype("datetime64[ns]").astype("int64")
    Xall = Xall[(Xall.t0 >= s0) & (Xall.t0 < s1)]
    di = ((Xall.t0 - s0) // DAY).astype(int).to_numpy()
    for lab, m in (("ETS_days", mE), ("other_days", mO)):
        mm = Xall.mag.to_numpy()[m[di]]
        bs = [maxc(rng.choice(mm, len(mm))) for _ in range(1000)]
        small = ((mm >= Mc[r] - 0.5) & (mm < Mc[r])).sum(); big = (mm >= Mc[r]).sum()
        det.append(dict(region=r, days=lab, n_events=len(mm), Mc=maxc(mm), Mc_q05=np.nanquantile(bs, 0.05),
                        Mc_q95=np.nanquantile(bs, 0.95), masking_index=small / max(big, 1), n_ge_Mc=int(big),
                        n_days=int(m.sum())))
    return pd.DataFrame(out), pd.DataFrame(det), B


# ------------------------------------------------------------------------------------------- main
def main():
    A = pd.read_parquet(f"{T}/eq_local.parquet")
    A["t"] = A.t.astype("datetime64[ns]")
    Mc = completeness(A)
    summ = dict(Mc=Mc)
    st_rows, dl_rows, dose_rows, q2_rows, det_rows, sh_rows = [], [], [], [], [], []
    for reg in ("SVI", "WA", "WA_pre"):
        r = "WA" if reg.startswith("WA") else "SVI"
        K, W, days, dP, ep = load_lfe(reg)
        X = q1_trials(A, reg, Mc)
        big = X.t0[X.mag >= 3.0].to_numpy()
        TT, rnd, loc = q1_hours(reg, X, K, W, days, big)
        keep = [c for c in TT.columns if c not in ("t",)]
        TT[keep].to_parquet(f"{T}/q1_trials_{reg}.parquet", index=False)
        btw = TT.status.eq("between"); ets = TT.status.eq("ets")
        sets = {"all": btw}
        for c in CLASSES + ["overriding_deep"]:
            sets[c] = btw & TT.cls.eq(c)
            sets[f"{c}_mccrory"] = btw & TT.cls_mccrory.eq(c)
        for lab, m in (("M<2.5", TT.mag < 2.5), ("M2.5-3.5", (TT.mag >= 2.5) & (TT.mag < 3.5)), ("M>=3.5", TT.mag >= 3.5),
                       ("d<30", TT.dist_km < 30), ("d30-75", (TT.dist_km >= 30) & (TT.dist_km < 75)),
                       ("d75-150", TT.dist_km >= 75), ("dyn>=1kPa", TT.dyn_kpa >= 1),
                       ("dyn0.1-1kPa", (TT.dyn_kpa >= 0.1) & (TT.dyn_kpa < 1)), ("dyn<0.1kPa", TT.dyn_kpa < 0.1),
                       ("static>=1Pa", TT.static_kpa >= 1e-3), ("beneath_intraslab", TT.beneath), ("near30", TT.near)):
            sets[lab] = btw & m
        for z in ("updip", "lfe_zone", "downdip", "outside_strike"):
            sets[f"zone_{z}"] = btw & TT.q3.eq(z)
            for c in ("crustal", "intraslab"):
                sets[f"zone_{z}_{c}"] = btw & TT.q3.eq(z) & TT.cls.eq(c)
        for g in sorted(TT.grp.dropna().unique()):
            sets[f"grp_{g}"] = btw & TT.grp.eq(g)
        for nm, sel in sets.items():
            for wn in (("H6", "H1", "H24") if nm in ("all", "crustal", "intraslab", "near_interface") else ("H6",)):
                st_rows.append(dict(region=reg, set=nm, window=wn, headline=bool(wn == "H6" and nm in ["all"] + CLASSES),
                                    **stack(TT, rnd, sel, wn)))
        for nm, sel in (("ets_all", ets), ("ets_dyn>=1kPa", ets & (TT.dyn_kpa >= 1)), ("ets_near30", ets & TT.near)):
            for wn in ("H1", "H6"):
                st_rows.append(dict(region=reg, set=nm, window=f"{wn}_vs_local", headline=False, **stack(TT, loc, sel, wn, ets=True)))
        for nm in ("all", "crustal", "intraslab", "near_interface"):
            sh_rows.append(dict(region=reg, set=nm, **shift_check(TT, reg, W, days, sets[nm])))
        for col in ("dyn_kpa", "static_kpa"):
            dose_rows.append(dict(region=reg, test="H6_between", **dose(TT[btw], col)))
        # delayed (G test 2) with P-weighted daily sums (headline) and counts (column)
        Tr = X.assign(inst_s=X.t0, sigma_kpa=X.dyn_kpa)
        for wtag, dtab in (("P", dP), ("count", days)):
            for nm in ["all"] + CLASSES + ["beneath_intraslab", "near30", "zone_updip", "zone_lfe_zone", "zone_downdip"]:
                if nm == "all":
                    sel = Tr.mag > -9
                elif nm in CLASSES:
                    sel = Tr.cls.eq(nm)
                elif nm == "beneath_intraslab":
                    sel = Tr.beneath
                elif nm == "near30":
                    sel = Tr.near
                else:
                    sel = Tr.q3.eq(nm[5:])
                for win in (3, 10):
                    if wtag == "count" and (win == 10 or nm not in ["all"] + CLASSES):
                        continue
                    GT.delayed.__defaults__[-1].clear()
                    Rr, res, _ = GT.delayed(reg, Tr, dtab, big, sel, win)
                    if len(Rr):
                        Rr.assign(set=nm, win=win, weight=wtag).to_csv(f"{T}/q1_delayed_{reg}_{nm}_{win}_{wtag}.csv", index=False)
                    dl_rows.append(dict(region=reg, set=nm, days=f"+1..+{win}", weight=wtag,
                                        headline=bool(win == 3 and wtag == "P" and nm in ["all"] + CLASSES), **res))
        summ[reg] = dict(n_trials=len(TT), status=TT.status.value_counts().to_dict(),
                         named=TT.loc[TT.named_abcd, ["trial", "t_PT", "place", "mag", "cls", "dist_km", "status", "H6_data",
                                                      "mu_data"]].to_dict("records"))
        print(reg, json.dumps(summ[reg], default=str)[:800], flush=True)
        if reg != "WA_pre":
            dP2 = dP.copy()
            dP2["n_mz_cnt"] = days.n_mz.to_numpy(); dP2["n_data_cnt"] = days.n_data.to_numpy()
            q, d, _ = q2(A, r, Mc, dP2, ep, K)
            q2_rows.append(q); det_rows.append(d)
            print(reg, "Q2 done", flush=True)
    pd.DataFrame(st_rows).to_csv(f"{T}/q1_hours_stacks.csv", index=False)
    pd.DataFrame(sh_rows).to_csv(f"{T}/q1_shift_check.csv", index=False)
    pd.DataFrame(dl_rows).to_csv(f"{T}/q1_delayed_stacks.csv", index=False)
    pd.DataFrame(dose_rows).to_csv(f"{T}/q1_dose.csv", index=False)
    pd.concat(q2_rows).to_csv(f"{T}/q2_rates.csv", index=False)
    pd.concat(det_rows).to_csv(f"{T}/q2_detection.csv", index=False)
    json.dump(summ, open(f"{T}/g2_summary.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
