#!/usr/bin/env python
"""s10_measure -- S10 pilot step 3 (PLAN §2-6): per (family, station, epoch) series: grand stack, S re-pick, cc
realign, P-weighted day stacks (+ P>=0.75 and scrambled-control stacks), S-anchored stretch vs the leave-day-out epoch
reference with 200-draw bootstrap errors (2-4 and 3-5 s), mirror (noise) stretch, per-episode stacks, pairwise doublet
graph, differential step injection, precision vs number of repeats. eps = -dv/v.

    python s10_measure.py [--stations NET.STA,...] [--workers 8]
"""
from __future__ import annotations
import os, sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S

EPS = np.linspace(-0.01, 0.01, 201)
DE = EPS[1] - EPS[0]
WINS = {"24": (2.0, 4.0), "35": (3.0, 5.0)}
MIR0 = -12.5                       # mirror start (inside the stored span, clear of the 5 % taper)
NBOOT = 200
INJ = (0.001, -0.001, 0.003, -0.003)
NS_PREC = (10, 20, 50, 100, 200)
RNG_SEED = 20261011
RAIL = 0.0099                      # |eps| at the +-1 % grid edge = no correlation peak in range -> not a measurement (Rev. 2)


# ------------------------------------------------------------------------------------------ stretch core
class Ref:
    """reference trace (on S.TAU) prepared for the stretch bank of one window."""
    def __init__(self, r, tS, win, mirror_end=None):
        from scipy.interpolate import CubicSpline
        w0, w1 = win
        if mirror_end is None:
            tr = S.TAU - tS
            self.idx = np.flatnonzero((tr >= w0) & (tr <= w1))
            x = tr[self.idx]
            sp = CubicSpline(tr, r)
        else:
            # time-reversed noise segment of the same length, ending at mirror_end, placed at lapse w0..w1
            iend = int(round((tS + mirror_end - S.T0W) * S.FS))
            K = int(round((w1 - w0) * S.FS)) + 1
            M = 30
            self.idx = iend - np.arange(K)
            kk = np.arange(-M, K + M)
            ii = np.clip(iend - kk, 0, S.NW - 1)
            sp = CubicSpline(w0 + kk / S.FS, r[ii])
            x = w0 + np.arange(K) / S.FS
        B = sp(x[None, :] / (1.0 + EPS[:, None]))       # ref(tau/(1+eps)): cur delayed by (1+eps) <=> eps = -dv/v
        B = B - B.mean(1, keepdims=True)
        self.B = B / np.linalg.norm(B, axis=1, keepdims=True)

    def measure(self, C):
        """C (m, NW) current traces -> eps (m,), cc (m,)"""
        X = np.atleast_2d(C)[:, self.idx].astype(np.float64)
        X = X - X.mean(1, keepdims=True)
        nx = np.linalg.norm(X, axis=1, keepdims=True); nx[nx == 0] = np.nan
        cc = (X / nx) @ self.B.T                                  # (m, n_eps)
        i = np.nanargmax(np.nan_to_num(cc, nan=-9), axis=1)
        ii = np.clip(i, 1, len(EPS) - 2)
        r = np.arange(len(i))
        y0, y1, y2 = cc[r, ii - 1], cc[r, ii], cc[r, ii + 1]
        d2 = y0 - 2 * y1 + y2
        delta = np.where((d2 < 0) & (i == ii), 0.5 * (y0 - y2) / np.where(d2 == 0, -1, d2), 0.0)
        eps = EPS[i] + np.clip(delta, -1, 1) * DE
        return eps, cc[r, i]


def stretch_trace(x, tS, a):
    """inject: cur(tau) = x(tau/(1+a)) about S -> reads eps = +a"""
    from scipy.interpolate import CubicSpline
    tr = S.TAU - tS
    return CubicSpline(tr, x)(tr / (1.0 + a))


# ------------------------------------------------------------------------------------------ helpers
def env(x):
    from scipy.signal import hilbert
    e = np.abs(hilbert(x))
    k = int(0.25 * S.FS)
    return np.convolve(e, np.ones(k) / k, "same")


def pick_S(G, approx=0.0):
    e = np.zeros(S.NW)
    for c in range(3):
        if np.isfinite(G[c]).all():
            e += env(G[c]) ** 2
    m = np.flatnonzero((S.TAU >= approx - 2.0) & (S.TAU <= approx + 2.0))
    return float(S.TAU[m[np.argmax(e[m])]])


def wmean(W, w):
    """W (n, NW) with possible NaN rows; weights w -> weighted mean over finite rows"""
    ok = np.isfinite(W[:, 0])
    if not ok.any():
        return np.full(S.NW, np.nan), 0.0
    return (w[ok, None] * W[ok]).sum(0) / w[ok].sum(), float(w[ok].sum())


def realign(W, G, tS, lagmax=0.5):
    """shift each window by the 3-C cc lag (+-lagmax) of tau in [tS-1, tS+2] against G; returns shifted W, lags, cc"""
    n = W.shape[0]
    L = int(round(lagmax * S.FS))
    seg = np.flatnonzero((S.TAU >= tS - 1.0) & (S.TAU <= tS + 2.0))
    lags = np.arange(-L, L + 1)
    num = np.zeros((n, len(lags))); den_w = np.zeros((n, len(lags))); den_g = 0.0
    for c in range(3):
        if not np.isfinite(G[c]).all():
            continue
        g = G[c, seg]; g = g - g.mean(); den_g += g @ g
        Wc = np.nan_to_num(W[:, c])
        for k, l in enumerate(lags):
            x = Wc[:, seg + l]
            x = x - x.mean(1, keepdims=True)
            num[:, k] += x @ g; den_w[:, k] += (x * x).sum(1)
    cc = num / np.sqrt(np.maximum(den_w, 1e-30) * max(den_g, 1e-30))
    i = np.argmax(cc, 1); ii = np.clip(i, 1, len(lags) - 2); r = np.arange(n)
    y0, y1, y2 = cc[r, ii - 1], cc[r, ii], cc[r, ii + 1]; d2 = y0 - 2 * y1 + y2
    delta = np.where((d2 < 0) & (i == ii), 0.5 * (y0 - y2) / np.where(d2 == 0, -1, d2), 0.0)
    lag = (lags[i] + np.clip(delta, -1, 1)) / S.FS            # window content at tau+lag matches G at tau
    f = np.fft.rfftfreq(S.NW, 1 / S.FS)
    ph = np.exp(2j * np.pi * f[None, :] * lag[:, None])        # x(tau + lag)
    out = W.copy()
    for c in range(3):
        ok = np.isfinite(W[:, c, 0])
        if ok.any():
            out[ok, c] = np.fft.irfft(np.fft.rfft(W[ok, c].astype(np.float64), axis=1) * ph[ok], n=S.NW, axis=1)
    return out.astype(np.float32), lag, cc[r, i]


def boot_stacks(Wd, w, nb, rng):
    """Wd (n, NW) finite rows, weights w -> (nb, NW) bootstrap weighted stacks"""
    n = Wd.shape[0]
    cnt = np.stack([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(nb)]).astype(np.float64)
    ww = cnt * w[None, :]
    return (ww @ Wd) / ww.sum(1, keepdims=True)


def clusters(dates, eps_tab):
    """episode id per date (r23 episodes +-2 d), others clustered by 7-d gaps"""
    d = pd.to_datetime(pd.Series(dates)).to_numpy()
    lab = np.full(len(d), "", dtype=object)
    for r in eps_tab.itertuples():
        m = (d >= np.datetime64(r.start) - np.timedelta64(2, "D")) & (d <= np.datetime64(r.end) + np.timedelta64(2, "D"))
        lab[m & (lab == "")] = r.episode
    o = np.argsort(d); prev = None; cid = 0
    for j in o:
        if lab[j] == "":
            if prev is None or (d[j] - prev) > np.timedelta64(7, "D"):
                cid += 1
            lab[j] = f"c{pd.Timestamp(d[j]).strftime('%Y-%m')}_{cid}"
            prev = d[j]
    return lab


# ------------------------------------------------------------------------------------------ graph
def graph_solve(n, pi, pj, d, sig):
    """v_j - v_i = d, weights 1/sig^2, 4 Huber IRLS passes, zero mean; -> v, sigma_v, chi2_red, n_pairs"""
    if len(pi) == 0:
        return np.full(n, np.nan), np.full(n, np.nan), np.nan
    A = np.zeros((len(pi), n)); A[np.arange(len(pi)), pi] = -1; A[np.arange(len(pi)), pj] = 1
    w0 = 1.0 / np.maximum(sig, 1e-6) ** 2
    w = w0.copy(); v = np.zeros(n)
    for it in range(5):
        N = A.T @ (A * w[:, None]); b = A.T @ (w * d)
        Ni = np.linalg.pinv(N, rcond=1e-10)
        v = Ni @ b
        r = (A @ v - d) / np.maximum(sig, 1e-6)
        c = 1.345
        w = w0 * np.where(np.abs(r) <= c, 1.0, c / np.maximum(np.abs(r), 1e-12))
    deg = np.bincount(np.r_[pi, pj], minlength=n)
    chi2 = float(np.sum(w * (A @ v - d) ** 2) / max(len(pi) - n + 1, 1))
    sv = np.sqrt(np.clip(np.diag(Ni), 0, None)) * np.sqrt(max(chi2, 1.0))
    v = v - v[deg > 0].mean()
    v[deg == 0] = np.nan; sv[deg == 0] = np.nan
    return v, sv, chi2


def pairs_measure(St, tS, win, sig):
    """all pairs i<j of day stacks St (n, NW): eps_ij with day i as reference"""
    n = St.shape[0]
    pi, pj, dd, ss, cc = [], [], [], [], []
    for i in range(n - 1):
        R = Ref(St[i], tS, win)
        e, c = R.measure(St[i + 1:])
        for k, j in enumerate(range(i + 1, n)):
            pi.append(i); pj.append(j); dd.append(e[k]); cc.append(c[k]); ss.append(np.hypot(sig[i], sig[j]))
    return np.array(pi, int), np.array(pj, int), np.array(dd), np.array(ss), np.array(cc)


# ------------------------------------------------------------------------------------------ one series
def run_series(fam, ns, ep, X, W, sp, X2, W2, EPI):
    import zlib
    rng = np.random.default_rng(RNG_SEED + fam * 1000 + zlib.crc32(f'{ns}|{ep}'.encode()) % 997)
    data = (X.series == "data").to_numpy()
    okZ = np.isfinite(W[:, 0, 0])
    Xd, Wd = X[data & okZ].reset_index(drop=True), W[data & okZ]
    Xc, Wc = X[~data & okZ].reset_index(drop=True), W[~data & okZ]
    out = dict(fam=fam, netsta=ns, epoch=ep, n_data=len(Xd), n_c1=len(Xc))
    if len(Xd) < 6 or Xd.date.nunique() < 3:
        return out, None
    P = Xd.P.to_numpy()
    G0 = np.stack([wmean(Wd[:, c], P)[0] for c in range(3)])
    tS0 = pick_S(G0, 0.0)
    Wd, lag_d, cc_d = realign(Wd, G0, tS0)
    if len(Xc):
        Wc, lag_c, cc_c = realign(Wc, G0, tS0)
    ok2 = np.isfinite(W2[:, 0, 0]) if len(X2) else np.zeros(0, bool)
    X2, W2 = X2[ok2].reset_index(drop=True), W2[ok2]
    if len(X2):
        W2, lag_2, cc_2 = realign(W2, G0, tS0)
    G1 = np.stack([wmean(Wd[:, c], P)[0] for c in range(3)])
    tS = pick_S(G1, tS0)
    SP = sp if np.isfinite(sp) else 6.0
    mend = -(SP + 1.0)
    tr = S.TAU - tS
    coda = {k: np.sqrt(np.mean(G1[0][(tr >= a) & (tr <= b)] ** 2)) for k, (a, b) in WINS.items()}
    mirr = np.sqrt(np.mean(G1[0][(tr >= MIR0) & (S.TAU >= -13.5) & (tr <= mend)] ** 2))
    Spk = np.sqrt(np.mean(G1[0][(tr >= -0.5) & (tr <= 1.0)] ** 2))
    out.update(tS_pick0=tS0, tS_pick=tS, sp=SP, rms_coda24=coda["24"], rms_coda35=coda["35"], rms_mirror=mirr,
               rms_S=Spk, causal24=coda["24"] / mirr, causal35=coda["35"] / mirr, lag_med_abs=float(np.median(np.abs(lag_d))),
               realign_cc_med=float(np.median(cc_d)), realign_cc_ctrl_med=float(np.median(cc_c)) if len(Xc) else np.nan,
               realign_cc_ctrl2_med=float(np.median(cc_2)) if len(X2) else np.nan,
               n_days=int(Xd.date.nunique()))
    # ---------------- day stacks
    days = np.array(sorted(Xd.date.unique()))
    Gw = (P[:, None] * Wd[:, 0]).sum(0); Ptot = P.sum()
    hi = Xd.hiP.to_numpy()
    rows, St, Wsum = [], [], []
    for d in days:
        m = (Xd.date == d).to_numpy()
        w = P[m]; Zd = Wd[m, 0].astype(np.float64)
        Sd = (w @ Zd) / w.sum()
        ref = (Gw - w @ Zd) / (Ptot - w.sum())
        R = {k: Ref(ref, tS, win) for k, win in WINS.items()}
        RM = {k: Ref(ref, tS, win, mirror_end=mend) for k, win in WINS.items()}
        r = dict(fam=fam, netsta=ns, epoch=ep, date=d, n=int(m.sum()), sumP=float(w.sum()), n_hi=int((m & hi).sum()))
        Bs = boot_stacks(Zd, w, NBOOT, rng) if m.sum() >= 2 else None
        for k in WINS:
            e, c = R[k].measure(Sd)
            r[f"eps{k}"], r[f"cc{k}"] = float(e[0]), float(c[0])
            r[f"sig{k}"] = float(np.std(R[k].measure(Bs)[0], ddof=1)) if Bs is not None else np.nan
            em, cm = RM[k].measure(Sd)
            r[f"mir{k}"], r[f"mircc{k}"] = float(em[0]), float(cm[0])
        # P >= 0.75 (era 0.90) equal-weight
        mh = m & hi
        if mh.sum() >= 3:
            Zh = Wd[mh, 0].astype(np.float64); Sh = Zh.mean(0); Bh = boot_stacks(Zh, np.ones(mh.sum()), NBOOT, rng)
            for k in WINS:
                e, c = R[k].measure(Sh); r[f"hi_eps{k}"] = float(e[0]); r[f"hi_sig{k}"] = float(np.std(R[k].measure(Bh)[0], ddof=1))
        # scrambled control (same day, c1 rows, equal weight, <= n data)
        if len(Xc):
            mc = np.flatnonzero((Xc.date == d).to_numpy())
            if len(mc) > m.sum():
                mc = np.sort(rng.choice(mc, m.sum(), replace=False))
            r["n_ctrl"] = int(len(mc))
            if len(mc) >= 3:
                Zc = Wc[mc, 0].astype(np.float64); Sc = Zc.mean(0); Bc = boot_stacks(Zc, np.ones(len(mc)), NBOOT, rng)
                for k in WINS:
                    e, c = R[k].measure(Sc); r[f"c_eps{k}"], r[f"c_cc{k}"] = float(e[0]), float(c[0])
                    r[f"c_sig{k}"] = float(np.std(R[k].measure(Bc)[0], ddof=1))
        else:
            r["n_ctrl"] = 0
        # Revision 1: frame-time scrambled control (same c1 rows read at t + offset, no site shift)
        if len(X2):
            mc = np.flatnonzero((X2.date == d).to_numpy())
            if len(mc) > m.sum():
                mc = np.sort(rng.choice(mc, m.sum(), replace=False))
            r["n_ctrl2"] = int(len(mc))
            if len(mc) >= 3:
                Zc = W2[mc, 0].astype(np.float64); Sc = Zc.mean(0); Bc = boot_stacks(Zc, np.ones(len(mc)), NBOOT, rng)
                for k in WINS:
                    e, c = R[k].measure(Sc); r[f"c2_eps{k}"], r[f"c2_cc{k}"] = float(e[0]), float(c[0])
                    r[f"c2_sig{k}"] = float(np.std(R[k].measure(Bc)[0], ddof=1))
        # horizontals (2-4 s only, 100 draws)
        for c, nm in ((1, "N"), (2, "E")):
            okc = m & np.isfinite(Wd[:, c, 0])
            if okc.sum() >= 3 and np.isfinite(G1[c]).all():
                okall = np.isfinite(Wd[:, c, 0])
                Gc = (P[okall, None] * Wd[okall, c]).sum(0); refc = (Gc - P[okc] @ Wd[okc, c]) / (P[okall].sum() - P[okc].sum())
                Rc = Ref(refc, tS, WINS["24"]); Hc = Wd[okc, c].astype(np.float64)
                e, cc_ = Rc.measure((P[okc] @ Hc) / P[okc].sum())
                r[f"{nm}_eps24"], r[f"{nm}_cc24"] = float(e[0]), float(cc_[0])
                r[f"{nm}_sig24"] = float(np.std(Rc.measure(boot_stacks(Hc, P[okc], 100, rng))[0], ddof=1))
        rows.append(r); St.append(Sd); Wsum.append(w.sum())
    D = pd.DataFrame(rows)
    St = np.array(St); Wsum = np.array(Wsum)
    D["cluster"] = clusters(D.date, EPI)
    # Revision 2: railed values (no cc peak inside +-1 %) are kept as *_raw and set to NaN
    for pre in ("", "hi_", "c_", "c2_", "N_", "E_"):
        for k in WINS:
            e = f"{pre}eps{k}"
            if e in D:
                D[f"{e}_raw"] = D[e]; rl = D[e].abs() >= RAIL
                D[f"{pre}rail{k}"] = rl.where(D[e].notna())
                D.loc[rl, e] = np.nan
                sg = f"{pre}sig{k}"
                if sg in D:
                    D.loc[rl, sg] = np.nan
    for k in WINS:
        D[f"mirrail{k}"] = D[f"mir{k}"].abs() >= RAIL
    # ---------------- per-episode stacks (leave-episode-out reference)
    erows = []
    Xd_cl = dict(zip(D.date, D.cluster))
    cl_row = np.array([Xd_cl[d] for d in Xd.date])
    for cl in pd.unique(D.cluster):
        m = cl_row == cl
        if m.sum() < 3:
            continue
        w = P[m]; Zd = Wd[m, 0].astype(np.float64)
        Sd = (w @ Zd) / w.sum(); ref = (Gw - w @ Zd) / max(Ptot - w.sum(), 1e-9)
        if Ptot - w.sum() <= 0:
            continue
        r = dict(fam=fam, netsta=ns, epoch=ep, cluster=cl, date_mid=sorted(Xd.date[m])[int(m.sum() // 2)],
                 n=int(m.sum()), sumP=float(w.sum()))
        Bs = boot_stacks(Zd, w, NBOOT, rng)
        for k, win in WINS.items():
            R = Ref(ref, tS, win); e, c = R.measure(Sd)
            r[f"eps{k}"], r[f"cc{k}"], r[f"sig{k}"] = float(e[0]), float(c[0]), float(np.std(R.measure(Bs)[0], ddof=1))
        erows.append(r)
    E = pd.DataFrame(erows)
    for k in WINS:
        if len(E):
            E[f"eps{k}_raw"] = E[f"eps{k}"]; rl = E[f"eps{k}"].abs() >= RAIL; E[f"rail{k}"] = rl
            E.loc[rl, [f"eps{k}", f"sig{k}"]] = np.nan
    # ---------------- pairwise graph (Z, both windows)
    for k, win in WINS.items():
        good = D[f"eps{k}"].notna().to_numpy()
        sig = D[f"sig{k}"].fillna(D[f"sig{k}"].median()).to_numpy()
        pi, pj, dd, ss, cc = pairs_measure(St, tS, win, sig)
        keep = good[pi] & good[pj] & (np.abs(dd) < RAIL)
        out[f"graph_pairs_railed_{k}"] = int((good[pi] & good[pj] & (np.abs(dd) >= RAIL)).sum())
        pi, pj, dd, ss = pi[keep], pj[keep], dd[keep], ss[keep]
        v, sv, chi2 = graph_solve(len(D), pi, pj, dd, ss)
        D[f"g_eps{k}"], D[f"g_sig{k}"] = v, sv
        out[f"graph_chi2_{k}"] = chi2; out[f"graph_pairs_{k}"] = int(len(pi))
    # ---------------- differential step injection (natural vs injected, Z)
    n = len(D)
    inj = []
    if n >= 6:
        b = n // 2
        pre = np.arange(n) < b
        cl = D.cluster.to_numpy()
        last_pre = pre & (cl == cl[b - 1])
        Gpre = (Wsum[pre, None] * St[pre]).sum(0); Gpost = (Wsum[~pre, None] * St[~pre]).sum(0)
        for k, win in WINS.items():
            e_, c_ = Ref(Gpre, tS, win).measure(Gpost); out[f"prepost_cc{k}"] = float(c_[0]); out[f"prepost_eps{k}"] = float(e_[0])
        early_pre = pre & ~last_pre
        for k, win in WINS.items():
            good = D[f"eps{k}"].notna().to_numpy()
            sig = D[f"sig{k}"].fillna(D[f"sig{k}"].median()).to_numpy(); wts = np.where(good, 1 / np.maximum(sig, 1e-6) ** 2, 0.0)

            def daily(Stx):
                G = (Wsum[:, None] * Stx).sum(0); out_ = np.zeros(n)
                for i in range(n):
                    ref = (G - Wsum[i] * Stx[i]) / (Wsum.sum() - Wsum[i])
                    out_[i] = Ref(ref, tS, win).measure(Stx[i])[0][0]
                return out_

            def contrasts(v):
                wm = lambda m: np.sum(wts[m] * v[m]) / np.sum(wts[m]) if m.any() else np.nan
                return wm(~pre) - wm(pre), (wm(last_pre) - wm(early_pre)) if early_pre.any() else np.nan
            nat_d = daily(St)
            pi, pj, dd, ss, cc = pairs_measure(St, tS, win, sig)
            kp = good[pi] & good[pj] & (np.abs(dd) < RAIL)
            nat_g = graph_solve(n, pi[kp], pj[kp], dd[kp], ss[kp])[0]
            for a in INJ:
                Sti = St.copy()
                for i in np.flatnonzero(~pre):
                    Sti[i] = stretch_trace(St[i], tS, a)
                inj_d = daily(Sti)
                pi, pj, dd, ss, cc = pairs_measure(Sti, tS, win, sig)
                inj_g = graph_solve(n, pi[kp], pj[kp], dd[kp], ss[kp])[0]
                for est, vn, vi in (("daily", nat_d, inj_d), ("graph", nat_g, inj_g)):
                    cn, ln = contrasts(np.nan_to_num(vn)); ci, li = contrasts(np.nan_to_num(vi))
                    inj.append(dict(fam=fam, netsta=ns, epoch=ep, win=k, est=est, a=a, boundary=D.date[b],
                                    recovery=(ci - cn) / a, leakage=(li - ln) / a if np.isfinite(li) else np.nan))
    INJT = pd.DataFrame(inj)
    # ---------------- precision vs number of repeats (busiest days, disjoint halves, equal weights)
    prec = []
    cnt = D.set_index("date").n.sort_values(ascending=False)
    for d in cnt.index[:3]:
        cl = D.set_index("date").cluster[d]
        pool_days = [d] + [x for x in D.date[D.cluster == cl] if x != d]
        for nn in NS_PREC:
            m = np.zeros(len(Xd), bool); used = []
            for x in pool_days:
                if m.sum() >= 2 * nn:
                    break
                m |= (Xd.date == x).to_numpy(); used.append(x)
            if m.sum() < 2 * nn:
                continue
            idx = np.flatnonzero(m)
            Zp = Wd[idx, 0].astype(np.float64)
            ref = (Gw - P[idx] @ Zp) / (Ptot - P[idx].sum())
            if Ptot - P[idx].sum() <= 0:
                continue
            for k, win in WINS.items():
                R = Ref(ref, tS, win); diffs = []; railed = []
                for _ in range(50):
                    pp = rng.permutation(len(idx))
                    A_, B_ = Zp[pp[:nn]].mean(0), Zp[pp[nn:2 * nn]].mean(0)
                    e = R.measure(np.stack([A_, B_]))[0]; diffs.append(e[0] - e[1]); railed.append(bool(np.any(np.abs(e) >= RAIL)))
                diffs = np.array(diffs)
                prec.append(dict(fam=fam, netsta=ns, epoch=ep, day=d, n_pool_days=len(used), n=nn, win=k,
                                 sigma=float(np.std(diffs, ddof=1) / np.sqrt(2)), rail_frac=float(np.mean(railed)),
                                 sigma_unrailed=float(np.std(diffs[~np.array(railed)], ddof=1) / np.sqrt(2))
                                 if (~np.array(railed)).sum() > 5 else np.nan))
    PREC = pd.DataFrame(prec)
    stacks = dict(G1=G1, tS=tS, sp=SP, n=len(Xd))
    return out, (D, E, INJT, PREC, stacks)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--stations", default=""); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--tag", default="")
    ap.add_argument("--rs", default="", help="alternative row_station table (KELB+KLNB merged series)")
    a = ap.parse_args()
    S.ensure(); log = S.Log(f"{S.LOGD}/measure{a.tag}.log")
    X = pd.read_parquet(f"{S.TB}/rows.parquet")
    RS = pd.read_parquet(a.rs or f"{S.TB}/row_station.parquet")
    O = pd.read_csv(f"{S.TB}/offsets_sp.csv").set_index(["fam", "netsta"])
    EPI = pd.read_csv(f"{S.ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")
    stas = a.stations.split(",") if a.stations else sorted(RS.netsta.unique())
    jobs = []
    for ns in stas:
        p = f"{S.WIN}/{ns}.npz"
        if not os.path.exists(p):
            log(f"{ns}: no windows"); continue
        jobs.append(ns)
    import multiprocessing as mp
    res_out, res_D, res_E, res_I, res_P, res_S = [], [], [], [], [], {}

    def series_iter():
        for ns in jobs:
            z = np.load(f"{S.WIN}/{ns}.npz")
            pos = {r: i for i, r in enumerate(z["row_id"])}
            p2 = f"{S.WIN}/{ns}__c1f.npz"
            z2 = np.load(p2) if os.path.exists(p2) else None
            pos2 = {r: i for i, r in enumerate(z2["row_id"])} if z2 is not None else {}
            R = RS[RS.netsta == ns].merge(X, on=["row_id", "date"])
            for (fam, ep), g in R.groupby(["fam", "epoch"]):
                g = g[g.row_id.isin(pos)].reset_index(drop=True)
                W = z["W"][[pos[r] for r in g.row_id]]
                g2 = g[g.series == "c1"].reset_index(drop=True)
                if z2 is not None:
                    g2 = g2[g2.row_id.isin(pos2)].reset_index(drop=True)
                    W2 = z2["W"][[pos2[r] for r in g2.row_id]]
                else:
                    g2 = g2.iloc[:0]; W2 = np.zeros((0, 3, S.NW), np.float32)
                yield (int(fam), ns, ep, g, W, float(O.sp.get((fam, ns), np.nan)), g2, W2, EPI)

    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers) as p:
        for out, res in p.imap_unordered(_run, series_iter()):
            res_out.append(out)
            if res is not None:
                D, E, I, Pr, st = res
                res_D.append(D); res_E.append(E); res_I.append(I); res_P.append(Pr)
                res_S[f"{out['fam']}|{out['netsta']}|{out['epoch']}"] = st
            log(f"  {out['fam']} {out['netsta']} {out['epoch']}: n_data {out['n_data']} days {out.get('n_days', 0)} "
                f"causal24 {out.get('causal24', np.nan):.2f} {time.time() - t0:.0f}s")
    sfx = a.tag
    pd.DataFrame(res_out).to_csv(f"{S.TB}/series{sfx}.csv", index=False)
    pd.concat(res_D, ignore_index=True).to_parquet(f"{S.TB}/daily{sfx}.parquet", index=False)
    pd.concat(res_E, ignore_index=True).to_parquet(f"{S.TB}/episodes{sfx}.parquet", index=False)
    pd.concat(res_I, ignore_index=True).to_csv(f"{S.TB}/injection{sfx}.csv", index=False)
    pd.concat(res_P, ignore_index=True).to_csv(f"{S.TB}/precision{sfx}.csv", index=False)
    np.savez(f"{S.TB}/grandstacks{sfx}.npz", keys=np.array(list(res_S)), G1=np.array([v["G1"] for v in res_S.values()]),
             tS=np.array([v["tS"] for v in res_S.values()]), sp=np.array([v["sp"] for v in res_S.values()]),
             n=np.array([v["n"] for v in res_S.values()]), tau=S.TAU)
    log(f"measure done {time.time() - t0:.0f}s")
    open(f"{S.OUT}/MEASURE{sfx}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


def _run(args):
    try:
        return run_series(*args)
    except Exception as e:
        import traceback
        return dict(fam=args[0], netsta=args[1], epoch=args[2], n_data=-1, error=f"{e} {traceback.format_exc()[-600:]}"), None


if __name__ == "__main__":
    main()
