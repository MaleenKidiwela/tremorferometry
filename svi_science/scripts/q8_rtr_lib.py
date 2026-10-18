"""Q8 RTR detector library (rules frozen in svi_science/results/q8_rtr/PLAN.md §2).

Per episode: events (t hours since window start, x along-strike km). Candidate lines x = x0 + v (t - t0) over
[t0, t0+D]; corridor half-width W. Score S = -log10 Poisson P(>= N | mu), mu = space/time-independence expectation
inside the line's bounding rectangle. Streak shape: N >= 5, >= 3 of 4 time quarters occupied, Spearman >= 0.5.
Greedy overlap removal (30 %). Direction class from a 4-km onset map.
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import maximum_filter

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q8_rtr"

DT = 1.0 / 6.0          # h
DX = 1.0                # km
W = 3.0                 # km corridor half-width
SPEEDS = [1, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20]
DURS = [1.5, 2, 3, 4, 5, 6, 8]
SPAN_MIN, SPAN_MAX = 8.0, 80.0
S_KEEP = 2.0
NMIN = 5
RHO_MIN = 0.5
OVERLAP = 0.30
# ---- debug iteration 1 (2026-10-08, after step-1 FAIL; PLAN.md §5 D1)
BG = "smooth"          # was "independence"
BG_HALF_H = 24.0
NFAM_MIN = 5
FAM_RHO_MIN = 0.6
FAM_MAX_SHARE = 0.5
COVER_MIN = 0.5
BOX_NMS = True
ONSET_BIN = 4.0
ONSET_MIN_EV = 5
BEHIND_MIN = 0.6
BACK_H = 12.0


def _sat(H):
    S = np.zeros((H.shape[0] + 1, H.shape[1] + 1))
    S[1:, 1:] = H.cumsum(0).cumsum(1)
    return S


def scan(t, x, T_end, xlo_q=None, xhi_q=None, bg=None):
    """Return DataFrame of local-maximum candidate lines with S >= S_KEEP (no shape test yet)."""
    if xlo_q is None:
        xlo_q, xhi_q = np.percentile(x, [2, 98])
    xg0 = np.floor(xlo_q - W - 2)
    nx = int(np.ceil((xhi_q + W + 2 - xg0) / DX)) + 1
    nt = int(np.ceil(T_end / DT)) + 1
    it = np.clip((t / DT).astype(int), 0, nt - 1)
    ix = np.clip(((x - xg0) / DX).astype(int), 0, nx - 1)
    H = np.zeros((nt, nx)); np.add.at(H, (it, ix), 1.0)
    S2 = _sat(H)                                   # (nt+1, nx+1)
    R = np.zeros((nt, nx + 1)); R[:, 1:] = H.cumsum(1)
    if bg is None:
        bg = BG
    if bg == "smooth":
        from scipy.ndimage import uniform_filter
        G = uniform_filter(H, size=(int(round(2 * BG_HALF_H / DT)) + 1, 3), mode="nearest")
        SG = _sat(G)
        RG = np.zeros((nt, nx + 1)); RG[:, 1:] = G.cumsum(1)
    wb = int(round(2 * W / DX))
    x0s = np.arange(np.ceil(xlo_q), np.floor(xhi_q) + 1, 1.0)
    out = []
    for D in DURS:
        nD = int(round(D / DT))
        i0 = np.arange(0, nt - nD)
        if len(i0) == 0:
            continue
        for sp in SPEEDS:
            span = sp * D
            if span < SPAN_MIN or span > SPAN_MAX:
                continue
            for sg in (1, -1):
                v = sg * sp
                xe = x0s + v * D
                ok = (xe >= xlo_q) & (xe <= xhi_q)
                if not ok.any():
                    continue
                X0 = x0s[ok]
                # rectangle bins
                A = np.floor((np.minimum(X0, X0 + v * D) - W - xg0) / DX).astype(int)
                B = np.floor((np.maximum(X0, X0 + v * D) + W - xg0) / DX).astype(int) + 1
                A = np.clip(A, 0, nx); B = np.clip(B, 0, nx)
                I0 = i0[:, None]; I1 = I0 + nD
                Nrect = S2[I1, B[None, :]] - S2[I0, B[None, :]] - S2[I1, A[None, :]] + S2[I0, A[None, :]]
                Npar = np.zeros((len(i0), len(X0))); mu = np.zeros_like(Npar)
                for k in range(nD):
                    xc = X0 + v * (k + 0.5) * DT
                    a = np.clip(np.floor((xc - W - xg0) / DX).astype(int), 0, nx)
                    b = np.clip(a + wb + 1, 0, nx)
                    j = I0 + k
                    Npar += R[j, b[None, :]] - R[j, a[None, :]]
                    nt_j = R[j, B[None, :]] - R[j, A[None, :]]
                    box = S2[I1, b[None, :]] - S2[I0, b[None, :]] - S2[I1, a[None, :]] + S2[I0, a[None, :]]
                    mu += nt_j * box
                mu = np.where(Nrect > 0, mu / np.maximum(Nrect, 1), 0.0)
                if bg == "smooth":
                    mu = np.zeros_like(Npar)
                    for k in range(nD):
                        xc = X0 + v * (k + 0.5) * DT
                        a = np.clip(np.floor((xc - W - xg0) / DX).astype(int), 0, nx)
                        b = np.clip(a + wb + 1, 0, nx)
                        mu += RG[I0 + k, b[None, :]] - RG[I0 + k, a[None, :]]
                Sc = np.where(Npar >= NMIN, -stats.poisson.logsf(Npar - 1, np.maximum(mu, 1e-9)) / np.log(10), 0.0)
                Sc = np.where(np.isfinite(Sc), Sc, 999.0)
                mx = maximum_filter(Sc, size=(int(round(2 / DT)) + 1, 7), mode="constant")
                sel = (Sc >= S_KEEP) & (Sc == mx)
                ii, jj = np.nonzero(sel)
                for a_, b_ in zip(ii, jj):
                    out.append((i0[a_] * DT, D, v, X0[b_], Npar[a_, b_], mu[a_, b_], Nrect[a_, b_], Sc[a_, b_]))
    return pd.DataFrame(out, columns=["t0", "D", "v", "x0", "N", "mu", "Nrect", "S"])


def onset_map(t, x):
    b = np.floor(x / ONSET_BIN).astype(int)
    df = pd.DataFrame(dict(b=b, t=t))
    g = df.groupby("b").t
    on = g.quantile(0.10)[g.size() >= ONSET_MIN_EV]
    return on.to_dict()


def _path_onsets(onset, xs):
    return np.array([onset.get(int(np.floor(z / ONSET_BIN)), np.nan) for z in xs])


def refine(cands, t, x, onset, unit=None):
    """Shape test + greedy overlap removal + direction class. t sorted ascending required."""
    order = np.argsort(t, kind="stable"); assert np.all(order == np.arange(len(t)))
    cands = cands.sort_values("S", ascending=False).reset_index(drop=True)
    owner = np.full(len(t), -1)
    kept = []
    for r in cands.itertuples(index=False):
        lo = np.searchsorted(t, r.t0); hi = np.searchsorted(t, r.t0 + r.D)
        tt = t[lo:hi]; xx = x[lo:hi]
        xc = r.x0 + r.v * (tt - r.t0)
        m = np.abs(xx - xc) <= W
        idx = np.arange(lo, hi)[m]
        n = len(idx)
        if n < NMIN:
            continue
        q = np.floor((t[idx] - r.t0) / r.D * 4).clip(0, 3).astype(int)
        nq = len(np.unique(q))
        if nq < 3:
            continue
        sgn = np.sign(r.v)
        rho = stats.spearmanr(t[idx], sgn * x[idx])[0] if np.ptp(x[idx]) > 0 else 0.0
        if not (rho >= RHO_MIN):
            continue
        if unit is not None and NFAM_MIN:
            uu = pd.Series(unit[idx])
            vc = uu.value_counts()
            if len(vc) < NFAM_MIN or vc.iloc[0] / n > FAM_MAX_SHARE:
                continue
            fm = pd.DataFrame(dict(u=unit[idx], t=t[idx], x=x[idx])).groupby("u").agg(t=("t", "median"), x=("x", "median"))
            frho = stats.spearmanr(fm.t, sgn * fm.x)[0]
            if not (frho >= FAM_RHO_MIN):
                continue
            pb = np.floor((x[idx] - min(r.x0, r.x0 + r.v * r.D)) / 2.0).astype(int)
            nb = int(np.ceil(abs(r.v) * r.D / 2.0))
            cover = len(np.unique(pb[(pb >= 0) & (pb < nb)])) / max(nb, 1)
            if cover < COVER_MIN:
                continue
        else:
            frho, cover, vc = np.nan, np.nan, pd.Series([n])
        if BOX_NMS and kept:
            lo_x, hi_x = sorted([r.x0, r.x0 + r.v * r.D])
            clash = False
            for kk in kept:
                ov_t = min(r.t0 + r.D, kk["t0"] + kk["D"]) - max(r.t0, kk["t0"])
                if ov_t <= 0.5 * min(r.D, kk["D"]):
                    continue
                a2, b2 = sorted([kk["x_line_start"], kk["x_line_end"]])
                ov_x = min(hi_x, b2) - max(lo_x, a2)
                if ov_x > 0.5 * min(hi_x - lo_x, b2 - a2):
                    clash = True; break
            if clash:
                continue
        own = owner[idx]
        if (own >= 0).any():
            shared = pd.Series(own[own >= 0]).value_counts()
            sizes = np.array([kept[k]["N_ev"] for k in shared.index])
            if np.max(shared.values / np.minimum(n, sizes)) > OVERLAP:
                continue
        k = len(kept)
        owner[idx[owner[idx] < 0]] = k
        te = t[idx]; xe = x[idx]
        # direction class
        xs = r.x0 + r.v * np.linspace(0, r.D, 41)
        po = _path_onsets(onset, xs)
        have = ~np.isnan(po)
        behind = float((po[have] < r.t0).mean()) if have.sum() else np.nan
        q1 = np.nanmedian(po[:11]) if np.isfinite(po[:11]).any() else np.nan
        q4 = np.nanmedian(po[-11:]) if np.isfinite(po[-11:]).any() else np.nan
        dq = q1 - q4      # >0: moving toward earlier-activated region
        if np.isfinite(behind) and behind >= BEHIND_MIN and np.isfinite(dq) and dq >= BACK_H:
            cls = "RTR"
        elif np.isfinite(behind) and behind >= BEHIND_MIN and np.isfinite(dq) and dq <= -BACK_H:
            cls = "forward"
        else:
            cls = "other"
        nf = len(idx) // 4 if len(idx) >= 8 else 2
        ts_slope = stats.theilslopes(xe, te)[0] if np.ptp(te) > 0 else np.nan
        rec = dict(t0=r.t0, D=r.D, v=r.v, x0=r.x0, N=r.N, mu=r.mu, Nrect=r.Nrect, S=r.S, N_ev=n, n_quart=nq, rho=rho,
                   t_first=te.min(), t_last=te.max(), x_start=float(np.median(xe[:nf])), x_end=float(np.median(xe[-nf:])),
                   span_ev=float(np.ptp(xe)), v_ev=ts_slope, behind=behind, onset_dq_h=dq, cls=cls,
                   x_line_start=r.x0, x_line_end=r.x0 + r.v * r.D, n_fam=len(vc), fam_rho=frho, cover=cover,
                   top_fam_share=float(vc.iloc[0] / n))
        rec["idx"] = idx
        kept.append(rec)
    return kept


def detect(t, x, T_end, xlo_q=None, xhi_q=None, onset=None, unit=None):
    """t (h, sorted), x (km). Returns list of kept streak dicts (with event index arrays)."""
    if onset is None:
        onset = onset_map(t, x)
    c = scan(t, x, T_end, xlo_q, xhi_q)
    if len(c) == 0:
        return []
    return refine(c, t, x, onset, unit)


# ------------------------------------------------------------------ controls
def control_c1(t, day, rng):
    """time redrawn uniformly within its own family-day (day = integer UTC day index relative to window, hours)."""
    return day * 24.0 + rng.uniform(0, 24.0, size=len(t))


def control_c2(x, unit, day, rng):
    """within each UTC day, family-day blocks swap positions (times kept)."""
    df = pd.DataFrame(dict(u=unit, d=day, x=x))
    pos = df.groupby("u").x.first()
    newx = np.array(x, float).copy()
    for d, g in df.groupby("d"):
        us = g.u.unique()
        perm = rng.permutation(us)
        mp = dict(zip(us, pos.loc[perm].values))
        newx[g.index.values] = g.u.map(mp).values
    return newx
