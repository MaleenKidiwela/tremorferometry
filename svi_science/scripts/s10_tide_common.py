"""s10_tide_common -- shared machinery for the S10 tide / sub-daily follow-up (rules: results/s10_dvv/tide_subday/PLAN.md).

Series preparation reuses s10_measure (G0 -> pick_S -> realign -> G1 -> pick_S) unchanged and caches the realigned Z
windows. The stretch measurement is the s10_measure.Ref estimator; `bank_u` is Ref.__init__ without the final row
normalisation (so references that are weighted sums of day stacks can be built exactly by linearity), and
`measure_bank` is Ref.measure vectorised over many references at once. Both are checked against Ref in selftest().
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_measure as M

OUT = f"{S.ROOT}/svi_science/results/s10_dvv/tide_subday"
TB, FIG, LOGD, CACHE = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/logs", f"{OUT}/cache"
PILOT_TB = f"{S.ROOT}/svi_science/results/s10_dvv/pilot/tables"
PILOT_WIN = f"{S.ROOT}/svi_science/results/s10_dvv/pilot/win"
WINS = M.WINS
EPS = M.EPS
DE = M.DE
RAIL = M.RAIL
SITE_OF = {"PB.B009": "PGC", "PB.B010": "PGC", "PB.B011": "PGC", "CN.PGC": "PGC"}
FAMS = (10, 87, 17)


def site(ns):
    return SITE_OF.get(ns, ns)


def ensure():
    for d in (TB, FIG, LOGD, CACHE):
        os.makedirs(d, exist_ok=True)


# ------------------------------------------------------------------ estimator (exactly Ref, vectorised)
def bank_u(r, tS, win):
    """unnormalised, row-demeaned stretch bank of reference r (Ref.__init__ before the row normalisation)."""
    from scipy.interpolate import CubicSpline
    w0, w1 = win
    tr = S.TAU - tS
    idx = np.flatnonzero((tr >= w0) & (tr <= w1))
    x = tr[idx]
    B = CubicSpline(tr, r)(x[None, :] / (1.0 + EPS[:, None]))
    return B - B.mean(1, keepdims=True), idx


def seg_idx(tS, win):
    tr = S.TAU - tS
    return np.flatnonzero((tr >= win[0]) & (tr <= win[1]))


def norm_bank(Bu):
    """row-normalise a bank (..., n_eps, K)."""
    return Bu / np.linalg.norm(Bu, axis=-1, keepdims=True)


def measure_bank(X, Bn):
    """X (..., m, K) segments of current stacks; Bn (..., n_eps, K) normalised banks (broadcast over leading dims).
    -> eps (..., m), cc (..., m); NaN where the stack is empty. Same arithmetic as s10_measure.Ref.measure."""
    X = X - X.mean(-1, keepdims=True)
    nx = np.linalg.norm(X, axis=-1, keepdims=True)
    bad = (nx[..., 0] == 0) | ~np.isfinite(nx[..., 0])
    nx = np.where(nx == 0, np.nan, nx)
    cc = np.matmul(X / nx, np.swapaxes(Bn, -1, -2))            # (..., m, n_eps)
    cc0 = np.nan_to_num(cc, nan=-9)
    i = np.argmax(cc0, axis=-1)
    ii = np.clip(i, 1, len(EPS) - 2)
    y0 = np.take_along_axis(cc, (ii - 1)[..., None], -1)[..., 0]
    y1 = np.take_along_axis(cc, ii[..., None], -1)[..., 0]
    y2 = np.take_along_axis(cc, (ii + 1)[..., None], -1)[..., 0]
    d2 = y0 - 2 * y1 + y2
    delta = np.where((d2 < 0) & (i == ii), 0.5 * (y0 - y2) / np.where(d2 == 0, -1, d2), 0.0)
    eps = EPS[i] + np.clip(delta, -1, 1) * DE
    ccm = np.take_along_axis(cc, i[..., None], -1)[..., 0]
    eps = np.where(bad, np.nan, eps); ccm = np.where(bad, np.nan, ccm)
    return eps, ccm


def unrail(eps, cc):
    """pilot Revision 2: |eps| >= 0.99 % -> not a measurement."""
    rl = np.abs(eps) >= RAIL
    return np.where(rl, np.nan, eps), np.where(rl | ~np.isfinite(eps), np.nan, cc)


def selftest():
    rng = np.random.default_rng(1)
    r = np.convolve(rng.standard_normal(S.NW), np.ones(5) / 5, "same")
    c = np.stack([M.stretch_trace(r, 0.3, a) + 0.3 * rng.standard_normal(S.NW) for a in (0.002, -0.004, 0.0)])
    out = []
    for k, win in WINS.items():
        R = M.Ref(r, 0.3, win)
        e0, c0 = R.measure(c)
        Bu, idx = bank_u(r, 0.3, win)
        e1, c1 = measure_bank(c[:, idx], norm_bank(Bu))
        out.append(float(np.max(np.abs(e0 - e1))) + float(np.max(np.abs(c0 - c1))))
    return max(out)


# ------------------------------------------------------------------ series preparation (s10_measure.run_series, first part)
def series_list():
    RS = pd.read_parquet(f"{PILOT_TB}/row_station.parquet")
    g = RS.groupby(["netsta", "epoch"]).size().reset_index()
    return [(ns, ep) for ns, ep in zip(g.netsta, g.epoch)]


def prepare_station(ns, X, RS, force=False):
    """all (family, epoch) series of one station -> cache npz per series. Returns list of cache paths."""
    out = []
    p = f"{PILOT_WIN}/{ns}.npz"
    if not os.path.exists(p):
        return out
    z = np.load(p)
    pos = {r: i for i, r in enumerate(z["row_id"])}
    p2 = f"{PILOT_WIN}/{ns}__c1f.npz"
    z2 = np.load(p2) if os.path.exists(p2) else None
    pos2 = {r: i for i, r in enumerate(z2["row_id"])} if z2 is not None else {}
    R = RS[RS.netsta == ns].merge(X, on=["row_id", "date"])
    for (fam, ep), g in R.groupby(["fam", "epoch"]):
        fn = f"{CACHE}/{int(fam)}_{ns}_{ep}.npz"
        if os.path.exists(fn) and not force:
            out.append(fn); continue
        g = g[g.row_id.isin(pos)].reset_index(drop=True)
        gd = g[g.series == "data"].reset_index(drop=True)
        W = z["W"][[pos[r] for r in gd.row_id]]
        rms = z["rms_raw"][[pos[r] for r in gd.row_id]]
        ok = np.isfinite(W[:, 0, 0])
        gd, W, rms = gd[ok].reset_index(drop=True), W[ok], rms[ok]
        if len(gd) < 6 or gd.date.nunique() < 3:
            continue
        g2 = g[g.series == "c1"].reset_index(drop=True)
        if z2 is not None:
            g2 = g2[g2.row_id.isin(pos2)].reset_index(drop=True)
            W2 = z2["W"][[pos2[r] for r in g2.row_id]]
            ok2 = np.isfinite(W2[:, 0, 0]); g2, W2 = g2[ok2].reset_index(drop=True), W2[ok2]
        else:
            g2 = g2.iloc[:0]; W2 = np.zeros((0, 3, S.NW), np.float32)
        P = gd.P.to_numpy()
        G0 = np.stack([M.wmean(W[:, c], P)[0] for c in range(3)])
        tS0 = M.pick_S(G0, 0.0)
        Wr, lag, rcc = M.realign(W, G0, tS0)
        if len(g2):
            W2r, lag2, rcc2 = M.realign(W2, G0, tS0)
        else:
            W2r, lag2, rcc2 = W2, np.zeros(0), np.zeros(0)
        G1 = np.stack([M.wmean(Wr[:, c], P)[0] for c in range(3)])
        tS = M.pick_S(G1, tS0)
        tr = S.TAU - tS
        sw = (tr >= -0.5) & (tr <= 1.0)
        logamp = np.log10(np.sqrt(np.mean(Wr[:, 0][:, sw].astype(np.float64) ** 2, 1)) * rms[:, 0] + 1e-30)
        np.savez(fn, fam=int(fam), netsta=ns, epoch=ep, tS=tS, row_id=gd.row_id.to_numpy(), date=gd.date.to_numpy().astype(str),
                 P=P, Z=Wr[:, 0].astype(np.float32), lag=lag, rcc=rcc, logamp=logamp,
                 row_id2=g2.row_id.to_numpy(), date2=g2.date.to_numpy().astype(str), Z2=W2r[:, 0].astype(np.float32),
                 G1=G1[0])
        out.append(fn)
    return out
