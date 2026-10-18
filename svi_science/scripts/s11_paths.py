"""s11_paths -- the two S11 data paths rebuilt from the day-stack store (tables/store/<NET.STA>.npz), with optional
waveform-level injection (each day stack stretched about tS by a day-specific eps before any stacking / SVD).

raw path (Merlin D3): episode stack = Ws-weighted mean of the day stacks (= P-weighted mean of the windows), leave-episode-
  out reference, joint 3C S-anchored stretch (s11_measure.Ref3). No bootstrap here (sigma comes from s11_measure).
svd path (user recipe): days with sum P >= 3, SVD-Wiener per component (s10_svd_measure.svdw, noise tau < tS - 8 s), joint
  3C stretch against the realigned epoch grand stack, days CC > 0.8 and |eps| < 0.99 %, CC-weighted episode mean;
  sigma = pooled within-episode day scatter / sqrt(n_eff).
eps = -dv/v."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s11_common as C
from s10_measure import stretch_trace
from s10_svd_measure import svdw
from s11_measure import Ref3

DAY_MIN_P = 3.0
CC_DAY = 0.8


def load_store(ns, tag=""):
    z = np.load(f"{C.TB}/store{tag}/{ns}.npz")
    keys = z["keys"]; nd = z["n_days"]; n2 = z["n_c2"]
    o1 = np.r_[0, np.cumsum(nd)]; o2 = np.r_[0, np.cumsum(n2)]
    St, Ws, days, Sc, Cn, c2d = z["St"], z["Ws"], z["days"], z["Sc"], z["Cn"], z["c2_day"]
    out = {}
    for i, k in enumerate(keys):
        out[str(k)] = dict(days=days[o1[i]:o1[i + 1]], St=St[o1[i]:o1[i + 1]], Ws=Ws[o1[i]:o1[i + 1]],
                           c2_day=c2d[o2[i]:o2[i + 1]], Sc=Sc[o2[i]:o2[i + 1]], Cn=Cn[o2[i]:o2[i + 1]],
                           G1=z["G1"][i].astype(np.float64), tS=float(z["tS"][i]), sp=float(z["sp"][i]))
    return out


def unit_map():
    A = pd.read_csv(f"{C.MTB}/active_days.csv")
    m = dict(zip(A.date, A.episode))
    pe = f"{C.TB}/ext/active_days.csv"                # Revision 2: dates of the added families (superset, same labels)
    if os.path.exists(pe):
        E = pd.read_csv(pe)
        for d, e in zip(E.date, E.episode):
            m.setdefault(d, e)
    return m


def inject(St, tS, amp):
    """St (n,3,NW) float, amp (n,) eps-units stretch per day -> stretched copy"""
    if amp is None or not np.any(amp):
        return St.astype(np.float64)
    X = St.astype(np.float64).copy()
    for k in np.flatnonzero(amp != 0):
        for c in range(3):
            if np.any(X[k, c]):
                X[k, c] = stretch_trace(X[k, c], tS, float(amp[k]))
    return X


def _comps(G1):
    return tuple(c for c in range(3) if np.isfinite(G1[c]).all())


def raw_cells(st, unit, amp=None, wins=None, control=False):
    """episode values of the raw path from the store (natural or injected)."""
    wins = wins or list(C.WINS)
    tS, G1 = st["tS"], st["G1"]; cs = _comps(G1)
    if control:
        days = st["days"][st["c2_day"]]; X = inject(st["Sc"], tS, None if amp is None else amp[st["c2_day"]]); Ws = st["Cn"]
        Xd = inject(st["St"], tS, amp); Wd = st["Ws"]
    else:
        days = st["days"]; X = inject(st["St"], tS, amp); Ws = st["Ws"]; Xd, Wd = X, Ws
    ep = np.array([unit.get(d, None) for d in days], dtype=object)
    epd = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
    Gs = np.einsum("nc,nct->ct", Wd, Xd); Wt = Wd.sum(0)
    rows = []
    for e in pd.unique(ep):
        if e is None:
            continue
        m = ep == e; md = epd == e
        we = Wd[md].sum(0)
        if np.any((Wt - we)[list(cs)] <= 0) or Ws[m][:, 0].sum() <= 0:
            continue
        Sd = np.full((3, S.NW), np.nan); ref = np.full((3, S.NW), np.nan)
        for c in cs:
            w = Ws[m, c]
            if w.sum() <= 0:
                continue
            Sd[c] = (w @ X[m, c]) / w.sum()
            ref[c] = (Gs[c] - Wd[md, c] @ Xd[md, c]) / (Wt[c] - Wd[md, c].sum())
        r = dict(episode=e)
        for wn in wins:
            e_, c_ = Ref3(ref, tS, C.WINS[wn]).measure(Sd[None], cs)
            r[f"eps3c_{wn}"], r[f"cc3c_{wn}"] = float(e_[0]), float(c_[0])
        rows.append(r)
    D = pd.DataFrame(rows, columns=["episode"] + [f"{p}3c_{wn}" for wn in wins for p in ("eps", "cc")])
    for wn in wins:
        if len(D):
            D.loc[D[f"eps3c_{wn}"].abs() >= C.RAIL, f"eps3c_{wn}"] = np.nan
    return D


def svd_days(st, unit, amp=None, wins=None, control=False):
    """day values of the SVD path (natural or injected) -> DataFrame(date, episode, eps/cc per window)"""
    wins = wins or list(C.WINS)
    tS, G1 = st["tS"], st["G1"]; cs = _comps(G1)
    if control:
        days = st["days"][st["c2_day"]]; Ws = st["Cn"]; X = inject(st["Sc"], tS, None if amp is None else amp[st["c2_day"]])
        keep = np.all(Ws[:, list(cs)] >= 3, 1)
        ref = G1
    else:
        days = st["days"]; Ws = st["Ws"]; X = inject(st["St"], tS, amp)
        keep = (Ws[:, 0] >= DAY_MIN_P) & np.all(Ws[:, list(cs)] > 0, 1)
        if amp is not None and np.any(amp):
            # grand stack rebuilt from the injected day stacks (it is the P-weighted mean of all realigned windows)
            Wt = Ws.sum(0); X0 = st["St"].astype(np.float64)
            ref = G1.copy()
            for c in cs:
                ref[c] = G1[c] + (Ws[:, c] @ (X[:, c] - X0[:, c])) / Wt[c]
        else:
            ref = G1
    if keep.sum() < 3:
        return pd.DataFrame()
    days = days[keep]; X = X[keep]
    noise = np.flatnonzero(S.TAU < tS - 8)
    Xf = np.full_like(X, np.nan)
    for c in cs:
        Xf[:, c] = svdw(X[:, c], noise)
    D = pd.DataFrame(dict(date=days, episode=[unit.get(d, None) for d in days]))
    for wn in wins:
        e_, c_ = Ref3(ref, tS, C.WINS[wn]).measure(Xf, cs)
        D[f"eps3c_{wn}"], D[f"cc3c_{wn}"] = e_, c_
    return D


def svd_cells(D, wins=None, control=False):
    """episode aggregation of SVD day values: days CC > 0.8 (data) / any cc (control), non-railed, CC-weighted mean;
    sigma = pooled within-episode weighted scatter / sqrt(n_eff)"""
    wins = wins or list(C.WINS)
    if not len(D):
        return pd.DataFrame()
    D = D[D.episode.notna()]
    out = None
    for wn in wins:
        e, c = f"eps3c_{wn}", f"cc3c_{wn}"
        K = D[np.isfinite(D[e]) & (D[e].abs() < C.RAIL) & ((D[c] > CC_DAY) if not control else np.isfinite(D[c]))].copy()
        K["w"] = K[c].clip(lower=0.01)
        if not len(K):
            continue
        g = K.groupby("episode")
        mean = g.apply(lambda q: np.average(q[e], weights=q.w), include_groups=False)
        K["res"] = K[e] - K.episode.map(mean)
        neff = g.apply(lambda q: q.w.sum() ** 2 / (q.w ** 2).sum(), include_groups=False)
        nd = g.size()
        dof = (nd - 1).clip(lower=0).sum()
        if dof >= 3:
            sday = float(np.sqrt(np.sum(K.w * K.res ** 2) / K.w.sum() * len(K) / max(dof, 1)))
        elif len(K) >= 3:      # too few repeat days within episodes: scatter of all days about the series mean
            mu = np.average(K[e], weights=K.w)
            sday = float(np.sqrt(np.sum(K.w * (K[e] - mu) ** 2) / K.w.sum() * len(K) / (len(K) - 1)))
        else:
            sday = np.nan
        r = pd.DataFrame({f"eps3c_{wn}": mean, f"cc3c_{wn}": g[c].mean(), f"sig3c_{wn}": sday / np.sqrt(neff),
                          f"ndays_{wn}": nd})
        out = r if out is None else out.join(r, how="outer")
    if out is None:
        return pd.DataFrame()
    out = out.reset_index()
    for wn in wins:
        for p in ("eps", "cc", "sig", "ndays"):
            if f"{p}3c_{wn}" not in out and not (p == "ndays"):
                out[f"{p}3c_{wn}"] = np.nan
            if p == "ndays" and f"ndays_{wn}" not in out:
                out[f"ndays_{wn}"] = 0
    return out
