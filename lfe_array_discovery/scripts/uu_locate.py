#!/usr/bin/env python
"""The rebuilt locator.  Prefix `u_`.  One composable engine, run as a LADDER of variants so that
every change is measured on its own against Bostock's published family coordinates.

A variant is a dict:
  obs     't_env' (the current locator's envelope maximum) | 't_onS' (AIC onset) | 'xc' (station-
          pair envelope cross-correlation differential times) | 'bp' (envelope back-projection,
          no picks at all)
  robust  False -> the current rms of the demeaned residual; True -> the weighted mean |residual|
          about its median (L1), which removes the leverage of a single cycle-skipped pick
  weights None | 'q'  (quality weights from the stack peak height in MAD and the onset sharpness)
  sp      None | 'onP' | 'xc'   -- add the S-P observable (Chestler's second class)
  grid    'base' (the current narrow box) | 'wide'
  refine  nested local search after the coarse minimum (kills the grid-node quantisation)
  terms   'pub' (the a5/intraslab station terms in use now) | 'none' | 'loo' (terms estimated
          from the calibration families LEAVE-ONE-FAMILY-OUT, so no family sees its own term)
  depth   'free' | 'slab2' (depth fixed to the Slab2 interface at each node)
  vpvs    Vp/Vs of the Crosson model

Every family keeps its row in every variant's output; nothing is filtered (memory
`rank-dont-discard`).  Quality appears as columns (n_sta, cost, edge, ctrl_pass, boot_km).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uu_common as U                                                        # noqa: E402

MIN_STA = 5
MIN_PEAK = 4.0
BASE_GRID = dict(lat=(47.90, 49.31, 0.03), lon=(-124.90, -122.79, 0.04), dep=(18.0, 56.1, 2.0))
WIDE_GRID = dict(lat=(47.30, 49.80, 0.03), lon=(-125.80, -121.90, 0.04), dep=(6.0, 76.0, 2.0))
WA_BASE_GRID = dict(lat=(46.50, 48.80, 0.03), lon=(-125.20, -121.80, 0.04), dep=(15.0, 60.0, 2.0))
WA_WIDE_GRID = dict(lat=(46.30, 49.00, 0.03), lon=(-125.60, -121.60, 0.04), dep=(6.0, 76.0, 2.0))

LADDER = [
    ("L0_base_env_L2",  dict(obs="t_env", hcol="h_env", robust=False, wpow=0, grid="base")),
    ("L1_wide",         dict(obs="t_env", hcol="h_env", robust=False, wpow=0, grid="wide")),
    ("L2_robust",       dict(obs="t_env", hcol="h_env", robust=True,  wpow=0, grid="wide")),
    ("L3_snrw",         dict(obs="t_env", hcol="h_env", robust=True,  wpow=1, grid="wide")),
    ("L4_snrw2",        dict(obs="t_env", hcol="h_env", robust=True,  wpow=2, grid="wide")),
    ("L5_horiz",        dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide")),
    ("L6_terms",        dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo")),
    ("L7_refine",       dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True)),
    ("L8_bp",           dict(obs="bp",    hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo")),
    ("L8b_irls",        dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, irls=2)),
    ("L8c_irls_adapt",  dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, irls=2, adaptive=True)),
    ("L8d_irls_zprior", dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, irls=2, adaptive=True, zprior=8.0)),
    ("L9_xcdiff",       dict(obs="xc",    hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo")),
    ("L10_xc_refine",   dict(obs="xc",    hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True)),
    ("L11_xc_zprior",   dict(obs="xc",    hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, zprior=8.0)),
    ("L12_comb",        dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, irls=2, xc=True, a_pick=1.0, a_xc=1.0)),
    ("L13_comb_zprior", dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, irls=2, xc=True, a_pick=1.0, a_xc=1.0,
                             zprior=8.0)),
]
ABLATE = [
    ("A_onset",         dict(obs="t_onS", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True)),
    ("A_sp_onP",        dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, sp="onP")),
    ("A_sp_xc",         dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, sp="xc")),
    ("A_no_terms",      dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             refine=True)),
    ("A_L2_at_L7",      dict(obs="t_envH", hcol="h_envH", robust=False, wpow=2, grid="wide",
                             terms="loo", refine=True)),
    ("A_minpeak3",      dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, minpeak=3.0)),
    ("A_minpeak5",      dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, minpeak=5.0)),
    ("A_slab2_depth",   dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, depth="slab2")),
    ("A_wpow3",         dict(obs="t_envH", hcol="h_envH", robust=True, wpow=3, grid="wide",
                             terms="loo", refine=True)),
    ("A_wpow4",         dict(obs="t_envH", hcol="h_envH", robust=True, wpow=4, grid="wide",
                             terms="loo", refine=True)),
    ("A_adaptive",      dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, adaptive=True)),
    ("A_zprior8",       dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, zprior=8.0)),
    ("A_comb_xc2",      dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, irls=2, xc=True, a_pick=0.5, a_xc=1.5,
                             zprior=8.0)),
    ("A_xc_cc50",       dict(obs="xc",    hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, zprior=8.0, cc_min=0.5)),
    ("A_xc_noterms",    dict(obs="xc",    hcol="h_envH", robust=True, wpow=2, grid="wide",
                             refine=True, zprior=8.0)),
    ("A_vpvs_172",      dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, vpvs=1.72)),
    ("A_vpvs_185",      dict(obs="t_envH", hcol="h_envH", robust=True, wpow=2, grid="wide",
                             terms="loo", refine=True, vpvs=1.85)),
]


# ------------------------------------------------------------------ observables ---------------
def quality_weight(row, p=2):
    """Weight from MEASURED pick quality only (no truth anywhere in the formula).

    Calibration (results/u_locate/north/pick_residuals.csv): the long-stack envelope peak height
    in MAD predicts the pick residual almost deterministically -- of picks with h < 4 MAD only
    ~20 % land within 1 s of the time predicted by Bostock's published coordinates, at h >= 6 it
    is 81 % and at h >= 8 it is 90 % with a 0.15 s median.  The weight is therefore (h - 3)^p,
    p = 2 being the inverse-variance-like form; p = 0 reproduces the unweighted current locator."""
    h = np.clip(np.asarray(row["h_envH"], float), 0.0, 40.0)
    if p == 0:
        return np.ones(h.shape)
    return np.clip(h - 3.0, 0.05, 37.0) ** p


def family_obs(P, spec):
    """{grp: (stas, tS, tP, w)} from the pick table.  Nothing is dropped from the OUTPUT -- a
    family with too few admitted stations still gets a row, flagged `located=False`."""
    key = spec.get("obs", "t_env")
    tcol = "t_env" if key in ("bp", "xc") else key
    hcol = spec.get("hcol", "h_env")
    mp = float(spec.get("minpeak", MIN_PEAK))
    out = {}
    for g, gg in P.groupby("grp"):
        gg = gg[np.isfinite(gg[tcol])]
        if spec.get("adaptive"):
            # prefer the clean picks: use h >= 6 MAD if that still leaves 5 stations, else 5,
            # else the 4-MAD admission the current locator uses.  A cascade, not a cut: a family
            # that cannot reach 5 stations at any level keeps its row with located=False.
            for lev in (6.0, 5.0, mp):
                if (gg[hcol] >= lev).sum() >= MIN_STA:
                    gg = gg[gg[hcol] >= lev]
                    break
            else:
                gg = gg[gg[hcol] >= mp]
        else:
            gg = gg[gg[hcol] >= mp]
        tP = np.full(len(gg), np.nan)
        if len(gg):
            if spec.get("sp") == "onP":
                tP = gg.t_onP.to_numpy(float)
            elif spec.get("sp") == "xc":
                sp = gg.sp_xc.to_numpy(float)
                cc = gg.sp_cc.to_numpy(float)
                tP = np.where(np.isfinite(sp) & (cc > 0.5), gg[tcol].to_numpy(float) - sp, np.nan)
        w = quality_weight(gg, spec.get("wpow", 0)) if len(gg) else np.ones(0)
        out[int(g)] = (list(gg.sta), gg[tcol].to_numpy(float), tP, np.asarray(w, float))
    return out


# ------------------------------------------------------------------ station terms --------------
def loo_terms(P, truth, ST, vpvs=1.78, tcol="t_envH", hcol="h_envH", minpeak=4.0):
    """Per-station S term, estimated LEAVE-ONE-FAMILY-OUT from the calibration families only.

    term_s = median over calibration families f (excluding f itself) of the station's demeaned
    residual against the PUBLISHED coordinates of f.  The map returned is {(bkey, sta): term};
    a family not in the calibration set gets the all-family term."""
    xs, zs, TP, TS = U.tt_table(vpvs=vpvs)
    rows = []
    for r in truth.itertuples():
        gg = P[P.grp == r.grp]
        gg = gg[np.isfinite(gg[tcol]) & (gg[hcol] >= 6.0)]
        if len(gg) < 4:
            continue
        d = U.epi_km(r.t_lat, r.t_lon, ST.lat[gg.sta].to_numpy(float),
                     ST.lon[gg.sta].to_numpy(float))
        pred = U.interp2(xs, zs, TS, d, np.full(len(gg), r.t_dep))
        res = gg[tcol].to_numpy(float) - pred
        res = res - np.median(res)
        for s, rr in zip(gg.sta, res):
            rows.append(dict(bkey=int(r.bkey), sta=s, res=float(rr)))
    R = pd.DataFrame(rows)
    if R.empty:
        return {}, {}
    glob = R.groupby("sta").res.apply(lambda x: float(np.median(x[np.abs(x) < 3.0])
                                                      if (np.abs(x) < 3.0).any()
                                                      else np.median(x))).to_dict()
    loo = {}
    for b in R.bkey.unique():
        sub = R[R.bkey != b]
        loo[int(b)] = sub.groupby("sta").res.apply(
            lambda x: float(np.median(x[np.abs(x) < 3.0]) if (np.abs(x) < 3.0).any()
                            else np.median(x))).to_dict()
    return glob, loo


# ------------------------------------------------------------------ the search -----------------
def make_grid(stas, ST, gspec, vpvs, termS, termP, depth="free", zprior=0.0):
    G = U.Grid(stas, ST, gspec["lat"], gspec["lon"], gspec["dep"], vpvs=vpvs,
               termS=termS, termP=termP)
    G.mask = None
    G.sz = None
    if depth == "slab2":
        sz = U.slab2(G.la, G.lo)
        G.mask = np.abs(G.dp - sz) <= G.deps[1] - G.deps[0]
    elif zprior:
        G.sz = U.slab2(G.la, G.lo).astype(np.float32)
        G.zprior = float(zprior)
    return G


def cost_at(G, idx, tS, tP, w, spec):
    sel = np.isfinite(tS)
    TTs = G.TS[:, idx[sel]]
    c = U.misfit(TTs, tS[sel], w[sel], spec.get("robust", False))
    m = np.isfinite(tP) & sel
    if m.any():
        sp = (G.TS[:, idx[m]] - G.TP[:, idx[m]]) - (tS[m] - tP[m])[None, :]
        add = np.abs(sp).mean(1) if spec.get("robust") else np.sqrt((sp ** 2).mean(1))
        c = c + 0.6 * add
    if getattr(G, "sz", None) is not None:
        # soft Slab2 depth prior, expressed in seconds at 0.3 s per `zprior` km of departure.
        # NOTE (measured): Bostock's own published family depths sit a MEDIAN 9.7 km BELOW this
        # project's Slab2 field, so the prior is deliberately weak and is reported as a variant.
        c = c + 0.3 * np.abs(G.dp - G.sz) / G.zprior
    if getattr(G, "mask", None) is not None:
        c = np.where(G.mask, c, np.inf)
    return c


def refine_generic(ST, sel_stas, la0, lo0, dp0, vpvs, termS, termP, evalcost,
                   levels=((0.03, 0.04, 2.0), (0.006, 0.008, 0.4), (0.0012, 0.0016, 0.08))):
    """Nested local search around the coarse minimum: three levels, each a 7x7x7 box one coarse
    cell wide, so the final node spacing is ~0.13 km horizontally and 0.08 km in depth.  This is
    what removes the grid-node quantisation that put 75 % of the baseline solutions on a face."""
    xs, zs, TPt, TSt = U.tt_table(vpvs=vpvs)
    slat = ST.lat[sel_stas].to_numpy(float)
    slon = ST.lon[sel_stas].to_numpy(float)
    tsc = np.array([termS.get(s, 0.0) for s in sel_stas], float)
    tpc = np.array([termP.get(s, 0.0) for s in sel_stas], float)
    la, lo, dp, best = la0, lo0, dp0, np.inf
    for (dla, dlo, ddp) in levels:
        A = np.arange(-3, 4) * (dla / 3.0)
        Bv = np.arange(-3, 4) * (dlo / 3.0)
        C = np.arange(-3, 4) * (ddp / 3.0)
        LA, LO, DP = np.meshgrid(la + A, lo + Bv, np.maximum(dp + C, 4.0), indexing="ij")
        LA, LO, DP = LA.ravel(), LO.ravel(), DP.ravel()
        TS = np.empty((LA.size, len(sel_stas)), np.float32)
        TPp = np.empty_like(TS)
        for k in range(len(sel_stas)):
            d = U.epi_km(LA, LO, slat[k], slon[k])
            TS[:, k] = U.interp2(xs, zs, TSt, d, DP) + tsc[k]
            TPp[:, k] = U.interp2(xs, zs, TPt, d, DP) + tpc[k]
        c = evalcost(TS, TPp, LA, LO, DP)
        i = int(np.argmin(c))
        la, lo, dp, best = float(LA[i]), float(LO[i]), float(DP[i]), float(c[i])
    return la, lo, dp, best


def make_evalcost(spec, tS, tP, w, pairs, dt, cc, a_pick=1.0, a_xc=1.0):
    """Closure giving the cost of a set of candidate nodes from their (nnode, nsta) S and P
    tables.  Handles the pick term, the S-P term, the station-pair XC differential term and the
    soft Slab2 depth prior, so the coarse grid and the nested refinement use IDENTICAL code."""
    rb = bool(spec.get("robust", False))
    zp = float(spec.get("zprior", 0.0) or 0.0)
    use_pick = spec.get("obs") != "xc"
    selS = np.isfinite(tS) if use_pick else None
    mP = (np.isfinite(tP) & selS) if use_pick else None

    def f(TS, TPp, LA=None, LO=None, DP=None):
        c = np.zeros(TS.shape[0], np.float32)
        if use_pick and selS.any():
            c = c + U.misfit(TS[:, selS], tS[selS], w[selS], rb)
            if mP is not None and mP.any():
                sp = (TS[:, mP] - TPp[:, mP]) - (tS[mP] - tP[mP])[None, :]
                c = c + 0.6 * (np.abs(sp).mean(1) if rb else np.sqrt((sp ** 2).mean(1)))
            c = c * a_pick
        if pairs is not None and len(pairs):
            D = TS[:, pairs[:, 0]] - TS[:, pairs[:, 1]]
            r = D - dt[None, :]
            wp = cc / cc.sum()
            c = c + a_xc * ((np.abs(r) * wp[None, :]).sum(1) if rb
                            else np.sqrt(((r ** 2) * wp[None, :]).sum(1)))
        if zp and DP is not None:
            c = c + 0.3 * np.abs(DP - U.slab2(LA, LO)) / zp
        return c
    return f


# ------------------------------------------------------------------ XC differential ------------
def xc_differential(stacks, grp, stas, fs=U.FS, lmax=6.0):
    """Station-pair differential S times from ENVELOPE cross-correlation of the long stacks
    (Chestler & Creager cross-correlate stacked templates ACROSS stations to form differential
    times; between stations tens of km apart the waveforms are not similar but the S envelopes
    are, so the correlation is done on the envelope).  FFT-based, all pairs at once.
    Returns (pairs (m,2) indices into `stas`, dt (m,), cc (m,))."""
    import uu_picks as PK
    E, keep = [], []
    for i, s in enumerate(stas):
        k = f"{grp}_{s}"
        if k not in stacks:
            continue
        e = PK.envelope(np.asarray(stacks[k], float), fs)
        e = e - e.mean()
        n = np.sqrt((e ** 2).sum())
        E.append(e / max(n, 1e-30))
        keep.append(i)
    if len(E) < 2:
        return np.zeros((0, 2), int), np.zeros(0), np.zeros(0)
    E = np.asarray(E)
    N = E.shape[1]
    nl = int(round(lmax * fs))
    nf = 1
    while nf < 2 * N:
        nf *= 2
    F = np.fft.rfft(E, n=nf, axis=1)
    ps, dts, ccs = [], [], []
    for a_ in range(len(keep)):
        C = np.fft.irfft(F[a_][None, :] * np.conj(F[a_ + 1:]), n=nf, axis=1)
        if not C.size:
            continue
        # lag L>0 means station a_ is LATER than the other by L samples
        win = np.concatenate([C[:, -nl:], C[:, :nl + 1]], axis=1)
        k = np.argmax(win, axis=1)
        for m, b_ in enumerate(range(a_ + 1, len(keep))):
            ps.append((keep[a_], keep[b_]))
            dts.append((int(k[m]) - nl) / fs)
            ccs.append(float(win[m, k[m]]))
    return np.array(ps, int), np.array(dts, float), np.array(ccs, float)


# ------------------------------------------------------------------ back-projection ------------
def bp_locate(G, idx, ENV, w, t0_lo, t0_hi, fs=U.FS, batch=16384):
    """Envelope back-projection: max over (node, origin time) of the weighted envelope stack.
    ENV (k, N) already normalised per station.  Returns (node index, t0 seconds, power)."""
    try:
        import torch
        dev = "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        return None
    N = ENV.shape[1]
    TT = G.TS[:, idx]
    sh = np.rint((TT - TT.mean(1, keepdims=True)) * fs).astype(np.int32)
    k0 = max(int(round((t0_lo + U.LONG_PRE) * fs)), int(-sh.min()) + 1)
    k1 = min(int(round((t0_hi + U.LONG_PRE) * fs)), N - int(sh.max()) - 1)
    if k1 - k0 < 5:
        return None
    kg = torch.arange(k0, k1, device=dev)
    E = torch.from_numpy(np.ascontiguousarray(ENV, np.float32)).to(dev)
    wt = torch.from_numpy(np.asarray(w, np.float32)).to(dev)
    wt = wt / wt.sum()
    bestv, bestn, bestk = -1e30, 0, k0
    for b0 in range(0, sh.shape[0], batch):
        S = torch.from_numpy(sh[b0:b0 + batch]).to(dev).long()
        acc = torch.zeros((S.shape[0], kg.numel()), device=dev)
        for i in range(S.shape[1]):
            ix = (kg[None, :] + S[:, i:i + 1]).clamp_(0, N - 1)
            acc += wt[i] * E[i][ix]
        v, ki = acc.max(1)
        vv, ni = v.max(0)
        if float(vv) > bestv:
            bestv = float(vv)
            bestn = b0 + int(ni)
            bestk = k0 + int(ki[int(ni)])
        del acc, S
    return bestn, (bestk / fs - U.LONG_PRE), bestv


# ------------------------------------------------------------------ driver ---------------------
def run_variant(name, spec, P, stacks, ST, stas_all, gspecs, termmap, boot=0, nshuf=0,
                slab=True):
    gs = gspecs["wide" if spec.get("grid", "base") == "wide" else "base"]
    vpvs = spec.get("vpvs", 1.78)
    termS_glob, termP_glob, loo = termmap[spec.get("terms", "pub")]
    G = make_grid(stas_all, ST, gs, vpvs, termS_glob, termP_glob,
                  spec.get("depth", "free"), spec.get("zprior", 0.0))
    si = {s: i for i, s in enumerate(stas_all)}
    OBS = family_obs(P, spec)
    HQ = {int(g): gg[spec.get("hcol", "h_env")].to_numpy(float) for g, gg in P.groupby("grp")}
    a_pick = float(spec.get("a_pick", 1.0))
    a_xc = float(spec.get("a_xc", 1.0))
    use_xc = spec.get("obs") == "xc" or bool(spec.get("xc"))
    rng = np.random.default_rng(7)
    rows = []
    for g, (ss, tS, tP, w) in OBS.items():
        ss2 = [s for s in ss if s in si]
        keep = np.array([s in si for s in ss], bool)
        hq = HQ.get(g, np.array([0.0]))
        if len(ss2) < MIN_STA:
            rows.append(dict(grp=g, n_sta=len(ss2), located=False))
            continue
        tS, tP, w = tS[keep], tP[keep], w[keep]
        idx = np.array([si[s] for s in ss2])
        termS = dict(termS_glob)
        termP = dict(termP_glob)
        if spec.get("terms") == "loo" and g in loo:
            # LEAVE-ONE-FAMILY-OUT: the grid carries the all-family term, so a family's own
            # contribution is removed by shifting its OBSERVATIONS (exactly equivalent in a
            # demeaned differential misfit, and it avoids rebuilding the grid per family).
            termS = {s: float(loo[g].get(s, termS_glob.get(s, 0.0))) for s in ss2}
            termP = {s: termS[s] / vpvs for s in ss2}
            adj = np.array([termS[s] - float(termS_glob.get(s, 0.0)) for s in ss2])
            tS = tS - adj
            tP = tP - adj / vpvs
        else:
            termS = {s: float(termS_glob.get(s, 0.0)) for s in ss2}
            termP = {s: float(termP_glob.get(s, 0.0)) for s in ss2}

        # ---- back-projection is its own estimator (no picks at all)
        if spec.get("obs") == "bp":
            import uu_picks as PK
            ENV = []
            for s in ss2:
                e = PK.envelope(np.asarray(stacks[f"{g}_{s}"], float))
                e = e - np.median(e)
                sd = 1.4826 * np.median(np.abs(e - np.median(e)))
                ENV.append(e / max(sd, 1e-30))
            r = bp_locate(G, idx, np.asarray(ENV, np.float32), w, -6.0, 26.0)
            if r is None:
                rows.append(dict(grp=g, n_sta=len(ss2), located=False))
                continue
            n, t0, pw = r
            la, lo, dp, cost = float(G.la[n]), float(G.lo[n]), float(G.dp[n]), -float(pw)
            rows.append(dict(grp=g, n_sta=len(ss2), located=True, lat=la, lon=lo, dep=dp,
                             cost=float(cost), edge=G.edge(la, lo, dp), n_downweighted=0,
                             h_med=float(np.median(hq)), h_max=float(np.max(hq)),
                             slab2=float(U.slab2(la, lo)[0]),
                             dz=dp - float(U.slab2(la, lo)[0])))
            continue

        # ---- station-pair differential times (Chestler's differential-time class)
        pairs = dt = cc = None
        n_pair = 0
        if use_xc:
            ps, d_, c_ = xc_differential(stacks, g, ss2)
            if len(ps):
                ok = c_ >= float(spec.get("cc_min", 0.30))
                ps, d_, c_ = ps[ok], d_[ok], c_[ok]
            if len(ps) >= int(spec.get("min_pairs", 6)):
                pairs, dt, cc = ps, d_, c_
                n_pair = len(ps)
            elif spec.get("obs") == "xc":
                rows.append(dict(grp=g, n_sta=len(ss2), located=False, n_pair=0))
                continue

        def ev(TS, TPp, LA=None, LO=None, DP=None, ww=w, pp=pairs, dd=dt, ccv=cc):
            return make_evalcost(spec, tS, tP, ww, pp, dd, ccv, a_pick, a_xc)(
                TS, TPp, LA, LO, DP)

        cost = ev(G.TS[:, idx], G.TP[:, idx], G.la, G.lo, G.dp)
        if getattr(G, "mask", None) is not None:
            cost = np.where(G.mask, cost, np.inf)
        la, lo, dp, cost, n = U.best_node(G, cost)
        if spec.get("refine"):
            la, lo, dp, cost = refine_generic(ST, ss2, la, lo, dp, vpvs, termS, termP, ev)
        nrej = 0
        if spec.get("irls"):
            xs_, zs_, TPt_, TSt_ = U.tt_table(vpvs=vpvs)
            w2 = w.copy()
            for _ in range(int(spec["irls"])):
                dd_ = U.epi_km(la, lo, ST.lat[ss2].to_numpy(float), ST.lon[ss2].to_numpy(float))
                pr = U.interp2(xs_, zs_, TSt_, dd_, np.full(len(ss2), dp)) + np.array(
                    [termS.get(s_, 0.0) for s_ in ss2])
                r = tS - pr
                r = r - np.median(r[np.isfinite(r)]) if np.isfinite(r).any() else r
                m = max(float(1.4826 * np.median(np.abs(r[np.isfinite(r)]))), 0.20) \
                    if np.isfinite(r).any() else 1.0
                w2 = w * np.minimum(1.0, (2.5 * m / np.maximum(np.abs(r), 1e-6)) ** 2)
                nrej = int(np.nansum(np.abs(r) > 2.5 * m))

                def ev2(TS, TPp, LA=None, LO=None, DP=None, ww=w2):
                    return make_evalcost(spec, tS, tP, ww, pairs, dt, cc, a_pick, a_xc)(
                        TS, TPp, LA, LO, DP)
                c2 = ev2(G.TS[:, idx], G.TP[:, idx], G.la, G.lo, G.dp)
                if getattr(G, "mask", None) is not None:
                    c2 = np.where(G.mask, c2, np.inf)
                la, lo, dp, cost, n = U.best_node(G, c2)
                if spec.get("refine"):
                    la, lo, dp, cost = refine_generic(ST, ss2, la, lo, dp, vpvs, termS, termP,
                                                      ev2)
            w = w2
        d = dict(grp=g, n_sta=len(ss2), n_pair=n_pair, located=True, lat=la, lon=lo, dep=dp,
                 cost=float(cost), edge=G.edge(la, lo, dp), n_downweighted=nrej,
                 h_med=float(np.median(hq)), h_max=float(np.max(hq)))
        if nshuf:
            sh = []
            for _ in range(nshuf):
                pi = rng.permutation(idx.size)
                pp = None if pairs is None else np.column_stack([pi[pairs[:, 0]],
                                                                 pi[pairs[:, 1]]])
                cs = make_evalcost(spec, tS, tP, w, pp, dt, cc, a_pick, a_xc)(
                    G.TS[:, idx[pi]], G.TP[:, idx[pi]], G.la, G.lo, G.dp)
                sh.append(float(np.nanmin(cs)))
            d["cost_shuf_p05"] = float(np.quantile(sh, 0.05))
            d["ctrl_pass"] = bool(cost < d["cost_shuf_p05"])
        if boot:
            bl, bo, bd = [], [], []
            for _ in range(boot):
                pi = rng.integers(0, idx.size, idx.size)
                if len(set(pi.tolist())) < 3:
                    continue
                if pairs is not None and len(pairs):
                    inv = {}
                    for a_, v in enumerate(pi):
                        inv.setdefault(int(v), []).append(a_)
                    sel = [k for k, (i1, j1) in enumerate(pairs)
                           if i1 in inv and j1 in inv]
                    if len(sel) < 4:
                        continue
                    pp = np.array([[inv[pairs[k][0]][0], inv[pairs[k][1]][0]] for k in sel])
                    dd2, cc2 = dt[sel], cc[sel]
                else:
                    pp, dd2, cc2 = None, None, None
                cs = make_evalcost(spec, tS[pi], tP[pi], w[pi], pp, dd2, cc2, a_pick, a_xc)(
                    G.TS[:, idx[pi]], G.TP[:, idx[pi]], G.la, G.lo, G.dp)
                i2 = int(np.nanargmin(cs))
                bl.append(float(G.la[i2]))
                bo.append(float(G.lo[i2]))
                bd.append(float(G.dp[i2]))
            if len(bl) >= 10:
                bl, bo, bd = np.array(bl), np.array(bo), np.array(bd)
                d["boot_km"] = float(np.quantile(U.epi_km(bl, bo, la, lo), 0.68))
                d["boot_dz"] = float(np.quantile(np.abs(bd - dp), 0.68))
                d["boot_n"] = len(bl)
        if slab:
            sz = float(U.slab2(la, lo)[0])
            d["slab2"] = sz
            d["dz"] = dp - sz
        rows.append(d)
    R = pd.DataFrame(rows)
    R["variant"] = name
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default="north")
    ap.add_argument("--picks", required=True)
    ap.add_argument("--stacks", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--variants", default="ladder")
    ap.add_argument("--boot", type=int, default=0)
    ap.add_argument("--nshuf", type=int, default=0)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    B = U.bed(a.bed)
    if a.bed != "north":
        B.m.activate()
    ST = B.station_table()
    P = pd.read_parquet(a.picks)
    P = P[P.sta.isin(ST.index)]
    stas_all = sorted(P.sta.unique())
    Z = np.load(a.stacks, allow_pickle=True)
    stacks = {k: Z[k] for k in Z.files if k != "meta"} if a.variants != "ladder-nostack" else {}
    gspecs = (dict(base=BASE_GRID, wide=WIDE_GRID) if a.bed == "north"
              else dict(base=WA_BASE_GRID, wide=WA_WIDE_GRID))
    cS, cP, flag = U.station_terms_published(B, stas_all)
    U.log(f"{len(stas_all)} stations; published terms non-zero for {len(stas_all)-len(flag)}")
    termmap = {"pub": (cS, cP, {}), "none": ({}, {}, {})}
    if a.bed == "north":
        T = U.north_truth()
        T = T[T.set == "dev"]
        glob, loo = loo_terms(P, T, ST)
        fam_of = {int(r.grp): int(r.bkey) for r in T.itertuples()}
        loo_g = {g: loo.get(fam_of.get(g, -1), glob) for g in P.grp.unique()}
        termmap["loo"] = (glob, {k: v / 1.78 for k, v in glob.items()}, loo_g)
        json.dump(dict(global_terms=glob, published=cS, n_families=int(T.bkey.nunique())),
                  open(f"{a.outdir}/station_terms.json", "w"), indent=1)
    else:
        src = json.load(open(os.environ.get("U_TERMS",
                                            f"{U.OUT}/north/station_terms.json")))
        glob = {s: float(v) for s, v in src["global_terms"].items() if s in set(stas_all)}
        termmap["loo"] = (cS, cP, {})       # WA has no published truth -> published terms only
        U.log(f"WA: loo terms unavailable (no published WA truth); using published terms "
              f"({len(glob)} north terms carried only as a diagnostic)")
    names = [n for n, _ in LADDER] + [n for n, _ in ABLATE]
    todo = LADDER + ABLATE
    if a.only:
        want = set(a.only.split(","))
        todo = [(n, s) for n, s in todo if n in want]
    for name, spec in todo:
        out = f"{a.outdir}/loc_{name}.csv"
        if os.path.exists(out) and os.path.getsize(out) > 0:
            U.log(f"skip {name} (exists)")
            continue
        U.log(f"== {name} {spec}")
        R = run_variant(name, spec, P, stacks, ST, stas_all, gspecs, termmap,
                        boot=a.boot, nshuf=a.nshuf)
        R.to_csv(out, index=False)
        L = R[R.located == True]                                            # noqa: E712
        U.log(f"   {len(R)} families; located {len(L)}; "
              f"edge {int(L.edge.sum()) if len(L) else 0}")
    U.log("done " + ",".join(names))


if __name__ == "__main__":
    main()
