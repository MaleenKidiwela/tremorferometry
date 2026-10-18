"""s11_ball3_lib -- shared helpers for the S11 ballistic round-3 analysis (ballistic/round3/PLAN.md).
Group clock terms (generalised round-2 quartet terms), corrected pair tables, oriented reference frames, per-series
contrasts, random-effects pooling, conversions."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from s11_ball_analysis import demean
import s11_ball2_clock as CK

ROOT3 = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic/round3"
VPVS = 1.75
SPLIT = 2015
PKEY = ["fam", "pair", "epA", "epB", "bin"]
GROUPS = {"SVI": {"PGCsite": ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]},
          "WA_NW": {"G1": ["PB.B005", "PB.B006", "PB.B007", "UW.LRIV"], "G2": ["PB.B013", "PB.B943"]}}
REFSETS = {"SVI": {"BH": ["PB.B009", "PB.B010", "PB.B011"], "B011": ["PB.B011"], "B0910": ["PB.B009", "PB.B010"],
                   "PGC": ["CN.PGC"]},
           "WA_NW": {"G1BH": ["PB.B005", "PB.B006", "PB.B007"], "G2": ["PB.B013", "PB.B943"], "LRIV": ["UW.LRIV"],
                     "G1all": ["PB.B005", "PB.B006", "PB.B007", "UW.LRIV"]}}
HEAD_REFS = {"SVI": ("BH", "B011", "B0910"), "WA_NW": ("G1BH", "G2")}


def dirs(region):
    reg = region.upper(); o = f"{ROOT3}/{reg.lower()}"
    return reg, f"{o}/tables", f"{ROOT3}/tables", f"{ROOT3}/figures"


def load(region, name, tag=""):
    reg, TB, _, _ = dirs(region)
    p = f"{TB}/{name}{tag}.parquet"
    return pd.read_parquet(p) if os.path.exists(p) else pd.read_csv(f"{TB}/{name}{tag}.csv")


def node_group(region):
    return {s: g for g, mem in GROUPS[region].items() for s in mem}


def prep_pairs(pr, val="dS", sig="sig_dS", rail="railS", bins=None):
    d = pr[pr.kind == "data"].copy()
    if bins is not None:
        d = d[d.bin.isin(bins)]
    d.loc[d[rail], val] = np.nan
    d["episode"] = d["lab"]
    out = [demean(g, val, PKEY) for _, g in d.groupby("bin")]
    D = pd.concat(out, ignore_index=True)
    D["sig"] = np.maximum(D[sig], 0.5)
    D["nodeA"] = D.A + "|" + D.epA; D["nodeB"] = D.B + "|" + D.epB
    return D


def group_terms(D, region):
    """per (fam, bin, lab, group): LS node lags from within-group pairs, median gauge -> rows fam, bin, lab, node, k, var"""
    ng = node_group(region)
    D = D.assign(gA=D.A.map(ng), gB=D.B.map(ng))
    Q = D[D.gA.notna() & (D.gA == D.gB)]
    rows = []
    for (f, bk, lab, gname), g in Q.groupby(["fam", "bin", "lab", "gA"]):
        nodes = sorted(set(g.nodeA) | set(g.nodeB))
        if len(nodes) < 2:
            continue
        ix = {n: i for i, n in enumerate(nodes)}
        Am = np.zeros((len(g), len(nodes))); w = 1 / g.sig.to_numpy()
        Am[np.arange(len(g)), [ix[n] for n in g.nodeB]] = 1; Am[np.arange(len(g)), [ix[n] for n in g.nodeA]] = -1
        k = np.linalg.lstsq(Am * w[:, None], g.r.to_numpy() * w, rcond=None)[0]
        k = k - np.median(k)
        prec = np.abs(Am).T @ (w ** 2)
        for n in nodes:
            rows.append(dict(fam=f, bin=bk, lab=lab, node=n, group=gname, k=k[ix[n]], var=1 / max(prec[ix[n]], 1e-9),
                             date=g.date.min(), cc=float(np.nanmedian(np.r_[g.ccS_A, g.ccS_B]))))
    return pd.DataFrame(rows)


def correct(D, T, region):
    ng = node_group(region)
    lfo = T.set_index(["fam", "bin", "lab", "node"]).K_lfo.to_dict()
    kA = np.array([lfo.get((f, b, l, n), np.nan) for f, b, l, n in zip(D.fam, D.bin, D.lab, D.nodeA)])
    kB = np.array([lfo.get((f, b, l, n), np.nan) for f, b, l, n in zip(D.fam, D.bin, D.lab, D.nodeB)])
    inA = D.A.map(ng).notna().to_numpy(); inB = D.B.map(ng).notna().to_numpy()
    D = D.copy()
    D["kcorr"] = np.where(inA, np.nan_to_num(kA), 0.0) - np.where(inB, np.nan_to_num(kB), 0.0)
    D["has_corr"] = (inA | inB) & (~inA | np.isfinite(kA)) & (~inB | np.isfinite(kB))
    D["rc0"] = D.r + D.kcorr
    D["rc"] = D.rc0 - D.groupby(PKEY).rc0.transform("mean")
    return D


def oriented(D, cols=("rc", "r")):
    """long table: X (target), M (reference member), v = value of X - M"""
    keep = ["fam", "bin", "lab", "date", "sig", "n"] + [c for c in ("tmed",) if c in D]
    a = D[keep].copy(); a["X"], a["epX"], a["M"], a["epM"] = D.B, D.epB, D.A, D.epA
    b = D[keep].copy(); b["X"], b["epX"], b["M"], b["epM"] = D.A, D.epA, D.B, D.epB
    for c in cols:
        a[c] = D[c].to_numpy(); b[c] = -D[c].to_numpy()
    if "ccS_A" in D:
        a["ccX"], a["ccM"] = D.ccS_B.to_numpy(), D.ccS_A.to_numpy(); b["ccX"], b["ccM"] = D.ccS_A.to_numpy(), D.ccS_B.to_numpy()
    if "ccP3_A" in D:
        a["ccPX"], b["ccPX"] = D.ccP3_B.to_numpy(), D.ccP3_A.to_numpy()
    return pd.concat([a, b], ignore_index=True)


def ref_frame(L, refset, col):
    F = L[L.M.isin(refset) & ~L.X.isin(refset)]
    g = F.groupby(["fam", "X", "epX", "bin", "lab"])
    agg = dict(v=(col, "mean"), sig=("sig", lambda s: np.sqrt((s ** 2).mean() / len(s))), date=("date", "min"),
               nM=("M", "size"))
    if "ccX" in F:
        agg["ccX"] = ("ccX", "median")
    if "ccPX" in F:
        agg["ccPX"] = ("ccPX", "median")
    return g.agg(**agg).reset_index()


def sp_table(region, tag=""):
    geo = load(region, "series_geom", tag)
    SP = geo.groupby(["fam", "netsta"]).SmP_pick.median()
    return SP


def TS(sp):
    return VPVS / (VPVS - 1) * sp


def TP(sp):
    return sp / (VPVS - 1)


def series_contrast(d, late_mask_fn=None, minside=2):
    """per series (sid) decade contrast: mean(late) - mean(early) with var from sig; d has sid, r, sig, t"""
    L = d.t.to_numpy() >= SPLIT if late_mask_fn is None else late_mask_fn(d)
    g = pd.DataFrame(dict(sid=d.sid.to_numpy(), r=d.r.to_numpy(), s2=d.sig.to_numpy() ** 2, L=L))
    rows = []
    for sid, x in g.groupby("sid"):
        l, e = x[x.L], x[~x.L]
        if len(l) < minside or len(e) < minside:
            continue
        rows.append(dict(sid=sid, C=l.r.mean() - e.r.mean(), v=l.s2.sum() / len(l) ** 2 + e.s2.sum() / len(e) ** 2,
                         nl=len(l), ne=len(e)))
    return pd.DataFrame(rows)


def pool_re(c, s, groups=None):
    """inverse-variance mean of values c with sigma s; formal sigma inflated by sqrt(chi2/dof) (>=1);
    random-effects sigma from the scatter of per-group means (groups: labels, e.g. station) -> dict"""
    c = np.asarray(c, float); s = np.asarray(s, float)
    ok = np.isfinite(c) & np.isfinite(s) & (s > 0)
    c, s = c[ok], s[ok]
    if len(c) == 0:
        return dict(mean=np.nan, sig=np.nan, n=0)
    w = 1 / s ** 2; m = (w * c).sum() / w.sum(); sf = 1 / np.sqrt(w.sum())
    chi2 = float((w * (c - m) ** 2).sum()); dof = max(len(c) - 1, 1)
    out = dict(mean=m, sig_formal=sf, chi2dof=chi2 / dof, sig_infl=sf * max(1, np.sqrt(chi2 / dof)), n=int(len(c)))
    sig_re = np.nan; ng = 0
    if groups is not None:
        G = pd.DataFrame(dict(c=c, w=w, g=np.asarray(groups)[ok]))
        gm = G.groupby("g").apply(lambda x: pd.Series(dict(m=(x.w * x.c).sum() / x.w.sum(), W=x.w.sum())))
        ng = len(gm)
        if ng >= 3:
            Wg = gm.W.to_numpy(); mg = gm.m.to_numpy(); mm = (Wg * mg).sum() / Wg.sum()
            sig_re = float(np.sqrt((Wg * (mg - mm) ** 2).sum() / ((ng - 1) * Wg.sum())))
    out.update(sig_re=sig_re, n_groups=ng, sig=float(np.nanmax([out["sig_infl"], sig_re])))
    return out
