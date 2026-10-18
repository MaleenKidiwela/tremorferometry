"""s11_r3_common -- shared helpers for S11 round 3 (PLAN.md §6). Region via S11_REGION (SVI default, WA_NW).
Reads the existing S11 tables (never writes them); every output goes to s11_depth/round3/<svi|wa_nw>/.
eps = -dv/v; cell d = dv/v in %."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I

BASE = f"{C.ROOT}/svi_science/results/s11_depth"
SUB = "svi" if C.REGION == "SVI" else C.REGION.lower()
R3 = f"{BASE}/round3/{SUB}"
TB3, LOG3, FIG3 = f"{R3}/tables", f"{R3}/logs", f"{BASE}/round3/figures"
for _d in (TB3, LOG3, FIG3):
    os.makedirs(_d, exist_ok=True)
TAG = "_all" if C.REGION == "SVI" else ""          # SVI headline = 115 families (main + ext), PLAN §6 common
NMIN, CCMIN = 20, 0.35
V4 = ["ell20", "ell40", "ell80", "two5_40"]
V6 = ["ell10", "ell20", "ell40", "ell80", "two5_40", "two3_20"]
WINS3 = (C.HEAD, "2035", "3560")


def log(name):
    return C.Log(f"{LOG3}/{name}.log")


# ------------------------------------------------------------------------------ kernels: all variants
def kernel_table():
    """long table series x variant x win (3c only): a_s, a_m, a_d, a_b (existing + round-3 variants)"""
    parts = [pd.read_parquet(f"{C.TB}/kernel_coef.parquet")]
    if C.REGION == "SVI" and os.path.exists(f"{C.TB}/kernel_coef_ext.parquet"):
        parts.append(pd.read_parquet(f"{C.TB}/kernel_coef_ext.parquet"))
    for f in ("kernel_coef_r3.parquet", "kernel_coef_r3_ext.parquet"):
        p = f"{TB3}/{f}"
        if os.path.exists(p):
            k = pd.read_parquet(p)
            parts.append(k[k.variant.isin(["ell10", "two3_20"])])     # ell40 re-run is only a reproduction check
    K = pd.concat(parts, ignore_index=True)
    K = K[K.comp == "3c"][["series", "fam", "netsta", "variant", "win", "a_s", "a_m", "a_d", "a_b"]]
    return K.drop_duplicates(["series", "variant", "win"])


_KT = {}


def coefs(variant, win):
    """series -> (a_s, a_d) for one variant and window; series missing -> mean over the same (family, station)"""
    if "K" not in _KT:
        _KT["K"] = kernel_table()
    K = _KT["K"]
    k = K[(K.variant == variant) & (K.win == win)]
    return k.set_index("series")[["fam", "netsta", "a_s", "a_d"]]


def setup_coefs(setup, win, cells):
    """setup: 'ell40' | ('avg', [variants]) | ('w', {v: weight}) | ('site', {site: variant}) -> arrays a_s, a_d for cells"""
    key = (cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch).to_numpy()
    fam, ns = cells.fam.to_numpy(), cells.netsta.to_numpy()

    def one(v):
        k = coefs(v, win)
        a_s = pd.Series(key).map(k.a_s).to_numpy(float); a_d = pd.Series(key).map(k.a_d).to_numpy(float)
        miss = ~np.isfinite(a_d)
        if miss.any():
            alt = k.groupby(["fam", "netsta"])[["a_s", "a_d"]].mean()
            for i in np.flatnonzero(miss):
                if (fam[i], ns[i]) in alt.index:
                    a_s[i], a_d[i] = alt.loc[(fam[i], ns[i])]
        return a_s, a_d
    if isinstance(setup, str):
        return one(setup)
    kind, arg = setup
    if kind == "avg":
        rs = [one(v) for v in arg]
        return np.nanmean([r[0] for r in rs], 0), np.nanmean([r[1] for r in rs], 0)
    if kind == "w":                                # nan-aware weighted mean over variants
        num_s = np.zeros(len(cells)); num_d = np.zeros(len(cells)); den = np.zeros(len(cells))
        for v, w in arg.items():
            s_, d_ = one(v); ok = np.isfinite(s_) & np.isfinite(d_)
            num_s += np.where(ok, w * s_, 0); num_d += np.where(ok, w * d_, 0); den += np.where(ok, w, 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den > 0, num_s / den, np.nan), np.where(den > 0, num_d / den, np.nan)
    if kind == "site":
        a_s = np.full(len(cells), np.nan); a_d = np.full(len(cells), np.nan); site = cells.site.to_numpy()
        for v in set(arg.values()):
            s_, d_ = one(v); m = np.array([arg.get(x) == v for x in site])
            a_s[m], a_d[m] = s_[m], d_[m]
        return a_s, a_d
    raise ValueError(setup)


def with_setup(cells, setup, win=C.HEAD):
    c = cells.copy()
    c["a_s"], c["a_d"] = setup_coefs(setup, win, c)
    return c


def setup_name(setup):
    if isinstance(setup, str):
        return setup
    kind, arg = setup
    if kind == "avg":
        return "avg(" + ",".join(arg) + ")"
    return f"{kind}(" + ",".join(f"{k}:{v if isinstance(v, str) else round(v, 3)}" for k, v in arg.items()) + ")"


# ------------------------------------------------------------------------------ cells
def head_cells(win=C.HEAD, nmin=NMIN, control=False):
    """headline cell rule (PLAN §6 common), ell40 coefficients attached (setups replace them)"""
    if control:
        c = I.load_cells("raw", win, TAG, control=True, nmin=nmin, ccmin=CCMIN)
        return I.prepare("raw", win, TAG, "ell40", cells=c)
    return I.prepare("raw", win, TAG, "ell40", nmin=nmin, ccmin=CCMIN)


def lapse_cells():
    """joint-lapse cell set: 2-3.5 s and 3.5-6 s cells, own offsets per window (epoch suffix), column win"""
    out = []
    for w in ("2035", "3560"):
        c = head_cells(win=w); c["win"] = w; c["epoch_orig"] = c.epoch; c["epoch"] = c.epoch + f"_w{w}"; out.append(c)
    return pd.concat(out, ignore_index=True)


def lapse_with_setup(CL, setup):
    parts = []
    for w in ("2035", "3560"):
        q = CL[CL.win == w].copy(); q2 = q.assign(epoch=q.epoch_orig)
        q["a_s"], q["a_d"] = setup_coefs(setup, w, q2); parts.append(q)
    return pd.concat(parts, ignore_index=True)


def majors():
    import s11_q3 as Q3
    return list(Q3.MAJ15)


def dump(obj, path):
    json.dump(obj, open(path, "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))


def pt_now():
    return C.pt_now()
