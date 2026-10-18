"""s11_r3_inject -- waveform-level injection engine for S11 round 3 (PLAN §6). Day stacks of each series are stretched
about tS by a day-specific eps (s11_paths.inject) before any stacking, then re-measured by an estimator:
  raw   = the S11 raw path (episode stack vs leave-episode-out reference, joint 3C 2-5 s)  [s11_paths.raw_cells]
  graph = all-pairs graph of episode-node stacks (nodes fixed per series), REDPy solve (Weaver weights, Huber)
eps = -dv/v. Responses are returned per (series, run) as DataFrame(episode, eps)."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_paths as P
from s11_measure import Ref3

RAIL = C.RAIL
WIN = C.WINS[C.HEAD]
_G = {}


def load_store_file(path):
    z = np.load(path)
    keys = z["keys"]; nd = z["n_days"]; n2 = z["n_c2"]
    o1 = np.r_[0, np.cumsum(nd)]; o2 = np.r_[0, np.cumsum(n2)]
    St, Ws, days, Sc, Cn, c2d = z["St"], z["Ws"], z["days"], z["Sc"], z["Cn"], z["c2_day"]
    Nd = z["Nd"] if "Nd" in z else None
    out = {}
    for i, k in enumerate(keys):
        out[str(k)] = dict(days=days[o1[i]:o1[i + 1]], St=St[o1[i]:o1[i + 1]], Ws=Ws[o1[i]:o1[i + 1]],
                           c2_day=c2d[o2[i]:o2[i + 1]], Sc=Sc[o2[i]:o2[i + 1]], Cn=Cn[o2[i]:o2[i + 1]],
                           G1=z["G1"][i].astype(np.float64), tS=float(z["tS"][i]), sp=float(z["sp"][i]),
                           Nd=None if Nd is None else Nd[o1[i]:o1[i + 1]])
    return out


# ------------------------------------------------------------------------------ graph
def graph_pairs(St, tS, win=WIN):
    n = len(St); pi, pj, dd, cc = [], [], [], []
    for i in range(n - 1):
        R = Ref3(St[i], tS, win)
        e, c = R.measure(St[i + 1:], (0, 1, 2))
        pi += [i] * (n - 1 - i); pj += list(range(i + 1, n)); dd += list(e); cc += list(c)
    return np.array(pi, int), np.array(pj, int), np.array(dd, float), np.array(cc, float)


def graph_solve(n, pi, pj, d, cc):
    """REDPy solve (s11_redpy_common.solve): v_j - v_i = d -> v (zero mean), misclosure sigma, chi2"""
    import s11_redpy_common as RP
    ok = np.isfinite(d) & np.isfinite(cc)
    return RP.solve(n, pi[ok], pj[ok], d[ok], np.nan_to_num(cc[ok]))


def graph_from_stacks(St, tS):
    pi, pj, d, cc = graph_pairs(St, tS)
    v, sv, chi2, scale, npair = graph_solve(len(St), pi, pj, d, cc)
    mpc = np.array([np.nanmedian(cc[(pi == i) | (pj == i)]) if np.any((pi == i) | (pj == i)) else np.nan for i in range(len(St))])
    return v, sv, chi2, mpc


def episode_stacks_from_days(st, unit, amp, eps_list):
    tS = st["tS"]; X = P.inject(st["St"], tS, amp); Ws = st["Ws"]
    ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
    cs = [c for c in range(3) if np.isfinite(st["G1"][c]).all()]
    S = np.full((len(eps_list), 3, X.shape[-1]), np.nan)
    for k, e in enumerate(eps_list):
        m = ep == e
        for c in cs:
            w = Ws[m, c]
            if w.sum() > 0:
                S[k, c] = (w @ X[m, c]) / w.sum()
    return S


def est_raw(st, unit, amp, key):
    D = P.raw_cells(st, unit, amp=amp, wins=[C.HEAD])
    return D[["episode", f"eps3c_{C.HEAD}"]].rename(columns={f"eps3c_{C.HEAD}": "eps"})


def est_graph(st, unit, amp, key):
    nodes = _G["nodes"].get(key, [])
    if len(nodes) < 3:
        return pd.DataFrame(columns=["episode", "eps"])
    S = episode_stacks_from_days(st, unit, amp, nodes)
    v = graph_from_stacks(S, st["tS"])[0]
    return pd.DataFrame(dict(episode=nodes, eps=v))


# ------------------------------------------------------------------------------ engine
def _worker(key):
    try:
        st = _G["store"][key]
        runs = _G["runs_fn"](key, st, _G["unit"])
        out = {}
        for name, amp in runs.items():
            out[name] = _G["est"](st, _G["unit"], amp, key)
        return key, out
    except Exception as e:
        import traceback
        print("ERROR", key, e, traceback.format_exc()[-400:], flush=True)
        return key, {}


def run(keys, runs_fn, est, store_paths, workers, lg, nodes=None, unit=None):
    """keys: set of series keys; store_paths: list of folders holding <NET.STA>.npz; runs_fn(key, st, unit) -> {run: amp|None}"""
    import multiprocessing as mp
    unit = unit or P.unit_map()
    vals = {}
    for ns in C.STATIONS:
        store = {}
        for sp in store_paths:
            f = f"{sp}/{ns}.npz"
            if os.path.exists(f):
                store.update(load_store_file(f))
        ks = [k for k in store if k in keys]
        if not ks:
            continue
        _G.clear(); _G.update(store=store, unit=unit, runs_fn=runs_fn, est=est, nodes=nodes or {})
        with mp.get_context("fork").Pool(workers) as pool:
            for k, o in pool.imap_unordered(_worker, ks, chunksize=1):
                vals[k] = o
        lg(f"  inject {ns}: {len(ks)} series")
    return vals


def cell_vector(vals, cells, run, fallback_nat=True):
    """-> d (dv/v %) per cell for one run (NaN where missing); runs absent for a series fall back to 'nat'"""
    key = (cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch).to_numpy()
    lut = {}
    d = np.full(len(cells), np.nan)
    for i, (k, e) in enumerate(zip(key, cells.episode)):
        o = vals.get(k, {})
        q = o.get(run)
        if q is None and fallback_nat:
            q = o.get("nat")
        if q is None or not len(q):
            continue
        if (k, run) not in lut:
            lut[(k, run)] = dict(zip(q.episode, q.eps))
        v = lut[(k, run)].get(e)
        if v is not None and np.isfinite(v):
            d[i] = -100 * float(v)
    return d


def store_paths_P():
    """existing P-weighted stores of the region (SVI: main + ext)"""
    out = [f"{C.TB}/store"]
    if C.REGION == "SVI" and os.path.isdir(f"{C.TB}/store_ext"):
        out.append(f"{C.TB}/store_ext")
    return out
