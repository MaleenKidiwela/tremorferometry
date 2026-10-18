#!/usr/bin/env python
"""s11_r4_dive -- S11 round 4, task 2 (PLAN.md §8.2): the SVI 2019-11 candidate.
(a) carriers (per-cell contributions to DeltaD(2019-11) in the headline combined model H); (b) subsets near / far / borehole /
PGC site / non-PGC; (c) own scrambled controls (B2 >= 10-event control nodes, B1, B0) + look-elsewhere; (d) leave-one-family-
out, leave-one-station-out; (e) lapse (graph re-measured at 2-3.5 and 3.5-6 s, apparent model, day bootstrap); (f) ballistic
arrival times; (g) instrument / timing (clock terms, realignment lags, shift-free estimator, bias per ms, metadata, data days);
(h) tremor and earthquakes; (i) the early burst (q8_rtr streaks, carrier-family activity, sub-episode split). Verdict by the
pre-registered rule. SVI only.
    python s11_r4_dive.py [--workers 12]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json, io, urllib.request, zlib
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert os.environ.get("S11_REGION", "SVI").upper() == "SVI"
import s10_common as S
import s11_common as C
import s11_dist_common as D
import s11_r4_common as R4
import s11_r3_common as R
import s11_r3_fit as LF
import s11_r3_inject as J
import s11_r4_combined as CB
from s11_measure import Ref3

E19 = "2019-11"
BASE = R4.BASE
_G = {}


def dD(F, e=E19):
    g = LF.dvec(F, e)
    if g is None:
        return np.nan, np.nan, None
    return float(g @ F["th"]), float(np.sqrt(g @ F["C"] @ g)), g


def fit_dD(cells, label, rows, lg, **kw):
    cells, _ = R4.identifiable(cells.reset_index(drop=True))
    if E19 not in set(cells.episode) or cells.episode.nunique() < 3:
        rows.append(dict(item=label, dD=np.nan, sig=np.nan, n_cells=len(cells), n_cells_e=0)); return np.nan, np.nan
    F = R4.Comb(cells, **kw).solve(); v, s, _ = dD(F)
    rows.append(dict(item=label, dD=v, sig=s, z=v / s if s > 0 else np.nan, n_cells=len(cells), n_cells_e=int((cells.episode == E19).sum()),
                     n_stations_e=int(cells[cells.episode == E19].netsta.nunique())))
    lg(f"  {label}: ΔD(2019-11) {v:+.3f} ± {s:.3f} % ({len(cells)} cells, {int((cells.episode == E19).sum())} in 2019-11)")
    return v, s


# ------------------------------------------------------------------------------ (e) graph at other windows, day bootstrap
def _win_worker(key):
    try:
        st = _G["store"][key]; nodes = _G["nodes"].get(key, [])
        if len(nodes) < 3:
            return key, None
        unit = _G["unit"]; rng = np.random.default_rng(zlib.crc32(key.encode()))
        ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
        X = st["St"].astype(np.float64); Ws = st["Ws"]
        out = {}
        for wn in ("25", "2035", "3560"):
            win = C.WINS[wn]
            def stacks(sel_days):
                Sx = np.full((len(nodes), 3, X.shape[-1]), np.nan)
                for k, e in enumerate(nodes):
                    m = sel_days[k]
                    for c in range(3):
                        if not np.isfinite(st["G1"][c]).all():
                            continue
                        w = Ws[m, c]
                        if w.sum() > 0:
                            Sx[k, c] = (w @ X[m, c]) / w.sum()
                return Sx
            base = [np.flatnonzero(ep == e) for e in nodes]
            pi, pj, dd, cc = J.graph_pairs(stacks(base), st["tS"], win)
            v = J.graph_solve(len(nodes), pi, pj, dd, cc)[0]
            vb = []
            for b in range(50):
                sel = [rng.choice(ix, len(ix), replace=True) for ix in base]
                pi, pj, dd, cc = J.graph_pairs(stacks(sel), st["tS"], win)
                vb.append(J.graph_solve(len(nodes), pi, pj, dd, cc)[0])
            out[wn] = (np.asarray(v), np.nanstd(np.array(vb), 0, ddof=1))
        return key, (list(nodes), out)
    except Exception as e:
        import traceback
        print("ERROR", key, e, traceback.format_exc()[-400:], flush=True); return key, None


def lapse_graph(H, workers, lg):
    import multiprocessing as mp
    import s11_paths as P
    unit = P.unit_map()
    nr = H[~H.far]
    E = R4.r4_episode(False)
    nodes = {f"{f}|{n}|{e}": list(g.episode) for (f, n, e), g in E[(E.n >= R4.NODE_MIN) & E.b1_cc.notna()].groupby(["fam", "netsta", "epoch"])}
    D.set_grid(D.NEAR_GRID)
    vals = {}
    for ns in sorted(nr.netsta.unique()):
        store = {}
        for sp in (f"{R.TB3}/store_ccw", f"{R.TB3}/store_ccw_ext"):
            if os.path.exists(f"{sp}/{ns}.npz"):
                store.update(J.load_store_file(f"{sp}/{ns}.npz"))
        keys = [k for k in store if k in set(nr.series)]
        _G.clear(); _G.update(store=store, unit=unit, nodes=nodes)
        with mp.get_context("fork").Pool(workers) as pool:
            for k, o in pool.imap_unordered(_win_worker, keys, chunksize=1):
                if o is not None:
                    vals[k] = o
        lg(f"  lapse graph {ns}: {len(keys)} series")
    return vals


# ------------------------------------------------------------------------------ (i) sub-episode split
def _split_worker(key):
    try:
        st = _G["store"][key]; unit = _G["unit"]; rng = np.random.default_rng(zlib.crc32(key.encode()))
        days = st["days"]; ep = np.array([unit.get(d, None) for d in days], dtype=object)
        m19 = ep == E19
        if not m19.any():
            return key, None
        X = st["St"].astype(np.float64); Ws = st["Ws"]; tS = st["tS"]
        cs = [c for c in range(3) if np.isfinite(st["G1"][c]).all()]
        def stack(ix):
            Sx = np.full((3, X.shape[-1]), np.nan)
            for c in cs:
                w = Ws[ix, c]
                if w.sum() > 0:
                    Sx[c] = (w @ X[ix, c]) / w.sum()
            return Sx
        ref = stack(np.flatnonzero(~m19)); Rf = Ref3(ref, tS, C.WINS["25"])
        mid19 = R4.U.mid[E19]
        nb = [e for e in pd.unique(ep[~m19]) if isinstance(e, str) and e in R4.U.index and abs((R4.U.mid[e] - mid19).days) <= 3 * 365.25]
        vnb = [Rf.measure(stack(np.flatnonzero(ep == e))[None], cs)[0][0] for e in nb]
        vnb = [v for v in vnb if np.isfinite(v) and abs(v) < C.RAIL]
        if len(vnb) < 2:
            return key, None
        base = float(np.mean(vnb))
        out = {}
        early = m19 & (days <= "2019-11-10"); late = m19 & (days >= "2019-11-11")
        for name, mm in (("all", m19), ("early", early), ("late", late)):
            ix = np.flatnonzero(mm)
            if len(ix) == 0:
                continue
            v = Rf.measure(stack(ix)[None], cs)[0][0]
            vb = [Rf.measure(stack(rng.choice(ix, len(ix), replace=True))[None], cs)[0][0] for _ in range(100)]
            vb = np.array([x for x in vb if np.isfinite(x) and abs(x) < C.RAIL])
            ne = float(st["Nd"][ix].sum()) if st.get("Nd") is not None else np.nan
            out[name] = (-100 * (v - base), 100 * np.std(vb, ddof=1) if len(vb) > 10 else np.nan, int(len(ix)), ne)
        return key, out
    except Exception as e:
        import traceback
        print("ERROR", key, e, traceback.format_exc()[-400:], flush=True); return key, None


def split_episode(H, workers, lg):
    import multiprocessing as mp
    import s11_paths as P
    unit = P.unit_map(); nr = H[~H.far]
    D.set_grid(D.NEAR_GRID)
    res = {}
    for ns in sorted(nr.netsta.unique()):
        store = {}
        for sp in (f"{R.TB3}/store_ccw", f"{R.TB3}/store_ccw_ext"):
            if os.path.exists(f"{sp}/{ns}.npz"):
                store.update(J.load_store_file(f"{sp}/{ns}.npz"))
        keys = [k for k in store if k in set(nr.series)]
        _G.clear(); _G.update(store=store, unit=unit)
        with mp.get_context("fork").Pool(workers) as pool:
            for k, o in pool.imap_unordered(_split_worker, keys, chunksize=1):
                if o:
                    res[k] = o
    rows = []
    for k, o in res.items():
        for name, (v, s, nd, ne) in o.items():
            rows.append(dict(series=k, netsta=k.split("|")[1], fam=int(k.split("|")[0]), part=name, dvv=v, sig=s, n_days=nd, n_events=ne))
    T = pd.DataFrame(rows)
    # Amendment 2: sigma = event-bootstrap sigma of the 2019-11 cell x sqrt(n_events(2019-11) / n_events(part)); the day bootstrap
    # (kept as sig_dayboot) is degenerate for parts with 2-3 days
    ce = H[(H.episode == E19) & ~H.far].set_index("series").sig
    T["sig_dayboot"] = T.sig
    nall = T[T.part == "all"].set_index("series").n_events
    T["sig"] = [ce.get(k, np.nan) * np.sqrt(nall.get(k, np.nan) / ne) if ne and ne > 0 else np.nan for k, ne in zip(T.series, T.n_events)]
    T.to_csv(f"{R4.TB4}/dive_split_cells.csv", index=False)
    out = []
    for name, g in T[np.isfinite(T.sig) & (T.sig > 0) & np.isfinite(T.dvv)].groupby("part"):
        w = 1 / g.sig ** 2; m = float(np.sum(w * g.dvv) / w.sum()); s = float(1 / np.sqrt(w.sum()))
        chi = float(np.sum(w * (g.dvv - m) ** 2) / max(len(g) - 1, 1)); s *= np.sqrt(max(chi, 1.0))
        out.append(dict(part=name, dvv=m, sig=s, n_series=len(g), chi2=chi))
    O = pd.DataFrame(out); O.to_csv(f"{R4.TB4}/dive_split.csv", index=False)
    lg("sub-episode split (dv/v of the part minus the mean of neighbouring episodes, pooled; %): " + "; ".join(
        f"{r.part} {r.dvv:+.3f} ± {r.sig:.3f} (n {r.n_series})" for r in O.itertuples()))
    return O


# ------------------------------------------------------------------------------ helpers
def fdsn_epochs(ns, t0="2019-07-01", t1="2020-03-31"):
    net, sta = ns.split(".")
    url = (f"https://service.earthscope.org/fdsnws/station/1/query?net={net}&sta={sta}&level=channel&format=text"
           f"&starttime={t0}&endtime={t1}")
    try:
        txt = urllib.request.urlopen(url, timeout=40).read().decode()
        Dd = pd.read_csv(io.StringIO(txt), sep="|"); Dd.columns = [c.strip().lstrip("#").strip() for c in Dd.columns]
        Dd = Dd[Dd.Channel.str[:2].isin(["HH", "EH", "BH"])]
        return Dd[["Network", "Station", "Location", "Channel", "StartTime", "EndTime", "SampleRate", "SensorDescription"]].assign(netsta=ns)
    except Exception as e:
        return pd.DataFrame([dict(netsta=ns, error=str(e)[:200])])


def km(la1, lo1, la2, lo2):
    return D.km(la1, lo1, la2, lo2)


# ------------------------------------------------------------------------------ (g7) node-level realignment (Amendment 2)
def _ra_worker(key):
    try:
        from s10_measure import realign
        st = _G["store"][key]; nodes = _G["nodes"].get(key, [])
        if len(nodes) < 3:
            return key, None
        unit = _G["unit"]
        ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
        X = st["St"].astype(np.float64); Ws = st["Ws"]; tS = st["tS"]; G1 = st["G1"]
        Sx = np.full((len(nodes), 3, X.shape[-1]), np.nan)
        for k, e in enumerate(nodes):
            m = ep == e
            for c in range(3):
                if not np.isfinite(G1[c]).all():
                    continue
                w = Ws[m, c]
                if w.sum() > 0:
                    Sx[k, c] = (w @ X[m, c]) / w.sum()
        v0 = J.graph_from_stacks(Sx, tS)[0]
        Sr, lag, cc = realign(Sx.astype(np.float32), G1, tS)
        v1 = J.graph_from_stacks(Sr.astype(np.float64), tS)[0]
        return key, (list(nodes), np.asarray(v0), np.asarray(v1), 1e3 * np.asarray(lag), np.asarray(cc))
    except Exception as e:
        import traceback
        print("ERROR", key, e, traceback.format_exc()[-400:], flush=True); return key, None


def realign_graph(H, workers, lg):
    import multiprocessing as mp
    import s11_paths as P
    unit = P.unit_map(); vals = {}
    for far in (False, True):
        q = H[H.far == far]
        if not len(q):
            continue
        E = R4.r4_episode(far)
        nodes = {f"{f}|{n}|{e}": list(g.episode) for (f, n, e), g in E[(E.n >= R4.NODE_MIN) & E.b1_cc.notna()].groupby(["fam", "netsta", "epoch"])}
        D.set_grid(D.FAR_GRID if far else D.NEAR_GRID)
        paths = [f"{R4.TB4}/store_ccw_far"] if far else [f"{R.TB3}/store_ccw", f"{R.TB3}/store_ccw_ext"]
        for ns in sorted(q.netsta.unique()):
            store = {}
            for sp in paths:
                if os.path.exists(f"{sp}/{ns}.npz"):
                    store.update(J.load_store_file(f"{sp}/{ns}.npz"))
            keys = [k for k in store if k in set(q.series)]
            _G.clear(); _G.update(store=store, unit=unit, nodes=nodes)
            with mp.get_context("fork").Pool(workers) as pool:
                for k, o in pool.imap_unordered(_ra_worker, keys, chunksize=1):
                    if o is not None:
                        vals[k] = o
        lg(f"  node realignment {'far' if far else 'near'}: {sum(1 for k in vals if k in set(q.series))} series")
    D.set_grid(D.NEAR_GRID)
    return vals


def core(H, kw, tag, lg, majors):
    """sections (a)-(d) of PLAN §8.2 for one model form (kw passed to R4.Comb); tables get the suffix tag"""
    H = H.copy(); V = {}
    L = R4.Comb(H, **kw); F = L.solve()
    v0, s0, g = dD(F)
    lg(f"[{tag or 'H'}] ΔD(2019-11) {v0:+.3f} ± {s0:.3f} % (z {v0 / s0:.2f})")
    V["H"] = dict(dD=v0, sig=s0)
    w_op = g @ L.Op
    H["contrib"] = w_op * H.d.to_numpy(); H["op_w"] = w_op
    H.to_csv(f"{R4.TB4}/dive_cells_contrib{tag}.csv", index=False)
    car = []
    for col in ("fam", "netsta", "site", "far"):
        t = H.groupby(col).contrib.sum().sort_values(ascending=False)
        for k, x in t.items():
            car.append(dict(by=col, key=str(k), contrib=float(x), share=float(x / v0)))
    CA = pd.DataFrame(car); CA.to_csv(f"{R4.TB4}/dive_carriers{tag}.csv", index=False)
    st_share = CA[CA.by == "netsta"].set_index("key").share
    lg("carriers by station (share of ΔD): " + ", ".join(f"{k} {x:+.2f}" for k, x in st_share.items()))
    lg("carriers by site: " + ", ".join(f"{r.key} {r.share:+.2f}" for r in CA[CA.by == "site"].itertuples()))
    lg("carriers near/far: " + ", ".join(f"{r.key} {r.share:+.2f}" for r in CA[CA.by == "far"].itertuples()))
    fam_c = CA[CA.by == "fam"].copy(); fam_c["key"] = fam_c.key.astype(int)
    pos = fam_c[fam_c.share > 0].sort_values("share", ascending=False)          # same sign as DeltaD
    cs_ = pos.share.cumsum() / pos.share.sum()
    carriers = list(pos.key[cs_ <= 0.8]) + list(pos.key[cs_ > 0.8][:1])
    lg(f"carrier families (80 % of the same-sign contribution): {carriers}; top-10 family shares {fam_c.sort_values('share', ascending=False).head(10).set_index('key').share.round(2).to_dict()}")
    e19 = H[H.episode == E19]
    lg(f"2019-11 cells: {len(e19)} at {e19.netsta.nunique()} stations, {e19.fam.nunique()} families; near {int((~e19.far).sum())}, far {int(e19.far.sum())}")
    rows = []
    fit_dD(H, f"model{tag}", rows, lg, **kw)
    fit_dD(H[~H.far], f"near only{tag}", rows, lg, **kw)
    fit_dD(H[H.far], f"far only{tag}", rows, lg, **kw)
    fit_dD(H[H.netsta.str.startswith("PB.")], f"borehole only (PB, near + far){tag}", rows, lg, **kw)
    fit_dD(H[H.netsta.str.startswith("PB.") & ~H.far], f"near boreholes only{tag}", rows, lg, **kw)
    fit_dD(H[H.site == "PGC"], f"PGC site only{tag}", rows, lg, **kw)
    fit_dD(H[H.site != "PGC"], f"without the PGC site{tag}", rows, lg, **kw)
    fit_dD(H[H.netsta != "CN.PGC"], f"without CN.PGC (vault){tag}", rows, lg, **kw)
    cc = R4.control_cells(H); fit_dD(cc, f"scrambled control (B2 control nodes >= 10 events){tag}", rows, lg, **kw)
    cb = R4.control_cells(CB._ctrl_b1(H, "ell40")); fit_dD(cb, f"scrambled control (B1 reference path){tag}", rows, lg, **kw)
    T, comp = R4.q3(F, majors)
    T["z"] = T.dD / T.sig
    zobs = abs(v0 / s0)
    from scipy import stats
    nge = int((np.abs(T.z) >= zobs).sum()); expct = float(len(T) * 2 * stats.norm.sf(zobs))
    V["look_elsewhere"] = dict(n_episodes=len(T), n_ge=nge, expected_gaussian=expct, chi2_all=float(np.mean(T.z ** 2)),
                               p_any=float(1 - (1 - 2 * stats.norm.sf(zobs)) ** len(T)))
    lg(f"look-elsewhere: {nge} of {len(T)} episodes with |z| ≥ {zobs:.2f} (Gaussian expectation {expct:.2f}); mean z² {np.mean(T.z ** 2):.2f}")
    T.to_csv(f"{R4.TB4}/dive_all_episodes_z{tag}.csv", index=False)
    lf = []
    for f in sorted(e19.fam.unique()):
        x = []; fit_dD(H[H.fam != f], f"without family {f}", x, lambda m: None, **kw)
        lf.append(dict(kind="family", left_out=int(f), **{k: x[0][k] for k in ("dD", "sig")}))
    for s_ in sorted(e19.netsta.unique()):
        x = []; fit_dD(H[H.netsta != s_], f"without {s_}", x, lambda m: None, **kw)
        lf.append(dict(kind="station", left_out=s_, **{k: x[0][k] for k in ("dD", "sig")}))
    LO = pd.DataFrame(lf); LO["z"] = LO.dD / LO.sig; LO.to_csv(f"{R4.TB4}/dive_leaveout{tag}.csv", index=False)
    lofo = LO[LO.kind == "family"]; sg = np.sign(v0)
    V["lofo"] = dict(min=float(lofo.dD.min()), max=float(lofo.dD.max()), min_z=float((sg * lofo.z).min()),
                     worst=int(lofo.left_out[(sg * lofo.z).idxmin()]) if len(lofo) else None)
    los = LO[LO.kind == "station"]
    V["loso"] = dict(min=float(los.dD.min()), max=float(los.dD.max()), min_z=float((sg * los.z).min()))
    lg(f"leave-one-family-out: ΔD {lofo.dD.min():+.3f} … {lofo.dD.max():+.3f} (min signed z {V['lofo']['min_z']:.2f}, family {V['lofo']['worst']}); "
       f"leave-one-station-out {los.dD.min():+.3f} … {los.dD.max():+.3f} (min signed z {V['loso']['min_z']:.2f})")
    return dict(V=V, rows=rows, st_share=st_share, carriers=carriers, e19=e19, v0=v0, s0=s0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--test", action="store_true"); a = ap.parse_args()
    lg = R4.log("dive"); t0 = time.time()
    majors = R4.majors()
    H = CB.load_H(lg)
    CORE = {tag: core(H, kw, tag, lg, majors) for tag, kw in (("", {}), ("_noSc", dict(Sc=False)))}
    V = dict(CORE[""]["V"]); V["noSc"] = CORE["_noSc"]["V"]
    v0, s0 = CORE[""]["v0"], CORE[""]["s0"]
    rows = CORE[""]["rows"] + CORE["_noSc"]["rows"]
    st_share = CORE[""]["st_share"]
    # carriers for (f)-(i): from the model form whose DeltaD is positive (the candidate's sign; Amendment 1), else H
    cand = "" if v0 > 0 or CORE["_noSc"]["v0"] <= 0 else "_noSc"
    carriers = CORE[cand]["carriers"]; e19 = CORE[cand]["e19"]; st_share_c = CORE[cand]["st_share"]
    kw_c = {} if cand == "" else dict(Sc=False); vc = CORE[cand]["v0"]
    lg(f"carrier set for (f)-(i) taken from model form '{cand or 'H'}' (ΔD {CORE[cand]['v0']:+.3f})")
    fit_dD(H, "H without receiver term", rows, lg, Rr=False)
    fit_dD(H, "H without Sc and receiver term", rows, lg, Rr=False, Sc=False)
    rt = LF.temporal(R4.identifiable(H[~H.far])[0]).solve(); vv, ss, _ = dD(rt)
    rows.append(dict(item="round-3 temporal model on near B2 (r4 inflation)", dD=vv, sig=ss, z=vv / ss))
    lg(f"  round-3 temporal model, near B2 cells: {vv:+.3f} ± {ss:.3f}")
    try:
        b0c = R.head_cells(control=True); Fb0c = LF.temporal(b0c).solve(); vv, ss, _ = dD(Fb0c)
        rows.append(dict(item="scrambled control (B0, round-3 temporal model, near)", dD=vv, sig=ss, z=vv / ss if ss > 0 else np.nan))
        lg(f"  scrambled control B0: {vv:+.3f} ± {ss:.3f}")
    except Exception as e:
        lg(f"B0 control failed: {e}")
    # (e) lapse
    vals = lapse_graph(H, a.workers, lg)
    nr = H[~H.far].reset_index(drop=True)
    lap = {}
    for wn in ("25", "2035", "3560"):
        d = np.full(len(nr), np.nan); s_ = np.full(len(nr), np.nan)
        for i, (k, e) in enumerate(zip(nr.series, nr.episode)):
            if k in vals and e in vals[k][0]:
                j = vals[k][0].index(e); v, sb = vals[k][1][wn]
                if np.isfinite(v[j]) and abs(v[j]) < C.RAIL:
                    d[i] = -100 * v[j]; s_[i] = 100 * sb[j]
        c = nr.copy(); c["d"] = d; c["sig"] = s_ * c.infl
        c = c[np.isfinite(c.d) & np.isfinite(c.sig) & (c.sig > 0)]
        key = c.series; c = c[key.map(key.value_counts()) >= 3].reset_index(drop=True)
        c, _ = R4.identifiable(c)
        Fa = R4.Comb(c, apparent=True, Sc=False, Rr=False, MB=False, Ssite=False).solve(); va, sa, _ = dD(Fa)   # Amendment 2
        cw = R4.with_setup(c, "ell40", win=wn if wn != "25" else C.HEAD)
        cw = cw[np.isfinite(cw.a_d) & np.isfinite(cw.a_s)].reset_index(drop=True)
        Fk = R4.Comb(cw, Sc=False, Rr=False, MB=False).solve(); vk, sk_, _ = dD(Fk)
        ce = cw[cw.episode == E19]; ww = 1 / ce.sig ** 2
        lap[wn] = dict(apparent=va, apparent_sig=sa, kernel=vk, kernel_sig=sk_, n_cells=len(c), a_d=float(np.average(ce.a_d, weights=ww)),
                       a_s=float(np.average(ce.a_s, weights=ww)))
        lg(f"  lapse {wn}: apparent ΔE(2019-11) {va:+.3f} ± {sa:.3f}; with own kernel ΔD {vk:+.3f} ± {sk_:.3f}; 2019-11 cells mean a_d {lap[wn]['a_d']:.3f} a_s {lap[wn]['a_s']:.3f}")
    r_ = lap["3560"]["apparent"] / lap["2035"]["apparent"]
    sr = abs(r_) * np.hypot(lap["3560"]["apparent_sig"] / lap["3560"]["apparent"], lap["2035"]["apparent_sig"] / lap["2035"]["apparent"])
    pdp = lap["3560"]["a_d"] / lap["2035"]["a_d"]; psp = lap["3560"]["a_s"] / lap["2035"]["a_s"]
    zd, zs = (r_ - pdp) / sr, (r_ - psp) / sr
    call = "deep-like" if abs(zs) - abs(zd) >= 1 else ("shallow-like" if abs(zd) - abs(zs) >= 1 else "undecided")
    V["lapse"] = dict(windows=lap, ratio=r_, ratio_sig=sr, pred_deep=pdp, pred_shallow=psp, z_deep=zd, z_shallow=zs, call=call)
    lg(f"lapse ratio 3.5-6/2-3.5 {r_:.2f} ± {sr:.2f}; deep {pdp:.2f} (z {zd:+.2f}), shallow {psp:.2f} (z {zs:+.2f}) -> {call}")
    try:
        b0 = {}
        for wn in ("2035", "3560", "25"):
            c = R.head_cells(win=wn); c, _ = R4.identifiable(c)
            Fa = R4.Comb(c, apparent=True, Sc=False, Rr=False, MB=False, Ssite=False).solve(); va, sa, _ = dD(Fa); b0[wn] = (va, sa)
        V["lapse_B0"] = {k: dict(apparent=x[0], sig=x[1]) for k, x in b0.items()}
        lg(f"  lapse B0 path apparent: " + ", ".join(f"{k} {x[0]:+.3f} ± {x[1]:.3f}" for k, x in b0.items()))
    except Exception as e:
        lg(f"B0 lapse failed: {e}")
    # (f) ballistic arrival times
    BL = f"{BASE}/ballistic"
    ck3 = pd.read_csv(f"{BL}/round3/svi/tables/clock_terms.csv"); ck3 = ck3[(ck3.bin == "EP") & (ck3.lab == E19)]
    ck2 = pd.read_csv(f"{BL}/round2/tables/clock_terms.csv"); ck2 = ck2[(ck2.bin == "EP") & (ck2.lab == E19)]
    cm = pd.read_csv(f"{BL}/round3/tables/common_majors.csv"); cm = cm[(cm.region == "SVI") & (cm.episode == E19)]
    ec = pd.read_csv(f"{BL}/tables/episode_common.csv"); ec = ec[(ec.episode == E19) & (ec.kind == "data") & (ec.subset == "all")]
    ls_ = pd.read_parquet(f"{BL}/round3/svi/tables/lags_station.parquet")
    ls_ = ls_[(ls_.kind == "data") & (ls_.bin == "EP") & ~ls_.railS]
    ls_["r"] = ls_.lagS - ls_.groupby(["fam", "netsta", "epoch"]).lagS.transform("median")
    lag19 = ls_[ls_.lab == E19].groupby("netsta").r.agg(["median", "size"]).rename(columns={"median": "dS_ms_vs_series_median", "size": "n"})
    ck3.to_csv(f"{R4.TB4}/dive_clock_terms_2019_11.csv", index=False); lag19.to_csv(f"{R4.TB4}/dive_ballistic_S_2019_11.csv")
    lg("ballistic clock terms 2019-11 (round 3): " + ", ".join(f"{r.station} {r.K:+.1f} ± {r.sigK:.1f} ms" for r in ck3.itertuples()))
    lg("ballistic S lag 2019-11 vs series median (ms): " + ", ".join(f"{k} {r.dS_ms_vs_series_median:+.1f} (n {r.n})" for k, r in lag19.iterrows()))
    lg("ballistic round-1 episode-common terms 2019-11: " + "; ".join(f"{r.item} {r.common:+.1f} ms" for r in ec.itertuples()))
    if len(cm):
        lg("ballistic region S−P 2019-11: " + "; ".join(f"{r.kind} {r.dvv_constVpVs_pct:+.3f} ± {r.sig_pct:.3f} %" for r in cm.itertuples()))
    # expected S-P dv/v if the deep-weighted change were real: fraction of vertical S time in 25-45 km
    P_ = pd.read_csv(D.VEL)
    Fam = D.near_tables()["F"].set_index("fam")
    zf = float(np.average(Fam.depth_km.reindex(carriers).dropna())) if carriers else 35.0
    top, bot = P_.ztop_km.to_numpy(float), np.where(np.isfinite(P_.zbot_km.to_numpy(float)), P_.zbot_km.to_numpy(float), 200.0)
    vs = P_.vs.to_numpy(float)
    tt = lambda z0, z1: float(np.sum(np.clip(np.minimum(bot, z1) - np.maximum(top, z0), 0, None) / vs))
    frac = tt(25.0, zf) / tt(0.0, zf)
    V["ballistic"] = dict(clock_terms_ms={r.station: float(r.K) for r in ck3.itertuples()}, clock_sig_ms={r.station: float(r.sigK) for r in ck3.itertuples()},
                          S_lag_ms=lag19.dS_ms_vs_series_median.to_dict(),
                          SP_common_pct=float(cm[cm.kind == "data"].dvv_constVpVs_pct.iloc[0]) if len(cm) else np.nan,
                          SP_common_sig=float(cm[cm.kind == "data"].sig_pct.iloc[0]) if len(cm) else np.nan,
                          frac_S_time_25_45=frac, carrier_depth_km=zf, expected_SP_pct=vc * frac)
    lg(f"expected S−P dv/v if ΔD were deep: {vc:+.3f} × {frac:.2f} (S time fraction 25–{zf:.0f} km) = {vc * frac:+.3f} %")
    # (g) instrument / timing
    clockst = [s for s, k in V["ballistic"]["clock_terms_ms"].items() if abs(k) >= 3.0]
    x = []; v_nc, s_nc = fit_dD(H[~H.netsta.isin(clockst)], f"without stations with |clock| ≥ 3 ms ({clockst})", x, lg, **kw_c)
    rows += x
    E4 = R4.r4_episode(False)
    E4["lag_rel"] = E4.lag_med_ms - E4.groupby(["fam", "netsta", "epoch"]).lag_med_ms.transform("median")
    lagr = E4[E4.episode == E19].groupby("netsta").lag_rel.agg(["median", "mean", "size"])
    lagr.to_csv(f"{R4.TB4}/dive_realign_lag_2019_11.csv")
    lg("realignment lag of 2019-11 events minus the series median (ms): " + ", ".join(f"{k} {r['median']:+.2f} (n {int(r['size'])})" for k, r in lagr.iterrows()))
    allep = E4.groupby(["netsta", "episode"]).lag_rel.median().reset_index()
    allep.to_csv(f"{R4.TB4}/dive_realign_lag_all.csv", index=False)
    # shift-free estimator (B0 path): eps2p vs eps3c, same cells, round-3 temporal model
    Eb = pd.read_parquet(f"{C.TB}/episode_raw_all.parquet")
    c3 = R.head_cells()
    m = c3.merge(Eb[["fam", "netsta", "epoch", "episode", "eps2p_25", "shift2p_25"]], on=["fam", "netsta", "epoch", "episode"], how="left")
    F3 = LF.temporal(m).solve(); v3, s3, _ = dD(F3)
    m2 = m[np.isfinite(m.eps2p_25)].copy(); m2["d"] = -100 * m2.eps2p_25
    key = m2.series; m2 = m2[key.map(key.value_counts()) >= 3].reset_index(drop=True)
    F2 = LF.temporal(m2).solve(); v2, s2, _ = dD(F2)
    m3 = m2.copy(); m3["d"] = m.set_index(["series", "episode"]).d.reindex(list(zip(m3.series, m3.episode))).to_numpy()
    F3b = LF.temporal(m3).solve(); v3b, s3b, _ = dD(F3b)
    m["shift_rel_ms"] = 1e3 * (m.shift2p_25 - m.groupby("series").shift2p_25.transform("median"))
    sh = m[m.episode == E19].groupby("netsta").shift_rel_ms.agg(["median", "size"])
    V["shift_free"] = dict(B0_eps3c=v3, B0_eps3c_sig=s3, B0_eps3c_same_cells=v3b, B0_eps3c_same_cells_sig=s3b, B0_eps2p=v2, B0_eps2p_sig=s2,
                           shift_2019_11_ms=sh["median"].to_dict())
    lg(f"shift-free estimator (B0, round-3 temporal model): eps3c ΔD {v3:+.3f} ± {s3:.3f}; same cells {v3b:+.3f} ± {s3b:.3f}; eps2p {v2:+.3f} ± {s2:.3f}")
    lg("2019-11 residual shift (2-parameter fit) minus series median (ms): " + ", ".join(f"{k} {r['median']:+.1f} (n {int(r['size'])})" for k, r in sh.iterrows()))
    # (g7) node-level realignment before the graph stretch (Amendment 2; post code-test)
    RA = realign_graph(H, a.workers, lg)
    d0 = np.full(len(H), np.nan); d1 = np.full(len(H), np.nan); lg_ = np.full(len(H), np.nan)
    for i, (k, e) in enumerate(zip(H.series, H.episode)):
        if k in RA and e in RA[k][0]:
            j = RA[k][0].index(e); d0[i] = -100 * RA[k][1][j]; d1[i] = -100 * RA[k][2][j]; lg_[i] = RA[k][3][j]
    ok = np.isfinite(d0) & np.isfinite(d1)
    lg(f"node realignment: natural graph reproduces the cell values to {np.nanmax(np.abs(d0[ok] - H.d.to_numpy()[ok])):.1e} % ({int(ok.sum())} of {len(H)} cells)")
    Hn = H[ok].copy(); Hn["node_lag_ms"] = lg_[ok]
    Hn["node_lag_rel"] = Hn.node_lag_ms - Hn.groupby("series").node_lag_ms.transform("median")
    nl = Hn[Hn.episode == E19].groupby("netsta").node_lag_rel.agg(["median", "size"])
    nl.to_csv(f"{R4.TB4}/dive_node_lag_2019_11.csv")
    Hn.groupby(["netsta", "episode"]).node_lag_rel.median().reset_index().to_csv(f"{R4.TB4}/dive_node_lag_all.csv", index=False)
    lg("2019-11 node-stack lag vs the series median (ms, before realignment): " + ", ".join(f"{k} {r['median']:+.2f} (n {int(r['size'])})" for k, r in nl.iterrows()))
    Hr = Hn.copy(); Hr["d"] = d1[ok]
    ra = {}
    for tag, kw in (("", {}), ("_noSc", dict(Sc=False))):
        x = []; a0_, b0_ = fit_dD(Hn, f"cells with a realigned value{tag} (not realigned)", x, lg, **kw)
        a1_, b1_ = fit_dD(Hr, f"node-realigned graph{tag}", x, lg, **kw); rows += x
        ra[tag] = dict(before=a0_, before_sig=b0_, after=a1_, after_sig=b1_)
    V["node_realign"] = dict(results=ra, node_lag_2019_11_ms=nl["median"].to_dict())
    # stretch bias per ms of residual shift (grand stacks of the near store, 2-5 s, joint 3C)
    import s11_r3_inject as J_
    bias = []
    D.set_grid(D.NEAR_GRID)
    for ns in sorted(H[~H.far].netsta.unique()):
        st_ = J_.load_store_file(f"{R.TB3}/store_ccw/{ns}.npz")
        for k, st in list(st_.items())[:40]:
            G1 = st["G1"]; cs = tuple(c for c in range(3) if np.isfinite(G1[c]).all())
            if not cs:
                continue
            f = np.fft.rfftfreq(S.NW, 1 / S.FS); ph = np.exp(-2j * np.pi * f * 0.001)
            Gs = np.full_like(G1, np.nan)
            for c in cs:
                Gs[c] = np.fft.irfft(np.fft.rfft(G1[c]) * ph, n=S.NW)
            e_, _ = Ref3(G1, st["tS"], C.WINS["25"]).measure(Gs[None], cs)
            bias.append(dict(netsta=ns, series=k, eps_per_ms=float(e_[0])))
    BI = pd.DataFrame(bias); BI.to_csv(f"{R4.TB4}/dive_shift_bias.csv", index=False)
    bpm = float(np.nanmedian(BI.eps_per_ms)) * 100
    V["bias_dvv_pct_per_ms_delay"] = -bpm
    lg(f"stretch bias of a 1 ms residual delay in 2-5 s: eps {bpm:+.3f} % → dv/v {-bpm:+.3f} % per ms (median of {len(BI)} series)")
    # metadata
    stations = sorted(set(st_share_c.index[st_share_c.abs() >= 0.05]) | set(e19.netsta[~e19.far]))
    meta = []
    for ns in stations:
        for e, a_, b_ in (S.EPOCHS.get(ns, []) or D.epochs_table(ns)):
            for dd in (a_, b_):
                if abs((pd.Timestamp(dd) - pd.Timestamp("2019-11-15")).days) <= 120:
                    meta.append(dict(netsta=ns, source="s10 epoch boundary", date=dd, detail=e))
        fd = fdsn_epochs(ns)
        if "error" in fd:
            meta.append(dict(netsta=ns, source="FDSN", date="", detail=f"query failed: {fd.error.iloc[0]}"))
        else:
            for r in fd.itertuples():
                for dd, lab in ((r.StartTime, "channel start"), (r.EndTime, "channel end")):
                    if isinstance(dd, str) and dd.strip() and abs((pd.Timestamp(dd[:10]) - pd.Timestamp("2019-11-15")).days) <= 120:
                        meta.append(dict(netsta=ns, source="FDSN", date=dd[:10], detail=f"{lab} {r.Location}.{r.Channel} {r.SensorDescription}"))
    ME = pd.DataFrame(meta) if meta else pd.DataFrame(columns=["netsta", "source", "date", "detail"])
    ME.to_csv(f"{R4.TB4}/dive_metadata.csv", index=False)
    lg(f"metadata events within ±120 d of 2019-11-15 at {stations}: {len(ME)} rows: " + "; ".join(f"{r.netsta} {r.date} {r.detail}" for r in ME.itertuples()))
    dd = []
    for ns in sorted(H.netsta.unique()):
        p = f"{R.TB3}/store_ccw/{ns}.npz" if ns in set(H.netsta[~H.far]) else f"{R4.TB4}/store_ccw_far/{ns}.npz"
        if os.path.exists(p):
            z = np.load(p); days = pd.to_datetime(pd.Series(z["days"]))
            for mo in ("2019-09", "2019-10", "2019-11", "2019-12", "2020-01"):
                dd.append(dict(netsta=ns, month=mo, n_day_stacks=int((days.dt.strftime("%Y-%m") == mo).sum())))
    pd.DataFrame(dd).pivot(index="netsta", columns="month", values="n_day_stacks").to_csv(f"{R4.TB4}/dive_data_days.csv")
    # (h) tremor + earthquakes
    Fc = Fam.reindex(carriers).dropna(subset=["lat"])
    la0, lo0 = float(Fc.lat.mean()), float(Fc.lon.mean())
    tr = pd.read_csv(f"{C.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["t"] = pd.to_datetime(tr.time)
    dmin = np.min([km(tr.lat.to_numpy(), tr.lon.to_numpy(), r.lat, r.lon) for r in Fc.itertuples()], 0)
    tr = tr[dmin <= 30.0]
    w = tr[(tr.t >= "2019-10-15") & (tr.t < "2019-12-16")]
    daily = w.groupby(w.t.dt.strftime("%Y-%m-%d")).size()
    daily.to_csv(f"{R4.TB4}/dive_tremor_daily.csv")
    U = R4.U
    tm = []
    for e in majors:
        if e not in U.index:
            continue
        s_, e_ = pd.Timestamp(U.start[e]), pd.Timestamp(U.end[e]) + pd.Timedelta(days=1)
        q = tr[(tr.t >= s_) & (tr.t < e_)]
        tm.append(dict(episode=e, start=str(U.start[e]), end=str(U.end[e]), n_tremor=len(q), n_days=int(q.t.dt.date.nunique())))
    TM = pd.DataFrame(tm); TM.to_csv(f"{R4.TB4}/dive_tremor_majors.csv", index=False)
    lg(f"PNSN tremor within 30 km of the carrier families: 2019-10-15…12-15 {len(w)} windows on {daily.size} days (peak {daily.idxmax() if len(daily) else ''} {daily.max() if len(daily) else 0}); "
       f"per major: " + ", ".join(f"{r.episode} {r.n_tremor}" for r in TM.itertuples()))
    eq = pd.read_parquet(f"{C.ROOT}/svi_science/results/g2_local/tables/eq_local.parquet")
    eq["tt"] = pd.to_datetime(eq.t)
    eq = eq[(eq.tt >= "2019-10-01") & (eq.tt < "2020-01-01") & (eq.type == "earthquake") & ~eq.removed.fillna(False).astype(bool)]
    eq["dist"] = km(eq.lat.to_numpy(), eq.lon.to_numpy(), la0, lo0)
    sel = eq[((eq.dist <= 100) & (eq.mag >= 2.0)) | ((eq.dist <= 150) & (eq.mag >= 3.0))].sort_values("tt")
    sel[["id", "t", "t_PT", "lat", "lon", "depth", "mag", "dist"]].to_csv(f"{R4.TB4}/dive_earthquakes.csv", index=False)
    lg(f"earthquakes 2019-10…12 (M ≥ 2 ≤ 100 km, M ≥ 3 ≤ 150 km of {la0:.2f}, {lo0:.2f}): {len(sel)}; largest "
       + ("; ".join(f"M{r.mag:.1f} {r.t_PT} {r.dist:.0f} km z {r.depth:.0f}" for r in sel.nlargest(5, "mag").itertuples()) if len(sel) else "none"))
    V["seismicity"] = dict(n_tremor_window=len(w), n_eq=len(sel), max_mag=float(sel.mag.max()) if len(sel) else None,
                           tremor_2019_11=int(TM.set_index("episode").n_tremor.get(E19, -1)) if len(TM) else None,
                           tremor_median_other_majors=float(TM[TM.episode != E19].n_tremor.median()) if len(TM) else None)
    # (i) burst
    st_rows = []
    import glob
    for f in sorted(glob.glob(f"{C.ROOT}/svi_science/results/q8_rtr/tables/streaks_all_abs*.parquet")):
        d_ = pd.read_parquet(f); d_ = d_[(d_.ep == E19) & (d_["mode"] == "data")]
        if not len(d_):
            continue
        d_ = d_.assign(st_utc=pd.to_datetime(d_.start)); d_ = d_[(d_.st_utc >= "2019-11-08") & (d_.st_utc < "2019-11-12")]
        d_["start_PT"] = (d_.st_utc - pd.Timedelta(hours=8)).dt.strftime("%Y-%m-%d %H:%M PST")
        d_["rate_per_h"] = d_.N_ev / d_.dur_h
        d_["table"] = os.path.basename(f)
        st_rows.append(d_[["table", "uid", "start", "start_PT", "dur_h", "N_ev", "n_fam", "rate_per_h", "v_ev", "x_start", "x_end", "cls"]])
    BU = pd.concat(st_rows, ignore_index=True) if st_rows else pd.DataFrame()
    BU.to_csv(f"{R4.TB4}/dive_bursts.csv", index=False)
    if len(BU):
        top = BU.sort_values("N_ev", ascending=False).iloc[0]
        lg(f"largest early burst: {top.start_PT} ({top.table}), {top.N_ev} events in {top.dur_h:.1f} h ({top.rate_per_h:.0f}/h), {top.n_fam} families, x {top.x_start:.0f}→{top.x_end:.0f} km, {top.cls}")
    X = pd.read_parquet(f"{C.MTB}/rows.parquet")
    pe = f"{C.TB}/ext/rows.parquet"
    if os.path.exists(pe):
        X = pd.concat([X, pd.read_parquet(pe)], ignore_index=True)
    X = X[(X.series == "data") & (X.date >= "2019-11-05") & (X.date <= "2019-12-02")]
    act = X.assign(carrier=X.fam.isin(carriers)).groupby(["date", "carrier"]).P.sum().unstack(fill_value=0)
    act.to_csv(f"{R4.TB4}/dive_family_activity.csv")
    lg("ΣP per UTC day 2019-11 (carrier / other families): " + ", ".join(f"{d_[5:]} {r.get(True, 0):.0f}/{r.get(False, 0):.0f}" for d_, r in act.iterrows()))
    SPL = split_episode(H, a.workers, lg)
    V["split"] = SPL.set_index("part")[["dvv", "sig", "n_series"]].to_dict("index") if len(SPL) else {}
    # ---------------- verdict
    RW = pd.DataFrame(rows); RW.to_csv(f"{R4.TB4}/dive_subsets.csv", index=False)
    get = lambda lab: RW.set_index("item").loc[lab] if lab in set(RW.item) else None
    vj = f"{R4.TB4}/r4_verdict.json"
    VJ = json.load(open(vj))["verdict"] if os.path.exists(vj) else {}
    car_st = list(st_share_c[st_share_c > 0].sort_values(ascending=False).index)
    cum = st_share_c[car_st].cumsum() / st_share_c[car_st].sum(); top50 = list(cum.index[: int(np.searchsorted(cum.to_numpy(), 0.5)) + 1])
    shv = {k: v for k, v in V["shift_free"]["shift_2019_11_ms"].items() if k in top50 and np.isfinite(v)}
    timing_i = bool(np.isfinite(v_nc) and abs(v_nc) < 2 * s_nc and abs(v_nc) < s_nc) if clockst else False
    timing_ii = bool(abs(v2) < s2 and v3b >= 3 * s3b and len(shv) >= 1 and all(abs(v) >= 3.0 for v in shv.values()))
    meta_hit = bool(len(ME) and ME[ME.netsta.isin(top50) & (ME.date.astype(str).str.len() >= 10)].date.apply(
        lambda x: abs((pd.Timestamp(x) - pd.Timestamp("2019-11-15")).days) <= 30).any())
    rc_ = V["node_realign"]["results"][cand]
    timing_iii = bool(np.isfinite(rc_["after"]) and abs(rc_["after"]) < rc_["after_sig"] and abs(rc_["before"]) >= 2 * rc_["before_sig"])
    artefact = bool(timing_i or timing_ii or timing_iii or meta_hit)
    for tag, vkey in (("", "H_ell40"), ("_noSc", "HnoSc_ell40")):            # pre-registered on H; H without Sc = Amendment 1 column
        vt, st_ = CORE[tag]["v0"], CORE[tag]["s0"]; VV = CORE[tag]["V"]
        gates_ok = bool(VJ.get(vkey, {}).get("wording_all_pass")) if VJ else None
        ctrl = get(f"scrambled control (B2 control nodes >= 10 events){tag}")
        near, far = get(f"near only{tag}"), get(f"far only{tag}")
        c_ok = bool(ctrl is not None and np.isfinite(ctrl.dD) and abs(ctrl.dD) < 2 * ctrl.sig)
        nf_ok = bool(near is not None and far is not None and np.sign(near.dD) == np.sign(far.dD) == np.sign(vt)
                     and abs(near.dD / near.sig) >= 1.5 and abs(far.dD / far.sig) >= 1.5)
        real = bool(abs(vt / st_) >= 3 and gates_ok and c_ok and VV["lofo"]["min_z"] >= 2 and nf_ok and call != "shallow-like"
                    and not artefact)
        blocking = [k for k, ok in (("|ΔD| ≥ 3σ", abs(vt / st_) >= 3), ("model gates/rules", bool(gates_ok)), ("own control", c_ok),
                                    ("leave-one-family-out ≥ 2σ", VV["lofo"]["min_z"] >= 2), ("near and far same sign ≥ 1.5σ", nf_ok),
                                    ("lapse not shallow-like", call != "shallow-like"), ("no timing artefact", not (timing_i or timing_ii or timing_iii)),
                                    ("no metadata event", not meta_hit)) if not ok]
        V["verdict" + tag] = dict(call="real deep-weighted change" if real else ("instrument or timing artefact" if artefact else "undecided"),
                                  dD=vt, sig=st_, blocking=blocking, timing_rule_i=timing_i, timing_rule_ii=timing_ii, timing_rule_iii=timing_iii, metadata_hit=meta_hit,
                                  stations_top50=top50, clock_stations=clockst, gates_ok=gates_ok, control_ok=c_ok, near_far_ok=nf_ok)
        lg(f"VERDICT 2019-11 {tag or 'H'}: {V['verdict' + tag]}")
    V["carriers"] = dict(families=[int(x) for x in carriers], station_share=st_share_c.round(3).to_dict(), from_model=cand or "H")
    R4.dump(V, f"{R4.TB4}/dive_verdict.json")
    lg(f"dive done {time.time() - t0:.0f}s")
    if not a.test:
        open(f"{R4.R4}/DIVE_DONE", "w").write(R4.pt_now() + "\n")


if __name__ == "__main__":
    main()
