#!/usr/bin/env python
"""s11_r5_pgc -- S11 round 5, task 4 (PLAN.md §9.4): SVI 2019-11 per instrument at the PGC site.
CN.PGC (surface vault) and PB.B009 / B010 / B011 (boreholes ~190-220 m) are co-located: a change in all of them is the
site (or real); a change in one is that instrument.
  (A) kernels at the profile windows 2-3, 3-4, 4-5, 5-6 s for the PGC-site pairs with a 2019-11 cell (s11_r3_kernel physics
      and seeds; the stored 2-5 / 2-3.5 / 3.5-6 s coefficients must be reproduced within 1e-6);
  (B) graph re-measured on the round-3 cc-weighted stores at 2-5 s and the profile windows (same nodes, 50-draw event bootstrap);
  (C) per-instrument apparent DeltaE(2019-11) (coefficient 1, offsets, sum-to-zero), control, common-value call;
  (D) per-instrument model DeltaD (round-5 selected form; all non-PGC-site cells + that instrument);
  (E) sub-episode parts per instrument (whole, <= 11-10, 11-11..21, 11-22..25, >= 11-26, UTC days);
  (F) metadata: FDSN channel epochs (EarthScope / NRCan), s10 epochs, data days 2019-10..12;
  (G) lapse profile vs deep / shallow / uniform / timing shapes.
    python s11_r5_pgc.py [--workers 12] [--npart 100000]
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
import s11_r3_common as R
import s11_r3_fit as LF
import s11_r3_inject as J
import s11_r4_common as R4
import s11_r5_common as R5
from scipy import stats

E19 = "2019-11"
INSTR = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
PW = {"25": (2.0, 5.0), "2035": (2.0, 3.5), "3560": (3.5, 6.0), "23": (2.0, 3.0), "34": (3.0, 4.0), "45": (4.0, 5.0), "56": (5.0, 6.0)}
PROF = ["23", "34", "45", "56"]
PARTS = [("all", None, None), ("le_1110", None, "2019-11-10"), ("1111_1121", "2019-11-11", "2019-11-21"),
         ("1122_1125", "2019-11-22", "2019-11-25"), ("ge_1126", "2019-11-26", None)]
_G = {}


# ------------------------------------------------------------------------------ (A) kernels at the profile windows
def _kpair(args):
    import s11_r3_kernel as RK
    i, fam, ns, src, rec, npart = args
    out = {}
    for var in RK.VAR3:
        seed = C.SEED + 1000 * i + (500000 if RK.VAR3[var][2] else 0)
        out[var] = RK.simulate(src, rec, rec[:2], var, npart, seed)
    return i, fam, ns, out


def kernels(series, workers, npart, lg):
    import multiprocessing as mp
    import s11_kernel_mc as K
    import s11_r3_kernel as RK
    st, c0, _, _ = RK.strike_axis(); RK.AX["v"], RK.AX["c"] = st, c0
    Fm, Sm, _ = K.geometry(None); Fe, Se, _ = K.geometry(f"{C.TB}/ext/families.csv")
    want = {(int(k.split("|")[0]), k.split("|")[1]) for k in series}
    jobs = []
    i = 0
    for fam in sorted(Fm):
        for ns in C.STATIONS:
            if (fam, ns) in want:
                jobs.append((i, fam, ns, Fm[fam], Sm[ns], npart, "main"))
            i += 1
    i = 600; got = {(j[1], j[2]) for j in jobs}
    for fam in sorted(Fe):
        for ns in C.STATIONS:
            if (fam, ns) in want and (fam, ns) not in got:
                jobs.append((i, fam, ns, Fe[fam], Se[ns], npart, "ext"))
            i += 1
    lg(f"kernels: {len(jobs)} pairs ({sum(j[6] == 'ext' for j in jobs)} ext) for {len(series)} series, {npart} particles")
    src = {(j[1], j[2]): j[6] for j in jobs}
    res = {}
    with mp.get_context("fork").Pool(workers) as p:
        for i_, fam, ns, out in p.imap_unordered(_kpair, [j[:6] for j in jobs]):
            res[(fam, ns)] = out
    Wm = K.lapse_weights(None); We = K.lapse_weights(f"{C.TB}/store_ext")
    rows = []
    for key in series:
        fam, ns, ep = key.split("|"); fam = int(fam)
        if (fam, ns) not in res:
            continue
        Eb = (We if src[(fam, ns)] == "ext" else Wm).get(key)
        if Eb is None:
            Eb = (Wm if src[(fam, ns)] == "ext" else We).get(key)
        if Eb is None:
            continue
        for v, k in res[(fam, ns)].items():
            for wn, win in PW.items():
                wb = K.window_weights(Eb, win, "3c")
                if wb.sum() <= 0:
                    continue
                av = K.coef(k, wb)
                rows.append(dict(series=key, fam=fam, netsta=ns, epoch=ep, variant=v, win=wn, a_s=float(C.group_coef(av, "s")),
                                 a_m=float(C.group_coef(av, "m")), a_d=float(C.group_coef(av, "d")), a_b=float(C.group_coef(av, "b")),
                                 unit_sum=float(np.nansum(av))))
    KC = pd.DataFrame(rows)
    old = R.kernel_table()
    m = KC[KC.win.isin(["25", "2035", "3560"])].merge(old, on=["series", "variant", "win"], suffixes=("", "_old"))
    dmax = float(np.nanmax(np.abs(m[["a_s", "a_d"]].to_numpy() - m[["a_s_old", "a_d_old"]].to_numpy()))) if len(m) else np.nan
    lg(f"kernel reproduction (2-5, 2-3.5, 3.5-6 s; {len(m)} rows): max |Δa| {dmax:.1e} -> {'OK' if dmax <= 1e-6 else 'MISMATCH'}; "
       f"unit sum max |Σa−1| {np.nanmax(np.abs(KC.unit_sum - 1)):.1e}")
    KC.to_csv(f"{R5.TB5}/pgc_kernel_profile.csv", index=False)
    return KC, dmax


# ------------------------------------------------------------------------------ (B) graph at the profile windows
def _win_worker(key):
    try:
        st = _G["store"][key]; nodes = _G["nodes"].get(key, [])
        if len(nodes) < 3:
            return key, None
        unit = _G["unit"]; rng = np.random.default_rng(zlib.crc32(key.encode()))
        ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
        X = st["St"].astype(np.float64); Ws = st["Ws"]
        base = [np.flatnonzero(ep == e) for e in nodes]

        def stacks(sel):
            Sx = np.full((len(nodes), 3, X.shape[-1]), np.nan)
            for k in range(len(nodes)):
                m = sel[k]
                for c in range(3):
                    if not np.isfinite(st["G1"][c]).all():
                        continue
                    w = Ws[m, c]
                    if w.sum() > 0:
                        Sx[k, c] = (w @ X[m, c]) / w.sum()
            return Sx
        S0 = stacks(base)
        boots = [stacks([rng.choice(ix, len(ix), replace=True) for ix in base]) for _ in range(50)]
        out = {}
        for wn, win in PW.items():
            pi, pj, dd, cc = J.graph_pairs(S0, st["tS"], win)
            v = J.graph_solve(len(nodes), pi, pj, dd, cc)[0]
            vb = []
            for Sb in boots:
                pi, pj, dd, cc = J.graph_pairs(Sb, st["tS"], win)
                vb.append(J.graph_solve(len(nodes), pi, pj, dd, cc)[0])
            out[wn] = (np.asarray(v), np.nanstd(np.array(vb), 0, ddof=1))
        return key, (list(nodes), out)
    except Exception as e:
        import traceback
        print("ERROR", key, e, traceback.format_exc()[-400:], flush=True); return key, None


def profile_graph(P, workers, lg):
    import multiprocessing as mp
    import s11_paths as Pp
    unit = Pp.unit_map()
    E = R4.r4_episode(False)
    nodes = {f"{f}|{n}|{e}": list(g.episode) for (f, n, e), g in E[(E.n >= R4.NODE_MIN) & E.b1_cc.notna()].groupby(["fam", "netsta", "epoch"])}
    D.set_grid(D.NEAR_GRID)
    vals = {}
    for ns in INSTR:
        store = {}
        for sp in (f"{R.TB3}/store_ccw", f"{R.TB3}/store_ccw_ext"):
            if os.path.exists(f"{sp}/{ns}.npz"):
                store.update(J.load_store_file(f"{sp}/{ns}.npz"))
        keys = [k for k in store if k in set(P.series)]
        _G.clear(); _G.update(store=store, unit=unit, nodes=nodes)
        t0 = time.time()
        with mp.get_context("fork").Pool(workers) as pool:
            for k, o in pool.imap_unordered(_win_worker, keys, chunksize=1):
                if o is not None:
                    vals[k] = o
        lg(f"  profile graph {ns}: {len(keys)} series {time.time() - t0:.0f}s")
    rows = []
    for i, (k, e) in enumerate(zip(P.series, P.episode)):
        if k not in vals or e not in vals[k][0]:
            continue
        j = vals[k][0].index(e)
        for wn in PW:
            v, sb = vals[k][1][wn]
            if np.isfinite(v[j]) and abs(v[j]) < C.RAIL and np.isfinite(sb[j]) and sb[j] > 0:
                rows.append(dict(series=k, episode=e, win=wn, d=-100 * v[j], sig=100 * sb[j] * P.infl.iloc[i]))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ helpers
def app_dE(c, e=E19):
    """apparent episode model on a cell set (coefficient 1, offsets, sum-to-zero) -> (dE, sig)"""
    c = c.reset_index(drop=True)
    key = c.series; c = c[key.map(key.value_counts()) >= 3].reset_index(drop=True)
    if e not in set(c.episode) or c.episode.nunique() < 3:
        return np.nan, np.nan, 0
    F = R5.Model5(c, kind="episode", unit=None, MB=False, apparent=True).solve()
    g = LF.dvec(F, e)
    if g is None:
        return np.nan, np.nan, 0
    return float(g @ F["th"]), float(np.sqrt(g @ F["C"] @ g)), int((c.episode == e).sum())


def common(x, s):
    w = 1 / s ** 2; m = float(np.sum(w * x) / w.sum()); sm = float(1 / np.sqrt(w.sum())); chi = float(np.sum(w * (x - m) ** 2))
    return m, sm, chi, float(stats.chi2.sf(chi, len(x) - 1)) if len(x) > 1 else np.nan


def fdsn(ns, t0="2019-07-18", t1="2020-03-14"):
    net, sta = ns.split(".")
    base = ("https://www.earthquakescanada.nrcan.gc.ca/fdsnws/station/1/query" if net == "CN"
            else "https://service.earthscope.org/fdsnws/station/1/query")
    url = f"{base}?net={net}&sta={sta}&level=channel&format=text&starttime=2015-01-01&endtime=2026-01-01"
    try:
        txt = urllib.request.urlopen(url, timeout=60).read().decode()
        Dd = pd.read_csv(io.StringIO(txt), sep="|"); Dd.columns = [c.strip().lstrip("#").strip() for c in Dd.columns]
        Dd["StartTime"] = Dd.StartTime.astype(str).str[:19]; Dd["EndTime"] = Dd.EndTime.fillna("").astype(str).str[:19]
        Dd["netsta"] = ns; Dd["source"] = base.split("/")[2]
        lo, hi = pd.Timestamp(t0), pd.Timestamp(t1)
        st_ = pd.to_datetime(Dd.StartTime, errors="coerce"); en = pd.to_datetime(Dd.EndTime.replace("", np.nan), errors="coerce")
        Dd["starts_in_window"] = (st_ >= lo) & (st_ <= hi); Dd["ends_in_window"] = (en >= lo) & (en <= hi)
        Dd["active_2019_11_15"] = (st_ <= "2019-11-15") & (en.isna() | (en >= "2019-11-15"))
        return Dd
    except Exception as e:
        return pd.DataFrame([dict(netsta=ns, error=str(e)[:200])])


# ------------------------------------------------------------------------------ (E) sub-episode parts per instrument
def _split_worker(key):
    from s11_measure import Ref3
    try:
        st = _G["store"][key]; unit = _G["unit"]
        days = st["days"].astype(str); ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
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
        Rf = Ref3(stack(np.flatnonzero(~m19)), tS, C.WINS["25"])
        mid19 = R4.U.mid[E19]
        nb = [e for e in pd.unique(ep[~m19]) if isinstance(e, str) and e in R4.U.index and abs((R4.U.mid[e] - mid19).days) <= 3 * 365.25]
        vnb = [Rf.measure(stack(np.flatnonzero(ep == e))[None], cs)[0][0] for e in nb]
        vnb = [v for v in vnb if np.isfinite(v) and abs(v) < C.RAIL]
        if len(vnb) < 2:
            return key, None
        base = float(np.mean(vnb)); out = {}
        for name, a, b in PARTS:
            mm = m19.copy()
            if a:
                mm &= days >= a
            if b:
                mm &= days <= b
            ix = np.flatnonzero(mm)
            if len(ix) == 0:
                continue
            v = Rf.measure(stack(ix)[None], cs)[0][0]
            ne = float(st["Nd"][ix].sum()) if st.get("Nd") is not None else np.nan
            out[name] = (-100 * (v - base), int(len(ix)), ne)
        return key, out
    except Exception as e:
        import traceback
        print("ERROR", key, e, traceback.format_exc()[-400:], flush=True); return key, None


def split_parts(H, workers, lg):
    import multiprocessing as mp
    import s11_paths as Pp
    unit = Pp.unit_map(); D.set_grid(D.NEAR_GRID); res = {}
    nr = H[~H.far & H.netsta.isin(INSTR)]
    for ns in INSTR:
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
    ce = nr[nr.episode == E19].set_index("series").sig
    for k, o in res.items():
        nall = o.get("all", (np.nan, 0, np.nan))[2]
        for name, (v, nd, ne) in o.items():
            s = ce.get(k, np.nan) * np.sqrt(nall / ne) if ne and ne > 0 and np.isfinite(nall) else np.nan
            rows.append(dict(series=k, netsta=k.split("|")[1], fam=int(k.split("|")[0]), part=name, dvv=v, sig=s, n_days=nd, n_events=ne))
    T = pd.DataFrame(rows); T.to_csv(f"{R5.TB5}/pgc_split_cells.csv", index=False)
    out = []
    for (ns, part), g in T[np.isfinite(T.sig) & (T.sig > 0) & np.isfinite(T.dvv)].groupby(["netsta", "part"]):
        w = 1 / g.sig ** 2; m = float(np.sum(w * g.dvv) / w.sum()); s = float(1 / np.sqrt(w.sum()))
        chi = float(np.sum(w * (g.dvv - m) ** 2) / max(len(g) - 1, 1)); s *= np.sqrt(max(chi, 1.0))
        out.append(dict(netsta=ns, part=part, dvv=m, sig=s, n_series=len(g), chi2=chi, n_events=float(g.n_events.sum())))
    O = pd.DataFrame(out); O.to_csv(f"{R5.TB5}/pgc_split.csv", index=False)
    for r in O.itertuples():
        lg(f"  split {r.netsta} {r.part}: {r.dvv:+.3f} ± {r.sig:.3f} % (series {r.n_series}, events {r.n_events:.0f})")
    return O


# ------------------------------------------------------------------------------ (G) lapse profile shapes
def shape_fit(x, s, shapes):
    out = {}
    for name, f in shapes.items():
        w = 1 / s ** 2; k = float(np.sum(w * f * x) / np.sum(w * f * f)); chi = float(np.sum(w * (x - k * f) ** 2))
        out[name] = dict(k=k, k_sig=float(1 / np.sqrt(np.sum(w * f * f))), chi2=chi)
    best = min(out, key=lambda n: out[n]["chi2"])
    dec = all(out[n]["chi2"] - out[best]["chi2"] >= 4 for n in out if n != best)
    return out, best, dec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--npart", type=int, default=100000)
    ap.add_argument("--test", action="store_true")
    a = ap.parse_args()
    if a.test:                                   # code test: one instrument, few particles, outputs in tables/test, no marker
        R5.TB5 = f"{R5.TB5}/test"; os.makedirs(R5.TB5, exist_ok=True); INSTR[:] = ["PB.B010"]; a.npart = 3000
    lg = R5.log("pgc" + ("_test" if a.test else "")); t0 = time.time()
    H, _ = R5.load_H(lg)
    P = H[H.netsta.isin(INSTR) & ~H.far].reset_index(drop=True)
    V = dict(written=R5.pt_now())
    lg(f"PGC-site near cells {len(P)} ({P.groupby('netsta').size().to_dict()}); 2019-11 cells {P[P.episode == E19].groupby('netsta').size().to_dict()}")
    # (C) per-instrument apparent DeltaE at 2-5 s (H cells), control, call
    rows = []
    for ns in INSTR + ["ALL"]:
        c = P if ns == "ALL" else P[P.netsta == ns]
        v, s, n = app_dE(c)
        cc = R4.control_cells(c); vc, sc, nc = app_dE(cc)
        rows.append(dict(netsta=ns, dE=v, sig=s, z=v / s if s > 0 else np.nan, n_cells_e=n, n_cells=len(c), ctrl_dE=vc, ctrl_sig=sc, ctrl_n=nc))
        lg(f"apparent ΔE(2019-11) {ns}: {v:+.3f} ± {s:.3f} % (z {v / s if s > 0 else np.nan:.2f}; {n} cells in 2019-11); scrambled control {vc:+.3f} ± {sc:.3f}")
    A = pd.DataFrame(rows); A.to_csv(f"{R5.TB5}/pgc_instrument_apparent.csv", index=False)
    a4 = A[A.netsta.isin(INSTR) & np.isfinite(A.dE) & (A.sig > 0)]
    m, sm, chi, p = common(a4.dE.to_numpy(), a4.sig.to_numpy())
    call = None; why = ""
    if p >= 0.05 and m >= 2.5 * sm and (a4.dE > 0).sum() >= 3:
        call = "SITE-WIDE"
    else:
        for r in a4.itertuples():
            o = a4[a4.netsta != r.netsta]
            mo, so, _, _ = common(o.dE.to_numpy(), o.sig.to_numpy())
            if r.dE >= 2.5 * r.sig and abs(mo) <= 1.5 * so and abs(r.dE - mo) >= 2.5 * np.hypot(r.sig, so):
                call = f"ONE INSTRUMENT ({r.netsta})"; why = f"others' common value {mo:+.3f} ± {so:.3f}"
        if call is None:
            call = "MIXED/UNDECIDED"; why = "instruments ≥ 2σ: " + ", ".join(r.netsta for r in a4.itertuples() if r.dE >= 2 * r.sig)
    vb = a4[a4.netsta != "CN.PGC"]; mb, sb, _, _ = common(vb.dE.to_numpy(), vb.sig.to_numpy())
    pv = a4[a4.netsta == "CN.PGC"]
    contrast = (float(pv.dE.iloc[0] - mb), float(np.hypot(pv.sig.iloc[0], sb))) if len(pv) else (np.nan, np.nan)
    V["instrument"] = dict(common=m, common_sig=sm, chi2=chi, p=p, call=call, why=why, vault_minus_borehole=contrast[0],
                           vault_minus_borehole_sig=contrast[1], borehole_mean=mb, borehole_mean_sig=sb)
    lg(f"common value over the 4 instruments {m:+.3f} ± {sm:.3f} (χ² {chi:.2f}, p {p:.3f}) -> {call} {why}; vault − borehole mean "
       f"{contrast[0]:+.3f} ± {contrast[1]:.3f}")
    # (D) model DeltaD per instrument: selected form on non-PGC-site cells + the instrument
    sel = json.load(open(f"{R5.R5}/selected.json"))["selected"] if os.path.exists(f"{R5.R5}/selected.json") else "E2"
    mrows = []
    for setup in ("two3_20", "ell40"):
        for ns in INSTR + ["none", "ALL"]:
            sub = H[H.site != "PGC"] if ns == "none" else (H if ns == "ALL" else pd.concat([H[H.site != "PGC"], H[H.netsta == ns]]))
            cs, good = R5.setup_cells(sub.reset_index(drop=True), setup); cs = cs[good].reset_index(drop=True)
            cs, _ = R4.identifiable(cs)
            F = R5.factory(R5.FORMS[sel])(cs).solve(); g = LF.dvec(F, E19)
            v, s = (float(g @ F["th"]), float(np.sqrt(g @ F["C"] @ g))) if g is not None else (np.nan, np.nan)
            mrows.append(dict(form=sel, setup=setup, with_instrument=ns, dD=v, sig=s, z=v / s if s > 0 else np.nan, n_cells=len(cs)))
            lg(f"model {sel}|{setup} non-PGC-site + {ns}: ΔD(2019-11) {v:+.3f} ± {s:.3f}")
    pd.DataFrame(mrows).to_csv(f"{R5.TB5}/pgc_instrument_model.csv", index=False)
    # (F) metadata
    meta = pd.concat([fdsn(ns) for ns in INSTR], ignore_index=True)
    meta.to_csv(f"{R5.TB5}/pgc_metadata_channels.csv", index=False)
    ev = []
    for ns in INSTR:
        for e, a_, b_ in S.EPOCHS.get(ns, []):
            for dd in (a_, b_):
                if abs((pd.Timestamp(dd) - pd.Timestamp("2019-11-15")).days) <= 120:
                    ev.append(dict(netsta=ns, source="s10 epoch", date=dd, detail=e))
        mm = meta[(meta.netsta == ns)]
        if "error" in mm and mm.error.notna().any():
            ev.append(dict(netsta=ns, source="FDSN", date="", detail="query failed: " + str(mm.error.iloc[0])))
            continue
        for r in mm[(mm.starts_in_window) | (mm.ends_in_window)].itertuples():
            ev.append(dict(netsta=ns, source=r.source, date=r.StartTime if r.starts_in_window else r.EndTime,
                           detail=f"{r.Channel} {'start' if r.starts_in_window else 'end'} {r.SensorDescription} {r.SampleRate}"))
    EV = pd.DataFrame(ev, columns=["netsta", "source", "date", "detail"]); EV.to_csv(f"{R5.TB5}/pgc_metadata_events.csv", index=False)
    near30 = EV[EV.date.astype(str).str.len() >= 10]
    near30 = near30[[abs((pd.Timestamp(str(d)[:10]) - pd.Timestamp("2019-11-15")).days) <= 30 for d in near30.date]] if len(near30) else near30
    V["metadata"] = dict(events_120d=EV.to_dict("records"), events_30d=near30.to_dict("records"))
    lg(f"metadata events within ±120 d: {len(EV)}; within ±30 d: {len(near30)}: " + "; ".join(f"{r.netsta} {r.date} {r.detail}" for r in EV.itertuples()))
    if call.startswith("ONE INSTRUMENT"):
        inst = call.split("(")[1].rstrip(")")
        V["instrument"]["artefact"] = bool((near30.netsta == inst).any()) if len(near30) else False
    # data days / events per instrument, 2019-10..12 (store days)
    import s11_paths as Pp
    unit = Pp.unit_map(); dd = []
    for ns in INSTR:
        store = {}
        for sp in (f"{R.TB3}/store_ccw", f"{R.TB3}/store_ccw_ext"):
            if os.path.exists(f"{sp}/{ns}.npz"):
                store.update(J.load_store_file(f"{sp}/{ns}.npz"))
        days = sorted({str(d) for st in store.values() for d in st["days"] if "2019-10-01" <= str(d) <= "2019-12-31"})
        e19 = [d for d in days if unit.get(d) == E19]
        dd.append(dict(netsta=ns, n_days_2019_10_12=len(days), n_days_2019_11_episode=len(e19), first=days[0] if days else "", last=days[-1] if days else ""))
    pd.DataFrame(dd).to_csv(f"{R5.TB5}/pgc_data_days.csv", index=False)
    lg("data days 2019-10..12 (stores): " + ", ".join(f"{r['netsta']} {r['n_days_2019_10_12']} ({r['n_days_2019_11_episode']} in the episode)" for r in dd))
    # (E) parts
    split_parts(H, a.workers, lg)
    # (A) kernels for the 2019-11 PGC-site series
    s19 = sorted(set(P[P.episode == E19].series))
    KC, dmax = kernels(s19, a.workers, a.npart, lg)
    V["kernel_repro_max"] = dmax
    # (B) graph at the profile windows
    PG = profile_graph(P, a.workers, lg); PG.to_csv(f"{R5.TB5}/pgc_profile_cells.csv", index=False)
    chk = PG[PG.win == "25"].merge(P[["series", "episode", "d"]], on=["series", "episode"], suffixes=("", "_H"))
    lg(f"profile graph 2-5 s vs H cells: max |Δd| {np.nanmax(np.abs(chk.d - chk.d_H)):.1e} % ({len(chk)} cells)")
    prof = []
    for ns in INSTR + ["ALL"]:
        for wn in PW:
            c = PG[PG.win == wn].merge(P[["series", "episode", "fam", "netsta", "epoch", "site"]], on=["series", "episode"])
            if ns != "ALL":
                c = c[c.netsta == ns]
            v, s, n = app_dE(c)
            row = dict(netsta=ns, win=wn, t0=PW[wn][0], t1=PW[wn][1], dE=v, sig=s, n_cells_e=n)
            ce = c[c.episode == E19].series.unique()
            for var in ("two3_20", "ell40", "ell10"):
                k = KC[(KC.variant == var) & (KC.win == wn) & KC.series.isin(ce)]
                row[f"a_d_{var}"] = float(k.a_d.mean()) if len(k) else np.nan; row[f"a_s_{var}"] = float(k.a_s.mean()) if len(k) else np.nan
            tt = np.linspace(PW[wn][0], PW[wn][1], 200); row["timing_f"] = float(tt.mean() / np.mean(tt ** 2))
            prof.append(row)
    PR = pd.DataFrame(prof); PR.to_csv(f"{R5.TB5}/pgc_profile.csv", index=False)
    shapes_out = {}
    for ns in INSTR + ["ALL"]:
        q = PR[(PR.netsta == ns) & PR.win.isin(PROF) & np.isfinite(PR.dE) & (PR.sig > 0)]
        if len(q) < 3:
            continue
        for var in ("two3_20", "ell40"):
            sh = {"deep": q[f"a_d_{var}"].to_numpy(), "shallow": q[f"a_s_{var}"].to_numpy(), "uniform": np.ones(len(q)), "timing": q.timing_f.to_numpy()}
            if not all(np.isfinite(v).all() for v in sh.values()):
                continue
            o, best, dec = shape_fit(q.dE.to_numpy(), q.sig.to_numpy(), sh)
            shapes_out[f"{ns}|{var}"] = dict(fits=o, best=best, decisive=dec, values=q.dE.round(4).tolist(), sigmas=q.sig.round(4).tolist())
            lg(f"profile {ns} [{var}]: ΔE 2-3/3-4/4-5/5-6 s = {', '.join(f'{x:+.3f}±{y:.3f}' for x, y in zip(q.dE, q.sig))}; χ² "
               + ", ".join(f"{n} {o[n]['chi2']:.2f}" for n in o) + f" -> best {best} ({'decisive' if dec else 'not decisive'})")
    V["profile"] = shapes_out
    R5.dump(V, f"{R5.TB5}/pgc_verdict.json")
    lg(f"pgc done {time.time() - t0:.0f}s")
    if not a.test:
        open(f"{R5.R5}/PGC_DONE", "w").write(R5.pt_now() + "\n")


if __name__ == "__main__":
    main()
