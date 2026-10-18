#!/usr/bin/env python
"""s11_v2 -- S11 Revision 1 (PLAN.md §3, MERLIN_REVIEW_1.md): cell rule v2, smooth parameterisation, Q1 (deep annual
amplitude) and Q2 (decade contrast of deep), injections A/B incl. coherent kernel error, Q1 gates G1-G5, SNB diagnostic,
figure. Raw path, no new measurement. dv/v in %; eps = -dv/v.
    python s11_v2.py [--workers 12] [--nperm 500]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I
import s11_paths as P

CC2 = 0.35
KNOTS = [2007, 2011, 2015, 2019, 2023]
PRIOR = 1.0
DEC1, DEC2 = C.DEC1, ("2015-01-01", "2025-12-31")          # C.DEC1: SVI 2004-2014, WA 2007-2014 (PLAN §5)
STEP_DATE, PEAK = "2015-01-01", 233.0
U = C.episodes()


# ------------------------------------------------------------------------------ basis and model
def tdec(dates):
    d = pd.to_datetime(pd.Series(dates)).reset_index(drop=True)
    return (d.dt.year + (d.dt.dayofyear - 1) / 365.25).to_numpy(), d.dt.dayofyear.to_numpy()


def basis(t, doy, terms):
    ph = 2 * np.pi * (doy - 1) / 365.25
    out = {}
    for tm in terms:
        if tm == "lin":
            out[tm] = t - 2015.0
        elif tm.startswith("h"):
            out[tm] = np.maximum(0.0, t - float(tm[1:]))
        elif tm == "c1":
            out[tm] = np.cos(ph)
        elif tm == "s1":
            out[tm] = np.sin(ph)
        elif tm == "c2":
            out[tm] = np.cos(2 * ph)
        elif tm == "s2":
            out[tm] = np.sin(2 * ph)
    return out


SLOW = ["lin"] + [f"h{k}" for k in KNOTS]
MODELS = {"smooth": (["c1", "s1", "c2", "s2"] + SLOW, SLOW + ["c1", "s1"]),
          "slow+annual": (["c1", "s1"] + SLOW, SLOW + ["c1", "s1"])}


def fit(cells, model="smooth", prior=PRIOR, d=None, mid=None):
    Sterms, Dterms = MODELS[model]
    m = cells.episode.map(U.mid) if mid is None else mid
    t, doy = tdec(m)
    bS, bD = basis(t, doy, Sterms), basis(t, doy, Dterms)
    names = [("D", "", tm) for tm in Dterms]
    sites = sorted(cells.site.unique())
    names += [("S", s, tm) for s in sites for tm in Sterms]
    idx = {nm: i for i, nm in enumerate(names)}
    G = np.zeros((len(cells), len(names)))
    ad, as_ = cells.a_d.to_numpy(), cells.a_s.to_numpy(); site = cells.site.to_numpy()
    for tm in Dterms:
        G[:, idx[("D", "", tm)]] = ad * bD[tm]
    for s in sites:
        ms = site == s
        for tm in Sterms:
            G[ms, idx[("S", s, tm)]] = as_[ms] * bS[tm][ms]
    dd = cells.d.to_numpy() if d is None else d
    w = 1 / cells.sig.to_numpy() ** 2
    Gt, dt, nser = I.demean(cells, G, dd, w)
    N = Gt.T @ (Gt * w[:, None]); b = Gt.T @ (w * dt)
    Ci = np.linalg.inv(N + np.eye(len(names)) / prior ** 2)
    th = Ci @ b
    r = dt - Gt @ th
    peff = float(np.trace(Ci @ N)); chi2 = float(np.sum(w * r ** 2) / max(len(dd) - nser - peff, 1))
    Cs = Ci * max(chi2, 1.0)
    return dict(th=th, C=Cs, idx=idx, names=names, chi2=chi2, n=len(dd), nser=nser, N=N, b=b)


def amp(F, kind="D", site="", c="c1", s="s1"):
    i, j = F["idx"][(kind, site, c)], F["idx"][(kind, site, s)]
    a, b = F["th"][i], F["th"][j]; A = float(np.hypot(a, b))
    J = np.zeros(len(F["th"])); J[i], J[j] = a / max(A, 1e-12), b / max(A, 1e-12)
    pk = (np.degrees(np.arctan2(b, a)) % 360) / 360 * 365.25 + 1
    return A, float(np.sqrt(J @ F["C"] @ J)), float(pk), (a, b)


def q2vec(F):
    g = np.zeros(len(F["th"]))
    for (a0, a1), sg in ((DEC2, 1.0), (DEC1, -1.0)):
        grid = pd.date_range(a0, a1, freq="MS")
        t, doy = tdec(grid)
        bb = basis(t, doy, SLOW)
        for tm in SLOW:
            g[F["idx"][("D", "", tm)]] += sg * bb[tm].mean()
    return g


def Q(F):
    A, sA, pk, cs = amp(F)
    g = q2vec(F)
    return dict(Q1=A, Q1_sig=sA, Q1_peak=pk, Q1_c=cs[0], Q1_s=cs[1], Q2=float(g @ F["th"]), Q2_sig=float(np.sqrt(g @ F["C"] @ g)),
                chi2=F["chi2"], n_cells=F["n"], n_series=F["nser"])


def cells_v2(win=C.HEAD, variant="ell40", ccmin=CC2, control=False):
    if control:
        c = I.load_cells("raw", win, control=True, ccmin=ccmin)
        return I.prepare("raw", win, variant=variant, cells=c)
    return I.prepare("raw", win, variant=variant, ccmin=ccmin)


def site_amp_curve(F, site):
    x = np.arange(1, 366); ph = 2 * np.pi * (x - 1) / 365.25
    v = np.zeros(len(x))
    for tm, f in (("c1", np.cos(ph)), ("s1", np.sin(ph)), ("c2", np.cos(2 * ph)), ("s2", np.sin(2 * ph))):
        k = ("S", site, tm) if site else ("D", "", tm)
        if k in F["idx"]:
            v += F["th"][F["idx"][k]] * f
    return x, v


# ------------------------------------------------------------------------------ kernel-free in-phase amplitudes
def demeaned(c):
    w = 1 / c.sig ** 2; key = c.fam.astype(str) + "|" + c.netsta + "|" + c.epoch
    return c.d - (c.d * w).groupby(key).transform("sum") / w.groupby(key).transform("sum")


def inphase(c, peak):
    c = c.copy(); c["dm"] = demeaned(c)
    doy = c.episode.map(U.mid).dt.dayofyear.to_numpy()
    x = np.cos(2 * np.pi * (doy - peak) / 365.25); w = 1 / c.sig.to_numpy() ** 2
    X = np.column_stack([np.ones(len(x)), x]); N = X.T @ (X * w[:, None]); m = np.linalg.solve(N, X.T @ (w * c.dm.to_numpy()))
    r = c.dm.to_numpy() - X @ m; chi2 = float(np.sum(w * r ** 2) / max(len(x) - 2, 1))
    return float(m[1]), float(np.sqrt(np.linalg.inv(N)[1, 1] * max(chi2, 1.0)))


# ------------------------------------------------------------------------------ injections
_G = {}


def _inj_series(key):
    try:
        st = _G["store"][key]; unit = _G["unit"]
        out = {"nat": P.raw_cells(st, unit, wins=[C.HEAD])[["episode", f"eps3c_{C.HEAD}"]]}
        t = pd.to_datetime(pd.Series(st["days"]))
        post = (t >= STEP_DATE).to_numpy(); ann = np.cos(2 * np.pi * (t.dt.dayofyear.to_numpy() - PEAK) / 365.25)
        for var in ("ell40", "two5_40", "ell80"):
            K = _G["K"][var]
            if key not in K.index:
                continue
            a_s, a_d = float(K.loc[key, "a_s"]), float(K.loc[key, "a_d"])
            for run, dvv in (("A", 0.003 * a_d * post), ("B", 0.001 * a_s * ann)):
                out[f"{run}_{var}"] = P.raw_cells(st, unit, amp=-dvv, wins=[C.HEAD])[["episode", f"eps3c_{C.HEAD}"]]
        return key, out
    except Exception as e:
        print("ERROR", key, e, flush=True); return key, {}


def run_injections(cells, workers, log):
    import multiprocessing as mp
    unit = P.unit_map()
    Kall = pd.read_parquet(f"{C.TB}/kernel_coef.parquet")
    Kall = Kall[(Kall.comp == "3c") & (Kall.win == C.HEAD)]
    K = {v: Kall[Kall.variant == v].set_index("series") for v in ("ell40", "two5_40", "ell80")}
    need = set(cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch)
    vals = {}
    for ns in C.STATIONS:
        if not os.path.exists(f"{C.TB}/store/{ns}.npz"):        # region runs: stations without any measured series
            continue
        _G.clear(); _G.update(store=P.load_store(ns), unit=unit, K=K)
        keys = [k for k in _G["store"] if k in need]
        with mp.get_context("fork").Pool(workers) as pool:
            for k, o in pool.imap_unordered(_inj_series, keys):
                vals[k] = o
        log(f"  injections {ns}: {len(keys)} series")
    ckey = (cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch).to_numpy()

    def dvec(run):
        d = np.full(len(cells), np.nan)
        for i, (k, e) in enumerate(zip(ckey, cells.episode)):
            q = vals.get(k, {}).get(run)
            if q is None or not len(q):
                continue
            m = q.episode == e
            if m.any():
                d[i] = -100 * float(q[f"eps3c_{C.HEAD}"][m].iloc[0])
        return d
    dn = dvec("nat")
    log(f"  natural-from-store vs measured: max |diff| {np.nanmax(np.abs(dn - cells.d.to_numpy())):.1e} %")
    rows = []
    for var in ("ell40", "two5_40", "ell80"):
        for run in ("A", "B"):
            di = dvec(f"{run}_{var}"); ok = np.isfinite(di) & np.isfinite(dn)
            cc = cells[ok].reset_index(drop=True)
            Fn, Fi = fit(cc, d=dn[ok]), fit(cc, d=di[ok])
            dth = Fi["th"] - Fn["th"]
            Fd = dict(Fn); Fd["th"] = dth
            A1, _, pk, _ = amp(Fd); q2 = float(q2vec(Fd) @ dth)
            r = dict(run=run, generated_with=var, inverted_with="ell40", n_cells=int(ok.sum()), dQ1=A1, dQ1_peak=pk, dQ2=q2)
            if run == "A":
                r["Q2_recovery"] = q2 / 0.3; r["pass"] = bool(0.7 <= r["Q2_recovery"] <= 1.3 and A1 < 0.03)
            else:
                r["Q1_leak"] = A1 / 0.1
                for s in C.SITES:
                    if ("S", s, "c1") in Fd["idx"]:
                        r[f"S_{s}_annual_rec"] = amp(Fd, "S", s)[0] / 0.1
                r["pass"] = bool(r["Q1_leak"] < 0.25)
            rows.append(r)
            log("  inj " + ", ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--nperm", type=int, default=500)
    a = ap.parse_args()
    log = C.Log(f"{C.LOGD}/v2.log"); t0 = time.time()
    cells = cells_v2()
    log(f"v2 cells: {len(cells)} in {cells.fam.astype(str).add(cells.netsta + cells.epoch).nunique()} series, "
        f"{cells.episode.nunique()} episodes; by site {cells.site.value_counts().to_dict()}")
    cells.to_csv(f"{C.TB}/v2_cells.csv", index=False)
    # ---------------- models
    rows = []
    Fh = fit(cells, "smooth"); qh = Q(Fh)
    for model, cc_, label in (("smooth", cells, "head"), ("slow+annual", cells, "slow+annual")):
        F = fit(cc_, model); q = Q(F)
        for f_ in (2.0, 0.5):
            qq = Q(fit(cc_, model, prior=PRIOR * f_))
            q[f"prior_x{f_:g}_dQ1_over_sig"] = (qq["Q1"] - q["Q1"]) / q["Q1_sig"]
            q[f"prior_x{f_:g}_dQ2_over_sig"] = (qq["Q2"] - q["Q2"]) / q["Q2_sig"]
        rows.append(dict(model=model, label=label, variant="ell40", **q))
    for var in ("ell20", "ell80", "two5_40", "ell40_s2"):
        rows.append(dict(model="smooth", label=f"kernel_{var}", variant=var, **Q(fit(cells_v2(variant=var), "smooth"))))
    c08 = cells_v2(ccmin=0.8)
    rows.append(dict(model="smooth", label="cc0.8", variant="ell40", **Q(fit(c08, "smooth"))))
    noSNB = cells[cells.site != "SNB"].reset_index(drop=True)
    # joint lapse
    cl = []
    for w in ("2035", "3560"):
        c = cells_v2(win=w); c["epoch"] = c.epoch + f"_w{w}"; cl.append(c)
    CL = pd.concat(cl, ignore_index=True)
    rows.append(dict(model="smooth", label="joint-lapse", variant="ell40", **Q(fit(CL, "smooth"))))
    # episode model -> per-episode bound
    par, summ = I.invert(cells, mode="episode")
    D = par[par.kind == "D"]
    rows.append(dict(model="episode", label="episode", variant="ell40", n_cells=summ["n_cells"], n_series=summ["n_series"],
                     chi2=summ["chi2red"], episode_D_sig_median=float(D.sig.median()),
                     episode_D_2sig_bound=float(2 * D.sig.median()), prior_x2_change=summ["prior_x2_change"],
                     prior_x05_change=summ["prior_x0.5_change"]))
    par.to_csv(f"{C.TB}/v2_inv_episode.csv", index=False)
    M = pd.DataFrame(rows); M.to_csv(f"{C.TB}/v2_models.csv", index=False)
    for r in rows:
        log("model " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    # ---------------- Q1 gates
    G = []
    pk0 = qh["Q1_peak"]
    dpk = lambda p: abs((p - pk0 + 182.6) % 365.25 - 182.6)
    ok1 = True
    for s in C.SITES:
        q = Q(fit(cells[cells.site != s].reset_index(drop=True)))
        within = abs(q["Q1"] - qh["Q1"]) < 2 * q["Q1_sig"]
        pk_ok = dpk(q["Q1_peak"]) <= 60 if q["Q1"] / q["Q1_sig"] >= 1 else True
        ok1 &= within and pk_ok
        G.append(dict(gate="G1", stat=f"without_{s}", Q1=q["Q1"], Q1_sig=q["Q1_sig"], peak=q["Q1_peak"], Q2=q["Q2"], Q2_sig=q["Q2_sig"],
                      passed=bool(within and pk_ok)))
    G.append(dict(gate="G1", stat="all_sites", passed=bool(ok1)))
    qo = Q(fit(cells[cells.fam % 2 == 1].reset_index(drop=True))); qe = Q(fit(cells[cells.fam % 2 == 0].reset_index(drop=True)))
    w2 = abs(qo["Q1"] - qe["Q1"]) < 2 * np.hypot(qo["Q1_sig"], qe["Q1_sig"])
    pko = abs((qo["Q1_peak"] - qe["Q1_peak"] + 182.6) % 365.25 - 182.6) <= 60 if min(qo["Q1"] / qo["Q1_sig"], qe["Q1"] / qe["Q1_sig"]) >= 1 else True
    G.append(dict(gate="G2", stat="odd", Q1=qo["Q1"], Q1_sig=qo["Q1_sig"], peak=qo["Q1_peak"], Q2=qo["Q2"], Q2_sig=qo["Q2_sig"]))
    G.append(dict(gate="G2", stat="even", Q1=qe["Q1"], Q1_sig=qe["Q1_sig"], peak=qe["Q1_peak"], Q2=qe["Q2"], Q2_sig=qe["Q2_sig"],
                  passed=bool(w2 and pko)))
    cc = cells_v2(control=True)
    if len(cc) >= 20:
        qc = Q(fit(cc))
        G.append(dict(gate="G3", stat="scrambled_control", Q1=qc["Q1"], Q1_sig=qc["Q1_sig"], peak=qc["Q1_peak"], Q2=qc["Q2"],
                      Q2_sig=qc["Q2_sig"], n=len(cc), passed=bool(qc["Q1"] < 2 * qc["Q1_sig"])))
    else:
        G.append(dict(gate="G3", stat="scrambled_control", n=len(cc), passed=None))
    rng = np.random.default_rng(C.SEED); eps_u = cells.episode.unique(); mids = U.mid.reindex(eps_u).to_numpy()
    perm = []
    for _ in range(a.nperm):
        mp_ = dict(zip(eps_u, rng.permutation(mids)))
        perm.append(Q(fit(cells, mid=cells.episode.map(mp_)))["Q1"])
    perm = np.array(perm)
    G.append(dict(gate="G3col", stat="episode_date_permutation_p", Q1=qh["Q1"], value=float((1 + np.sum(perm >= qh["Q1"])) / (1 + len(perm))),
                  perm_q95=float(np.quantile(perm, 0.95)), passed=None))
    # G4 station amplitudes (kernel-free in-phase at the pooled kernel-free peak)
    cc_all = cells.copy(); cc_all["dm"] = demeaned(cc_all)
    fpool = C.annual_fit(cc_all.episode.map(U.mid), cc_all.dm.to_numpy(), cc_all.sig.to_numpy())
    pkp = fpool["peak_doy"]
    st = []
    for ns in C.STATIONS:
        q = cells[cells.netsta == ns]
        if len(q) < 10:
            continue
        A_, s_ = inphase(q, pkp); w = 1 / q.sig ** 2
        st.append(dict(netsta=ns, site=C.SITE[ns], n=len(q), A=A_, A_sig=s_, a_d=float(np.average(q.a_d, weights=w)),
                       a_s=float(np.average(q.a_s, weights=w))))
    ST = pd.DataFrame(st)
    for col in ("a_d", "a_s"):
        x, y, s_ = ST[col].to_numpy(), ST.A.to_numpy(), ST.A_sig.to_numpy()
        k = np.sum(x * y / s_ ** 2) / np.sum(x ** 2 / s_ ** 2); chi = float(np.sum((y - k * x) ** 2 / s_ ** 2))
        ST[f"pred_{col}"] = k * x
        G.append(dict(gate="G4", stat=f"chi2_prop_{col}", value=chi, scale=float(k), dof=len(x) - 1,
                      p=float(stats.chi2.sf(chi, len(x) - 1))))
    chd = [g for g in G if g.get("stat") == "chi2_prop_a_d"][0]; chs = [g for g in G if g.get("stat") == "chi2_prop_a_s"][0]
    G.append(dict(gate="G4", stat="deep_pattern_supported", value=chs["value"] - chd["value"],
                  passed=bool(chd["p"] > 0.05 and chs["value"] - chd["value"] > 4)))
    ST.to_csv(f"{C.TB}/v2_station_amp.csv", index=False)
    # G5 lapse ratio
    cA, cB = cells_v2(win="2035"), cells_v2(win="3560")
    A1_, s1_ = inphase(cA, pkp); A2_, s2_ = inphase(cB, pkp)
    ratio = A2_ / A1_; sr = abs(ratio) * np.hypot(s1_ / A1_, s2_ / A2_)
    Kc = pd.read_parquet(f"{C.TB}/kernel_coef.parquet"); Kc = Kc[(Kc.variant == "ell40") & (Kc.comp == "3c")]
    kk = Kc.pivot_table(index="series", columns="win", values=["a_d", "a_s"])
    sers = set(cA.fam.astype(str) + "|" + cA.netsta + "|" + cA.epoch) | set(cB.fam.astype(str) + "|" + cB.netsta + "|" + cB.epoch)
    kk = kk[kk.index.isin(sers)]
    pd_ = float(kk[("a_d", "3560")].mean() / kk[("a_d", "2035")].mean()); ps_ = float(kk[("a_s", "3560")].mean() / kk[("a_s", "2035")].mean())
    zd, zs = (ratio - pd_) / sr, (ratio - ps_) / sr
    G.append(dict(gate="G5", stat="lapse_ratio_3560_over_2035", value=float(ratio), sig=float(sr), pred_deep=pd_, pred_shallow=ps_,
                  z_deep=float(zd), z_shallow=float(zs), A2035=A1_, A3560=A2_,
                  passed=bool(abs(zs) - abs(zd) >= 1)))
    # SNB diagnostic
    sn = []
    q = cells[cells.netsta == "CN.SNB"]
    for ep, g in q.groupby("epoch"):
        if len(g) >= 8:
            A_, s_ = inphase(g, pkp); sn.append(dict(epoch=ep, n=len(g), A=A_, A_sig=s_, a_d=g.a_d.mean()))
    fs = C.annual_fit(q.episode.map(U.mid), demeaned(q).to_numpy(), q.sig.to_numpy())
    sn.append(dict(epoch="all_free_phase", n=len(q), A=fs["amp"], A_sig=fs["amp_sig"], peak=fs["peak_doy"]))
    qn = Q(fit(noSNB)); sn.append(dict(epoch="Q1_without_SNB", A=qn["Q1"], A_sig=qn["Q1_sig"], peak=qn["Q1_peak"]))
    pd.DataFrame(sn).to_csv(f"{C.TB}/v2_snb.csv", index=False)
    for r in sn:
        log("SNB " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    GT = pd.DataFrame(G); GT.to_csv(f"{C.TB}/v2_gates.csv", index=False)
    for r in G:
        log("gate " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    # ---------------- injections
    INJ = run_injections(cells, a.workers, log); INJ.to_csv(f"{C.TB}/v2_injection.csv", index=False)
    # ---------------- figure
    try:
        fig_annual(Fh, qh, ST, pkp)
    except Exception as e:                           # a figure failure must never block V2_DONE
        log(f"FIGURE ERROR fig_annual (results unaffected): {e}")
    json.dump(dict(Q=qh, written=C.pt_now()), open(f"{C.TB}/v2_headline.json", "w"), indent=1)
    log(f"v2 done {time.time() - t0:.0f}s")
    open(f"{C.OUT}/V2_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


def fig_annual(F, q, ST, pkp):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6",
                         "grid.linewidth": 0.5, "font.size": 9})
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
    x, v = site_amp_curve(F, "")
    ax[0].plot(x, v, "-", color="k", lw=2.5, label=f"deep D (25–45 km): {q['Q1']:.3f} % amplitude, peak day {q['Q1_peak']:.0f}")
    if len(C.SITES) <= len(CAT):                     # SVI: unchanged
        col = {s: CAT[i] for i, s in enumerate(C.SITES)}
    else:                                           # region runs: top-4 sites by cells coloured, the rest a grey bundle
        rank = ST.groupby("site").n.sum().sort_values(ascending=False).index.tolist()
        col = {s: (CAT[rank.index(s)] if s in rank[:len(CAT)] else None) for s in C.SITES}
    grey = 0
    for i, s in enumerate(C.SITES):
        if ("S", s, "c1") in F["idx"]:
            x, v = site_amp_curve(F, s); A = amp(F, "S", s)[0]
            if col[s] is None:
                ax[0].plot(x, v, "-", color="#b8b8b4", lw=0.7, zorder=0,
                           label="other sites' shallow S (grey)" if grey == 0 else None); grey += 1
            else:
                ax[0].plot(x, v, "-", color=col[s], lw=1.5, label=f"shallow S {s} (0–5 km): 1-yr amplitude {A:.3f} %")
    ax[0].axhline(0, color="grey", lw=0.6, ls=":")
    ax[0].set(xlabel="day of year", ylabel="dv/v (%)", title="smooth model, v2 cells: annual terms (σ in tables, not drawn)")
    ax[0].legend(fontsize=7)
    wts = 1 / ST.A_sig ** 2; wts = wts / wts.max()
    sc = ax[1].scatter(ST.a_d, ST.A, c=wts, cmap="viridis", vmin=0, vmax=1, s=60, edgecolor="k", lw=0.4, zorder=3)
    for r in ST.itertuples():
        ax[1].annotate(r.netsta.split(".")[1], (r.a_d, r.A), xytext=(4, 3), textcoords="offset points", fontsize=8)
    o = np.argsort(ST.a_d.to_numpy())
    ax[1].plot(ST.a_d.to_numpy()[o], ST.pred_a_d.to_numpy()[o], "-", color="k", lw=1.2, label="best fit ∝ a_d (deep)")
    ax[1].plot(ST.a_d.to_numpy(), ST.pred_a_s.to_numpy(), "x", color="#eb6834", ms=8, mew=1.5, label="best fit ∝ a_s (shallow), per station")
    ax[1].set(xlabel="mean deep coefficient a_d of the station's cells", ylabel=f"in-phase annual dv/v amplitude (%), peak day {pkp:.0f}",
              title="per-station annual amplitude (kernel-free)")
    ax[1].legend(fontsize=7); fig.colorbar(sc, ax=ax[1], label="relative weight 1/σ²")
    fig.tight_layout(); fig.savefig(f"{C.FIG}/fig_v2_annual.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
