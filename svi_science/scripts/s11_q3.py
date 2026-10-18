#!/usr/bin/env python
"""s11_q3 -- S11 Revision 2 (PLAN §4): Q3 = DeltaD(e) under the temporal model (S per site smooth = 2 annual harmonics +
linear + 4-yr hinges; D free per episode; series offsets free), 15-major composite, configurations (a) 76 families
n >= 50, (b) 76 families n >= 20, (c) HEADLINE 115 families n >= 20 (cc >= 0.35, 1/sigma^2, control inflation with the
B010 fix); gates on (c): prior x2/x0.5, scrambled control, one-episode deep pulse, common shallow pulse, leave-one-site-
out; figure. Raw path. dv/v in %; eps = -dv/v.
    python s11_q3.py [--workers 12]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I
import s11_paths as P
import s11_v2 as V

MAJ15_SVI = ["2003-02", "2004-07", "2005-09", "2008-05", "2010-08b1", "2011-08", "2012-08", "2013-09", "2015-12", "2018-06",
         "2019-08", "2019-11", "2020-11", "2023-04", "2024-08"]
MAJ15 = C.MAJ_LIST if C.REGION != "SVI" else MAJ15_SVI     # region runs: list fixed by count (PLAN §5, tables/maj_list.csv)
WIN3 = 3 * 365.25
CC2 = 0.35
PRIOR = 1.0


def combine():
    """tag _all = S10-multi families + added families"""
    E = pd.read_parquet(f"{C.TB}/episode_raw.parquet"); S_ = pd.read_csv(f"{C.TB}/series.csv")
    pe = f"{C.TB}/episode_raw_ext.parquet"
    if os.path.exists(pe):
        E = pd.concat([E, pd.read_parquet(pe)], ignore_index=True); S_ = pd.concat([S_, pd.read_csv(f"{C.TB}/series_ext.csv")], ignore_index=True)
    E.to_parquet(f"{C.TB}/episode_raw_all.parquet", index=False); S_.to_csv(f"{C.TB}/series_all.csv", index=False)
    return os.path.exists(pe)


def cells_for(tag, nmin, control=False):
    if control:
        c = I.load_cells("raw", C.HEAD, tag, control=True, nmin=nmin, ccmin=CC2)
        return I.prepare("raw", C.HEAD, tag, cells=c)
    return I.prepare("raw", C.HEAD, tag, nmin=nmin, ccmin=CC2)


# ------------------------------------------------------------------------------ temporal model
def fit(cells, prior=PRIOR, d=None):
    U = C.episodes()
    m = cells.episode.map(U.mid); t, doy = V.tdec(m)
    Sterms = V.MODELS["smooth"][0]; bS = V.basis(t, doy, Sterms)
    eps_ = sorted(cells.episode.unique(), key=lambda e: U.mid[e])
    names = [("D", "", e) for e in eps_]
    sites = sorted(cells.site.unique())
    names += [("S", s, tm) for s in sites for tm in Sterms]
    idx = {n: i for i, n in enumerate(names)}
    G = np.zeros((len(cells), len(names)))
    G[np.arange(len(cells)), [idx[("D", "", e)] for e in cells.episode]] = cells.a_d.to_numpy()
    for s in sites:
        ms = (cells.site == s).to_numpy()
        for tm in Sterms:
            G[ms, idx[("S", s, tm)]] = cells.a_s.to_numpy()[ms] * bS[tm][ms]
    dd = cells.d.to_numpy() if d is None else d
    w = 1 / cells.sig.to_numpy() ** 2
    Gt, dt, ns = I.demean(cells, G, dd, w)
    N = Gt.T @ (Gt * w[:, None]); Ci = np.linalg.inv(N + np.eye(len(names)) / prior ** 2)
    th = Ci @ (Gt.T @ (w * dt)); r = dt - Gt @ th
    chi = float(np.sum(w * r ** 2) / max(len(dd) - ns - np.trace(Ci @ N), 1))
    return dict(th=th, C=Ci * max(chi, 1.0), idx=idx, eps=eps_, chi2=chi, n=len(dd), nser=ns,
                ncell=cells.groupby("episode").size(), U=U)


def dvec(F, e):
    U = F["U"]; g = np.zeros(len(F["th"]))
    oth = [x for x in F["eps"] if x != e and abs((U.mid[x] - U.mid[e]).days) <= WIN3]
    if not oth:
        return None
    g[F["idx"][("D", "", e)]] = 1.0
    for x in oth:
        g[F["idx"][("D", "", x)]] -= 1.0 / len(oth)
    return g


def q3(F):
    rows, gs = [], []
    for e in F["eps"]:
        g = dvec(F, e)
        if g is None:
            continue
        rows.append(dict(episode=e, mid=F["U"].mid[e], cls=F["U"].cls.get(e, ""), major15=e in MAJ15,
                         n_cells=int(F["ncell"].get(e, 0)), dD=float(g @ F["th"]), sig=float(np.sqrt(g @ F["C"] @ g))))
        if e in MAJ15:
            gs.append(g)
    T = pd.DataFrame(rows)
    comp = None
    if gs:
        gc = np.mean(gs, 0); comp = dict(composite=float(gc @ F["th"]), composite_sig=float(np.sqrt(gc @ F["C"] @ gc)), n_majors=len(gs))
    return T, comp


# ------------------------------------------------------------------------------ pulse injections (waveform level)
_G = {}


def _ser(key):
    try:
        st = _G["store"][key]; unit = _G["unit"]; K = _G["K"]
        out = {"nat": P.raw_cells(st, unit, wins=[C.HEAD])[["episode", f"eps3c_{C.HEAD}"]]}
        if key not in K.index:
            return key, out
        a_s, a_d = float(K.loc[key, "a_s"]), float(K.loc[key, "a_d"])
        ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
        for e in MAJ15:
            m = ep == e
            if not m.any():
                continue
            out[f"deep_{e}"] = P.raw_cells(st, unit, amp=-(0.003 * a_d) * m, wins=[C.HEAD])[["episode", f"eps3c_{C.HEAD}"]]
            out[f"shal_{e}"] = P.raw_cells(st, unit, amp=-(0.001 * a_s) * m, wins=[C.HEAD])[["episode", f"eps3c_{C.HEAD}"]]
        mall = np.isin(ep, MAJ15)
        if mall.any():
            out["shal_all15"] = P.raw_cells(st, unit, amp=-(0.001 * a_s) * mall, wins=[C.HEAD])[["episode", f"eps3c_{C.HEAD}"]]
        return key, out
    except Exception as ex:
        print("ERROR", key, ex, flush=True); return key, {}


def injections(cells, workers, log):
    import multiprocessing as mp
    unit = P.unit_map()
    K = I.kernels("ell40"); K = K[K.win == C.HEAD].set_index("series")
    need = set(cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch)
    vals = {}
    for ns in C.STATIONS:
        store = {}
        for tg in ("", "_ext"):
            if os.path.exists(f"{C.TB}/store{tg}/{ns}.npz"):
                store.update(P.load_store(ns, tg))
        _G.clear(); _G.update(store=store, unit=unit, K=K)
        keys = [k for k in store if k in need]
        with mp.get_context("fork").Pool(workers) as pool:
            for k, o in pool.imap_unordered(_ser, keys):
                vals[k] = o
        log(f"  pulses {ns}: {len(keys)} series")
    ckey = (cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch).to_numpy()

    def dv(run):
        d = np.full(len(cells), np.nan)
        for i, (k, e) in enumerate(zip(ckey, cells.episode)):
            o = vals.get(k, {}); q = o.get(run, o.get("nat"))
            if q is None or not len(q):
                continue
            mm = q.episode == e
            if mm.any():
                d[i] = -100 * float(q[f"eps3c_{C.HEAD}"][mm].iloc[0])
        return d
    dn = dv("nat")
    log(f"  natural-from-store vs measured: max |diff| {np.nanmax(np.abs(dn - cells.d.to_numpy())):.1e} %, missing {int(np.isnan(dn).sum())}")
    rows = []
    for e in MAJ15 + ["all15"]:
        for kind, a_ in (("deep", 0.3), ("shal", 0.1)):
            if e == "all15" and kind == "deep":
                continue
            di = dv(f"{kind}_{e}"); ok = np.isfinite(di) & np.isfinite(dn)
            cc = cells[ok].reset_index(drop=True)
            Fn, Fi = fit(cc, d=dn[ok]), fit(cc, d=di[ok])
            if e == "all15":
                gs = [dvec(Fn, x) for x in MAJ15 if ("D", "", x) in Fn["idx"]]
                gc = np.mean([g for g in gs if g is not None], 0); val = float(gc @ (Fi["th"] - Fn["th"]))
            else:
                if ("D", "", e) not in Fn["idx"]:
                    continue
                g = dvec(Fn, e); val = float(g @ (Fi["th"] - Fn["th"]))
            r = dict(episode=e, kind=kind, injected=a_, response=val, ratio=val / a_, n_cells=int(ok.sum()))
            rows.append(r)
            log(f"  pulse {kind} {e}: Δ(ΔD) {val:+.4f} % → ratio {val / a_:+.3f}")
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    log = C.Log(f"{C.LOGD}/q3.log"); t0 = time.time()
    have_ext = combine()
    configs = [("a_76fam_n50", "", 50), ("b_76fam_n20", "", 20)] + ([("c_115fam_n20", "_all", 20)] if have_ext else [])
    if C.REGION != "SVI":
        configs = [("a_n50", "", 50), ("b_n20", "", 20)]          # PLAN §5: headline = n >= 20 (last)
    summ, allT = [], []
    for name, tag, nmin in configs:
        cells = cells_for(tag, nmin)
        F = fit(cells); T, comp = q3(F)
        T["config"] = name; allT.append(T)
        mj = T[T.major15]
        r = dict(config=name, n_cells=F["n"], n_series=F["nser"], n_episodes=len(F["eps"]), chi2=F["chi2"],
                 sig_major_median=float(mj.sig.median()), sig_major_min=float(mj.sig.min()), sig_major_max=float(mj.sig.max()),
                 n_major_sig_le_0p1=int((mj.sig <= 0.10).sum()), **(comp or {}))
        r["acceptance"] = bool(r["n_major_sig_le_0p1"] >= 5)
        summ.append(r)
        log(f"{name}: " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
        cells.to_csv(f"{C.TB}/q3_cells_{name}.csv", index=False)
    pd.DataFrame(summ).to_csv(f"{C.TB}/q3_summary.csv", index=False)
    TT = pd.concat(allT, ignore_index=True); TT.to_csv(f"{C.TB}/q3_dD.csv", index=False)
    # ---------------- gates on the headline configuration (last one)
    name, tag, nmin = configs[-1]
    cells = cells_for(tag, nmin); F = fit(cells); T, comp = q3(F); T = T.set_index("episode")
    G = []
    mj = [e for e in MAJ15 if e in T.index]
    for f_ in (2.0, 0.5):
        Tf, cf = q3(fit(cells, prior=PRIOR * f_)); Tf = Tf.set_index("episode")
        mv = np.abs((Tf.dD[mj] - T.dD[mj]) / T.sig[mj])
        mvc = abs(cf["composite"] - comp["composite"]) / comp["composite_sig"]
        G.append(dict(gate="prior", stat=f"x{f_:g}", value=float(mv.max()), composite=float(mvc), passed=bool(mv.max() < 0.1 and mvc < 0.1)))
    cc = cells_for(tag, nmin, control=True)
    Fc = fit(cc); Tc, cc_comp = q3(Fc); Tc = Tc.set_index("episode")
    meas = [e for e in mj if e in Tc.index and Tc.n_cells[e] >= 5]
    frac = float(np.mean([abs(Tc.dD[e]) < 2 * Tc.sig[e] for e in meas])) if meas else np.nan
    Tc.to_csv(f"{C.TB}/q3_control_dD.csv")
    G.append(dict(gate="control", stat="frac_majors_within_2sig", value=frac, n=len(meas), n_cells=len(cc),
                  composite=cc_comp["composite"] if cc_comp else np.nan, composite_sig=cc_comp["composite_sig"] if cc_comp else np.nan,
                  passed=bool(frac >= 0.9 and (cc_comp is None or abs(cc_comp["composite"]) < 2 * cc_comp["composite_sig"]))
                  if meas else None))
    lo = []
    for s in C.SITES:
        q = cells[cells.site != s].reset_index(drop=True)
        Ts, cs = q3(fit(q)); Ts = Ts.set_index("episode")
        ee = [e for e in mj if e in Ts.index]
        ok = [abs(Ts.dD[e] - T.dD[e]) < 2 * Ts.sig[e] for e in ee]
        fr = float(np.mean(ok)) if ok else np.nan
        Ts.assign(left_out=s).to_csv(f"{C.TB}/q3_loso_{s}.csv")
        G.append(dict(gate="loso", stat=f"without_{s}", value=fr, n=len(ee), composite=cs["composite"] if cs else np.nan,
                      composite_sig=cs["composite_sig"] if cs else np.nan, passed=bool(fr >= 0.9) if ok else None))
        lo.append(Ts.assign(left_out=s))
    INJ = injections(cells, a.workers, log); INJ.to_csv(f"{C.TB}/q3_pulses.csv", index=False)
    dp = INJ[(INJ.kind == "deep") & (INJ.episode != "all15")].ratio
    sp = INJ[(INJ.kind == "shal") & (INJ.episode != "all15")].ratio.abs()
    G.append(dict(gate="deep_pulse", stat="median_recovery", value=float(dp.median()), lo=float(dp.min()), hi=float(dp.max()),
                  passed=bool(0.8 <= dp.median() <= 1.1)))
    G.append(dict(gate="shallow_pulse", stat="median_leak", value=float(sp.median()), lo=float(sp.min()), hi=float(sp.max()),
                  composite=float(INJ[(INJ.kind == "shal") & (INJ.episode == "all15")].ratio.iloc[0]) if (INJ.episode == "all15").any() else np.nan,
                  passed=bool(sp.median() <= 0.35)))
    GT = pd.DataFrame(G); GT.to_csv(f"{C.TB}/q3_gates.csv", index=False)
    for r in G:
        log("gate " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    T.reset_index().to_csv(f"{C.TB}/q3_dD_headline.csv", index=False)
    json.dump(dict(config=name, composite=comp, written=C.pt_now()), open(f"{C.TB}/q3_headline.json", "w"), indent=1, default=str)
    figure(TT, configs, Tc)
    log(f"q3 done {time.time() - t0:.0f}s")
    open(f"{C.OUT}/Q3_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


def figure(TT, configs, Tc):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6",
                         "grid.linewidth": 0.5, "font.size": 9})
    head = configs[-1][0]
    T = TT[TT.config == head].sort_values("mid"); T["mid"] = pd.to_datetime(T.mid)
    fig, ax = plt.subplots(2, 1, figsize=(15, 8.5), gridspec_kw=dict(height_ratios=[2.2, 1]))
    ax[0].plot(T.mid, T.dD, "-", color="#9a9a96", lw=0.8, zorder=1)
    sc = ax[0].scatter(T.mid, T.dD, c=np.log10(T.n_cells.clip(lower=1)), cmap="viridis", s=14, zorder=2, lw=0)
    M = T[T.major15]
    ax[0].scatter(M.mid, M.dD, c=np.log10(M.n_cells), cmap="viridis", s=70, edgecolor="k", lw=1.0, zorder=3)
    for r in M.itertuples():
        ax[0].annotate(r.episode, (r.mid, r.dD), xytext=(0, 7), textcoords="offset points", fontsize=7, ha="center")
    ax[0].axhline(0, color="grey", lw=0.6, ls=":")
    ax[0].set(ylabel="ΔD (%): deep-weighted dv/v of the episode\nminus neighbouring episodes (±3 yr)",
              title=f"Q3 = ΔD(e), temporal model, {head}; large outlined dots = the 15 majors; σ in tables (median for majors "
                    f"{M.sig.median():.3f} %), not drawn")
    cb = fig.colorbar(sc, ax=ax[0], pad=0.01); cb.set_label("log10 cells in the episode")
    # composite panel: majors' DeltaD per configuration + composite line
    xs = np.arange(len(MAJ15))
    cols = ["#2a78d6", "#eb6834", "#1baf7a"]
    for i, (name, _, _) in enumerate(configs):
        q = TT[(TT.config == name) & TT.major15].set_index("episode").reindex(MAJ15)
        ax[1].plot(xs, q.dD, "-o", ms=4, lw=1.2, color=cols[i % 3], label=f"{name}: composite {q.dD.mean():+.3f} %")
        ax[1].axhline(q.dD.mean(), color=cols[i % 3], lw=0.8, ls="--")
    if Tc is not None and len(Tc):
        qc = Tc.reindex(MAJ15)
        ax[1].plot(xs, qc.dD, ":", marker="x", ms=4, lw=1, color="#52514e", label="scrambled control (headline)")
    ax[1].set_xticks(xs); ax[1].set_xticklabels(MAJ15, rotation=35, fontsize=7); ax[1].axhline(0, color="grey", lw=0.6, ls=":")
    ax[1].set(ylabel="ΔD (%)", title="the 15 majors per configuration; dashed = 15-major composite (equal-weight mean)")
    ax[1].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(f"{C.FIG}/fig_q3_dD.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
