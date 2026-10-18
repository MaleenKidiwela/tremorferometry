#!/usr/bin/env python
"""s11_dist_invert -- S11 distant-station lever, step 3 (distant/PLAN.md §6-§7): natural fits for N0 / N1 / F1 / F0 / FO
(smooth: Q1, Q2, M/B annual; temporal: Q3 per major + composite), N0 reproduction check, waveform-level injections through
the full raw chain on the near AND far stores (A deep step, B shallow annual, Ps / Pd one-major pulses, Ps(all majors),
mid-crust annual; kernel-error generation with two5_40 and ell80), responses per configuration and kernel setup (NOM,
AVG4), lever verdicts. Writes distant/<region>/tables/{inv_*.csv, inv_summary.json}; marker INVERT_DONE.
    S11_REGION=SVI python s11_dist_invert.py --workers 10
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s11_common as C
import s11_dist_common as D
import s11_dist_model as Mo
import s11_invert as I
import s11_paths as P
import s11_v2 as V
import s11_q3 as Q3

GEN = ("ell40", "two5_40", "ell80")
HEAD = C.HEAD
MAJ = Mo.MAJ


def configs(cells):
    nr = cells[~cells.far].reset_index(drop=True); fr = cells[cells.far].reset_index(drop=True)
    out = {"N0": (nr, False), "N1": (nr, True)}
    if len(fr):
        out.update({"F1": (cells, True), "F0": (cells, False), "FO": (fr, True)})
    return out


# ------------------------------------------------------------------------------ injections
_G = {}


def _ser(key):
    try:
        st = _G["store"][key]; unit = _G["unit"]; K = _G["K"]
        rc = lambda amp=None: P.raw_cells(st, unit, amp=amp, wins=[HEAD])[["episode", f"eps3c_{HEAD}"]]
        out = {"nat": rc()}
        t = pd.to_datetime(pd.Series(st["days"]))
        post = (t >= "2015-01-01").to_numpy().astype(float)
        ann = np.cos(2 * np.pi * (t.dt.dayofyear.to_numpy() - 233.0) / 365.25)
        ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object)
        inmaj = np.isin(ep, MAJ).astype(float)
        for g in GEN:
            if key not in K[g].index:
                continue
            a_s, a_d = float(K[g].loc[key, "a_s"]), float(K[g].loc[key, "a_d"])
            out[f"A_{g}"] = rc(-(0.003 * a_d) * post)
            out[f"B_{g}"] = rc(-(0.001 * a_s) * ann)
            if inmaj.any():
                out[f"Psall_{g}"] = rc(-(0.001 * a_s) * inmaj)
        if key in K["ell40"].index:
            k = K["ell40"].loc[key]
            out["Mann"] = rc(-(0.001 * float(k.a_m)) * ann)
            for e in MAJ:
                m = (ep == e).astype(float)
                if not m.any():
                    continue
                out[f"Pd_{e}"] = rc(-(0.003 * float(k.a_d)) * m)
                out[f"Ps_{e}"] = rc(-(0.001 * float(k.a_s)) * m)
        return key, out
    except Exception as ex:
        import traceback
        print("ERROR", key, ex, traceback.format_exc()[-600:], flush=True); return key, {}


def run_injections(cells, workers, log):
    import multiprocessing as mp
    unit = P.unit_map()
    need = set(cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch)
    vals = {}
    jobs = []                                   # (grid, kernel tables, store loader)
    Kn = {g: Mo._kern(g, False) for g in GEN}; Kn = {g: k[k.win == HEAD].set_index("series") for g, k in Kn.items()}
    near_st = sorted(cells[~cells.far].netsta.unique()); far_st = sorted(cells[cells.far].netsta.unique())
    for far, stas in ((False, near_st), (True, far_st)):
        if far:
            Kf = {g: Mo._kern(g, True) for g in GEN}; K = {g: k[k.win == HEAD].set_index("series") for g, k in Kf.items()}
        else:
            K = Kn
        D.set_grid(D.FAR_GRID if far else D.NEAR_GRID)
        for ns in stas:
            store = {}
            if far:
                with Mo.tb(D.DTB):
                    store.update(P.load_store(ns, Mo.FAR_TAG))
            else:
                for tg in ("", "_ext"):
                    if os.path.exists(f"{C.TB}/store{tg}/{ns}.npz"):
                        store.update(P.load_store(ns, tg))
            keys = [k for k in store if k in need]
            if not keys:
                continue
            sp_ = f"{D.DTB}/store{Mo.FAR_TAG}/{ns}.npz" if far else f"{C.TB}/store/{ns}.npz"
            tau = np.load(sp_)["tau"] if os.path.exists(sp_) else None
            if tau is not None and len(tau) != S.NW:
                raise RuntimeError(f"grid mismatch at {ns}: store {len(tau)} vs grid {S.NW}")
            _G.clear(); _G.update(store=store, unit=unit, K=K)
            t0 = time.time()
            with mp.get_context("fork").Pool(workers) as pool:
                for k, o in pool.imap_unordered(_ser, keys):
                    vals[k] = o
            log(f"  injections {'far' if far else 'near'} {ns}: {len(keys)} series {time.time() - t0:.0f}s")
    D.set_grid(D.NEAR_GRID)
    ckey = (cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch).to_numpy()
    runs = sorted({r for o in vals.values() for r in o})
    out = {}
    for run in runs:
        d = np.full(len(cells), np.nan)
        for i, (k, e) in enumerate(zip(ckey, cells.episode)):
            o = vals.get(k, {}); q = o.get(run, o.get("nat") if run.startswith(("Pd_", "Ps_", "Psall")) else None)
            if q is None or not len(q):
                continue
            mm = (q.episode == e).to_numpy()
            if mm.any():
                d[i] = -100 * float(q[f"eps3c_{HEAD}"].to_numpy()[mm][0])
        out[run] = d
    dn = out["nat"]
    for far in (False, True):
        m = (cells.far == far).to_numpy() & np.isfinite(dn)
        if m.any():
            log(f"  natural-from-store vs measured ({'far' if far else 'near'}): max |diff| {np.nanmax(np.abs(dn[m] - cells.d.to_numpy()[m])):.1e} %, "
                f"missing {int(((cells.far == far).to_numpy() & ~np.isfinite(dn)).sum())}")
    return out


def responses(cfgs, inj, log):
    """injected - natural responses for every configuration (cells common to all configurations: finite natural values)"""
    rows = []
    dn_all = inj["nat"]
    for name, (cc, MB) in cfgs.items():
        sel = cc._row.to_numpy()
        for setup, ccs in (("NOM", cc), ("AVG4", cc.assign(**{k: cc[f"{k}_AVG4"] for k in ("a_s", "a_m", "a_d", "a_b")}))):
            dn = dn_all[sel]; okn = np.isfinite(dn)
            kinds = ("smooth", "temporal") + (("temporal_Sc",) if name in ("N1", "F1") and setup == "NOM" else ())
            for kind_ in kinds:
                kind, Sc = ("temporal", True) if kind_ == "temporal_Sc" else (kind_, False)
                setup_ = "NOM_Sc" if Sc else setup
                def resp(run):
                    di = inj.get(run, np.full(len(dn_all), np.nan))[sel]; ok = okn & np.isfinite(di)
                    if ok.sum() < 20:
                        return None, None
                    c_ = ccs[ok].reset_index(drop=True)
                    F = Mo.fit(c_, kind, MB=MB, d=(di - dn)[ok], Sc=Sc)
                    return F, ok
                if kind == "smooth":
                    for g in GEN:
                        if setup == "AVG4" and g == "ell40":
                            pass
                        F, ok = resp(f"A_{g}")
                        if F is not None:
                            A1, _, _, _ = V.amp(F); q2 = float(V.q2vec(F) @ F["th"])
                            rows.append(dict(config=name, setup=setup, run="A", gen=g, metric="Q2_recovery", value=q2 / 0.3, dQ1=A1, n=int(ok.sum())))
                        F, ok = resp(f"B_{g}")
                        if F is not None:
                            A1, _, pk, _ = V.amp(F)
                            rows.append(dict(config=name, setup=setup, run="B", gen=g, metric="Q1_leak", value=A1 / 0.1, peak=pk, n=int(ok.sum())))
                    F, ok = resp("Mann")
                    if F is not None:
                        rows.append(dict(config=name, setup=setup, run="Mann", gen="ell40", metric="Q1_leak", value=V.amp(F)[0] / 0.1, n=int(ok.sum())))
                else:
                    for g in GEN:
                        F, ok = resp(f"Psall_{g}")
                        if F is not None:
                            gc = Mo.comp_vec(F)
                            if gc is not None:
                                rows.append(dict(config=name, setup=setup_, run="Psall", gen=g, metric="composite_leak", value=float(gc @ F["th"]) / 0.1, n=int(ok.sum())))
                        F, ok = resp(f"B_{g}")
                        if F is not None:
                            gc = Mo.comp_vec(F)
                            if gc is not None:
                                rows.append(dict(config=name, setup=setup_, run="B", gen=g, metric="composite_leak", value=float(gc @ F["th"]) / 0.1, n=int(ok.sum())))
                    if setup_ == "AVG4":
                        continue
                    for e in MAJ:
                        for kk, amp_ in (("Pd", 0.3), ("Ps", 0.1)):
                            F, ok = resp(f"{kk}_{e}")
                            if F is None or ("D", "", e) not in F["idx"]:
                                continue
                            g = Q3.dvec(F, e)
                            if g is None:
                                continue
                            rows.append(dict(config=name, setup=setup_, run=kk, gen="ell40", major=e, metric="ratio", value=float(g @ F["th"]) / amp_, n=int(ok.sum())))
        log(f"  responses {name} done")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--no-inject", action="store_true")
    a = ap.parse_args()
    D.ensure(); log = C.Log(f"{D.DLOG}/invert{Mo.FAR_TAG}.log"); t0 = time.time()
    cells = Mo.all_cells()
    for v in Mo.VARS4:
        cv = pd.concat([Mo.with_variant(cells[~cells.far], v, False), Mo.with_variant(cells[cells.far], v, True)]) if cells.far.any() \
            else Mo.with_variant(cells, v, False)
        cv = cv.loc[cells.index]
        for k in ("a_s", "a_m", "a_d", "a_b"):
            cells[f"{k}_{v}"] = cv[k].to_numpy()
    for k in ("a_s", "a_m", "a_d", "a_b"):
        cells[f"{k}_AVG4"] = cells[[f"{k}_{v}" for v in Mo.VARS4]].mean(1)
    cells["_row"] = np.arange(len(cells)); cells["mid"] = cells.episode.map(Mo.U.mid)
    cells.to_csv(f"{D.DTB}/inv{Mo.FAR_TAG}_cells.csv", index=False)
    fr = cells[cells.far]
    log(f"cells: near {int((~cells.far).sum())}, far {len(fr)} in {fr.fam.astype(str).add(fr.netsta + fr.epoch).nunique()} far series, "
        f"{fr.netsta.nunique()} far stations, sites {sorted(fr.site.unique())}; far a_s/a_m/a_d/a_b median "
        f"{fr.a_s.median():.3f}/{fr.a_m.median():.3f}/{fr.a_d.median():.3f}/{fr.a_b.median():.3f}; near "
        f"{cells[~cells.far].a_s.median():.3f}/{cells[~cells.far].a_m.median():.3f}/{cells[~cells.far].a_d.median():.3f}/{cells[~cells.far].a_b.median():.3f}")
    w = 1 / cells.sig ** 2
    sw = w.groupby(cells.site).sum() / w.sum()
    log(f"weight share by site: {sw.sort_values(ascending=False).round(3).to_dict()}; far share {float(w[cells.far].sum() / w.sum()):.3f}")
    cfgs = configs(cells)
    # ---------------- natural fits
    rows, q3rows = [], []
    for name, (cc, MB) in cfgs.items():
        for setup in ("NOM", "AVG4"):
            ccs = cc if setup == "NOM" else cc.assign(**{k: cc[f"{k}_AVG4"] for k in ("a_s", "a_m", "a_d", "a_b")})
            Fs = Mo.fit(ccs, "smooth", MB=MB); q = Mo.Qs(Fs)
            Ft = Mo.fit(ccs, "temporal", MB=MB); T, comp = Mo.q3tab(Ft)
            mj = T[T.major15]
            r = dict(config=name, setup=setup, MB=MB, n_far=int(cc.far.sum()), **q,
                     q3_composite=comp.get("composite"), q3_composite_sig=comp.get("composite_sig"), n_majors=comp.get("n_majors"),
                     q3_sig_major_median=float(mj.sig.median()) if len(mj) else np.nan, q3_n_major_sig_le_0p1=int((mj.sig <= 0.1).sum()),
                     q3_chi2_majors=float(np.mean((mj.dD / mj.sig) ** 2)) if len(mj) else np.nan, q3_chi2fit=Ft["chi2"])
            rows.append(r)
            T["config"] = name; T["setup"] = setup; q3rows.append(T)
            log(f"{name}/{setup}: Q1 {q['Q1']:.4f} ± {q['Q1_sig']:.4f} (day {q['Q1_peak']:.0f}), Q2 {q['Q2']:+.4f} ± {q['Q2_sig']:.4f}, "
                f"Q3 comp {comp.get('composite', np.nan):+.4f} ± {comp.get('composite_sig', np.nan):.4f}, σ major median {r['q3_sig_major_median']:.4f} "
                f"({r['q3_n_major_sig_le_0p1']} ≤ 0.1), M {q.get('M_annual', np.nan):.3f} ± {q.get('M_annual_sig', np.nan):.3f}, "
                f"B {q.get('B_annual', np.nan):.3f} ± {q.get('B_annual_sig', np.nan):.3f}")
    for name in ("N1", "F1"):                       # Amendment 1 column: temporal + episode-locked regional shallow Sc(e)
        if name not in cfgs:
            continue
        cc, MB = cfgs[name]
        Ft = Mo.fit(cc, "temporal", MB=MB, Sc=True); T, comp = Mo.q3tab(Ft); mj = T[T.major15]
        rows.append(dict(config=name, setup="NOM_Sc", MB=MB, n_far=int(cc.far.sum()), n_cells=len(cc), q3_composite=comp.get("composite"),
                         q3_composite_sig=comp.get("composite_sig"), n_majors=comp.get("n_majors"),
                         q3_sig_major_median=float(mj.sig.median()) if len(mj) else np.nan, q3_n_major_sig_le_0p1=int((mj.sig <= 0.1).sum()),
                         q3_chi2_majors=float(np.mean((mj.dD / mj.sig) ** 2)) if len(mj) else np.nan, q3_chi2fit=Ft["chi2"]))
        T["config"] = name; T["setup"] = "NOM_Sc"; q3rows.append(T)
        log(f"{name}/T+Sc: Q3 comp {comp.get('composite', np.nan):+.4f} ± {comp.get('composite_sig', np.nan):.4f}, σ major median "
            f"{rows[-1]['q3_sig_major_median']:.4f} ({rows[-1]['q3_n_major_sig_le_0p1']} ≤ 0.1)")
    NM = pd.DataFrame(rows); NM.to_csv(f"{D.DTB}/inv{Mo.FAR_TAG}_models.csv", index=False)
    pd.concat(q3rows, ignore_index=True).to_csv(f"{D.DTB}/inv{Mo.FAR_TAG}_q3.csv", index=False)
    # ---------------- N0 reproduction (stored near results)
    n0 = NM[(NM.config == "N0") & (NM.setup == "NOM")].iloc[0]
    ref = pd.read_csv(f"{D.NEAR_TB}/q3_summary.csv") if os.path.exists(f"{D.NEAR_TB}/q3_summary.csv") else None
    if ref is not None:
        rr = ref.iloc[-1]
        dq = abs(float(rr.composite) - n0.q3_composite); ds = abs(float(rr.composite_sig) - n0.q3_composite_sig)
        log(f"N0 reproduction vs {D.NEAR_TB}/q3_summary.csv ({rr.config}): composite {rr.composite:+.4f} ± {rr.composite_sig:.4f} vs "
            f"{n0.q3_composite:+.4f} ± {n0.q3_composite_sig:.4f}; |Δ| {dq:.1e} / {ds:.1e} -> {'OK' if max(dq, ds) <= 0.002 else 'MISMATCH'}")
        if max(dq, ds) > 0.002:
            open(f"{D.DOUT}/INVERT{Mo.FAR_TAG}_STOP", "w").write(f"N0 reproduction mismatch {dq:.4f} {ds:.4f}\n"); return
    # prior test on F1 (and N1)
    pri = []
    for name in ("N1", "F1"):
        if name not in cfgs:
            continue
        cc, MB = cfgs[name]
        Fs = Mo.fit(cc, "smooth", MB=MB); q = Mo.Qs(Fs); Ft = Mo.fit(cc, "temporal", MB=MB); T, comp = Mo.q3tab(Ft)
        for f_ in (2.0, 0.5):
            q_ = Mo.Qs(Mo.fit(cc, "smooth", MB=MB, prior=f_)); T_, c_ = Mo.q3tab(Mo.fit(cc, "temporal", MB=MB, prior=f_))
            m = T.set_index("episode").join(T_.set_index("episode")[["dD"]], rsuffix="_p")
            mj = m[m.major15]
            Ts0, cs0 = Mo.q3tab(Mo.fit(cc, "temporal", MB=MB, Sc=True)); Ts1, cs1 = Mo.q3tab(Mo.fit(cc, "temporal", MB=MB, prior=f_, Sc=True))
            ms = Ts0.set_index("episode").join(Ts1.set_index("episode")[["dD"]], rsuffix="_p"); ms = ms[ms.major15]
            pri.append(dict(config=name, prior_factor=f_, Sc_dcomp_over_sig=(cs1["composite"] - cs0["composite"]) / cs0["composite_sig"],
                            Sc_max_major_over_sig=float(np.max(np.abs(ms.dD_p - ms.dD) / ms.sig)), dQ1_over_sig=(q_["Q1"] - q["Q1"]) / q["Q1_sig"], dQ2_over_sig=(q_["Q2"] - q["Q2"]) / q["Q2_sig"],
                            dcomp_over_sig=(c_["composite"] - comp["composite"]) / comp["composite_sig"],
                            max_major_over_sig=float(np.max(np.abs(mj.dD_p - mj.dD) / mj.sig)), M_annual=q_.get("M_annual"), B_annual=q_.get("B_annual")))
    PR = pd.DataFrame(pri); PR.to_csv(f"{D.DTB}/inv{Mo.FAR_TAG}_prior.csv", index=False)
    log("prior test:\n" + PR.round(3).to_string())
    if a.no_inject:
        return
    # ---------------- injections
    inj = run_injections(cells, a.workers, log)
    np.savez_compressed(f"{D.DTB}/inv{Mo.FAR_TAG}_injected.npz", runs=np.array(sorted(inj)), d=np.array([inj[k] for k in sorted(inj)]))
    R = responses(cfgs, inj, log); R.to_csv(f"{D.DTB}/inv{Mo.FAR_TAG}_responses.csv", index=False)
    # ---------------- summary
    summ = dict(region=D.REG, written=C.pt_now(), n_near=int((~cells.far).sum()), n_far=int(cells.far.sum()))
    for name in cfgs:
        g = R[(R.config == name) & (R.setup == "NOM")]
        pd_ = g[g.run == "Pd"].value; ps_ = g[g.run == "Ps"].value
        def one(run, gen, metric):
            x = g[(g.run == run) & (g.gen == gen) & (g.metric == metric)].value
            return float(x.iloc[0]) if len(x) else np.nan
        nm = NM[(NM.config == name) & (NM.setup == "NOM")].iloc[0]
        summ[name] = dict(deep_pulse_median=float(pd_.median()) if len(pd_) else np.nan,
                          shallow_pulse_median_abs=float(ps_.abs().median()) if len(ps_) else np.nan,
                          shallow_pulse_median=float(ps_.median()) if len(ps_) else np.nan,
                          Psall_composite_leak=one("Psall", "ell40", "composite_leak"),
                          A_Q2_recovery={k: one("A", k, "Q2_recovery") for k in GEN},
                          B_Q1_leak={k: one("B", k, "Q1_leak") for k in GEN},
                          B_composite_leak={k: one("B", k, "composite_leak") for k in GEN},
                          Psall_composite_leak_gen={k: one("Psall", k, "composite_leak") for k in GEN},
                          Mann_Q1_leak=one("Mann", "ell40", "Q1_leak"),
                          Q1=float(nm.Q1), Q1_sig=float(nm.Q1_sig), Q2=float(nm.Q2), Q2_sig=float(nm.Q2_sig),
                          q3_composite=float(nm.q3_composite), q3_composite_sig=float(nm.q3_composite_sig),
                          q3_sig_major_median=float(nm.q3_sig_major_median), q3_n_major_sig_le_0p1=int(nm.q3_n_major_sig_le_0p1))
    for name in ("N1", "F1"):
        g = R[(R.config == name) & (R.setup == "NOM_Sc")]
        if name in summ and len(g):
            nm = NM[(NM.config == name) & (NM.setup == "NOM_Sc")].iloc[0]
            summ[name]["Sc_model"] = dict(deep_pulse_median=float(g[g.run == "Pd"].value.median()), shallow_pulse_median_abs=float(g[g.run == "Ps"].value.abs().median()),
                                          Psall_composite_leak=float(g[(g.run == "Psall") & (g.gen == "ell40")].value.iloc[0]) if len(g[(g.run == "Psall") & (g.gen == "ell40")]) else np.nan,
                                          q3_composite=float(nm.q3_composite), q3_composite_sig=float(nm.q3_composite_sig),
                                          q3_sig_major_median=float(nm.q3_sig_major_median))
    if "F1" in summ and "N1" in summ:
        f, n = summ["F1"], summ["N1"]
        summ["lever_works"] = bool(f["shallow_pulse_median_abs"] <= 0.8 * n["shallow_pulse_median_abs"] and f["q3_composite_sig"] <= n["q3_composite_sig"])
        summ["leak_factor_F1_over_N1"] = f["shallow_pulse_median_abs"] / n["shallow_pulse_median_abs"]
        summ["sig_composite_factor_F1_over_N1"] = f["q3_composite_sig"] / n["q3_composite_sig"]
        summ["sig_major_factor_F1_over_N1"] = f["q3_sig_major_median"] / n["q3_sig_major_median"]
        bl = [f["B_Q1_leak"][k] <= 0.8 * n["B_Q1_leak"][k] for k in GEN]
        summ["deep_annual_ambiguity_reduced"] = bool(all(bl) and f["Q1_sig"] <= n["Q1_sig"])
        summ["B_leak_factor_F1_over_N1"] = {k: f["B_Q1_leak"][k] / n["B_Q1_leak"][k] for k in GEN}
        summ["Q1_sig_factor_F1_over_N1"] = f["Q1_sig"] / n["Q1_sig"]
    json.dump(summ, open(f"{D.DTB}/inv{Mo.FAR_TAG}_summary.json", "w"), indent=1, default=float)
    log("SUMMARY " + json.dumps(summ, default=float))
    open(f"{D.DOUT}/INVERT{Mo.FAR_TAG}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"invert done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
