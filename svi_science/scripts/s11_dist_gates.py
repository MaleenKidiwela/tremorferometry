#!/usr/bin/env python
"""s11_dist_gates -- S11 distant-station lever, step 4 (distant/PLAN.md §8) on F1 (near + far cells, M and B free, NOM):
scrambled control, leave-one-site-out (every site incl. far sites), family split, lapse test on the far cells, far/near
annual amplitude test, site dominance; verdict json. Writes distant/<region>/tables/{gates.csv, gates_summary.json};
marker GATES_DONE.
    S11_REGION=SVI python s11_dist_gates.py
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
import s11_dist_common as D
import s11_dist_model as Mo
import s11_v2 as V
import s11_q3 as Q3

MAJ = Mo.MAJ


def q3_within(Ta, Tb, sig_col="sig_b", frac=0.9):
    m = Ta.set_index("episode")[["dD", "sig", "major15"]].join(Tb.set_index("episode")[["dD", "sig"]], rsuffix="_b", how="inner")
    m = m[m.major15]
    if not len(m):
        return np.nan, 0
    ok = np.abs(m.dD_b - m.dD) < 2 * m[sig_col]
    return float(ok.mean()), len(m)


def ratio_err(a, sa, b, sb):
    r = a / b
    return r, abs(r) * np.hypot(sa / a, sb / b)


def main():
    ap = argparse.ArgumentParser(); ap.parse_args()
    D.ensure(); log = C.Log(f"{D.DLOG}/gates{Mo.FAR_TAG}.log"); t0 = time.time()
    cells = Mo.all_cells()
    if not cells.far.any():
        log("no far cells -> gates not run"); return
    G = []
    Fs = Mo.fit(cells, "smooth"); q = Mo.Qs(Fs)
    Ft = Mo.fit(cells, "temporal"); T, comp = Mo.q3tab(Ft)
    G.append(dict(gate="F1", stat="headline", Q1=q["Q1"], Q1_sig=q["Q1_sig"], peak=q["Q1_peak"], Q2=q["Q2"], Q2_sig=q["Q2_sig"],
                  composite=comp["composite"], composite_sig=comp["composite_sig"]))
    pk0 = q["Q1_peak"]; dpk = lambda p: abs((p - pk0 + 182.6) % 365.25 - 182.6)
    # ---------------- scrambled control
    cc = Mo.all_cells(control=True)
    qc = Mo.Qs(Mo.fit(cc, "smooth")); Tc, cpc = Mo.q3tab(Mo.fit(cc, "temporal"))
    ncell = cc.groupby("episode").size()
    mj = Tc[Tc.major15 & Tc.episode.map(ncell).ge(5)]
    fr = float((np.abs(mj.dD) < 2 * mj.sig).mean()) if len(mj) else np.nan
    G.append(dict(gate="control", stat="Q1", Q1=qc["Q1"], Q1_sig=qc["Q1_sig"], n=len(cc), n_far=int(cc.far.sum()), passed=bool(qc["Q1"] < 2 * qc["Q1_sig"])))
    G.append(dict(gate="control", stat="Q3", frac_within=fr, n_majors=len(mj), composite=cpc.get("composite"), composite_sig=cpc.get("composite_sig"),
                  passed=bool(fr >= 0.9 and abs(cpc.get("composite", 0)) < 2 * cpc.get("composite_sig", np.inf))))
    # far-only control (column)
    ccf = cc[cc.far].reset_index(drop=True)
    if len(ccf) >= 30:
        dm = V.demeaned(ccf)
        G.append(dict(gate="control", stat="far_only_chi2_about_zero", value=float(np.mean((ccf.d / ccf.sig) ** 2)), n=len(ccf)))
    # ---------------- leave-one-site-out
    ok1 = ok3 = True
    for s in sorted(cells.site.unique()):
        c_ = cells[cells.site != s].reset_index(drop=True)
        qs = Mo.Qs(Mo.fit(c_, "smooth")); Ts, _ = Mo.q3tab(Mo.fit(c_, "temporal"))
        w1 = abs(qs["Q1"] - q["Q1"]) < 2 * qs["Q1_sig"]
        p1 = dpk(qs["Q1_peak"]) <= 60 if qs["Q1"] / qs["Q1_sig"] >= 1 else True
        f3, n3 = q3_within(T, Ts)
        ok1 &= bool(w1 and p1); ok3 &= bool(f3 >= 0.9) if n3 else True
        G.append(dict(gate="LOSO", stat=f"without_{s}", far_site=bool(cells[cells.site == s].far.all()), Q1=qs["Q1"], Q1_sig=qs["Q1_sig"],
                      peak=qs["Q1_peak"], Q2=qs["Q2"], Q2_sig=qs["Q2_sig"], q3_frac_within=f3, n_majors=n3, passed=bool(w1 and p1 and (f3 >= 0.9 if n3 else True))))
    G.append(dict(gate="LOSO", stat="all_sites_Q1", passed=bool(ok1))); G.append(dict(gate="LOSO", stat="all_sites_Q3", passed=bool(ok3)))
    # ---------------- family split
    co = cells[cells.fam % 2 == 1].reset_index(drop=True); ce = cells[cells.fam % 2 == 0].reset_index(drop=True)
    qo, qe = Mo.Qs(Mo.fit(co, "smooth")), Mo.Qs(Mo.fit(ce, "smooth"))
    w2 = abs(qo["Q1"] - qe["Q1"]) < 2 * np.hypot(qo["Q1_sig"], qe["Q1_sig"])
    pko = abs((qo["Q1_peak"] - qe["Q1_peak"] + 182.6) % 365.25 - 182.6) <= 60 if min(qo["Q1"] / qo["Q1_sig"], qe["Q1"] / qe["Q1_sig"]) >= 1 else True
    To, co_ = Mo.q3tab(Mo.fit(co, "temporal")); Te, ce_ = Mo.q3tab(Mo.fit(ce, "temporal"))
    m = To.set_index("episode")[["dD", "sig", "major15"]].join(Te.set_index("episode")[["dD", "sig"]], rsuffix="_e", how="inner")
    m = m[m.major15]
    f_fs = float((np.abs(m.dD - m.dD_e) < 2 * np.hypot(m.sig, m.sig_e)).mean()) if len(m) else np.nan
    cok = abs(co_.get("composite", 0) - ce_.get("composite", 0)) < 2 * np.hypot(co_.get("composite_sig", np.inf), ce_.get("composite_sig", np.inf))
    G.append(dict(gate="family_split", stat="Q1", Q1_odd=qo["Q1"], Q1_odd_sig=qo["Q1_sig"], Q1_even=qe["Q1"], Q1_even_sig=qe["Q1_sig"], passed=bool(w2 and pko)))
    G.append(dict(gate="family_split", stat="Q3", composite_odd=co_.get("composite"), composite_odd_sig=co_.get("composite_sig"),
                  composite_even=ce_.get("composite"), composite_even_sig=ce_.get("composite_sig"), frac_within=f_fs, n_majors=len(m),
                  passed=bool(cok and f_fs >= 0.9)))
    # ---------------- kernel-free annual amplitudes: near pooled peak, far/near ratio, lapse
    near = cells[~cells.far].copy(); far = cells[cells.far].copy()
    near["dm"] = V.demeaned(near)
    pool = C.annual_fit(near.episode.map(Mo.U.mid), near.dm.to_numpy(), near.sig.to_numpy()); pkp = pool["peak_doy"]
    An, sAn = V.inphase(near, pkp); Af, sAf = V.inphase(far, pkp)
    ff = far.copy(); ff["dm"] = V.demeaned(ff)
    poolf = C.annual_fit(ff.episode.map(Mo.U.mid), ff.dm.to_numpy(), ff.sig.to_numpy())
    wn, wf = 1 / near.sig ** 2, 1 / far.sig ** 2
    rr, srr = ratio_err(Af, sAf, An, sAn)
    pdp = float(np.average(far.a_d, weights=wf) / np.average(near.a_d, weights=wn))
    psp = float(np.average(far.a_s, weights=wf) / np.average(near.a_s, weights=wn))
    zd, zs = (rr - pdp) / srr, (rr - psp) / srr
    G.append(dict(gate="far_near_amp", stat="ratio_far_over_near", value=rr, sig=srr, A_near=An, A_near_sig=sAn, A_far=Af, A_far_sig=sAf,
                  near_peak=pkp, far_free_amp=poolf["amp"], far_free_amp_sig=poolf["amp_sig"], far_free_peak=poolf["peak_doy"],
                  pred_deep=pdp, pred_shallow=psp, z_deep=zd, z_shallow=zs, passed=bool(abs(zs) - abs(zd) >= 1)))
    for grp, cc_ in (("far", True), ("near", False)):
        cA = Mo.all_cells(win="2035"); cB = Mo.all_cells(win="3560")
        cA = cA[cA.far == cc_]; cB = cB[cB.far == cc_]
        if len(cA) < 20 or len(cB) < 20:
            G.append(dict(gate="lapse", stat=grp, n_2035=len(cA), n_3560=len(cB), passed=None)); continue
        A1, s1 = V.inphase(cA, pkp); A2, s2 = V.inphase(cB, pkp)
        r_, sr_ = ratio_err(A2, s2, A1, s1)
        pdl = float(cB.a_d.mean() / cA.a_d.mean()); psl = float(cB.a_s.mean() / cA.a_s.mean())
        G.append(dict(gate="lapse", stat=grp, value=r_, sig=sr_, A2035=A1, A2035_sig=s1, A3560=A2, A3560_sig=s2, pred_deep=pdl, pred_shallow=psl,
                      z_deep=(r_ - pdl) / sr_, z_shallow=(r_ - psl) / sr_, n_2035=len(cA), n_3560=len(cB),
                      passed=bool(abs((r_ - psl) / sr_) - abs((r_ - pdl) / sr_) >= 1)))
    # ---------------- site dominance
    w = 1 / cells.sig ** 2; sw = (w.groupby(cells.site).sum() / w.sum()).sort_values(ascending=False)
    G.append(dict(gate="site_dominance", stat="max_share", value=float(sw.iloc[0]), site=sw.index[0], passed=bool(sw.iloc[0] <= 0.5)))
    GT = pd.DataFrame(G); GT.to_csv(f"{D.DTB}/gates{Mo.FAR_TAG}.csv", index=False)
    for r in G:
        log("gate " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    summ = {f"{r['gate']}:{r['stat']}": r.get("passed") for r in G if "passed" in r}
    summ["written"] = C.pt_now()
    json.dump(summ, open(f"{D.DTB}/gates{Mo.FAR_TAG}_summary.json", "w"), indent=1, default=str)
    open(f"{D.DOUT}/GATES{Mo.FAR_TAG}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"gates done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
