#!/usr/bin/env python
"""s11_dist_b2 -- post-run column (coordinator request, 2026-10-12 ~9:15 pm PT; logged as post hoc in distant/NOTES):
(1) the best distant-lever model with the round-3 adopted B2 graph estimator for the NEAR cells (>= 50-event nodes,
    cc weights, bootstrap sigma; built exactly as s11_r3_q3.build) + the far B0 cells; temporal and T+Sc, M and B free;
    sigma per major, composite, DESIGN-level leaks (noise-free linear injection through the same model: 0.1 % x a_s in one
    major -> shallow leak; 0.3 % x a_d -> deep recovery; the waveform-level over-read of the reference path is absent
    in B2 by construction, so the design leak is B2's leak; calibrated by computing it also for F1-B0 where the
    waveform-level value is known), prior x2/x0.5, leave-PGC-out;
(2) SVI 2019-11 across every configuration + the far cells in that episode;
(3) same-station near-vs-far annual amplitudes (kernel-free, free phase) at stations used both as near and far.
Writes distant/<region>/tables/b2_*.csv. Read-only on round3/.
    S11_REGION=SVI python s11_dist_b2.py
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_dist_common as D
import s11_dist_model as Mo
import s11_invert as I
import s11_v2 as V
import s11_q3 as Q3
import s11_r3_q3 as RQ

MAJ = Mo.MAJ


def b2_near():
    E = RQ.r3_tables()
    q, _ = RQ.build(E, "B2")
    base = q.drop(columns=[c for c in ("a_s", "a_m", "a_d", "a_b", "series") if c in q])
    k = I.attach_kernels(base, Mo._kern("ell40", False), C.HEAD)
    k = k[np.isfinite(k.a_d) & np.isfinite(k.a_s)].reset_index(drop=True)
    k["far"] = False
    return k


def design(cells, Sc):
    """noise-free linear injections through the temporal model (M, B free): median shallow leak, median deep recovery"""
    lk, rc = [], []
    for e in MAJ:
        m = (cells.episode == e).to_numpy()
        if not m.any():
            continue
        for kind, amp, col, out in (("Ps", 0.1, "a_s", lk), ("Pd", 0.3, "a_d", rc)):
            d = np.where(m, amp * cells[col].to_numpy(), 0.0)
            F = Mo.fit(cells, "temporal", MB=True, d=d, Sc=Sc)
            if ("D", "", e) not in F["idx"]:
                continue
            g = Q3.dvec(F, e)
            if g is not None:
                out.append(float(g @ F["th"]) / amp)
    return (float(np.median(np.abs(lk))) if lk else np.nan), (float(np.median(rc)) if rc else np.nan)


def summarize(name, cells, Sc, log):
    F = Mo.fit(cells, "temporal", MB=True, Sc=Sc); T, comp = Mo.q3tab(F); mj = T[T.major15]
    lk, rc = design(cells, Sc)
    pr = []
    for f_ in (2.0, 0.5):
        Tf, cf = Mo.q3tab(Mo.fit(cells, "temporal", MB=True, prior=f_, Sc=Sc))
        mm = T.set_index("episode").join(Tf.set_index("episode")[["dD"]], rsuffix="_p"); mm = mm[mm.major15]
        pr.append(float(np.max(np.abs(mm.dD_p - mm.dD) / mm.sig)))
    r = dict(config=name, Sc=Sc, n_cells=len(cells), n_far=int(cells.far.sum()), n_majors=int(len(mj)),
             sig_major_median=float(mj.sig.median()), n_major_sig_le_0p1=int((mj.sig <= 0.1).sum()),
             composite=comp.get("composite"), composite_sig=comp.get("composite_sig"),
             chi2_majors=float(np.mean((mj.dD / mj.sig) ** 2)), n_beyond_2sig=int((np.abs(mj.dD / mj.sig) > 2).sum()),
             n_beyond_3sig=int((np.abs(mj.dD / mj.sig) >= 3).sum()), design_shallow_leak=lk, design_deep_recovery=rc,
             prior_x2_max_major=pr[0], prior_x05_max_major=pr[1])
    log(", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    T = T.assign(config=name, Sc=Sc)
    return r, T


def main():
    D.ensure(); log = C.Log(f"{D.DLOG}/b2.log")
    far = Mo.far_cells()
    nB0 = Mo.near_cells()
    nB2 = b2_near()
    sm = Mo.site_map(sorted(far.netsta.unique()))
    log(f"near B0 {len(nB0)}, near B2 {len(nB2)} (episodes {nB2.episode.nunique()}), far B0 {len(far)}")
    cfg = {"N1_B0": nB0, "F1_B0": pd.concat([nB0, far], ignore_index=True), "N1_B2": nB2, "F1_B2": pd.concat([nB2, far], ignore_index=True),
           "FO_B0": far.reset_index(drop=True)}
    rows, Ts = [], []
    for name, cc in cfg.items():
        for Sc in (False, True):
            r, T = summarize(name, cc.reset_index(drop=True), Sc, log); rows.append(r); Ts.append(T)
    pd.DataFrame(rows).to_csv(f"{D.DTB}/b2_models.csv", index=False)
    TT = pd.concat(Ts, ignore_index=True); TT.to_csv(f"{D.DTB}/b2_q3.csv", index=False)
    # leave-PGC-out / leave-one-site-out for the B2 combined models (majors' dD within 2 sigma)
    lo = []
    for name in ("F1_B2",):
        cc = cfg[name]
        for Sc in (False, True):
            T0, _ = Mo.q3tab(Mo.fit(cc, "temporal", MB=True, Sc=Sc))
            for s in sorted(cc.site.unique()):
                c_ = cc[cc.site != s].reset_index(drop=True)
                Ts_, cs_ = Mo.q3tab(Mo.fit(c_, "temporal", MB=True, Sc=Sc))
                m = T0.set_index("episode")[["dD", "sig", "major15"]].join(Ts_.set_index("episode")[["dD", "sig"]], rsuffix="_l", how="inner")
                m = m[m.major15]
                e19 = m.loc["2019-11"] if "2019-11" in m.index else None
                lo.append(dict(config=name, Sc=Sc, without=s, frac_within=float((np.abs(m.dD_l - m.dD) < 2 * m.sig_l).mean()), n_majors=len(m),
                               composite=cs_.get("composite"), composite_sig=cs_.get("composite_sig"),
                               dD_2019_11=float(e19.dD_l) if e19 is not None else np.nan, sig_2019_11=float(e19.sig_l) if e19 is not None else np.nan))
    pd.DataFrame(lo).to_csv(f"{D.DTB}/b2_loso.csv", index=False)
    # SVI 2019-11 across configurations + far cells in it
    if D.REG == "SVI":
        e = "2019-11"
        x = TT[TT.episode == e][["config", "Sc", "dD", "sig", "n_cells"]]
        q = pd.read_csv(f"{D.DTB}/inv_q3.csv"); q = q[q.episode == e][["config", "setup", "dD", "sig", "n_cells"]]
        log(f"2019-11 b2 table:\n{x.round(4).to_string()}\ninv table:\n{q.round(4).to_string()}")
        x.to_csv(f"{D.DTB}/b2_2019_11.csv", index=False)
        fc = far[far.episode == e]
        fc[["fam", "netsta", "epoch", "n", "d", "sig", "cc", "a_s", "a_m", "a_d", "a_b", "site"]].to_csv(f"{D.DTB}/far_cells_2019_11.csv", index=False)
        w = 1 / fc.sig ** 2
        log(f"far cells in 2019-11: {len(fc)} at {sorted(fc.netsta.unique())}; raw weighted mean dv/v {np.sum(w * fc.d) / np.sum(w):+.3f} (not demeaned)")
    # same-station near vs far annual amplitude (kernel-free, free phase), using both regions' inv_cells
    rows = []
    for reg in ("svi", "wa_nw"):
        p = f"{D.DROOT}/{reg}/tables/inv_cells.csv"
        if not os.path.exists(p):
            continue
        c = pd.read_csv(p); w = 1 / c.sig ** 2; key = c.fam.astype(str) + "|" + c.netsta + "|" + c.epoch
        c["dm"] = c.d - (c.d * w).groupby(key).transform("sum") / w.groupby(key).transform("sum")
        for (ns, fr_), g in c.groupby(["netsta", "far"]):
            if len(g) < 30:
                continue
            f = C.annual_fit(pd.to_datetime(g["mid"]), g.dm.to_numpy(), g.sig.to_numpy())
            ww = 1 / g.sig ** 2
            rows.append(dict(region=reg, netsta=ns, far=fr_, n=len(g), amp=f["amp"], amp_sig=f["amp_sig"], peak=f["peak_doy"],
                             a_s=float(np.average(g.a_s, weights=ww)), a_d=float(np.average(g.a_d, weights=ww)),
                             a_b=float(np.average(g.a_b, weights=ww)), hypo=np.nan))
    A = pd.DataFrame(rows)
    both = A.groupby("netsta").far.nunique(); A["shared"] = A.netsta.map(both) == 2
    A.to_csv(f"{D.DROOT}/station_near_far_annual.csv", index=False)
    log("same-station annual amplitudes (stations used as near in one region and far in the other):\n" +
        A[A.shared].sort_values(["netsta", "far"]).round(3).to_string())


if __name__ == "__main__":
    main()
