#!/usr/bin/env python
"""s14_veto_impact -- impact of the 2026-10-14 instrument vetoes (svi_science/instrument_vetoes.csv: B013 2012-07-13..
2013-07-24, B009 2007-01-08..2009-10-01) on the cheap S11 headline numbers: Q3 (per-major Delta D, 15/32-major composite)
and Q1/Q2 (v2 smooth model). Linear refits of the existing cell tables only (no waveform work, no gates, no injections).
Also a sensitivity column for WA: B013 series split at 2013-07-24 (lasting EH1 change after the outage; kept, not vetoed).

    S11_REGION=SVI   python s14_veto_impact.py
    S11_REGION=WA_NW python s14_veto_impact.py
-> svi_science/results/s14_cc/tables/veto_impact_<region>.csv / _q3_<region>.csv
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I
import s11_q3 as Q3
import s11_v2 as V2

OUT = f"{C.ROOT}/svi_science/results/s14_cc/tables"
B013_SPLIT = pd.Timestamp("2013-07-24")


def q3_cfg():
    if C.REGION == "SVI":
        if not os.path.exists(f"{C.TB}/episode_raw_all.parquet"):
            Q3.combine()
        return "_all", 20, "c_115fam_n20"
    return "", 20, "b_n20"


def split_b013(cells):
    c = cells.copy()
    mid = c.episode.map(C.episodes().mid)
    m = (c.netsta == "PB.B013") & (mid >= B013_SPLIT)
    c.loc[m, "epoch"] = c.loc[m, "epoch"] + "_post1307"
    return c


def run(veto, split=False):
    I.APPLY_VETO = veto
    tag, nmin, name = q3_cfg()
    cq = Q3.cells_for(tag, nmin)
    if split:
        cq = split_b013(cq)
    F = Q3.fit(cq); T, comp = Q3.q3(F)
    cv = V2.cells_v2()
    if split:
        cv = split_b013(cv)
    q = V2.Q(V2.fit(cv, "smooth"))
    mj = T[T.major15]
    r = dict(region=C.REGION, veto=veto, b013_split=split, q3_config=name, q3_cells=F["n"], q3_composite=comp["composite"],
             q3_composite_sig=comp["composite_sig"], q3_n_majors=comp["n_majors"], q3_sig_major_median=float(mj.sig.median()),
             q3_n_major_sig_le_0p1=int((mj.sig <= 0.10).sum()), q3_chi2_majors=float(np.mean((mj.dD / mj.sig) ** 2)),
             q3_n_beyond_2sig=int((np.abs(mj.dD / mj.sig) > 2).sum()), v2_cells=q["n_cells"], Q1=q["Q1"], Q1_sig=q["Q1_sig"],
             Q1_peak=q["Q1_peak"], Q2=q["Q2"], Q2_sig=q["Q2_sig"])
    return r, T.assign(veto=veto, b013_split=split), cq, cv


def main():
    rows, TT = [], []
    r0, T0, cq0, cv0 = run(False)
    r1, T1, cq1, cv1 = run(True)
    for r, cq, cv in ((r0, cq0, cv0), (r1, cq1, cv1)):
        r["q3_cells_B009_B013_in_veto_windows"] = int(I.veto_mask(cq.netsta, cq.episode.map(C.episodes().mid)).sum())
    w = 1 / cq0.sig ** 2
    vm = I.veto_mask(cq0.netsta, cq0.episode.map(C.episodes().mid))
    r1["q3_weight_share_removed"] = float(w[vm].sum() / w.sum())
    r1["q3_episodes_touched"] = ";".join(sorted(cq0.episode[vm].unique()))
    rows += [r0, r1]; TT += [T0, T1]
    if C.REGION != "SVI":
        r2, T2, _, _ = run(True, split=True); rows.append(r2); TT.append(T2)
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/veto_impact_{C.REGION}.csv", index=False)
    pd.concat(TT, ignore_index=True).to_csv(f"{OUT}/veto_impact_q3_{C.REGION}.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.T.to_string())
    a = T0.set_index("episode"); b = T1.set_index("episode")
    j = a[["dD", "sig", "n_cells", "major15"]].join(b[["dD", "sig", "n_cells"]], rsuffix="_veto")
    j["shift_over_sig"] = (j.dD_veto - j.dD) / j.sig
    print(j[(j.n_cells != j.n_cells_veto) | (j.shift_over_sig.abs() > 0.3)].round(3).to_string())


if __name__ == "__main__":
    main()
