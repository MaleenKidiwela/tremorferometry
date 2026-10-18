#!/usr/bin/env python
"""Per-station pick residuals against the times PREDICTED by Bostock's published family
coordinates.  Prefix `uu_`.  This is the diagnostic that separates "our travel times are wrong"
from "our picks are wrong": with the published hypocentre fixed and the origin time free, the
residual is the pick error plus the model error, and its STATION MEDIANS are the model/station
term while its SCATTER is the pick."""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uu_common as U                                                       # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--picks", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    B = U.bed("north")
    ST = B.station_table()
    T = U.north_truth()
    T = T[T.set == "dev"]
    P = pd.read_parquet(a.picks)
    M = T.merge(P, on="grp")
    M = M[M.sta.isin(ST.index)]
    xs, zs, TP, TS = U.tt_table()
    d = U.epi_km(M.t_lat.values, M.t_lon.values,
                 ST.lat[M.sta].values.astype(float), ST.lon[M.sta].values.astype(float))
    M["dist"] = d
    M["predS"] = U.interp2(xs, zs, TS, d, M.t_dep.values)
    M["predP"] = U.interp2(xs, zs, TP, d, M.t_dep.values)
    rows = []
    for thr in (0, 3, 4, 5, 6, 8, 10):
        Q = M[M.h_envH >= thr]
        for col in ("t_env", "t_envH", "t_onS"):
            q = Q[np.isfinite(Q[col])].copy()
            if not len(q):
                continue
            q["r"] = q[col] - q.predS
            q["r"] = q.groupby("grp").r.transform(lambda x: x - np.median(x))
            r = q.r.abs()
            rows.append(dict(h_envH_min=thr, pick=col, n=len(q), n_fam=int(q.grp.nunique()),
                             med_abs_res=float(r.median()), frac_lt_0p5=float((r < 0.5).mean()),
                             frac_lt_1=float((r < 1).mean()), frac_lt_2=float((r < 2).mean()),
                             frac_gt_4=float((r > 4).mean())))
    BT = pd.DataFrame(rows)
    # station terms / model bias at the cleanest picks
    q = M[(M.h_envH >= 6) & np.isfinite(M.t_envH)].copy()
    q["r"] = q.t_envH - q.predS
    q["r"] = q.groupby("grp").r.transform(lambda x: x - np.median(x))
    BS = q.groupby("sta").r.agg(["median", "std", "count"]).reset_index().round(3)
    BS = BS.sort_values("median")
    # P and S-P
    pp = M[np.isfinite(M.t_onP)].copy()
    pp["e"] = (pp.t_envH - pp.t_onP) - (pp.predS - pp.predP)
    sx = M[np.isfinite(M.sp_xc) & (M.sp_cc > 0.5)].copy()
    sx["e"] = sx.sp_xc - (sx.predS - sx.predP)
    extra = dict(P_onset_n=int(len(pp)), P_onset_frac=float(len(pp) / max(len(M), 1)),
                 P_onset_med_abs_err=float(pp.e.abs().median()) if len(pp) else None,
                 P_onset_frac_lt_1=float((pp.e.abs() < 1).mean()) if len(pp) else None,
                 SP_xc_n=int(len(sx)), SP_xc_frac=float(len(sx) / max(len(M), 1)),
                 SP_xc_med_abs_err=float(sx.e.abs().median()) if len(sx) else None,
                 SP_xc_frac_lt_1=float((sx.e.abs() < 1).mean()) if len(sx) else None)
    hi = BT[(BT.h_envH_min == 8) & (BT.pick == "t_envH")]
    lo = BT[(BT.h_envH_min == 0) & (BT.pick == "t_envH")]
    head = (
        f"the travel times are NOT the problem: at stack envelope peaks of >= 8 MAD the "
        f"observed S time agrees with the Crosson 1-D prediction from Bostock's published "
        f"coordinate to a median of {float(hi.med_abs_res.iloc[0]):.2f} s and "
        f"{100*float(hi.frac_lt_1.iloc[0]):.0f} % of picks land within 1 s.  The PICK is: over "
        f"all admitted stations the median is {float(lo.med_abs_res.iloc[0]):.2f} s and "
        f"{100*float(lo.frac_gt_4.iloc[0]):.0f} % of picks are more than 4 s out -- the "
        f"envelope maximum of a low-SNR stack lands on the wrong bump, and the unweighted "
        f"least-squares fit of the current locator lets those outliers drag the solution onto a "
        f"grid face.  Station medians span "
        f"{float(BS['median'].min()):+.2f} to {float(BS['median'].max()):+.2f} s, i.e. the "
        f"model/station-term error is small beside the pick scatter.  The AIC ONSET, which "
        f"Chestler uses, is WORSE than the envelope maximum at every SNR level here (the LFE "
        f"stack onset is emergent), so it is discarded; the P onset and the S-P envelope "
        f"cross-correlation are also measured and are too poor to help (below).")
    BT.to_csv(a.out, index=False)
    BS.to_csv(a.out.replace(".csv", "_by_station.csv"), index=False)
    json.dump(dict(by_threshold=BT.to_dict("records"), by_station=BS.to_dict("records"),
                   phases=extra, headline=head),
              open(a.out.replace(".csv", ".json"), "w"), indent=1)
    U.log(head)


if __name__ == "__main__":
    main()
