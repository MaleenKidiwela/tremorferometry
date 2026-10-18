#!/usr/bin/env python
"""S11 tomo step 3 (PLAN.md Amendment 6): SVI P/S picks on stacked family templates.

SVI-1 (set S1): e5r 2012 stacks, 116 families x 15 stations; frame t = ttS - 20 s + i/FS (rel. catalogue origin)
SVI-2 (set S2): S10 pilot grand stacks, 3 families x <=21 stations; frame offset_cut + tau
S = frozen Royer picker (largest horizontal extremum of Hilbert-rotated 3C in 10-90 % of the window),
P = V3 polarisation-weighted picker (headline), v0 frozen P as column.  Output tables/picks_S1V3.csv, picks_S2V3.csv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import resample_poly

ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, f"{ROOT}/lfe_array_discovery/scripts")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cwd = os.getcwd(); os.chdir(f"{ROOT}/lfe_array_discovery")
import a4f_hypo as H  # noqa: E402
os.chdir(cwd)
from s11_tomo_repick import pick_p_v3  # noqa: E402

T = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"
FS = 40.0


def picks_for(ST, meta):
    """ST {fam: {sta: (3,n) at 40 Hz}}; returns frozen-picker rows + V3 P."""
    P = H.make_picks(ST)
    P["tP_v0"], P["rP_v0"] = P.tP, P.rP
    tp, rp, wp = [], [], []
    for r in P.itertuples():
        if not np.isfinite(r.tS):
            tp.append(np.nan); rp.append(np.nan); wp.append(np.nan); continue
        x, y, w = pick_p_v3(ST[r.grp][r.sta], r.tS)
        tp.append(x); rp.append(y); wp.append(w)
    P["tP"], P["rP"], P["wP"] = tp, rp, wp
    return P


def main():
    # ---------------- SVI-1
    z = np.load(f"{ROOT}/lfe_catalog_repair/results/e5r_templates_2012.npz", allow_pickle=True)
    Q = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_template_qc_2012.csv")
    Q = Q.set_index(["bkey", "sta"])
    Fm = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv").set_index("bkey")
    ST, off = {}, {}
    for k in z.files:
        if "|" not in k:
            continue
        sta, b = k.split("|"); b = int(b)
        x = resample_poly(z[k].astype(np.float64), 2, 1, axis=-1)
        ST.setdefault(b, {})[sta] = x
        q = Q.loc[(b, sta)]
        off[(b, sta)] = float(q.ttS) - 20.0 - float(q.align_shift_s)
    P = picks_for(ST, None)
    o = np.array([off[(g, s)] for g, s in zip(P.grp, P.sta)])
    for c in ("tP", "tS", "tP_v0"):
        P[c] = P[c] + o
    P["frame_offset"] = o
    P["fam"] = P.grp.astype(str)
    P["lat"], P["lon"], P["depth_km"] = P.grp.map(Fm.lat), P.grp.map(Fm.lon), P.grp.map(Fm.dep)
    P["n_detections"] = P.grp.map(Fm.n_all)
    P["epi_km"] = [Q.loc[(g, s), "epi_km"] for g, s in zip(P.grp, P.sta)]
    P["ttS_model"] = [Q.loc[(g, s), "ttS"] for g, s in zip(P.grp, P.sta)]
    P["ttP_model"] = [Q.loc[(g, s), "ttP"] for g, s in zip(P.grp, P.sta)]
    P["snr_db"] = [Q.loc[(g, s), "snr_db"] for g, s in zip(P.grp, P.sta)]
    P["src"] = "e5r_templates_2012"
    P.to_csv(f"{T}/picks_S1V3.csv", index=False)
    ok = (P.rP >= 1.2) & (P.rS >= 1.2)
    print(f"[S1V3] {P.fam.nunique()} families, {len(P)} rows, rP&rS>=1.2 {ok.sum()} "
          f"(v0 {((P.rP_v0 >= 1.2) & (P.rS >= 1.2)).sum()}); families with >=6: "
          f"{(P[ok].groupby('fam').size() >= 6).sum()}", flush=True)
    d = (P.tS - P.ttS_model)[P.rS >= 1.2]
    print(f"[S1V3] S pick minus model ttS (rS>=1.2): median {d.median():.2f} s, IQR {d.quantile(.25):.2f}..{d.quantile(.75):.2f}")
    # ---------------- SVI-2
    g = np.load(f"{ROOT}/svi_science/results/s10_dvv/pilot/tables/grandstacks.npz", allow_pickle=True)
    O = pd.read_csv(f"{ROOT}/svi_science/results/s10_dvv/pilot/tables/offsets_sp.csv")
    O = O.set_index(["fam", "netsta"])
    keys = [str(k) for k in g["keys"]]
    K = pd.DataFrame([k.split("|") for k in keys], columns=["fam", "netsta", "epoch"])
    K["i"], K["n"] = range(len(K)), g["n"]
    K = K.sort_values("n", ascending=False).drop_duplicates(["fam", "netsta"])
    tau0 = float(g["tau"][0])
    ST2, off2 = {}, {}
    for r in K.itertuples():
        ST2.setdefault(int(r.fam), {})[r.netsta] = g["G1"][r.i].astype(np.float64)
        off2[(int(r.fam), r.netsta)] = float(O.loc[(int(r.fam), r.netsta), "offset"]) + tau0   # s10_windows.py:205 cuts at "offset"
    P2 = picks_for(ST2, None)
    o = np.array([off2[(f, s)] for f, s in zip(P2.grp, P2.sta)])
    for c in ("tP", "tS", "tP_v0"):
        P2[c] = P2[c] + o
    P2["frame_offset"] = o
    P2["fam"] = P2.grp.astype(str)
    P2["epoch"] = [K.set_index(["fam", "netsta"]).loc[(str(f), s), "epoch"] for f, s in zip(P2.grp, P2.sta)]
    P2["n_rows"] = [K.set_index(["fam", "netsta"]).loc[(str(f), s), "n"] for f, s in zip(P2.grp, P2.sta)]
    P2["src"] = "s10 pilot grandstacks"
    P2["lat"], P2["lon"], P2["depth_km"], P2["n_detections"] = np.nan, np.nan, np.nan, P2.n_rows
    P2.to_csv(f"{T}/picks_S2V3.csv", index=False)
    ok = (P2.rP >= 1.2) & (P2.rS >= 1.2)
    print(f"[S2V3] {P2.fam.nunique()} families, {len(P2)} rows, rP&rS>=1.2 {ok.sum()}; per family "
          f"{P2[ok].groupby('fam').size().to_dict()}", flush=True)


if __name__ == "__main__":
    main()
