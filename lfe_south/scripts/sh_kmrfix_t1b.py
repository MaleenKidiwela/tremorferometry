#!/usr/bin/env python
"""SH_KMRFIX task 1b: approximate effect of removing NC.KMR channels from the network sum on 2025-07-25 .. 2025-08-25
(segment A), from saved values only.  For each >= 8x candidate with a cut window, the KMR per-channel cc at lag 0
(t1_kmr_candidate_cc.parquet, reconstructed from the raw cut windows; channel-mean reproduces the stored cc_NC.KMR at
r = 0.97) is subtracted from the stored summed cc.  The daily MAD of the network sum is reduced by that channel's noise
variance (Gaussian: sigma = MAD / 0.6745; channel noise sigma = robust spread of the KMR channel cc at pre-arrival lags
of the same windows), the daily median is kept.  Candidates without a window (and 6-8x candidates, which cannot reach
9x by this change) keep their stored values.  Variants: drop_h (live horizontals only) and drop_all (all KMR).
The S1 rules (sh_v3_common.apply_rules, k=9, 4-s declustering over all families) are re-applied to all mf_A
candidates of the period with the corrected values.  Output results/sh_kmrfix/t1b_*."""
import os, sys, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import sh_checks_kmr as K
from sh_v3_common import apply_rules
RES = K.RES; EXT = K.EXT; OUT = f"{RES}/sh_kmrfix"
P0, P1 = pd.Timestamp("2025-07-25"), pd.Timestamp("2025-08-26")
CH = ["HHE", "HHN", "HHZ"]


def main():
    D = pd.read_parquet(f"{OUT}/t1_kmr_candidate_cc.parquet")
    D = D[(D.t >= pd.Timestamp("2025-06-01")) & (D.t < P1)].copy()
    # channel noise sigma from pre-arrival lags of the same windows (recompute quickly from the saved per-candidate table is
    # not possible -> use the robust spread of single-event cc over the 'healthy' period minus nothing: we re-read windows)
    # robust noise estimate: cc of the template against pre-arrival slices (index 0..60 starts) per period and channel
    W_ok = D.ok.values
    print("[t1b] candidates", len(D), "with windows", int(W_ok.sum()))
    # noise sigma per channel and period, from the cc distribution of candidates' channel cc around their median is
    # signal-contaminated; use instead the per-family-day KMR station noise (results/sh_conf/feat/fwd/<date>_noise.parquet:
    # channel-MEAN cc at 1000 random origins, noise_mad = 1.4826 MAD).  Per channel sigma ~ noise_mad * sqrt(n_live_ch).
    rows = []
    for d in pd.date_range("2025-06-01", "2025-08-25"):
        p = f"{RES}/sh_conf/feat/fwd/{d.date()}_noise.parquet"
        if not os.path.exists(p): continue
        N = pd.read_parquet(p); N = N[N.netsta == "NC.KMR"]
        rows.append(N[["date", "fam", "noise_med", "noise_mad"]])
    NZ = pd.concat(rows); NZ["fam"] = NZ.fam.astype(str)
    G = pd.read_csv(f"{EXT}/mf_A_diag.csv"); G["fam"] = G.fam.astype(str); G["date"] = G.day.astype(str)
    A = pd.read_parquet(f"{EXT}/mf_A.parquet"); A["fam"] = A.fam.astype(str)
    A = A[(A.t >= P0) & (A.t < P1)].reset_index(drop=True); A["date"] = A.t.dt.strftime("%Y-%m-%d")
    A = A.merge(G[["fam", "date", "mad"]], on=["fam", "date"], how="left").merge(NZ, on=["fam", "date"], how="left")
    A["key"] = np.round(A.t_epoch * 20).astype(np.int64)
    X = D[["fam", "key", "cc_HHE", "cc_HHN", "cc_HHZ", "ok"]]
    A = A.merge(X, on=["fam", "key"], how="left")
    A["med"] = A.ccsum - A.net_mad * A.mad
    nlive_k = A[["cc_HHE", "cc_HHN", "cc_HHZ"]].notna().sum(1)
    out = {}
    for var, chans in [("orig", []), ("drop_h", ["HHE", "HHN"]), ("drop_all", CH)]:
        B = A.copy()
        if chans:
            has = B.ok.fillna(False).astype(bool)
            sub = B[[f"cc_{c}" for c in chans]]
            nrem = sub.notna().sum(1)
            # live KMR channels of the day (HHE masked dead from 2025-08-03): per-channel noise sigma = noise_mad*sqrt(nlive_k)
            sig_ch = B.noise_mad * np.sqrt(nlive_k.clip(lower=1))
            sigma = B.mad / 0.6745
            sig2 = np.sqrt(np.clip(sigma ** 2 - nrem * sig_ch ** 2, 1e-12, None))
            mad2 = 0.6745 * sig2
            cs2 = B.ccsum - sub.fillna(0).sum(1)
            med2 = B.med - nrem * B.noise_med.fillna(0)
            B.loc[has, "ccsum"] = cs2[has]
            B.loc[has, "nlive"] = (B.nlive - nrem)[has]
            B.loc[has, "meancc"] = (cs2 / (B.nlive - nrem).clip(lower=1))[has]
            B.loc[has, "net_mad"] = ((cs2 - med2) / mad2)[has]
        E = apply_rules(B[["t", "fam", "ccsum", "meancc", "net_mad", "nlive"]].reset_index(drop=True), 9)
        out[var] = E
    per = lambda t: np.where(t < pd.Timestamp("2025-08-03"), "drift 07-25..08-02", "HHE stuck 08-03..08-25")
    R = []
    for var, E in out.items():
        for p, g in pd.Series(per(E.t), index=E.index).groupby(per(E.t)).groups.items():
            R.append(dict(variant=var, period=p, n_S1=len(g)))
        R.append(dict(variant=var, period="all 07-25..08-25", n_S1=len(E)))
    R = pd.DataFrame(R).pivot(index="period", columns="variant", values="n_S1")
    R["drop_h_vs_orig"] = R.drop_h / R.orig; R["drop_all_vs_orig"] = R.drop_all / R.orig
    # per family
    F = pd.DataFrame({v: E.groupby("fam").size() for v, E in out.items()}).fillna(0)
    F["drop_h_vs_orig"] = F.drop_h / F.orig; F["drop_all_vs_orig"] = F.drop_all / F.orig
    R.to_csv(f"{OUT}/t1b_effect_S1.csv"); F.to_csv(f"{OUT}/t1b_effect_S1_by_family.csv")
    # orig must reproduce the stored catalogue for the period
    CAT = pd.read_parquet(f"{EXT}/cat_S1_9x_internal.parquet"); CAT = CAT[(CAT.t >= P0) & (CAT.t < P1)]
    print(R.round(3).to_string()); print("stored catalogue in period:", len(CAT), "(orig recompute", len(out["orig"]), ")")
    print(F.describe().round(3).to_string())
    print("mad shrink median drop_h:", float(np.nanmedian(0.6745 * np.sqrt(np.clip((A.mad / 0.6745) ** 2 - 1 * (A.noise_mad * np.sqrt(nlive_k.clip(lower=1))) ** 2, 1e-12, None)) / A.mad)))


if __name__ == "__main__":
    main()
