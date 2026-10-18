#!/usr/bin/env python
"""Track E — does the winter excess collapse under p_real, or does p_real itself lean on winter?

The full score's inter-ETS DJF/JJA ratio goes UP in all four network epochs, which is the signature
of a score-level seasonal bias, not of a change in the real fraction.  The per-feature seasonal
medians say where it comes from: `n_sites_obs` (+1 site in DJF) and `rel_amp` (winter is ~1.4x
louder and rel_amp is referenced to a FIXED 2012 quiet-day rms), while every feature that is
normalised by the day's OWN distribution (z_site_med, conf_z, z_site_min, burst_z) is LOWER in
winter — the same direction as Track D's withheld-station moveout-shuffled test.

This script refits the identical model on a SEASON-ROBUST feature subset (only quantities
normalised by that family-day's own background) and repeats the comparison.
"""
from __future__ import annotations
import os, json
import numpy as np
import pandas as pd
import joblib
from sklearn.model_selection import GroupKFold
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e_train import mk_model, cv_scores, reliability

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
RES = f"{EXT}/results"
OUT = f"{RES}/trackE"

ROBUST = ["net_mad", "net_cc", "frac_sites_above", "z_site_med", "z_site_min", "z_ex_top",
          "conf_z", "conf_ze", "ze_site_med", "n_same10", "n_same60", "n_other5", "burst_z"]


def main():
    M = joblib.load(f"{OUT}/model_v1.joblib")
    A = pd.read_parquet(f"{OUT}/feat_ref2012.parquet")
    B = pd.read_parquet(f"{OUT}/feat_ref2010.parquet")
    for D in (A, B):
        D["match_own"] = D.match_own.fillna(False).astype(bool)
    e12 = (int(A[(A.ctrl == 0) & A.match_own].doy.min()), int(A[(A.ctrl == 0) & A.match_own].doy.max()))
    e10 = (int(B[(B.ctrl == 0) & B.match_own].doy.min()), int(B[(B.ctrl == 0) & B.match_own].doy.max()))
    A12 = A[(A.ctrl == 0) & A.doy.between(*e12)]
    Q12 = A[(A.ctrl == 0) & (A.doy < e12[0])]
    B10 = B[(B.ctrl == 0) & B.doy.between(*e10)].reset_index(drop=True)
    TR = pd.concat([A12.assign(src="ets"), Q12.assign(src="quiet")], ignore_index=True)
    y = TR.match_own.astype(int).to_numpy(); g = TR.date.to_numpy()
    m_ets = (TR.src == "ets").to_numpy()
    oof = cv_scores(M["kind"], ROBUST, TR, y, g)
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0).fit(oof[m_ets], y[m_ets])
    mod = mk_model(M["kind"], ROBUST); mod.fit(TR[ROBUST].to_numpy(), y)
    y10 = B10.match_own.to_numpy().astype(int)
    p10 = np.clip(iso.predict(mod.predict_proba(B10[ROBUST].to_numpy())[:, 1]), 1e-4, 1 - 1e-4)
    rep = dict(features=ROBUST,
               cv2012_auc_ets=float(roc_auc_score(y[m_ets], np.clip(iso.predict(oof[m_ets]), 0, 1))),
               valid2010_auc=float(roc_auc_score(y10, p10)),
               valid2010_prec_ge_08=float(y10[p10 >= 0.8].mean()) if (p10 >= 0.8).any() else np.nan,
               valid2010_n_ge_08=int((p10 >= 0.8).sum()))

    X = pd.read_parquet(f"{OUT}/feat_prod.parquet")
    R = pd.read_parquet(f"{OUT}/realness_v1.parquet")
    pr = np.clip(iso.predict(mod.predict_proba(X[ROBUST].to_numpy())[:, 1]), 1e-4, 1 - 1e-4)
    X["p_rob"] = np.where(X.persite_bad.fillna(False).to_numpy().astype(bool), np.nan, pr)
    S = X[["fam", "t", "p_rob"]].rename(columns={"fam": "family_id"}).drop_duplicates(["family_id", "t"])
    R = R.merge(S, on=["family_id", "t"], how="left")
    R[["family_id", "t", "p_rob"]].to_parquet(f"{OUT}/p_real_robust.parquet", index=False)

    rows = []
    for tier, mt in (("core", R.cat_core), ("ext", R.cat_ext)):
        for ep in sorted(R.epoch_id.unique()):
            m = mt.to_numpy() & (R.epoch_id == ep).to_numpy() & (~R.tremor_flag).to_numpy()
            r = dict(tier=tier, epoch=ep)
            for se in ("DJF", "JJA"):
                d = R[m & (R.season == se).to_numpy()]
                r[f"n_{se}"] = int(len(d))
                r[f"hp_full_{se}"] = int((d.p_real >= 0.8).sum())
                r[f"hp_rob_{se}"] = int((d.p_rob >= 0.8).sum())
            r["ratio_raw"] = r["n_DJF"] / max(r["n_JJA"], 1)
            r["ratio_hp_full"] = r["hp_full_DJF"] / max(r["hp_full_JJA"], 1e-9)
            r["ratio_hp_robust"] = r["hp_rob_DJF"] / max(r["hp_rob_JJA"], 1e-9)
            rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(f"{OUT}/seasonal_hp_robust.csv", index=False)
    rep["table"] = T.to_dict("records")

    # per-feature seasonal medians on the inter-ETS core population (the diagnosis itself)
    D = R[(~R.tremor_flag) & R.cat_core & R.p_real.notna()]
    fr = []
    for f in M["feats"]:
        a = D[D.season == "DJF"][f].astype(float); b = D[D.season == "JJA"][f].astype(float)
        sd = max(D[f].astype(float).std(), 1e-9)
        fr.append(dict(feature=f, med_DJF=a.median(), med_JJA=b.median(),
                       shift_sd=(a.median() - b.median()) / sd,
                       in_robust_set=f in ROBUST))
    FS = pd.DataFrame(fr).sort_values("shift_sd", key=abs, ascending=False)
    FS.to_csv(f"{OUT}/seasonal_feature_shift.csv", index=False)
    json.dump(rep, open(f"{OUT}/season_check.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in rep.items() if k != "table"}, indent=1))
    print(T.to_string(index=False))
    print(FS.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
