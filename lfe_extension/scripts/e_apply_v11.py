#!/usr/bin/env python
"""Track E v1.1 — re-apply the FROZEN classifier (no retraining) after Track C's cc-floor fix.

Covers EVERY current Track C detection at or above the record threshold (8.0 MAD), not just the
catalog members, so Track D can attach the score to whatever v1.1 membership it settles on.

MERGE KEY (documented, three equivalent forms, all carried as columns):
  * `month` + `det_row`   — `det_row` is the 0-BASED ROW INDEX into
                            `results/trackC/detections_<month>.parquet` exactly as stored on disk
                            (the file is written sorted by `t`; the index is taken after
                            `pd.read_parquet` with no re-sorting and BEFORE any net_mad filter).
  * `family_id` + `t`     — the semantic key: Bostock family (`fam`) and the detection time in
                            float64 epoch seconds, copied bit-for-bit from the Track C parquet.
                            Unique over the whole record (0 duplicates).
  * `det_uid` = f"{month}:{det_row}" — a single string key for convenience.

The model, the isotonic calibration and the 2010 log-odds anchor are loaded unchanged from
results/trackE/model_v1.joblib.  Nothing is refitted.
"""
from __future__ import annotations
import os, glob, json, argparse
import numpy as np
import pandas as pd
import joblib

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
RES = f"{EXT}/results"
OUT = f"{RES}/trackE"
PRIOR = {("tremor", "core"): 0.779, ("tremor", "ext"): 0.694,
         ("inter", "core"): 0.225, ("inter", "ext"): 0.057}
EPOCHS = [("E0", "2014-01", "2017-02"), ("E1", "2017-03", "2018-12"),
          ("E2", "2019-01", "2021-03"), ("E3", "2021-04", "2099-12")]
SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
          6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}


def rebase(p, prior_train, prior_pop):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    lr = (p / (1 - p)) / (prior_train / (1 - prior_train))
    o = lr * (prior_pop / (1 - prior_pop))
    return o / (1 + o)


def epoch_of(ym):
    for n, a, b in EPOCHS:
        if a <= ym <= b:
            return n
    return "E?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feat", default=f"{OUT}/feat_prod_v1_1.parquet")
    ap.add_argument("--model", default=f"{OUT}/model_v1.joblib")
    ap.add_argument("--hp", type=float, default=0.8)
    a = ap.parse_args()
    M = joblib.load(a.model)

    # ---- the key: row index in each month's Track C parquet, as stored ----
    keys = []
    for p in sorted(glob.glob(f"{RES}/trackC/detections_*.parquet")):
        ym = os.path.basename(p).replace("detections_", "").replace(".parquet", "")
        C = pd.read_parquet(p, columns=["fam", "t", "net_mad", "date", "tier_ext", "tier_core",
                                        "cat_ext", "cat_core", "tremor_flag", "quiet_flag",
                                        "hf_pass", "T_core_used", "T_ext_used", "tier_epoch"])
        C = C.reset_index().rename(columns={"index": "det_row"})
        C["month"] = ym
        keys.append(C[C.net_mad >= 8.0])
    K = pd.concat(keys, ignore_index=True)
    K["family_id"] = K.fam.astype(int)
    K["det_uid"] = K.month + ":" + K.det_row.astype(str)
    assert K.duplicated(["family_id", "t"]).sum() == 0, "(family_id, t) is not unique"

    X = pd.read_parquet(a.feat)
    X["family_id"] = X.fam.astype(int)
    raw = M["model"].predict_proba(X[M["feats"]].to_numpy())[:, 1]
    p = np.clip(M["iso"].predict(raw), 1e-4, 1 - 1e-4)
    bad = X.persite_bad.fillna(False).to_numpy().astype(bool) | X.n_sites_obs.isna().to_numpy()
    X["p_real"] = np.where(bad, np.nan, p)
    X["p_raw"] = np.where(bad, np.nan, raw)
    lo = np.log(X.p_real / (1 - X.p_real))
    X["p_real_cons"] = 1 / (1 + np.exp(-(lo + float(M.get("shift_2010", 0.0)))))

    featcols = ["family_id", "t", "p_real", "p_raw", "p_real_cons", "persite_bad"] + \
        [c for c in M["feats"]] + ["n_sites_obs", "top_z", "conf_wsh", "burst_excess", "burst_mu"]
    # `net_mad` is carried by the Track C key frame; taking it from BOTH sides would collide into
    # net_mad_x / net_mad_y and silently disappear from the delivered columns.
    featcols = [c for c in dict.fromkeys(featcols)
                if c in X.columns and (c in ("family_id", "t") or c not in K.columns)]
    R = K.merge(X[featcols], on=["family_id", "t"], how="left", validate="one_to_one")
    assert not any(c.endswith(("_x", "_y")) for c in R.columns), "merge produced suffixed columns"
    n_unscored = int(R.p_real.isna().sum())

    R["epoch_id"] = [epoch_of(m) for m in R.month]
    R["season"] = [SEASON[int(m[5:7])] for m in R.month]
    pp = np.where(R.tremor_flag.to_numpy(),
                  np.where(R.cat_core.to_numpy(), PRIOR[("tremor", "core")], PRIOR[("tremor", "ext")]),
                  np.where(R.cat_core.to_numpy(), PRIOR[("inter", "core")], PRIOR[("inter", "ext")]))
    R["prior_pop"] = pp
    R["p_real_prior"] = rebase(R.p_real.to_numpy(), M["prior_train"], pp)
    R["hp"] = R.p_real >= a.hp
    R["hp_cons"] = R.p_real_cons >= a.hp
    cols = ["det_uid", "month", "det_row", "family_id", "t", "date", "epoch_id", "season",
            "net_mad", "T_core_used", "T_ext_used", "tier_epoch", "tier_core", "tier_ext",
            "cat_core", "cat_ext", "hf_pass", "tremor_flag", "quiet_flag",
            "p_real", "p_raw", "p_real_cons", "p_real_prior", "prior_pop", "hp", "hp_cons",
            "persite_bad"] + [c for c in featcols if c not in ("family_id", "t", "p_real", "p_raw",
                                                               "p_real_cons", "persite_bad")]
    cols = [c for c in dict.fromkeys(cols) if c in R.columns]
    R[cols].to_parquet(f"{OUT}/realness_v1_1.parquet", index=False)

    # ---- inter-ETS purity / retention and the seasonal hp ratios under the NEW tiers ----
    cur = []
    for tier, mt in (("core", R.cat_core), ("ext", R.cat_ext)):
        D = R[mt.to_numpy() & (~R.tremor_flag).to_numpy() & R.p_real.notna().to_numpy()]
        for thr in np.round(np.arange(0.0, 0.96, 0.05), 2):
            m = D.p_real >= thr
            cur.append(dict(tier=tier, thr=float(thr), n=int(m.sum()),
                            retained=float(m.mean()) if len(D) else np.nan,
                            purity=float(D.p_real_prior[m].mean()) if m.any() else np.nan))
    CV = pd.DataFrame(cur); CV.to_csv(f"{OUT}/interets_purity_curve_v1_1.csv", index=False)
    sea = []
    for tier, mt in (("core", R.cat_core), ("ext", R.cat_ext)):
        for ep in ["E0", "E1", "E2", "E3"]:
            m = mt.to_numpy() & (R.epoch_id == ep).to_numpy() & (~R.tremor_flag).to_numpy()
            r = dict(tier=tier, epoch=ep,
                     T_core=float(R[R.epoch_id == ep].T_core_used.iloc[0]),
                     T_ext=float(R[R.epoch_id == ep].T_ext_used.iloc[0]))
            for se in ("DJF", "JJA"):
                d = R[m & (R.season == se).to_numpy()]
                r[f"n_{se}"] = int(len(d)); r[f"hp_{se}"] = int((d.p_real >= a.hp).sum())
                r[f"hp_rate_{se}"] = float((d.p_real >= a.hp).mean()) if len(d) else np.nan
            r["ratio_raw"] = r["n_DJF"] / max(r["n_JJA"], 1)
            r["ratio_hp"] = r["hp_DJF"] / max(r["hp_JJA"], 1e-9)
            sea.append(r)
    SE = pd.DataFrame(sea); SE.to_csv(f"{OUT}/seasonal_hp_v1_1.csv", index=False)

    rep = dict(n_rows=int(len(R)), n_scored=int(len(R) - n_unscored), n_unscored=n_unscored,
               n_persite_bad=int(R.persite_bad.fillna(False).sum()),
               n_cat_core=int(R.cat_core.sum()), n_cat_ext=int(R.cat_ext.sum()),
               hp_n=int(R.hp.sum()), hp_frac=float(R.hp.mean()),
               hp_n_cat_ext=int((R.hp & R.cat_ext).sum()),
               hp_cons_n=int(R.hp_cons.sum()),
               tiers={e: [float(R[R.epoch_id == e].T_core_used.iloc[0]),
                          float(R[R.epoch_id == e].T_ext_used.iloc[0])] for e in ["E0","E1","E2","E3"]},
               model="results/trackE/model_v1.joblib (frozen, NOT refitted)",
               prior_train=float(M["prior_train"]), logodds_shift_2010=float(M.get("shift_2010", 0.0)),
               hp_interets={})
    for tier, mt in (("core", R.cat_core), ("ext", R.cat_ext)):
        D = R[mt.to_numpy() & (~R.tremor_flag).to_numpy() & R.p_real.notna().to_numpy()]
        rep["hp_interets"][tier + "_n_interets"] = int(len(D))
        for thr in (0.5, 0.8, 0.9):
            m = D.p_real >= thr
            rep["hp_interets"][f"{tier}_ge_{thr}"] = dict(
                n=int(m.sum()), retained=float(m.mean()),
                purity=float(D.p_real_prior[m].mean()) if m.any() else np.nan)
    json.dump(rep, open(f"{OUT}/apply_report_v1_1.json", "w"), indent=1, default=float)
    print(json.dumps(rep, indent=1, default=float))
    print(SE.to_string(index=False))


if __name__ == "__main__":
    main()
