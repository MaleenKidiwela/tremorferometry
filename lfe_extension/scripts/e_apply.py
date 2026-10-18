#!/usr/bin/env python
"""Track E step 4 — score the whole v1 catalog and produce the delivered files + the analysis.

  results/trackE/realness_v1.parquet                    key + p_real + hp + every feature
  results/trackE/lfe_extension_catalog_v1_scored.parquet  v1 columns + p_real + hp  (v1 NEVER touched)

p_real       = the calibrated, 2010-validated P(this detection is a Bostock-matchable own-family LFE)
               on the ETS-day prior of the 2012 training population.
p_real_prior = the same likelihood ratio re-based on an independently measured realness prior for the
               population the event belongs to:  inside tremor windows the frozen point's
               precision_any_ets (0.779 core / 0.694 ext); outside them Track A step 3b's measured
               inter-ETS real fraction (0.225 core / 0.057 ext).  Reported, not validated.
hp           = p_real >= 0.8 (the coordinator's definition).
"""
from __future__ import annotations
import os, json, argparse
import numpy as np
import pandas as pd
import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
RES = f"{EXT}/results"
OUT = f"{RES}/trackE"
FIG = f"{EXT}/figures"
# Track A measurements, read from the frozen point / TRACK_A_RESULTS.md step 3b
PRIOR = {("tremor", "core"): 0.779, ("tremor", "ext"): 0.694,
         ("inter", "core"): 0.225, ("inter", "ext"): 0.057}


def rebase(p, prior_train, prior_pop):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    lr = (p / (1 - p)) / (prior_train / (1 - prior_train))
    o = lr * (prior_pop / (1 - prior_pop))
    return o / (1 + o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feat", default=f"{OUT}/feat_prod.parquet")
    ap.add_argument("--model", default=f"{OUT}/model_v1.joblib")
    ap.add_argument("--catalog", default=f"{RES}/trackD/lfe_extension_catalog_v1.parquet")
    ap.add_argument("--hp", type=float, default=0.8)
    ap.add_argument("--suffix", default="", help="smoke-test suffix on every output name")
    a = ap.parse_args()
    M = joblib.load(a.model)
    X = pd.read_parquet(a.feat)
    X = X[X.ctrl == 0].reset_index(drop=True)
    raw = M["model"].predict_proba(X[M["feats"]].to_numpy())[:, 1]
    X["p_real"] = np.clip(M["iso"].predict(raw), 1e-4, 1 - 1e-4)
    # a detection whose per-site pass is missing has NO independent-site or confirmation feature:
    # the booster would happily impute it, so it is left UNSCORED instead.
    nops = X.n_sites_obs.isna().to_numpy()
    if "persite_bad" in X.columns:
        nops = nops | X.persite_bad.fillna(False).to_numpy().astype(bool)
    X.loc[nops, "p_real"] = np.nan
    raw = np.where(nops, np.nan, raw)
    X["p_raw"] = raw
    # conservative reading: the log-odds anchor measured on the held-out 2010 network (6 sites);
    # p_real is anchored on the 2012 network (9 sites).  Production epochs sit between the two, so
    # the pair brackets the honest absolute level.
    sh = float(M.get("shift_2010", 0.0))
    lo = np.log(X.p_real / (1 - X.p_real))
    X["p_real_cons"] = 1 / (1 + np.exp(-(lo + sh)))

    C = pd.read_parquet(a.catalog)
    key = ["family_id", "t"]
    S = X.rename(columns={"fam": "family_id"})
    S = S.drop_duplicates(key)
    # keep only what the catalog does NOT already carry (net_mad, rel_amp, hfr, ... are the same
    # numbers from the same Track C rows) so the merge never produces _x/_y duplicates
    drop = {"date", "doy", "epoch", "tremor_flag", "quiet_flag", "n_tremor_day", "eq_flag",
            "tier_ext", "tier_core", "cat_ext", "cat_core", "match_own", "match_any", "ctrl"}
    cols_feat = key + [c for c in S.columns
                       if c not in drop and c not in key and c not in set(C.columns)]
    Cm = C.merge(S[cols_feat], on=key, how="left", validate="one_to_one")
    n_miss = int(Cm.p_real.isna().sum())
    # context prior
    is_core = Cm.cat_core.to_numpy()
    in_tremor = Cm.tremor_flag.to_numpy()
    pp = np.where(in_tremor,
                  np.where(is_core, PRIOR[("tremor", "core")], PRIOR[("tremor", "ext")]),
                  np.where(is_core, PRIOR[("inter", "core")], PRIOR[("inter", "ext")]))
    Cm["prior_pop"] = pp
    Cm["p_real_prior"] = rebase(Cm.p_real.to_numpy(), M["prior_train"], pp)
    Cm["hp"] = Cm.p_real >= a.hp
    Cm["hp_cons"] = Cm.p_real_cons >= a.hp
    Cm["hp_prior"] = Cm.p_real_prior >= a.hp

    keep = key + ["time", "date", "month", "epoch_id", "season", "tremor_flag", "quiet_flag",
                  "cat_core", "cat_ext", "p_real", "p_raw", "p_real_cons", "p_real_prior",
                  "prior_pop", "hp", "hp_cons", "hp_prior"] + [c for c in M["feats"]] + ["n_sites_obs", "top_z", "conf_wsh",
                                                           "burst_excess", "burst_mu"]
    keep = [c for c in dict.fromkeys(keep) if c in Cm.columns]
    Cm[keep].to_parquet(f"{OUT}/realness_v1{a.suffix}.parquet", index=False)
    Cm[list(C.columns) + ["p_real", "p_real_cons", "p_real_prior", "hp"]].to_parquet(
        f"{OUT}/lfe_extension_catalog_v1_scored{a.suffix}.parquet", index=False)

    rep = dict(n_catalog=int(len(C)), n_scored=int(len(C) - n_miss), n_unscored=n_miss,
               n_unscored_persite_bad=int(Cm.persite_bad.fillna(False).sum())
               if "persite_bad" in Cm.columns else 0,
               prior_train=float(M["prior_train"]), hp_cut=a.hp,
               hp_frac_all=float(Cm.hp.mean()), hp_n_all=int(Cm.hp.sum()),
               hp_cons_frac_all=float(Cm.hp_cons.mean()), hp_cons_n_all=int(Cm.hp_cons.sum()),
               logodds_shift_2010=float(M.get("shift_2010", 0.0)))

    # ---- distributions in / out of tremor windows, by season and epoch ----
    def agg(D):
        return dict(n=int(len(D)), median_p=float(D.p_real.median()),
                    frac_hp=float((D.p_real >= a.hp).mean()),
                    frac_p05=float((D.p_real >= 0.5).mean()),
                    mean_p_prior=float(D.p_real_prior.mean()),
                    mean_p_prior_hp=float(D.p_real_prior[D.p_real >= a.hp].mean())
                    if (D.p_real >= a.hp).any() else np.nan)
    tab = []
    for tier, mt in (("core", Cm.cat_core), ("ext", Cm.cat_ext)):
        for ctx, mc in (("tremor", Cm.tremor_flag), ("inter", ~Cm.tremor_flag)):
            for ep in ["ALL"] + sorted(Cm.epoch_id.unique()):
                me = np.ones(len(Cm), bool) if ep == "ALL" else (Cm.epoch_id == ep).to_numpy()
                for se in ["ALL"] + sorted(Cm.season.unique()):
                    ms = np.ones(len(Cm), bool) if se == "ALL" else (Cm.season == se).to_numpy()
                    D = Cm[mt.to_numpy() & mc.to_numpy() & me & ms]
                    if len(D) < 20:
                        continue
                    tab.append(dict(tier=tier, ctx=ctx, epoch=ep, season=se, **agg(D)))
    T = pd.DataFrame(tab)
    T.to_csv(f"{OUT}/p_real_by_context{a.suffix}.csv", index=False)

    # ---- inter-ETS purity vs retained fraction ----
    curves = []
    for tier, mt in (("core", Cm.cat_core), ("ext", Cm.cat_ext)):
        D = Cm[mt.to_numpy() & (~Cm.tremor_flag).to_numpy()]
        if not len(D):
            continue
        for thr in np.round(np.arange(0.0, 0.96, 0.05), 2):
            m = D.p_real >= thr
            curves.append(dict(tier=tier, thr=float(thr), n=int(m.sum()),
                               retained=float(m.mean()),
                               purity=float(D.p_real_prior[m].mean()) if m.any() else np.nan))
    CV = pd.DataFrame(curves)
    CV.to_csv(f"{OUT}/interets_purity_curve{a.suffix}.csv", index=False)

    # ---- the seasonal question: DJF vs JJA hp rate, per epoch ----
    sea = []
    for tier, mt in (("core", Cm.cat_core), ("ext", Cm.cat_ext)):
        for ep in sorted(Cm.epoch_id.unique()):
            m = mt.to_numpy() & (Cm.epoch_id == ep).to_numpy() & (~Cm.tremor_flag).to_numpy()
            row = dict(tier=tier, epoch=ep)
            for se in ("DJF", "JJA"):
                d = Cm[m & (Cm.season == se).to_numpy()]
                row[f"n_{se}"] = int(len(d))
                row[f"hp_rate_{se}"] = float((d.p_real >= a.hp).mean()) if len(d) else np.nan
                row[f"median_p_{se}"] = float(d.p_real.median()) if len(d) else np.nan
            row["rate_ratio_raw"] = (row["n_DJF"] / max(row["n_JJA"], 1))
            nhpD = row["n_DJF"] * (row["hp_rate_DJF"] if np.isfinite(row["hp_rate_DJF"]) else 0)
            nhpJ = row["n_JJA"] * (row["hp_rate_JJA"] if np.isfinite(row["hp_rate_JJA"]) else 0)
            row["rate_ratio_hp"] = nhpD / max(nhpJ, 1e-9)
            sea.append(row)
    SE = pd.DataFrame(sea)
    SE.to_csv(f"{OUT}/seasonal_hp{a.suffix}.csv", index=False)

    # ---- figures ----
    CAT4 = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # validated categorical order
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for lab, m, c in (("in tremor window", Cm.tremor_flag, "#eb6834"),
                      ("inter-ETS", ~Cm.tremor_flag, "#2a78d6")):
        ax[0].hist(Cm.p_real[m & Cm.cat_core].dropna(), bins=40, histtype="step", density=True,
                   color=c, label=f"{lab} (cat_core)")
        ax[0].hist(Cm.p_real[m & Cm.cat_ext & ~Cm.cat_core].dropna(), bins=40, histtype="step", density=True,
                   color=c, ls=":", label=f"{lab} (ext only)")
    ax[0].set_xlabel("p_real"); ax[0].set_ylabel("density"); ax[0].legend(fontsize=7)
    ax[0].set_title("p_real inside vs outside tremor windows")
    for i_, ep in enumerate(sorted(Cm.epoch_id.unique())):
        d = Cm[(Cm.epoch_id == ep) & Cm.cat_core].p_real.dropna()
        if len(d):
            ax[1].hist(d, bins=40, histtype="step", density=True, color=CAT4[i_ % 4], lw=1.6,
                       label=f"{ep} (n={len(d):,})")
    ax[1].set_xlabel("p_real"); ax[1].legend(fontsize=7); ax[1].set_title("by network epoch (cat_core)")
    for i_, se in enumerate(("DJF", "MAM", "JJA", "SON")):
        d = Cm[(Cm.season == se) & Cm.cat_core & ~Cm.tremor_flag].p_real.dropna()
        if len(d):
            ax[2].hist(d, bins=40, histtype="step", density=True, color=CAT4[i_ % 4], lw=1.6,
                       label=f"{se} (n={len(d):,})")
    ax[2].set_xlabel("p_real"); ax[2].legend(fontsize=7)
    ax[2].set_title("inter-ETS, by season (cat_core)")
    for _a in ax:
        for s_ in ("top", "right"):
            _a.spines[s_].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{FIG}/e_preal_distributions{a.suffix}.png", dpi=140); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 5))
    for tier, c, dy in (("core", "#2a78d6", 8), ("ext", "#eb6834", -14)):
        d = CV[CV.tier == tier]
        ax.plot(d.retained.clip(lower=1e-4), d.purity, "o-", color=c, label=f"cat_{tier}",
                ms=5, lw=2)
        for _, r in d.iterrows():
            if r.thr in (0.5, 0.8, 0.9):
                ax.annotate(f"p≥{r.thr:g}  n={int(r.n):,}", (max(r.retained, 1e-4), r.purity),
                            fontsize=7.5, xytext=(6, dy), textcoords="offset points", color=c)
    ax.set_xscale("log")
    ax.axhline(0.8, color="#52514e", ls="--", lw=1, label="80 % purity target")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.grid(color="#e8e7e3", lw=0.5); ax.set_axisbelow(True)
    ax.set_xlabel("fraction of inter-ETS events retained")
    ax.set_ylabel("estimated purity (mean p_real_prior)")
    ax.set_title("inter-ETS: purity vs retained fraction")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(f"{FIG}/e_interets_purity{a.suffix}.png", dpi=140); plt.close(fig)

    rep["hp_interets"] = {}
    for tier, mt in (("core", Cm.cat_core), ("ext", Cm.cat_ext)):
        D = Cm[mt.to_numpy() & (~Cm.tremor_flag).to_numpy()]
        for thr in (0.5, 0.8, 0.9):
            m = D.p_real >= thr
            rep["hp_interets"][f"{tier}_ge_{thr}"] = dict(
                n=int(m.sum()), retained=float(m.mean()),
                purity=float(D.p_real_prior[m].mean()) if m.any() else np.nan)
    json.dump(rep, open(f"{OUT}/apply_report{a.suffix}.json", "w"), indent=1, default=float)
    print("[apply]", json.dumps(rep, indent=1, default=float), flush=True)
    print(SE.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
