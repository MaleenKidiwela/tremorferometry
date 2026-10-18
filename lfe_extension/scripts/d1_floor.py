#!/usr/bin/env python
"""Track D, D1(c) -- the FALSE-ALARM FLOOR lambda_f(T, epoch, season)  (Merlin Step 6, component 3).

For every family, tier, epoch and season: the rate of catalog detections per OBSERVABLE QUIET
family-day, corrected by the measured real-inter-ETS-repeat fraction so that it is a false-alarm rate
and not merely a quiet-day detection rate.

  quiet family-day   Track C `quiet_flag`: no PNSN tremor detection inside the family box that UTC day
                     AND >= 3 days from any day that has one (PLAN amendment + Merlin Step 6).
  real fraction p    Track A step 3b (withheld-station stack test on the quiet-day detections):
                     p = 0.225 at T_core (95 % UB 0.325), p = 0.057 at T_ext (95 % UB 0.140).
                     lambda = (quiet-day rate) x (1 - p); the 95 % upper bound on p gives the LOW bracket
                     of lambda, i.e. the most generous reading of the detector.
  shrinkage          most families have few quiet observable days in a given (epoch, season) cell, so the
                     reported lambda is the gamma-Poisson posterior mean
                     (n_f + k*lambda_pool) / (D_f + k) with k = 10 template-days; the raw rate and the
                     counts are kept in the same row.

Seasons are meteorological (DJF / MAM / JJA / SON) so that the pre-registered "winter quiet-day FA within
2x summer" check (Merlin step 6) can be read straight off the table."""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_tiers  # noqa: E402  (Track C's network epochs and per-epoch tier thresholds)

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TC = f"{EXT}/results/trackC"
TD = f"{EXT}/results/trackD"
P_REAL = {"T_core": 0.225, "T_ext": 0.057}
P_REAL_UB = {"T_core": 0.325, "T_ext": 0.140}

SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
          9: "SON", 10: "SON", 11: "SON"}


EPOCHS, TIERS, TSRC = d_tiers.load()


def epoch_of(dates):
    return d_tiers.epoch_of(dates, EPOCHS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=float, default=10.0, help="shrinkage strength in template-days")
    a = ap.parse_args()
    print(f"[tiers] {TSRC}")
    fdf = sorted(glob.glob(f"{TC}/family_days_*.parquet"))
    dtf = sorted(glob.glob(f"{TC}/detections_*.parquet"))
    if not fdf:
        raise SystemExit("[d1c] no Track C months on disk yet")
    FD = pd.concat([pd.read_parquet(p, columns=["fam", "date", "observable", "quiet_flag", "n_sites"])
                    for p in fdf], ignore_index=True)
    FD = FD.drop_duplicates(subset=["fam", "date"])
    DT = pd.concat([pd.read_parquet(p, columns=["fam", "date", "quiet_flag", "cat_ext", "cat_core"])
                    for p in dtf], ignore_index=True)
    FD["epoch_id"] = epoch_of(FD.date.to_numpy().astype(str))
    FD["season"] = pd.DatetimeIndex(FD.date).month.map(SEASON)
    DT["epoch_id"] = epoch_of(DT.date.to_numpy().astype(str))
    DT["season"] = pd.DatetimeIndex(DT.date).month.map(SEASON)
    Q = FD[FD.observable & FD.quiet_flag]
    den = Q.groupby(["fam", "epoch_id", "season"]).size().rename("n_quiet_days")
    rows = []
    for tier, col in (("T_ext", "cat_ext"), ("T_core", "cat_core")):
        d = DT[DT.quiet_flag & DT[col]]
        num = d.groupby(["fam", "epoch_id", "season"]).size().rename("n_quiet_det")
        X = pd.concat([den, num], axis=1).reset_index()
        X["n_quiet_det"] = X.n_quiet_det.fillna(0.0)
        X = X[X.n_quiet_days.notna() & (X.n_quiet_days > 0)]
        X["tier"] = tier
        X["rate_raw"] = X.n_quiet_det / X.n_quiet_days
        pool = X.groupby(["epoch_id", "season"]).apply(
            lambda g: g.n_quiet_det.sum() / max(g.n_quiet_days.sum(), 1), include_groups=False).rename("rate_pool")
        X = X.merge(pool, on=["epoch_id", "season"], how="left")
        X["rate_shrunk"] = (X.n_quiet_det + a.k * X.rate_pool) / (X.n_quiet_days + a.k)
        X["p_real"] = P_REAL[tier]; X["p_real_ub"] = P_REAL_UB[tier]
        X["lambda_fa"] = X.rate_shrunk * (1 - P_REAL[tier])
        X["lambda_fa_low"] = X.rate_shrunk * (1 - P_REAL_UB[tier])
        X["lambda_fa_raw"] = X.rate_raw * (1 - P_REAL[tier])
        rows.append(X)
    L = pd.concat(rows, ignore_index=True).rename(columns={"fam": "family_id"})
    # ---- complete the grid.  Some (family, epoch, season) cells have NO quiet observable family-day at
    # all (2018-06 and 2021-10, for instance, are entirely inside tremor windows).  Rather than leave an
    # event without a floor, fall back in a fixed hierarchy and record which level was used.
    fams = sorted(set(FD.fam.unique()))
    eps = [e["epoch_id"] for e in EPOCHS]
    sea = ["DJF", "MAM", "JJA", "SON"]
    grid = pd.MultiIndex.from_product([["T_ext", "T_core"], fams, eps, sea],
                                      names=["tier", "family_id", "epoch_id", "season"]).to_frame(index=False)
    G = grid.merge(L, on=["tier", "family_id", "epoch_id", "season"], how="left")
    G["lambda_src"] = np.where(G.lambda_fa.notna(), "family-epoch-season", "")
    fe = L.groupby(["tier", "family_id", "epoch_id"]).apply(
        lambda g: (g.n_quiet_det.sum() + a.k * (g.rate_pool * g.n_quiet_days).sum()
                   / max(g.n_quiet_days.sum(), 1)) / (g.n_quiet_days.sum() + a.k),
        include_groups=False).rename("l_fe").reset_index()
    es = L.groupby(["tier", "epoch_id", "season"]).apply(
        lambda g: g.n_quiet_det.sum() / max(g.n_quiet_days.sum(), 1), include_groups=False).rename("l_es").reset_index()
    ee = L.groupby(["tier", "epoch_id"]).apply(
        lambda g: g.n_quiet_det.sum() / max(g.n_quiet_days.sum(), 1), include_groups=False).rename("l_e").reset_index()
    al = L.groupby(["tier"]).apply(lambda g: g.n_quiet_det.sum() / max(g.n_quiet_days.sum(), 1),
                                   include_groups=False).rename("l_a").reset_index()
    G = G.merge(fe, on=["tier", "family_id", "epoch_id"], how="left").merge(
        es, on=["tier", "epoch_id", "season"], how="left").merge(
        ee, on=["tier", "epoch_id"], how="left").merge(al, on=["tier"], how="left")
    pr = G.tier.map(P_REAL); prub = G.tier.map(P_REAL_UB)
    for col, lab in (("l_fe", "family-epoch"), ("l_es", "epoch-season"), ("l_e", "epoch"), ("l_a", "all")):
        need = G.lambda_fa.isna() & G[col].notna()
        G.loc[need, "lambda_fa"] = G.loc[need, col] * (1 - pr[need])
        G.loc[need, "lambda_fa_low"] = G.loc[need, col] * (1 - prub[need])
        G.loc[need, "rate_shrunk"] = G.loc[need, col]
        G.loc[need, "lambda_src"] = lab
    G["n_quiet_days"] = G.n_quiet_days.fillna(0); G["n_quiet_det"] = G.n_quiet_det.fillna(0)
    G["p_real"] = pr; G["p_real_ub"] = prub
    L = G
    L = L[["family_id", "tier", "epoch_id", "season", "n_quiet_days", "n_quiet_det", "rate_raw",
           "rate_pool", "rate_shrunk", "p_real", "p_real_ub", "lambda_fa", "lambda_fa_low",
           "lambda_fa_raw", "lambda_src"]]
    os.makedirs(TD, exist_ok=True)
    L.to_csv(f"{TD}/floor_family_epoch_season.csv", index=False)
    S = L[L.lambda_src == "family-epoch-season"].groupby(["tier", "epoch_id", "season"], as_index=False).agg(
        n_fam=("family_id", "size"), quiet_family_days=("n_quiet_days", "sum"),
        quiet_det=("n_quiet_det", "sum"))
    S["rate_pooled"] = S.quiet_det / S.quiet_family_days.clip(lower=1)
    S["lambda_pooled"] = S.rate_pooled * (1 - S.tier.map(P_REAL))
    S.to_csv(f"{TD}/floor_epoch_season_summary.csv", index=False)
    pd.set_option("display.width", 200)
    print(S.round(5).to_string(index=False))
    # the pre-registered seasonal-stability check
    for tier in ("T_ext", "T_core"):
        s = S[S.tier == tier].set_index(["epoch_id", "season"]).rate_pooled
        for eid in s.index.get_level_values(0).unique():
            v = s.loc[eid]
            if "DJF" in v.index and "JJA" in v.index and v.get("JJA", 0) > 0:
                print(f"[seasonal bar] {tier} {eid}: winter/summer quiet-day rate = "
                      f"{v['DJF'] / v['JJA']:.2f} (bar <= 2)")
    print(f"[write] {TD}/floor_family_epoch_season.csv, floor_epoch_season_summary.csv")


if __name__ == "__main__":
    main()
