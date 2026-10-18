#!/usr/bin/env python
"""g_trigger_build (2026-10-11) -- Study G step 3: event / candidate / observable-day / episode tables for SVI and
WA-NW (PLAN.md §0).  Reads only; writes results/g_trigger/tables/{svi,wa}_*.parquet.  CPU <= 8, < 30 GB.

    python svi_science/scripts/g_trigger_build.py
"""
import os, glob, json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
G = f"{ROOT}/svi_science/results/g_trigger/tables"
FEAT = f"{ROOT}/lfe_extension/results/refine/feat/prod"
WP = f"{ROOT}/lfe_extension_wa/results/wa_prod"


def seg_of(d):
    d = pd.Series(pd.to_datetime(d))
    s = np.full(len(d), "FWD", dtype=object)
    for a, b, n in (("2003-01-01", "2003-01-31", "CN4f"), ("2003-02-01", "2006-06-30", "ERA"),
                    ("2006-07-01", "2007-11-12", "B2f"), ("2007-11-13", "2009-08-26", "B3f"),
                    ("2009-08-27", "2013-12-31", "S15f")):
        s[((d >= a) & (d <= b)).to_numpy()] = n
    return s


def wa_era(d):
    y = pd.to_datetime(pd.Series(d)).dt.year.to_numpy()
    return np.where(y <= 2011, "E1", np.where(y <= 2014, "E2", "E3"))


def _fz(date):
    p = f"{FEAT}/{date}.parquet"
    if not os.path.exists(p):
        return None
    return pd.read_parquet(p, columns=["series", "fam", "t", "n_z2"]).assign(date=date)


def svi():
    E = pd.read_parquet(f"{ROOT}/svi_science/results/r23/tables/events_2003_2026.parquet")
    fw = E.src == "fw"
    dates = sorted(E.loc[fw, "date"].unique())
    with ProcessPoolExecutor(8) as ex:
        Z = pd.concat([z for z in ex.map(_fz, dates, chunksize=20) if z is not None], ignore_index=True)
    Z["tkey"] = np.round(Z.t, 2)
    Z = Z.drop_duplicates(["date", "series", "fam", "tkey"])
    F = E[fw].copy()
    F["tkey"] = np.round(F.time_utc.astype("int64") / 1e9, 2)
    F = F.merge(Z[["date", "series", "fam", "tkey", "n_z2"]], on=["date", "series", "fam", "tkey"], how="left",
                suffixes=("_x", ""))
    miss = F.n_z2.isna().mean()
    E.loc[fw, "n_z2"] = F.n_z2.to_numpy()
    print(f"SVI events {len(E)}; forward n_z2 joined, missing {miss:.4f}")
    E[["time_utc", "date", "series", "fam", "P", "seg", "n_z2", "detector", "src", "burst_day"]].to_parquet(
        f"{G}/svi_events.parquet", index=False)
    # all-candidate artefact gauge
    A = pd.read_parquet(f"{ROOT}/svi_science/results/s0/analysis_table.parquet", columns=["series", "time_utc", "date"])
    A = A.assign(seg="FWD")
    B = pd.read_parquet(f"{ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet",
                        columns=["series", "t", "date", "primary"])
    B = B[B.primary]
    B["time_utc"] = pd.to_datetime(B.t, unit="s").astype("datetime64[ns]")
    B["seg"] = seg_of(B.date)
    C = pd.concat([B[["series", "time_utc", "date", "seg"]], A[["series", "time_utc", "date", "seg"]]], ignore_index=True)
    C["time_utc"] = C.time_utc.astype("datetime64[ns]")
    C.to_parquet(f"{G}/svi_cand.parquet", index=False)
    print("SVI candidates", C.groupby("series").size().to_dict())
    O = pd.read_parquet(f"{ROOT}/svi_science/results/r23/tables/famdays_observable.parquet", columns=["date", "fam", "seg"])
    O = O.groupby(["date", "seg"]).size().rename("n_obs").reset_index()
    O.to_parquet(f"{G}/svi_obsdays.parquet", index=False)
    ep = pd.read_csv(f"{ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")
    ep[["episode", "start", "end", "major", "big", "kind", "validated"]].to_csv(f"{G}/svi_episodes.csv", index=False)


def _wa_month(p):
    ym = os.path.basename(p)[2:9]
    cols = ["date", "series", "fid", "k", "t", "P_own_WAVE_PNET", "n_sites"]
    D = pd.read_parquet(p, columns=cols + ["net", "reg"], filters=[("net", "==", "core"), ("reg", "==", "NW")])
    cand = D[["series", "t", "date"]].copy()
    D = D[D.P_own_WAVE_PNET >= 0.75].drop(columns=["net", "reg"])
    if len(D):
        fl = pd.read_parquet(f"{WP}/post_v2/dedup/flags/flags_{ym}.parquet", columns=["date", "series", "fid", "k", "dup_D"])
        n0 = len(D)
        D = D.merge(fl, on=["date", "series", "fid", "k"], how="left")
        assert len(D) == n0
        z = []
        for d in D.date.unique():
            q = f"{WP}/mf/prod/det_{d}.parquet"
            if os.path.exists(q):
                X = pd.read_parquet(q, columns=["net", "series", "fid", "k", "n_z2"], filters=[("net", "==", "core")])
                z.append(X.assign(date=d).drop(columns="net"))
        if z:
            Z = pd.concat(z).drop_duplicates(["date", "series", "fid", "k"])
            D = D.merge(Z, on=["date", "series", "fid", "k"], how="left")
        else:
            D["n_z2"] = np.nan
        assert len(D) == n0
    return D, cand


def wa():
    files = sorted(glob.glob(f"{WP}/conf/P/prod/P_*.parquet"))
    with ProcessPoolExecutor(8) as ex:
        R = list(ex.map(_wa_month, files))
    E = pd.concat([r[0] for r in R if len(r[0])], ignore_index=True)
    E["time_utc"] = pd.to_datetime(E.t, unit="s").astype("datetime64[ns]")
    E["dup_D"] = E.dup_D.fillna(False).astype(bool)
    E["era"] = wa_era(E.date)
    print("WA P>=0.75 rows", E.groupby("series").size().to_dict(), "once", E[~E.dup_D].groupby("series").size().to_dict(),
          "n_z2 missing", round(E.n_z2.isna().mean(), 4))
    E = E[~E.dup_D]
    E[["time_utc", "date", "series", "fid", "P_own_WAVE_PNET", "n_sites", "n_z2", "era"]].to_parquet(f"{G}/wa_events.parquet", index=False)
    C = pd.concat([r[1] for r in R], ignore_index=True)
    C["time_utc"] = pd.to_datetime(C.t, unit="s").astype("datetime64[ns]")
    C["era"] = wa_era(C.date)
    C[["series", "time_utc", "date", "era"]].to_parquet(f"{G}/wa_cand.parquet", index=False)
    print("WA candidates", C.groupby("series").size().to_dict())
    r = pd.read_parquet(f"{WP}/post_v2/dedup/rates_daily_v2_dedup.parquet", columns=["unit_kind", "region", "date", "fid"],
                        filters=[("unit_kind", "==", "family"), ("region", "==", "NW")])
    O = r.groupby("date").size().rename("n_obs").reset_index()
    O["era"] = wa_era(O.date)
    O.to_parquet(f"{G}/wa_obsdays.parquet", index=False)
    ep = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/episode_rates/episodes_nw_final.csv")
    ep[["episode", "start", "end", "major", "kind"]].to_csv(f"{G}/wa_episodes.csv", index=False)


if __name__ == "__main__":
    svi()
    wa()
