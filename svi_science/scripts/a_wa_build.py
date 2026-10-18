#!/usr/bin/env python
"""a_wa_build (2026-10-07) -- Part 2 item A, step 1: build the WA-NW event / exposure / family tables for the
repeat of S7 (between-episode LFEs).  Rules in svi_science/results/a_wa_between/PLAN_A.md.  CPU only, < 15 GB.
Reads only; writes into svi_science/results/a_wa_between/.

    python svi_science/scripts/a_wa_build.py
"""
from __future__ import annotations
import os, sys, glob, json, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
WP = f"{ROOT}/lfe_extension_wa/results/wa_prod"
OUT = f"{ROOT}/svi_science/results/a_wa_between"
sys.path.insert(0, f"{ROOT}/svi_science/scripts")
import s3_build as S3B          # import only: slab_grid / _bilin / project (no side effects)

SPAN = json.load(open(f"{WP}/post_v2/catalogue_span.json"))
LOG = open(f"{OUT}/a_build.log", "w")


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def families():
    M = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/nw_families_meta.csv")[["fid", "family_uid"]]
    C = pd.read_csv(f"{ROOT}/data/catalogue/families_v1_3.csv")
    C = C[C.region == "Royer_NW_WA"].set_index("family_uid")
    F = M.copy()
    F["source_group"] = F.family_uid.map(C.nw_source_group)
    F["array_group"] = F.family_uid.map(C.nw_array_group)
    F["copy_of"] = F.family_uid.map(C.nw_copy_of)
    F["primary_uid"] = np.where(F.copy_of.notna(), F.copy_of, F.family_uid)
    for c_out, c_in in (("lat", "lat"), ("lon", "lon"), ("dep", "depth_km"), ("quality", "location_quality"), ("erh", "erh")):
        F[c_out] = F.primary_uid.map(C[c_in])
        F[c_out + "_own"] = F.family_uid.map(C[c_in])
    lons, lats, Z = S3B.slab_grid()
    F["slab2_dep"] = [-float(S3B._bilin(lons, lats, Z, lo, la)) if np.isfinite(lo) and np.isfinite(la) else np.nan
                      for lo, la in zip(F.lon, F.lat)]
    fr = json.load(open(f"{ROOT}/svi_science/results/s3/build_info.json"))["nw_frame"]
    a, d, _, _ = S3B.project(fr, F.lat.values, F.lon.values)
    F["along_km"] = a; F["downdip_km"] = d
    F["ab"] = F.quality.isin(["A", "B"])
    log(f"families {len(F)}; copies {int(F.copy_of.notna().sum())}; source groups {F.source_group.nunique()}; "
        f"primary quality {F.quality.value_counts().to_dict()}; missing dep {int(F.dep.isna().sum())}")
    return F


def events(F):
    cols = ["date", "net", "series", "fid", "k", "t", "reg", "role", "P_own_WAVE_PNET", "n_sta_day"]
    parts = []
    for p in sorted(glob.glob(f"{WP}/conf/P/prod/P_*.parquet")):
        ym = p.split("_")[-1][:7]
        if ym < "2009-08":
            continue
        D = pd.read_parquet(p, columns=cols, filters=[("net", "==", "core"), ("reg", "==", "NW")])
        D = D[D.P_own_WAVE_PNET >= 0.5]
        if not len(D):
            continue
        fl = pd.read_parquet(f"{WP}/post_v2/dedup/flags/flags_{ym}.parquet",
                             columns=["date", "series", "fid", "k", "det_id", "dup_D", "win_fid_D", "win_km_D", "dup_B"])
        n0 = len(D)
        D = D.merge(fl, on=["date", "series", "fid", "k"], how="left")
        assert len(D) == n0, ym
        parts.append(D.drop(columns=["net", "reg"]))
    E = pd.concat(parts, ignore_index=True)
    E["time_utc"] = pd.to_datetime(E.t, unit="s").astype("datetime64[ns]")
    X = pd.read_parquet(f"{WP}/post/tremor_family_day.parquet").set_index(["fid", "date"])
    Y = X.reindex(pd.MultiIndex.from_arrays([E.fid.to_numpy(), E.date.to_numpy()]))
    n30 = Y.n30.fillna(0).to_numpy(); n60 = Y.n60.fillna(0).to_numpy()
    cov = (E.date >= SPAN["covered_from"]) & (E.date <= SPAN["covered_to"])
    E["day_class"] = np.where(cov, np.where(n30 > 0, "tremor", "quiet"), "unknown")
    E["n_trem30"] = np.where(cov, n30, np.nan); E["n_trem60"] = np.where(cov, n60, np.nan)
    sg = F.set_index("fid").source_group
    E["dup_D"] = E.dup_D.fillna(False).astype(bool)
    E["once_D"] = ~E.dup_D
    same = E.win_fid_D.map(sg).to_numpy() == E.fid.map(sg).to_numpy()
    E["once_group"] = ~(E.dup_D & same)
    log(f"event rows (P>=0.5, core, NW, 2009-08..): {len(E)}; without flag row: {int(E.det_id.isna().sum())}")
    return E


def exposure():
    R = pd.read_parquet(f"{WP}/post_v2/dedup/rates_daily_v2_dedup.parquet",
                        columns=["unit_kind", "fid", "date", "region", "era", "day_class", "n_v2_data", "n_v2_data_D"])
    R = R[(R.region == "NW") & (R.unit_kind == "family")].drop(columns=["unit_kind", "region"])
    X = pd.read_parquet(f"{WP}/post/tremor_family_day.parquet")
    R = R.merge(X, on=["fid", "date"], how="left")
    R["n30"] = R.n30.fillna(0); R["n60"] = R.n60.fillna(0)
    R = R[(R.date >= SPAN["covered_from"])].reset_index(drop=True)
    R["fid"] = R.fid.astype(int)
    log(f"scanned NW family-days (covered span start ->): {len(R)}")
    return R


def live_core():
    rows = []
    for p in sorted(glob.glob(f"{WP}/mf/prod/run_*.json")):
        d = json.load(open(p))
        rows.append((d["date"], len(d["members"].get("core", []))))
    L = pd.DataFrame(rows, columns=["date", "n_core_live"])
    log(f"live-core days {len(L)}")
    return L


def main():
    os.makedirs(OUT, exist_ok=True)
    F = families(); F.to_csv(f"{OUT}/families_wa.csv", index=False)
    E = events(F)
    # acceptance: data P>=0.75 counts per era/class vs rates_daily_v2_dedup sums
    R = exposure()
    L = live_core(); L.to_csv(f"{OUT}/live_core_by_day.csv", index=False)
    R = R.merge(L, on="date", how="left")
    R.to_parquet(f"{OUT}/exposure_famday.parquet")
    E = E.merge(L, on="date", how="left")
    E.to_parquet(f"{OUT}/events_all.parquet")
    era = np.where(E.date < "2012-01-01", "E1", np.where(E.date < "2015-01-01", "E2", "E3"))
    h = E[(E.P_own_WAVE_PNET >= 0.75)].assign(era=era[E.P_own_WAVE_PNET >= 0.75])
    a = h.groupby(["era", "day_class", "series"]).agg(raw=("fid", "size"), once_D=("once_D", "sum")).unstack("series")
    r = R.groupby(["era", "day_class"])[["n_v2_data", "n_v2_data_D"]].sum()
    log("acceptance (events P>=0.75 from P files vs rates_daily_v2_dedup, 2009-08-06 ->):\n" + a.to_string() + "\n" + r.to_string())
    log("done")


if __name__ == "__main__":
    main()
