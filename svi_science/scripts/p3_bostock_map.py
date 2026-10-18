#!/usr/bin/env python
"""p3_bostock_map (2026-10-07) -- P3 step 1: Bostock (2015) 2003-2013 catalogue -> our 115 family ids.
Rules frozen in svi_science/results/p3_bostock/PLAN_P3.md.  CPU only.
    python svi_science/scripts/p3_bostock_map.py
"""
import os, sys, json
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/p3_bostock"
BOST = f"{ROOT}/lfe_catalog_repair/results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
FAMS12 = f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv"
S0 = f"{ROOT}/svi_science/results/s0"


def main():
    os.makedirs(OUT, exist_ok=True)
    R = {}
    B = pd.read_csv(BOST)
    FAM = pd.read_csv(FAMS12)
    F = pd.read_csv(f"{S0}/families.csv")
    R["n_rows"] = int(len(B)); R["n_keys_catalogue"] = int(B.groupby(["lat", "lon", "dep"]).ngroups)
    k = {(round(r.lat, 4), round(r.lon, 4), round(r.dep, 4)): int(r.bkey) for r in FAM.itertuples()}
    R["n_keys_e5r"] = len(k)
    B["fid"] = [k.get((round(a, 4), round(b, 4), round(c, 4)), -1) for a, b, c in zip(B.lat, B.lon, B.dep)]
    R["n_unmapped_rows"] = int((B.fid < 0).sum())
    B["time_utc"] = (pd.to_datetime(dict(year=B.year, month=B.month, day=B.day, hour=B.hour, minute=B["min"],
                                         second=B.sec.astype(int))) + pd.to_timedelta(B.sec % 1.0, unit="s")).astype("datetime64[ns]")
    B["date"] = B.time_utc.dt.strftime("%Y-%m-%d")
    # (a) production families mapped
    got = set(B.fid[B.fid >= 0])
    R["prod_fids_mapped"] = int(len(set(F.fid) & got)); R["prod_fids_total"] = int(len(F))
    R["prod_fids_missing"] = sorted(set(F.fid) - got)
    # (b) per-year totals
    raw_y = B.groupby("year").size(); map_y = B[B.fid >= 0].groupby("year").size()
    R["per_year_raw"] = {int(y): int(v) for y, v in raw_y.items()}
    R["per_year_mapped_130"] = {int(y): int(v) for y, v in map_y.items()}
    R["per_year_mapped_115"] = {int(y): int(v) for y, v in B[B.fid.isin(F.fid)].groupby("year").size().items()}
    R["check_b_per_year_equal"] = bool(raw_y.equals(map_y))
    # (c) vs svi_ref_common.bostock() for 2010-2012
    sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
    try:
        import svi_ref_common as SR
        G = SR.bostock(years=(2010, 2011, 2012))
        # svi_ref_common's 'year' is the year of t_origin (catalogue time + offset); compare on catalogue year
        a = B[B.year.isin([2010, 2011, 2012]) & (B.fid >= 0)].groupby(["fid"]).size()
        b = G.groupby("fam").size()
        diff = a.subtract(b, fill_value=0)
        R["check_c_vs_svi_ref_common"] = dict(n_ours=int(a.sum()), n_svi_ref=int(b.sum()),
                                              n_fam_diff=int((diff != 0).sum()), max_abs_diff=int(diff.abs().max()),
                                              note="per family, 2010-2012 pooled; svi_ref_common drops rows whose family "
                                                   "has no accepted 2012 offset (t_origin NaN)")
    except Exception as e:  # pragma: no cover
        R["check_c_vs_svi_ref_common"] = dict(error=repr(e))
    # (d) 2012 per family vs e5r n_all
    c12 = B[(B.year == 2012) & (B.fid >= 0)].groupby("fid").size()
    nall = FAM.set_index("bkey").n_all
    d = c12.reindex(nall.index, fill_value=0) - nall
    R["check_d_2012_vs_n_all"] = dict(n_fam_diff=int((d != 0).sum()), max_abs_diff=int(d.abs().max()))
    cc = R["check_c_vs_svi_ref_common"]
    R["PASS"] = bool(R["prod_fids_mapped"] >= 110 and R["check_b_per_year_equal"] and R["check_d_2012_vs_n_all"]["n_fam_diff"] == 0
                     and cc.get("n_fam_diff", 1) == 0)
    # coverage
    days = np.array(sorted(B.date.unique()), dtype="datetime64[D]")
    gaps = np.diff(days).astype(int)
    st = np.r_[days[:1], days[1:][gaps > 5]]; en = np.r_[days[:-1][gaps > 5], days[-1:]]
    win = pd.DataFrame(dict(first_day=st.astype(str), last_day=en.astype(str)))
    win["n_days_with_events"] = [int(((days >= s) & (days <= e)).sum()) for s, e in zip(st, en)]
    win["span_days"] = (en - st).astype(int) + 1
    win["n_events_115"] = [int(B[B.fid.isin(F.fid) & (B.date >= str(s)) & (B.date <= str(e))].shape[0]) for s, e in zip(st, en)]
    win.to_csv(f"{OUT}/bostock_scanned_windows.csv", index=False)
    R["n_days_with_events"] = int(days.size)
    R["n_days_2003_02_26_to_2013_10_13"] = int((days[-1] - days[0]).astype(int) + 1)
    R["share_of_days_with_events"] = float(days.size / R["n_days_2003_02_26_to_2013_10_13"])
    R["windows"] = win.to_dict("records")
    # mapped table, 115 families
    M = B[B.fid.isin(F.fid)][["fid", "time_utc", "date", "year", "mag", "dep"]].copy()
    Fi = F.set_index("fid")
    M["group"] = M.fid.map(Fi.group); M["downdip_km"] = M.fid.map(Fi.downdip_km); M["along_km"] = M.fid.map(Fi.along_strike_km)
    M = M.sort_values("time_utc").reset_index(drop=True)
    M.to_parquet(f"{OUT}/bostock_events_mapped.parquet", index=False)
    yg = M.groupby(["year", "group"]).size().unstack(fill_value=0)
    yg.to_csv(f"{OUT}/bostock_year_by_group.csv")
    json.dump(R, open(f"{OUT}/step1_checks.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in R.items() if k != "windows"}, indent=1, default=str))
    print(win.to_string(index=False))
    print(yg.to_string())


if __name__ == "__main__":
    main()
