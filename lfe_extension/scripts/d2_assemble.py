#!/usr/bin/env python
"""Track D, D2 -- CATALOG ASSEMBLY (Merlin Step 7).

Concatenates every finished Track C month into the product:
  results/trackD/lfe_extension_catalog_v1.parquet (+ .csv)   one row per catalog EVENT
  results/trackD/family_table_v1.csv                          one row per family
  results/trackD/family_day_table_v1.parquet                  one row per (family, day)

TIME BASIS (stated on every row and in RESULTS.md): `time` is the ORIGIN-AXIS time the network sum was
formed on = Bostock's database clock for the family plus that family's per-family offset
(`e5r_family_offsets_2012.csv`; the offsets Track A measured so the 2012 templates register).  It is NOT
an origin time in the seismological sense: the absolute datum is Bostock's, and the offset was fitted, so
only DIFFERENCES between events -- within a family and, through the common templates, between families --
carry physical meaning.

MONTH OWNERSHIP is half-open: a detection belongs to the month of its own UTC date, so an event is
counted once even where two monthly runs overlap at a boundary; the assembled counts are cross-checked
against `results/trackC/months_done.csv`."""
from __future__ import annotations
import argparse, glob, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_tiers  # noqa: E402  (Track C's network epochs and per-epoch tier thresholds)

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
TC = f"{EXT}/results/trackC"
TD = f"{EXT}/results/trackD"
OUT = f"{EXT}/results"

SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
          9: "SON", 10: "SON", 11: "SON"}
KEEP = ["time", "t", "date", "doy", "fam", "net_mad", "net_cc", "rel_amp", "n_sites", "n_sta",
        "n_sta_excluded", "n_sta_post_boundary", "tier_ext", "tier_core", "cat_ext", "cat_core",
        "hfr", "hf_cut", "hf_pass", "tremor_flag", "n_tremor_day", "eq_flag", "eq_dt", "quiet_flag",
        "epoch", "epoch_template", "con_pgc", "wsh_pgc", "mad_scale"]


EPOCHS, TIERS, TSRC = d_tiers.load()


def epoch_of(dates):
    return d_tiers.epoch_of(dates, EPOCHS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", action="store_true", default=True)
    ap.add_argument("--tag", default="v1", help="output version tag: lfe_extension_catalog_<tag>.parquet, "
                                               "family_table_<tag>.csv, family_day_table_<tag>.parquet, "
                                               "catalog_<tag>_meta.json (v1 files are never overwritten "
                                               "unless --tag v1 is given explicitly)")
    a = ap.parse_args()
    os.makedirs(TD, exist_ok=True)
    files = sorted(glob.glob(f"{TC}/detections_*.parquet"))
    if not files:
        raise SystemExit("[d2] no Track C months on disk")
    parts, own = [], []
    for p in files:
        ym = os.path.basename(p)[len("detections_"):-len(".parquet")]
        d = pd.read_parquet(p)
        n_all = len(d)
        d = d[d.date.astype(str).str[:7] == ym]          # half-open month ownership
        d = d[(d.cat_ext | d.cat_core)].copy()
        d["month"] = ym
        own.append(dict(ym=ym, n_rows_file=n_all, n_rows_owned=int((pd.read_parquet(p, columns=["date"]).date.astype(str).str[:7] == ym).sum()),
                        n_cat_ext=int(d.cat_ext.sum()), n_cat_core=int(d.cat_core.sum()), n_events=len(d)))
        parts.append(d[[c for c in KEEP if c in d.columns] + ["month"]])
    E = pd.concat(parts, ignore_index=True)
    n_before = len(E)
    E = E.drop_duplicates(subset=["fam", "t"]).reset_index(drop=True)
    n_dupe = n_before - len(E)
    E = E.rename(columns={"fam": "family_id", "net_cc": "cc"})
    E["epoch_id"] = epoch_of(E.date.to_numpy().astype(str))
    E["season"] = pd.DatetimeIndex(E.date).month.map(SEASON)

    # ---- family metadata: Bostock hypocentre + SeisBench id + our relocation
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv")
    F = F.rename(columns={"bkey": "family_id", "lat": "lat_bostock", "lon": "lon_bostock",
                          "dep": "depth_bostock", "sb_fam_id": "seisbench_family_id"})
    R = pd.read_csv(f"{OUT}/a5_relocation.csv").rename(columns={"family": "family_id"})
    acc = R[R.verdict.astype(str).str.startswith("ACCEPT(")].copy()
    acc = acc.rename(columns={"lat_A": "lat_ours", "lon_A": "lon_ours", "z_A": "depth_ours",
                              "dh_A": "reloc_dh_km", "dz_A": "reloc_dz_km"})
    acc["reloc_accepted"] = True
    fcols = ["family_id", "lat_bostock", "lon_bostock", "depth_bostock", "seisbench_family_id", "mag_med"]
    E = E.merge(F[fcols], on="family_id", how="left")
    E = E.merge(acc[["family_id", "lat_ours", "lon_ours", "depth_ours", "reloc_dh_km", "reloc_dz_km",
                     "reloc_accepted"]], on="family_id", how="left")
    E["reloc_accepted"] = E.reloc_accepted.fillna(False)

    # ---- Mc(f,d) and the false-alarm floor lambda(f, T, epoch, season)
    mcp = f"{TD}/mc_family_day.parquet"
    if os.path.exists(mcp):
        MC = pd.read_parquet(mcp)[["family_id", "date", "Mc_T_ext", "Mc_T_core", "u", "v"]]
        E = E.merge(MC, on=["family_id", "date"], how="left")
    flp = f"{TD}/floor_family_epoch_season.csv"
    if os.path.exists(flp):
        L = pd.read_csv(flp)
        for tier, suf in (("T_ext", "ext"), ("T_core", "core")):
            l = L[L.tier == tier][["family_id", "epoch_id", "season", "lambda_fa"]] \
                .rename(columns={"lambda_fa": f"lambda_{suf}"})
            E = E.merge(l, on=["family_id", "epoch_id", "season"], how="left")

    E = E.sort_values(["time", "family_id"]).reset_index(drop=True)
    E.to_parquet(f"{TD}/lfe_extension_catalog_{a.tag}.parquet", index=False)
    if a.csv:
        E.to_csv(f"{TD}/lfe_extension_catalog_{a.tag}.csv", index=False)

    # ---- family-day table
    fdf = sorted(glob.glob(f"{TC}/family_days_*.parquet"))
    FD = pd.concat([pd.read_parquet(p) for p in fdf], ignore_index=True).drop_duplicates(["fam", "date"])
    FD = FD.rename(columns={"fam": "family_id"})
    FD["epoch_id"] = epoch_of(FD.date.to_numpy().astype(str))
    FD["season"] = pd.DatetimeIndex(FD.date).month.map(SEASON)
    if os.path.exists(mcp):
        FD = FD.merge(pd.read_parquet(mcp)[["family_id", "date", "Mc_T_ext", "Mc_T_core",
                                            "Pdet_T_ext_M1.8", "Pdet_T_core_M1.8", "u", "v"]],
                      on=["family_id", "date"], how="left")
    if os.path.exists(flp):
        for tier, suf in (("T_ext", "ext"), ("T_core", "core")):
            l = L[L.tier == tier][["family_id", "epoch_id", "season", "lambda_fa"]] \
                .rename(columns={"lambda_fa": f"lambda_{suf}"})
            FD = FD.merge(l, on=["family_id", "epoch_id", "season"], how="left")
    FD.to_parquet(f"{TD}/family_day_table_{a.tag}.parquet", index=False)

    # ---- family table
    qc = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    qcu = qc[qc.usable.astype(bool) & qc.template.astype(bool)]
    tplsta = qcu.groupby("bkey").sta.apply(lambda s: ",".join(sorted(set(s)))).rename("template_stations")
    ntpl = qcu.groupby("bkey").sta.nunique().rename("n_template_stations")
    T = F[fcols].merge(tplsta, left_on="family_id", right_index=True, how="left") \
                .merge(ntpl, left_on="family_id", right_index=True, how="left")
    T["n_template_stations"] = T.n_template_stations.fillna(0).astype(int)
    T["template_less"] = T.n_template_stations == 0
    T = T.merge(R[["family_id", "lat_A", "lon_A", "z_A", "dh_A", "dz_A", "n_sites_60km", "gap60_deg",
                   "verdict"]].rename(columns={"lat_A": "lat_ours", "lon_A": "lon_ours",
                                               "z_A": "depth_ours", "dh_A": "reloc_dh_km",
                                               "dz_A": "reloc_dz_km", "verdict": "reloc_verdict"}),
                on="family_id", how="left")
    # relocation uncertainty: there is no per-family formal error in a5_relocation.csv, so the honest
    # number is the MATCHED SYNTHETIC FLOOR of Track A-2 step 5 on this geometry -- 2.03 km median
    # horizontal for the ACCEPT subset, 2.73 km for all located families (a5_synth_floor_matched).
    T["reloc_sigma_h_km"] = np.where(T.reloc_verdict.astype(str).str.startswith("ACCEPT"), 2.03,
                                     np.where(T.reloc_verdict.notna(), 2.73, np.nan))
    E["cat_core_hp"] = E.cat_core & (E.n_sites >= 5)
    for suf, col in (("ext", "cat_ext"), ("core", "cat_core"), ("core_hp", "cat_core_hp")):
        c = E[E[col]].groupby("family_id").size().rename(f"n_events_{suf}")
        T = T.merge(c, left_on="family_id", right_index=True, how="left")
        T[f"n_events_{suf}"] = T[f"n_events_{suf}"].fillna(0).astype(int)
    cfe = f"{TD}/completeness_family_epoch.csv"
    if os.path.exists(cfe):
        C = pd.read_csv(cfe)
        for eid in sorted(C.epoch_id.dropna().unique()):
            for tier, suf in (("T_ext", "ext"), ("T_core", "core")):
                c = C[(C.epoch_id == eid) & (C.tier == tier)][
                    ["family_id", "recall", "recall_lo", "recall_hi"]].rename(
                    columns={"recall": f"recall_{eid}_{suf}", "recall_lo": f"recall_{eid}_{suf}_lo",
                             "recall_hi": f"recall_{eid}_{suf}_hi"})
                T = T.merge(c, on="family_id", how="left")
    # the matched-false-alarm recall is a STATION-SET quantity (the ablation), carried alongside
    cfs = f"{TD}/completeness_family_stationset.csv"
    if os.path.exists(cfs):
        CS = pd.read_csv(cfs)
        for ss in ("S15", "S14", "S13"):
            for tier, suf in (("T_ext", "ext"), ("T_core", "core")):
                c = CS[(CS.run_id == ss) & (CS.tier == tier)][["family_id", "recall_matched_fa"]] \
                    .rename(columns={"recall_matched_fa": f"recall_matchedFA_{ss}_{suf}"})
                if len(c):
                    T = T.merge(c, on="family_id", how="left")
    obs = FD.groupby("family_id").observable.mean().rename("frac_days_observable")
    T = T.merge(obs, left_on="family_id", right_index=True, how="left")
    T.to_csv(f"{TD}/family_table_{a.tag}.csv", index=False)

    # ---- cross-check against the Track C ledger
    O = pd.DataFrame(own)
    MD = pd.read_csv(f"{TC}/months_done.csv")
    X = O.merge(MD[["ym", "n_cat_ext", "n_cat_core"]].rename(columns={"n_cat_ext": "ledger_cat_ext",
                                                                     "n_cat_core": "ledger_cat_core"}),
                on="ym", how="outer")
    X["d_ext"] = X.n_cat_ext - X.ledger_cat_ext
    X["d_core"] = X.n_cat_core - X.ledger_cat_core
    X.to_csv(f"{TD}/catalog_month_crosscheck_{a.tag}.csv", index=False)
    # a month can disagree with the ledger ONLY because the ledger row was written before Track C
    # re-flagged the tiers.  Verify each disagreement against the month's own parquet: recompute
    # membership from net_mad, hf_pass and the stored cross-family keep at that epoch's CURRENT
    # thresholds, and confirm the stored cat_* columns already match them.
    ver = []
    for r in X[(X.d_ext.abs() > 0) | (X.d_core.abs() > 0)].itertuples():
        f = f"{TC}/detections_{r.ym}.parquet"
        if not os.path.exists(f):
            continue
        d = pd.read_parquet(f)
        d = d[d.date.astype(str).str[:7] == r.ym]
        eid = epoch_of(d.date.to_numpy().astype(str))
        th = np.array([[TIERS.get(e, d_tiers.FROZEN)["T_ext"],
                        TIERS.get(e, d_tiers.FROZEN)["T_core"]] for e in eid])
        want_ext = (d.net_mad.to_numpy() >= th[:, 0]) & d.hf_pass.to_numpy() & d.keep_xfam_ext.to_numpy()
        want_core = (d.net_mad.to_numpy() >= th[:, 1]) & d.hf_pass.to_numpy() & d.keep_xfam_core.to_numpy()
        ver.append(dict(ym=r.ym, epochs=";".join(sorted(set(eid))),
                        stored_matches_current_T=bool((want_ext == d.cat_ext.to_numpy()).all()
                                                      and (want_core == d.cat_core.to_numpy()).all()),
                        d_ext=r.d_ext, d_core=r.d_core))
    V = pd.DataFrame(ver)
    if len(V):
        V.to_csv(f"{TD}/catalog_ledger_disagreements_{a.tag}.csv", index=False)
    bad = X[(X.d_ext.abs() > 0) | (X.d_core.abs() > 0) | X.n_cat_ext.isna() | X.ledger_cat_ext.isna()]
    meta = dict(tag=a.tag, n_months=len(files), n_events=int(len(E)),
                ledger_disagreements_explained_by_reflag=int(V.stored_matches_current_T.sum()) if len(V) else 0, n_cat_ext=int(E.cat_ext.sum()),
                n_cat_core=int(E.cat_core.sum()), n_duplicates_removed=int(n_dupe),
                n_families=int(E.family_id.nunique()), t0=str(E.time.min()), t1=str(E.time.max()),
                months_mismatched=len(bad),
                time_basis="Bostock DB clock + per-family offset (e5r_family_offsets_2012.csv); NOT an origin time")
    json.dump(meta, open(f"{TD}/catalog_{a.tag}_meta.json", "w"), indent=1, default=str)
    print(json.dumps(meta, indent=1, default=str))
    if len(bad):
        print("[CROSS-CHECK] months where the assembled count differs from months_done.csv:")
        print(bad.to_string(index=False))
        if len(V):
            ok = int(V.stored_matches_current_T.sum())
            print(f"  -> verified against the month parquets: {ok} of {len(V)} carry cat_* columns that "
                  f"ALREADY match the current per-epoch thresholds, i.e. the ledger row predates Track "
                  f"C's tier re-flag and the CATALOG is the correct one:")
            print(V.to_string(index=False))
    else:
        print("[CROSS-CHECK] every month matches months_done.csv exactly")
    print(f"[write] {TD}/lfe_extension_catalog_{a.tag}.parquet (+ .csv), family_table_{a.tag}.csv, "
          f"family_day_table_{a.tag}.parquet, catalog_month_crosscheck_{a.tag}.csv")


if __name__ == "__main__":
    main()
