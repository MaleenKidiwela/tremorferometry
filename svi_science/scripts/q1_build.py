"""Q1 build: patch x episode count tables (all episodes, MAJOR + minor) and episode covariates.

Rules: svi_science/results/q1_patch_budget/PLAN.md (pre-registered).
Writes results/q1_patch_budget/:
  counts_long.parquet   catalogue, variant, patchdef, patch, episode, n, t_med (median event time, UTC)
  patches.csv           catalogue, patchdef, patch, as_km, dd_km, lat, lon, n_fam
  episodes.csv          catalogue, episode, start, end, w0, w1, major, era
  tremor_ep.csv         region, episode, n_trem, trem_p05, trem_p95, trem_speed (km/d), trem_dir (+1 NW-ward / -1)
  det_patch_ep.csv      catalogue, patchdef, patch, episode, mean_sta, frac_obs
"""
import glob
import importlib.util
import json
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q1_patch_budget"
ER = f"{ROOT}/lfe_extension_wa/results/episode_rates"
spec = importlib.util.spec_from_file_location("b", f"{ROOT}/svi_science/scripts/s3_build.py")
b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
rng = np.random.default_rng(20261007)


def ts_slope(t, a):
    t = np.asarray(t, float); a = np.asarray(a, float)
    i, j = np.triu_indices(len(t), 1)
    dt = t[j] - t[i]; da = a[j] - a[i]
    m = (np.abs(dt) > 1e-9) & (np.abs(da) > 1e-6)
    s = da[m] / dt[m]
    return float(np.median(s)) if len(s) else np.nan


def bins10(x, lo, hi, shift=0.0):
    edges = np.arange(lo + shift, hi + shift + 10, 10)
    i = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    return i, edges


def count_table(ev, cat, variant, patchdef, col):
    g = ev.groupby([col, "ep"])
    n = g.size().rename("n")
    t = g.time.apply(lambda s: s.sort_values().iloc[len(s) // 2]).rename("t_med")
    d = pd.concat([n, t], axis=1).reset_index().rename(columns={col: "patch", "ep": "episode"})
    d["patch"] = d.patch.astype(str)
    d.insert(0, "patchdef", patchdef); d.insert(0, "variant", variant); d.insert(0, "catalogue", cat)
    return d


def era_of(s):
    s = pd.Timestamp(s)
    return "E1" if s < pd.Timestamp("2012-01-01") else ("E2" if s < pd.Timestamp("2015-01-01") else "E3")


def main():
    info = {}
    rows, prow, erow = [], [], []
    # ------------------------------------------------------------------ SVI
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    G = pd.read_csv(f"{ROOT}/svi_science/results/s0/groups.csv")
    lo, hi = np.floor(F.along_strike_km.min() / 10) * 10, np.ceil(F.along_strike_km.max() / 10) * 10
    bi, ed = bins10(F.along_strike_km.values, lo, hi)
    F["bin10"] = [f"{ed[i]:+.0f}..{ed[i + 1]:+.0f}" for i in bi]
    for _, r in G.iterrows():
        prow.append(dict(catalogue="SVI", patchdef="group", patch=str(int(r.group)), as_km=r.along_strike_km, dd_km=r.downdip_km,
                         lat=r.lat, lon=r.lon, n_fam=int(r.n_families)))
    for k, g in F.groupby("bin10"):
        prow.append(dict(catalogue="SVI", patchdef="bin10", patch=k, as_km=float(k.split("..")[0]) + 5, dd_km=g.downdip_km.mean(),
                         lat=g.lat.mean(), lon=g.lon.mean(), n_fam=len(g)))
    Esvi = b.padded_windows(pd.read_csv(f"{ER}/episodes_svi_final.csv"))
    for _, r in Esvi.iterrows():
        erow.append(dict(catalogue="SVI", episode=r.episode, start=r.start, end=r.end, w0=r.w0, w1=r.w1, major=bool(r.major), era="SVI"))
    for tab, tag in [("analysis_table.parquet", ""), ("analysis_table_S13.parquet", "S13_")]:
        T = pd.read_parquet(f"{ROOT}/svi_science/results/s0/{tab}", columns=["run", "series", "fam", "P_same", "keep_xfam8", "time_utc"],
                            filters=[("P_same", ">=", 0.5)])
        info[f"svi_{tag}runs"] = T.run.astype(str).value_counts().to_dict()
        T = T.rename(columns={"time_utc": "time"})
        T["ep"] = b.assign(T.time, Esvi)
        T = T[T.ep.notna()].merge(F[["fid", "group", "bin10"]].rename(columns={"fid": "fam"}), on="fam", how="left")
        T["group"] = T.group.astype(int)
        for series, sv in [("data", ""), ("c1", "c1_")]:
            d = T[T.series == series]
            sel = {"once075": d[(d.P_same >= 0.75) & d.keep_xfam8], "raw075": d[d.P_same >= 0.75], "once050": d[d.keep_xfam8]}
            for v, dv in sel.items():
                if series == "c1" and v == "raw075":
                    continue
                if tag and v == "raw075":
                    continue
                for pdf in ["group", "bin10"]:
                    if (tag or series == "c1") and pdf == "bin10":
                        continue
                    rows.append(count_table(dv, "SVI", f"{tag}{sv}{v}", pdf, pdf))
        del T
    # detectability (SVI)
    D = pd.read_parquet(f"{ROOT}/svi_science/results/s0/famday_detectability.parquet", columns=["fid", "date", "observable", "fam_n_sta"])
    D["time"] = pd.to_datetime(D.date) + pd.Timedelta(hours=12)
    D["ep"] = b.assign(D.time, Esvi); D = D[D.ep.notna()]
    D = D.merge(F[["fid", "group", "bin10"]], on="fid")
    det = []
    for pdf in ["group", "bin10"]:
        q = D.groupby([pdf, "ep"]).agg(mean_sta=("fam_n_sta", lambda s: float(np.mean(s))), frac_obs=("observable", "mean")).reset_index()
        q = q.rename(columns={pdf: "patch", "ep": "episode"}); q["patch"] = q.patch.astype(str)
        q.insert(0, "patchdef", pdf); q.insert(0, "catalogue", "SVI"); det.append(q)
    del D

    # ------------------------------------------------------------------ NW
    SG = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3.csv")
    AC = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3_array_confirmed.csv").rename(columns={"group": "source", "group_size": "source_size"})
    SG = SG.merge(AC, on="fid")
    src = pd.read_csv(f"{ROOT}/svi_science/results/s3/sources_nw.csv")
    ab = src[src.quality.isin(["A", "B"])]
    lo = np.floor(np.quantile(ab.as_km, 0.03) / 10) * 10; hi = np.ceil(np.quantile(ab.as_km, 0.97) / 10) * 10
    fam = SG[["fid", "source"]].merge(src[["source", "as_km", "dd_km", "lat", "lon"]], on="source")
    for nm, sh in [("bin10", 0.0), ("bin10s5", 5.0)]:
        i, ed = bins10(fam.as_km.values, lo, hi, sh)
        fam[nm] = [f"{ed[k]:+.0f}..{ed[k + 1]:+.0f}" for k in i]
        for k, g in fam.groupby(nm):
            prow.append(dict(catalogue="NW", patchdef=nm, patch=k, as_km=float(np.median(g.as_km)), dd_km=g.dd_km.mean(),
                             lat=g.lat.mean(), lon=g.lon.mean(), n_fam=len(g)))
    info["nw_bins"] = dict(lo=float(lo), hi=float(hi))
    Enw = b.padded_windows(pd.read_csv(f"{ER}/episodes_nw_final.csv"))
    for _, r in Enw.iterrows():
        erow.append(dict(catalogue="NW", episode=r.episode, start=r.start, end=r.end, w0=r.w0, w1=r.w1, major=bool(r.major), era=era_of(r.start)))
    C = pd.read_parquet(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet",
                        columns=["det_id", "time", "region", "family", "P_own_pnet", "selected_v2"],
                        filters=[("region", "==", "NW"), ("P_own_pnet", ">=", 0.5)])
    C["ep"] = b.assign(C.time, Enw); C = C[C.ep.notna()].copy()
    fl = []
    for f in sorted(glob.glob(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/dedup/flags/flags_20*.parquet")):
        q = pd.read_parquet(f, columns=["det_id", "series", "dup_D"])
        fl.append(q[(q.series == "data") & q.det_id.isin(C.det_id)].drop(columns="series"))
    fl = pd.concat(fl).drop_duplicates("det_id")
    C = C.merge(fl, on="det_id", how="left")
    info["nw_rows"] = int(len(C)); info["nw_rows_without_flag"] = int(C.dup_D.isna().sum())
    info["nw_selected_without_flag"] = int((C.selected_v2 & C.dup_D.isna()).sum())
    C = C.merge(fam[["fid", "bin10", "bin10s5"]].rename(columns={"fid": "family"}), on="family", how="left")
    info["nw_rows_without_patch"] = int(C.bin10.isna().sum())
    C = C[C.bin10.notna()]
    once = ~C.dup_D.fillna(False).astype(bool)
    sel = {"once075": C[C.selected_v2 & once], "raw075": C[C.selected_v2], "once050": C[once]}
    for v, dv in sel.items():
        for pdf in ["bin10", "bin10s5"]:
            if v != "once075" and pdf == "bin10s5":
                continue
            rows.append(count_table(dv, "NW", v, pdf, pdf))
    del C, sel
    # fixed-9 borehole run (WAVE-only, not de-duplicated)
    C9 = pd.read_parquet(f"{ROOT}/lfe_extension_wa/results/wa_fixed10/post_v2_core9/catalogue_v2.parquet",
                         columns=["time", "region", "family", "selected_v2"], filters=[("region", "==", "NW"), ("selected_v2", "==", True)])
    C9["ep"] = b.assign(C9.time, Enw); C9 = C9[C9.ep.notna()]
    C9 = C9.merge(fam[["fid", "bin10"]].rename(columns={"fid": "family"}), on="family", how="inner")
    info["nw_f9_rows"] = int(len(C9))
    rows.append(count_table(C9, "NW", "F9_raw075", "bin10", "bin10"))
    del C9
    # detectability (NW)
    W = pd.read_parquet(f"{ER}/wa_core_famday.parquet", columns=["fid", "date", "n_sta", "observable"])
    W["time"] = pd.to_datetime(W.date) + pd.Timedelta(hours=12)
    W["ep"] = b.assign(W.time, Enw); W = W[W.ep.notna()].merge(fam[["fid", "bin10", "bin10s5"]], on="fid")
    for pdf in ["bin10", "bin10s5"]:
        q = W.groupby([pdf, "ep"]).agg(mean_sta=("n_sta", "mean"), frac_obs=("observable", "mean")).reset_index()
        q = q.rename(columns={pdf: "patch", "ep": "episode"}); q.insert(0, "patchdef", pdf); q.insert(0, "catalogue", "NW"); det.append(q)
    del W

    # ------------------------------------------------------------------ Bostock 2003-2013 (S0 groups)
    B = pd.read_parquet(f"{ROOT}/svi_science/results/p3_bostock/bostock_events_mapped.parquet", columns=["fid", "time_utc", "group"])
    B = B.rename(columns={"time_utc": "time"})
    Eb = pd.read_csv(f"{ROOT}/svi_science/results/p3_bostock/episodes_bostock_final.csv")
    Eb["w0"] = pd.to_datetime(Eb.start); Eb["w1"] = pd.to_datetime(Eb.end) + pd.Timedelta(days=1)
    for _, r in Eb.iterrows():
        erow.append(dict(catalogue="BOSTOCK", episode=r.episode, start=r.start, end=r.end, w0=r.w0, w1=r.w1, major=True, era="B"))
    B["ep"] = b.assign(B.time, Eb)
    info["bostock_outside_episodes"] = int(B.ep.isna().sum())
    B = B[B.ep.notna()]
    rows.append(count_table(B, "BOSTOCK", "all", "group", "group"))
    # Bostock counted once: drop an event within +-2 s of an earlier kept event of another family (simple check variant)
    B = B.sort_values("time").reset_index(drop=True)
    t = B.time.values.astype("int64") / 1e9
    keep = np.ones(len(B), bool); last_t = -1e18
    for i in range(len(B)):
        if t[i] - last_t <= 2.0:
            keep[i] = False
        else:
            last_t = t[i]
    info["bostock_once_drop_frac"] = float(1 - keep.mean())
    rows.append(count_table(B[keep], "BOSTOCK", "once2s", "group", "group"))
    for _, r in G.iterrows():
        prow.append(dict(catalogue="BOSTOCK", patchdef="group", patch=str(int(r.group)), as_km=r.along_strike_km, dd_km=r.downdip_km,
                         lat=r.lat, lon=r.lon, n_fam=int(r.n_families)))

    CL = pd.concat(rows, ignore_index=True)
    CL.to_parquet(f"{OUT}/counts_long.parquet")
    pd.DataFrame(prow).to_csv(f"{OUT}/patches.csv", index=False)
    pd.DataFrame(erow).to_csv(f"{OUT}/episodes.csv", index=False)
    pd.concat(det, ignore_index=True).to_csv(f"{OUT}/det_patch_ep.csv", index=False)

    # ------------------------------------------------------------------ tremor per MAJOR episode
    P = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    P["time"] = pd.to_datetime(P.time).astype("datetime64[ns]")
    bi = json.load(open(f"{ROOT}/svi_science/results/s3/build_info.json"))
    trows = []
    for reg, E, flat, flon, fr in [("SVI", Esvi, F.lat.values, F.lon.values, bi["svi_frame"]),
                                   ("NW", Enw, ab.lat.values, ab.lon.values, bi["nw_frame"])]:
        m = b.near_mask(P.lat.values, P.lon.values, flat, flon, 10.0)
        Pr = P[m].copy(); Pr["ep"] = b.assign(Pr.time, E); Pr = Pr[Pr.ep.notna()]
        a, _, _, _ = b.project(fr, Pr.lat, Pr.lon); Pr["as_km"] = a
        for ep, g in Pr.groupby("ep"):
            g = g.sort_values("time")
            idx = np.sort(rng.choice(len(g), 1500, replace=False)) if len(g) > 1500 else np.arange(len(g))
            td = (g.time.values[idx] - g.time.values[0]) / np.timedelta64(1, "D")
            s = ts_slope(td, g.as_km.values[idx]) if len(g) >= 10 else np.nan
            trows.append(dict(region=reg, episode=ep, n_trem=len(g), trem_p05=float(np.quantile(g.as_km, 0.05)),
                              trem_p95=float(np.quantile(g.as_km, 0.95)), trem_speed=s,
                              trem_dir=(int(np.sign(s)) if np.isfinite(s) else 0)))
    pd.DataFrame(trows).to_csv(f"{OUT}/tremor_ep.csv", index=False)
    json.dump(info, open(f"{OUT}/build_info.json", "w"), indent=1, default=str)
    print(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
