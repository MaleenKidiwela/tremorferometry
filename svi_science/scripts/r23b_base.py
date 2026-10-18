#!/usr/bin/env python
"""r23b_base (2026-10-11) -- foundation for R23b (findings 1 and 6 on the 23-year SVI record).
Rules: svi_science/results/r23b/PLAN.md §1 (pre-registered). Read-only on every input; writes results/r23b/tables.

Outputs (tables/):
  windows.csv            145 r23 episodes with S3/Q1 padded windows (w0, w1), seg, seggroup, major, big
  fam_ep.parquet         family x episode: n_once, n_c1, R, R_kc, D (exposed days), X (sum n30), n30 env-only cols
  patch_ep.parquet       group x episode: same sums + C_eq (headline unit) + C_kc
  tremor_ep.csv          per episode: tremor p05/p95/speed/dir (source by date) + env-only columns
  det_patch_ep.csv       per group x episode: mean day-station count, fraction observable (D1 column)
  base_summary.json      acceptance numbers
"""
from __future__ import annotations
import importlib.util, json, os, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/r23b/tables"
R23 = f"{ROOT}/svi_science/results/r23/tables"
S0 = f"{ROOT}/svi_science/results/s0"
Q1 = f"{ROOT}/svi_science/results/q1_patch_budget"
BWR = f"{ROOT}/data/catalogue/svi_rates_v1_4_backward_famday.parquet"
FWR = f"{ROOT}/lfe_extension/results/refine/rates/famday.parquet"
PNSN_T0 = "2010-01-06"
spec = importlib.util.spec_from_file_location("s3b", f"{ROOT}/svi_science/scripts/s3_build.py")
s3b = importlib.util.module_from_spec(spec); spec.loader.exec_module(s3b)
SEGGROUP = {"CN4f": "ERA", "ERA": "ERA", "B2f": "B2B3", "B3f": "B2B3", "S15f": "S15f", "FWD": "FWD"}
RES = {}
LOG = open(f"{ROOT}/svi_science/results/r23b/logs/r23b_base.log", "w")


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"; print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def ts_slope(t, a):
    t = np.asarray(t, float); a = np.asarray(a, float)
    i, j = np.triu_indices(len(t), 1)
    dt = t[j] - t[i]; da = a[j] - a[i]
    m = (np.abs(dt) > 1e-9) & (np.abs(da) > 1e-6)
    s = da[m] / dt[m]
    return float(np.median(s)) if len(s) else np.nan


def main():
    rng = np.random.default_rng(20261007)
    F = pd.read_csv(f"{S0}/families.csv")
    fids = F.fid.astype(int).tolist()
    grp = dict(zip(F.fid, F.group.astype(int)))
    # ------------------------------------------------------------------ episodes + padded windows
    E = pd.read_csv(f"{R23}/episodes_2003_2026.csv")
    W = s3b.padded_windows(E[["episode", "start", "end", "major", "big", "seg", "kind", "S_rate", "n_once", "n_bostock",
                              "family_exposure_days", "gnss_D_mm", "big_S3"]].copy())
    W["seggroup"] = W.seg.map(SEGGROUP)
    W.to_csv(f"{OUT}/windows.csv", index=False)
    log(f"windows: {len(W)} episodes, {int(W.major.sum())} majors, {int(W.big.sum())} big")
    # ------------------------------------------------------------------ events
    V = pd.read_parquet(f"{R23}/events_2003_2026.parquet", columns=["series", "fam", "time_utc", "date", "seg"])
    V["ep"] = s3b.assign(V.time_utc, W)
    V = V[V.ep.notna()]
    n = V.groupby(["fam", "ep", "series"]).size().unstack(fill_value=0).reindex(columns=["data", "c1"], fill_value=0)
    n.columns = ["n_once", "n_c1"]
    # ------------------------------------------------------------------ rates per family-day
    rb = pd.read_parquet(BWR, columns=["fam", "date", "primary", "rate", "rate_kc", "n_sta_day", "sset"])
    rb = rb[rb.primary].drop(columns="primary")
    rf = pd.read_parquet(FWR, columns=["fam", "date", "rate", "n_sta_day", "sset"]).assign(rate_kc=lambda x: x.rate)
    R = pd.concat([rb, rf[rb.columns]], ignore_index=True)
    R = R[R.fam.isin(fids)]
    R["time"] = pd.to_datetime(R.date) + pd.Timedelta(hours=12)
    R["ep"] = s3b.assign(R.time, W)
    R = R[R.ep.notna()]
    rs = R.groupby(["fam", "ep"]).agg(R=("rate", "sum"), R_kc=("rate_kc", "sum"), n_days=("rate", "size"),
                                      frac_S15=("sset", lambda s: float((s == "S15").mean()))).reset_index()
    # observable fraction + day station count (D1 column)
    O = pd.read_parquet(f"{R23}/famdays_observable.parquet", columns=["fam", "date", "n_sta_day"])
    O["time"] = pd.to_datetime(O.date) + pd.Timedelta(hours=12); O["ep"] = s3b.assign(O.time, W); O = O[O.ep.notna()]
    nobs = O.groupby(["fam", "ep"]).size().rename("n_obs_days")
    # ------------------------------------------------------------------ tremor exposure per family-episode (r23 n30)
    X = pd.read_parquet(f"{R23}/famday_tremor.parquet")
    main = X[((X.src == "env") & (X.date < PNSN_T0)) | ((X.src == "pnsn") & (X.date >= PNSN_T0))].copy()
    main["time"] = pd.to_datetime(main.date) + pd.Timedelta(hours=12); main["ep"] = s3b.assign(main.time, W)
    main = main[main.ep.notna()]
    ex = main.groupby(["fam", "ep"]).agg(D=("n30", lambda s: int((s > 0).sum())), X=("n30", "sum"))
    # ------------------------------------------------------------------ family x episode table
    idx = pd.MultiIndex.from_product([fids, W.episode.tolist()], names=["fam", "ep"])
    FE = pd.DataFrame(index=idx).join(n).join(rs.set_index(["fam", "ep"])).join(nobs).join(ex)
    for c in ["n_once", "n_c1", "R", "R_kc", "n_days", "n_obs_days", "D", "X"]:
        FE[c] = FE[c].fillna(0)
    FE["R"] = FE.R.clip(lower=0); FE["R_kc"] = FE.R_kc.clip(lower=0)
    FE = FE.reset_index().rename(columns={"ep": "episode"})
    FE["group"] = FE.fam.map(grp)
    tot = FE.groupby("episode").agg(N_once=("n_once", "sum"), SR=("R", "sum"), SRkc=("R_kc", "sum"))
    FE = FE.join(tot, on="episode")
    FE["C_eq"] = np.where(FE.SR > 0, FE.N_once * FE.R / FE.SR, 0.0)
    FE["C_kc"] = np.where(FE.SRkc > 0, FE.N_once * FE.R_kc / FE.SRkc, 0.0)
    FE.to_parquet(f"{OUT}/fam_ep.parquet", index=False)
    PE = FE.groupby(["group", "episode"])[["n_once", "n_c1", "R", "R_kc", "C_eq", "C_kc"]].sum().reset_index()
    # median event time per group-episode (Q1 offset O, used for T1)
    Vd = V[V.series == "data"].copy(); Vd["group"] = Vd.fam.map(grp)
    tm = Vd.groupby(["group", "ep"]).time_utc.apply(lambda s: s.sort_values().iloc[len(s) // 2]).rename("t_med").reset_index()
    PE = PE.merge(tm.rename(columns={"ep": "episode"}), on=["group", "episode"], how="left")
    PE.to_parquet(f"{OUT}/patch_ep.parquet", index=False)
    # detectability column
    Rn = R.groupby(["fam", "ep"]).n_sta_day.mean().rename("msta").reset_index()
    Dt = FE[["fam", "episode", "n_days", "n_obs_days", "group"]].merge(Rn.rename(columns={"ep": "episode"}), on=["fam", "episode"], how="left")
    Dt["frac_obs"] = np.where(Dt.n_days > 0, Dt.n_obs_days / Dt.n_days.clip(lower=1), 0)
    DT = Dt.groupby(["group", "episode"]).agg(mean_sta=("msta", "mean"), frac_obs=("frac_obs", "mean")).reset_index()
    DT["mean_sta"] = DT.mean_sta.fillna(DT.mean_sta.median())
    DT.to_csv(f"{OUT}/det_patch_ep.csv", index=False)
    # ------------------------------------------------------------------ tremor per episode
    TP = pd.read_parquet(f"{R23}/tremor_points.parquet")
    TP = TP[TP.dmin_fam_km <= 10.0].copy()
    TP["time"] = pd.to_datetime(TP.time).astype("datetime64[ns]")
    TP["ep"] = s3b.assign(TP.time, W); TP = TP[TP.ep.notna()]
    fr = json.load(open(f"{ROOT}/svi_science/results/s3/build_info.json"))["svi_frame"]
    a, _, _, _ = s3b.project(fr, TP.lat, TP.lon); TP["as_km"] = np.asarray(a)
    rows = []
    for (ep, src), g in TP.groupby(["ep", "src"]):
        g = g.sort_values("time")
        ii = np.sort(rng.choice(len(g), 1500, replace=False)) if len(g) > 1500 else np.arange(len(g))
        td = (g.time.values[ii] - g.time.values[0]) / np.timedelta64(1, "D")
        s = ts_slope(td, g.as_km.values[ii]) if len(g) >= 10 else np.nan
        rows.append(dict(episode=ep, src=src, n_trem=len(g), trem_p05=float(np.quantile(g.as_km, 0.05)),
                         trem_p95=float(np.quantile(g.as_km, 0.95)), trem_speed=s, trem_dir=int(np.sign(s)) if np.isfinite(s) else 0))
    T = pd.DataFrame(rows)
    Wd = W.set_index("episode")
    T["date_src"] = [("env" if str(Wd.loc[e, "start"]) < PNSN_T0 else "pnsn") for e in T.episode]
    main_t = T[T.src == T.date_src].drop(columns=["src", "date_src"])
    env_t = T[T.src == "env"].drop(columns=["src", "date_src"]).add_suffix("_env").rename(columns={"episode_env": "episode"})
    TR = W[["episode"]].merge(main_t, on="episode", how="left").merge(env_t, on="episode", how="left")
    TR.insert(0, "region", "SVI")
    TR.to_csv(f"{OUT}/tremor_ep.csv", index=False)
    # ------------------------------------------------------------------ acceptance vs Q1
    CL = pd.read_parquet(f"{Q1}/counts_long.parquet")
    q = CL[(CL.catalogue == "SVI") & (CL.variant == "once075") & (CL.patchdef == "group")]
    s3big = json.load(open(f"{ROOT}/svi_science/results/s3/build_info.json"))["big"]["SVI"]
    q = q[q.episode.isin(s3big)].assign(group=lambda x: x.patch.astype(int))[["group", "episode", "n"]]
    mine = PE[PE.episode.isin(s3big)][["group", "episode", "n_once"]]
    m = mine.merge(q, on=["group", "episode"], how="outer").fillna(0)
    RES["acc_counts_diff_cells"] = int((m.n_once != m.n).sum()); RES["acc_counts_total"] = [int(m.n_once.sum()), int(m.n.sum())]
    tq = pd.read_csv(f"{Q1}/tremor_ep.csv"); tq = tq[(tq.region == "SVI") & tq.episode.isin(s3big)].set_index("episode")
    tm_ = TR.set_index("episode").reindex(tq.index)
    RES["acc_tremor_max_dp05_km"] = float(np.abs(tm_.trem_p05 - tq.trem_p05).max())
    RES["acc_tremor_max_dp95_km"] = float(np.abs(tm_.trem_p95 - tq.trem_p95).max())
    RES["acc_pass"] = bool(RES["acc_counts_diff_cells"] == 0 and RES["acc_tremor_max_dp05_km"] <= 1 and RES["acc_tremor_max_dp95_km"] <= 1)
    log(f"acceptance: count cells differing {RES['acc_counts_diff_cells']} totals {RES['acc_counts_total']}; tremor max dp05 "
        f"{RES['acc_tremor_max_dp05_km']:.2f} dp95 {RES['acc_tremor_max_dp95_km']:.2f} km -> {'PASS' if RES['acc_pass'] else 'FAIL'}")
    # ------------------------------------------------------------------ tremor validation env vs PNSN, big ETS 2010-13
    b = W[W.big & (W.start >= "2010-01-06") & (W.start < "2014-01-01")].episode
    tv = TR.set_index("episode").loc[b]
    RES["env_vs_pnsn_big_2010_13"] = dict(episodes=list(b), dp05=(tv.trem_p05_env - tv.trem_p05).round(1).tolist(),
                                          dp95=(tv.trem_p95_env - tv.trem_p95).round(1).tolist(),
                                          dir_pnsn=tv.trem_dir.tolist(), dir_env=tv.trem_dir_env.tolist(),
                                          n_pnsn=tv.n_trem.tolist(), n_env=tv.n_trem_env.tolist())
    med = max(float(np.median(np.abs(tv.trem_p05_env - tv.trem_p05))), float(np.median(np.abs(tv.trem_p95_env - tv.trem_p95))))
    RES["env_validated"] = bool(med <= 15.0); RES["env_vs_pnsn_median_abs_km"] = med
    log(f"env vs PNSN extents (2010-13 big): {RES['env_vs_pnsn_big_2010_13']} -> validated {RES['env_validated']} (median |d| {med:.1f} km)")
    # summaries
    big = W[W.big].set_index("episode")
    s = PE[PE.episode.isin(big.index)].groupby("episode")[["n_once", "n_c1", "R", "C_eq"]].sum()
    RES["big_totals"] = s.round(1).reset_index().to_dict("records")
    log("big ETS totals:\n" + s.round(1).to_string())
    json.dump(RES, open(f"{OUT}/base_summary.json", "w"), indent=1, default=str)
    log("done")


if __name__ == "__main__":
    main()
