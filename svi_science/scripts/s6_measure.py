"""S6 step 2: amplitudes, S-window check, best stations, relative size dM, station-epoch correction, completeness.
Rules: results/s6/PLAN_S6.md. Input: results/s6/env/env_*.npz (s6_extract.py) + targets.parquet + s0 analysis table.
Output: results/s6/sizes.parquet (one row per measured detection, all kept), station_steps.csv, family_stations.csv,
family_sizes.csv, measure_summary.json."""
import glob, json
import numpy as np, pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s6"
STAS = ["CN.NLLB", "CN.PFB", "CN.PGC", "CN.SNB", "CN.VGZ", "PB.B001", "PB.B003", "PB.B004", "PB.B005", "PB.B006", "PB.B007",
        "PB.B009", "PB.B010", "PB.B011", "PB.B926", "UW.LRIV"]
SB = slice(145, 185)      # S window -0.5..+3.5 s about expected S (S at 15.0 s, 0.1 s bins)
NBN = slice(5, 100)       # noise -14.5..-5 s
CHANGES = [("CN.PGC", "2015-01-12"), ("CN.NLLB", "2017-02-23"), ("CN.VGZ", "2017-03-03"), ("CN.SNB", "2017-06-08"),
           ("CN.PGC", "2017-08-10"), ("UW.LRIV", "2019-05-30")]


def load():
    T = pd.read_parquet(f"{OUT}/targets.parquet")
    A = pd.read_parquet(f"{ROOT}/svi_science/results/s0/analysis_table.parquet")
    A = A[A.series == "data"]
    rows = []; envs = []
    T["head"] = (T.P_same >= 0.75) & T.keep_xfam8 & (T.kind == "event")
    head = set(T.det_id[T["head"]])
    for f in sorted(glob.glob(f"{OUT}/env/env_*.npz")):
        z = np.load(f)
        e = z["env"]
        with np.errstate(invalid="ignore"):
            aS = np.sqrt(np.nanmean(e[:, SB], axis=1)); aN = np.sqrt(np.nanmean(e[:, NBN], axis=1))
            fS = np.isfinite(e[:, SB]).mean(axis=1)
        d = pd.DataFrame(dict(det_id=z["det_id"], sta=np.array(STAS)[z["sta"]], A_S=aS, A_N=aN, fracS=fS, ncomp=z["ncomp"],
                              offset_source=z["offset_source"], resp_epoch=z["resp_epoch"]))
        rows.append(d)
        h = np.isin(z["det_id"], list(head))
        if h.any():
            with np.errstate(invalid="ignore", divide="ignore"):
                envs.append((z["det_id"][h], np.array(STAS)[z["sta"][h]], (e[h] / np.nanmean(e[h][:, NBN], axis=1, keepdims=True)).astype(np.float32)))
    R = pd.concat(rows, ignore_index=True)
    R = R.merge(T[["det_id", "fam", "time_utc", "P_same", "keep_xfam8", "net_mad", "kind", "head"]], on="det_id")
    return T, R, envs


def window_check(envs, R):
    """median normalized envelope per family x station over headline detections; peak position"""
    det = np.concatenate([x[0] for x in envs]); sta = np.concatenate([x[1] for x in envs]); E = np.vstack([x[2] for x in envs])
    fam = pd.Series(R.drop_duplicates("det_id").set_index("det_id").fam).reindex(det).values
    df = pd.DataFrame(dict(fam=fam, sta=sta)); out = []; med = {}
    for (f, s), ix in df.groupby(["fam", "sta"]).groups.items():
        m = np.nanmedian(E[np.asarray(ix)], axis=0)
        seg = m[100:300]
        if not np.isfinite(seg).any(): continue
        pk = 100 + int(np.nanargmax(seg)); med[(f, s)] = m
        out.append(dict(fam=f, sta=s, n_env=len(ix), peak_s=(pk + 0.5) * 0.1 - 15.0, peak_ratio=float(m[pk]),
                        s_window_ratio=float(np.nanmean(m[SB]))))
    W = pd.DataFrame(out); W["window_ok"] = W.peak_s.between(-1.5, 4.5)
    np.savez_compressed(f"{OUT}/family_station_median_env.npz", keys=np.array([f"{f}|{s}" for f, s in med]),
                        env=np.vstack(list(med.values())))
    return W


def epoch_of(sta, t):
    ch = [pd.Timestamp(d) for s, d in CHANGES if s == sta]
    return int(np.sum([t >= c for c in ch]))


def pipeline(R, tag):
    """R: one row per detection x station with logS, logN, sn, fam, sta, time_utc, epoch, head, window_ok."""
    H = R[R["head"]]
    n_head_f = H.groupby("fam").det_id.nunique()
    FS = H.groupby(["fam", "sta"]).agg(med_sn=("sn", "median"), n=("logS", "count")).reset_index()
    FS["presence"] = FS.n / FS.fam.map(n_head_f)
    wok = R.drop_duplicates(["fam", "sta"]).set_index(["fam", "sta"]).window_ok
    FS["eligible"] = [bool(wok.get((f, s_), False)) for f, s_ in zip(FS.fam, FS.sta)]
    best = {}; weak = {}
    for f, g in FS.groupby("fam"):
        g = g[g.eligible].sort_values("med_sn", ascending=False)
        q = g[(g.med_sn >= 2) & (g.presence >= 0.5)]
        if len(q) >= 3: best[f] = q.sta.head(5).tolist(); weak[f] = False
        else:
            best[f] = g[g.presence >= 0.5].sta.head(3).tolist() or g.sta.head(3).tolist(); weak[f] = True
    FS["best"] = [s_ in best.get(f, []) for f, s_ in zip(FS.fam, FS.sta)]
    FS["weak_family"] = FS.fam.map(weak)
    R = R.copy(); R["best"] = [s_ in best.get(f, []) for f, s_ in zip(R.fam, R.sta)]
    B = R[R.best & R.logS.notna()].copy()

    def dm(B, col):
        ref = B[B["head"]].groupby(["fam", "sta"])[col].median().rename("ref")
        B = B.join(ref, on=["fam", "sta"]); B["r"] = B[col] - B["ref"]
        g = B.groupby("det_id").r
        out = pd.DataFrame(dict(dM=g.median(), n_used=g.count()))
        out.loc[out.n_used < 2, "dM"] = np.nan
        return B, out

    B1, D1 = dm(B, "logS")
    B1 = B1.sort_values("det_id", kind="stable")
    res = []
    for det, g in B1.groupby("det_id", sort=True).r:
        v = g.values
        res.append(np.full(len(v), np.nan) if len(v) < 3 else np.array([v[i] - np.median(np.delete(v, i)) for i in range(len(v))]))
    B1["resid"] = np.concatenate(res)
    Ep = B1[B1["head"]].groupby(["sta", "epoch"]).resid.agg(["median", "count"]).reset_index().rename(columns={"median": "e_res"})
    Ep.loc[Ep["count"] < 20, "e_res"] = 0.0                      # too few rows to estimate a station-epoch term
    steps = []
    for s_, d in CHANGES:
        c = pd.Timestamp(d); e_after = epoch_of(s_, c)
        allst = R[(R.sta == s_) & R["head"] & R.logS.notna()]
        refall = allst.groupby("fam").logS.median()
        a0 = allst[(allst.time_utc >= c - pd.Timedelta("30D")) & (allst.time_utc < c)]
        a1 = allst[(allst.time_utc >= c) & (allst.time_utc < c + pd.Timedelta("30D"))]
        r0 = (a0.logS - a0.fam.map(refall)).median(); r1 = (a1.logS - a1.fam.map(refall)).median()
        e0 = Ep[(Ep.sta == s_) & (Ep.epoch == e_after - 1)]; e1 = Ep[(Ep.sta == s_) & (Ep.epoch == e_after)]
        steps.append(dict(station=s_, change=d, n_month_before=len(a0), n_month_after=len(a1), month_median_before=r0,
                          month_median_after=r1, month_step=r1 - r0,
                          resid_epoch_before=float(e0.e_res.iloc[0]) if len(e0) else np.nan,
                          resid_epoch_after=float(e1.e_res.iloc[0]) if len(e1) else np.nan,
                          n_resid_before=int(e0["count"].iloc[0]) if len(e0) else 0, n_resid_after=int(e1["count"].iloc[0]) if len(e1) else 0))
    S = pd.DataFrame(steps); S["resid_step"] = S.resid_epoch_after - S.resid_epoch_before
    S["incomplete_gt_0p05"] = S.resid_step.abs() > 0.05
    S.to_csv(f"{OUT}/station_steps_{tag}.csv", index=False); Ep.to_csv(f"{OUT}/station_epoch_residuals_{tag}.csv", index=False)
    B = B.merge(Ep[["sta", "epoch", "e_res"]], on=["sta", "epoch"], how="left"); B["e_res"] = B.e_res.fillna(0.0)
    B["logS_c"] = B.logS - B.e_res; B["logN_c"] = B.logN - B.e_res
    B2, D2 = dm(B, "logS_c"); BN, DN = dm(B, "logN_c")
    D = D2.join(D1.dM.rename("dM_senscorr_only")).join(DN.dM.rename("dN")).join(B.groupby("det_id").sn.median().rename("sn_med"))
    D["low_sn"] = D.sn_med < 2
    FS.to_csv(f"{OUT}/family_stations_{tag}.csv", index=False)
    return D, weak, B, S


def completeness(out, tag=""):
    out["z"] = out.dM - out.dN
    H = out[(out.kind == "event") & (out.P_same >= 0.75) & out.keep_xfam8 & out.z.notna()]
    def maxc(z):
        h, ed = np.histogram(z, bins=np.arange(np.floor(z.min() * 10) / 10, z.max() + 0.2, 0.1))
        return ed[np.argmax(h)] + 0.05 + 0.1
    mcz = {f: maxc(g.z.values) for f, g in H.groupby("fam") if len(g) >= 30}
    mcg = {k: maxc(g.z.values) for k, g in H.groupby("group")}
    out["Mc_z"] = [mcz.get(f, mcg.get(gp, np.nan)) for f, gp in zip(out.fam, out.group)]
    nfm = out[out.kind == "event"].groupby(["fam", "ym"]).dN.median().rename("dN_fm")
    out = out.join(nfm, on=["fam", "ym"])
    out["Mc_fm"] = out.Mc_z + out.dN_fm
    out["complete"] = out.dM >= out.Mc_fm
    return out


def main():
    T, R, envs = load()
    R["epoch"] = [epoch_of(s_, t) for s_, t in zip(R.sta, R.time_utc)]
    W = window_check(envs, R); W.to_csv(f"{OUT}/window_check.csv", index=False)
    R = R.merge(W[["fam", "sta", "window_ok"]], on=["fam", "sta"], how="left"); R["window_ok"] = R.window_ok.fillna(False)
    # RMS estimator (pre-registered)
    Rr = R.copy()
    Rr["logS"] = np.log10(Rr.A_S); Rr["logN"] = np.log10(Rr.A_N); Rr["sn"] = Rr.A_S / Rr.A_N
    Rr.loc[~np.isfinite(Rr.logS) | (Rr.fracS < 0.5), ["logS", "logN"]] = np.nan
    Drms, weak_rms, Brms, Srms = pipeline(Rr, "rms")
    # matched-filter estimator (Deviation 1, headline)
    M = pd.read_parquet(f"{OUT}/mf_station.parquet")
    Rm = M.merge(R[["det_id", "sta", "fam", "time_utc", "head", "epoch", "window_ok"]], on=["det_id", "sta"], how="left")
    Rm["neg"] = ~(Rm.a_mf > 0)
    Rm["logN"] = np.log10(Rm.sig_noise)
    Rm["logS"] = np.where(Rm.neg, Rm.logN, np.log10(Rm.a_mf.clip(lower=1e-30)))
    Rm["sn"] = Rm.a_mf / Rm.sig_noise
    Rm.loc[~np.isfinite(Rm.logN), ["logS", "logN"]] = np.nan
    Dmf, weak, B, Smf = pipeline(Rm, "mf")
    Dmf = Dmf.join(Rm[Rm.det_id.isin(Dmf.index)].groupby("det_id").cc_mf.median().rename("cc_mf_med"))
    out = T.set_index("det_id").join(Dmf, how="left")
    out = out.join(Drms[["dM", "dN", "n_used", "sn_med", "low_sn"]].add_suffix("_rms"), how="left")
    out["n_used"] = out.n_used.fillna(0).astype(int)
    A = pd.read_parquet(f"{ROOT}/svi_science/results/s0/analysis_table.parquet")
    A = A[A.series == "data"][["fam", "time_utc", "episode", "major", "in_frozen_major", "in_episode", "quiet", "n_tremor_day",
                               "n_trem30", "n_sta_day", "u_mad", "t_since_fam_first_exposure_d", "t_since_ep_start_d",
                               "flag_ep_2017_03", "flag_ep_2025_10"]]
    out = out.reset_index().merge(A, on=["fam", "time_utc"], how="left")
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")[["fid", "dep", "mag_med", "group", "along_strike_km", "downdip_km", "lat", "lon"]]
    out = out.merge(F.rename(columns={"fid": "fam"}), on="fam", how="left")
    out["weak_family"] = out.fam.map(weak).fillna(True)
    out["M_approx"] = out.mag_med + out.dM
    out["ym"] = out.time_utc.dt.strftime("%Y-%m")
    out = completeness(out)
    # completeness for the RMS column
    tmp = completeness(out[["fam", "group", "kind", "P_same", "keep_xfam8", "ym"]].assign(dM=out.dM_rms, dN=out.dN_rms))
    out["complete_rms"] = tmp.complete.values
    out.to_parquet(f"{OUT}/sizes.parquet")
    fam_abs(B, F)
    hd = (out.kind == "event") & (out.P_same >= 0.75) & out.keep_xfam8
    summ = dict(n_targets=len(T), n_head=int(hd.sum()), n_head_dM=int((hd & out.dM.notna()).sum()),
                n_head_dM_rms=int((hd & out.dM_rms.notna()).sum()), n_head_complete=int((hd & out.complete).sum()),
                n_head_low_sn=int((hd & (out.low_sn == True)).sum()), n_head_low_sn_rms=int((hd & (out.low_sn_rms == True)).sum()),
                family_station_window_fail=int((~W.window_ok).sum()), family_station_n=len(W),
                weak_families_mf=int(sum(weak.values())), weak_families_rms=int(sum(weak_rms.values())), n_families=len(weak),
                median_station_sn_mf_head=float(Rm[Rm["head"]].sn.median()), median_station_sn_rms_head=float(Rr[Rr["head"]].sn.median()),
                frac_negative_mf_head=float(Rm[Rm["head"]].neg.mean()))
    json.dump(summ, open(f"{OUT}/measure_summary.json", "w"), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float)); print(Smf.round(3).to_string()); print(Srms.round(3).to_string())


def fam_abs(B, F):
    # absolute family size from the coherent stack amplitude (m/s) at each station (matched-filter ratios are relative)
    M = pd.read_csv(f"{OUT}/family_station_stacks.csv")
    M = M[(M.stack_src == "headline")].copy(); M["y"] = np.log10(M.stack_rms)
    M = M[M.set_index(["fam", "sta"]).index.isin(B.set_index(["fam", "sta"]).index)]   # the family's best stations only
    co = pd.read_csv(f"{ROOT}/lfe_extension/results/cut/svi_ext/channel_epochs.csv")  # no coords; use station table
    import sys; sys.path.insert(0, f"{ROOT}/lfe_extension/scripts"); import c_common as C
    ep = C.epochs(); xy = ep.groupby("netsta")[["lat", "lon"]].first()
    Fi = F.set_index("fid")
    R_ = []
    for r in M.itertuples():
        f = Fi.loc[r.fam]; s = xy.loc[r.sta]
        epi = C.epi_km(float(f.lat), float(f.lon), float(s.lat), float(s.lon)); R_.append(np.hypot(epi, float(f.dep)))
    M["R"] = R_
    fams = sorted(M.fam.unique()); stas = sorted(M.sta.unique())
    X = np.zeros((len(M), len(fams) + len(stas) + 1))
    for i, r in enumerate(M.itertuples()):
        X[i, fams.index(r.fam)] = 1; X[i, len(fams) + stas.index(r.sta)] = 1; X[i, -1] = -np.log10(r.R)
    # constraint: station terms sum to zero
    c = np.zeros((1, X.shape[1])); c[0, len(fams):len(fams) + len(stas)] = 1
    Xa = np.vstack([X, c * 100]); ya = np.r_[M.y.values, 0]
    b, *_ = np.linalg.lstsq(Xa, ya, rcond=None)
    out = pd.DataFrame(dict(fam=fams, M_f=b[:len(fams)])); out["spreading_n"] = b[-1]
    out = out.merge(F.rename(columns={"fid": "fam"}), on="fam")
    out["M_f_rel"] = out.M_f - out.M_f.median()
    out.to_csv(f"{OUT}/family_sizes.csv", index=False)
    print("spreading n", b[-1], "resid rms", np.sqrt(np.mean((X @ b - M.y.values) ** 2)))


if __name__ == "__main__":
    main()
