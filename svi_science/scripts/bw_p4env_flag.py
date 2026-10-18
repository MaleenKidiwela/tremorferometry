#!/usr/bin/env python
"""P4 quiet flag, attempt 2, step 2: envelope-CC tremor-day flag, calibration (once), flag file, validation (once),
sanity checks.  Rules: svi_science/results/backward/p4/envflag/ENVFLAG.md (written before computing).

  python bw_p4env_flag.py        (idempotent; calibration and validation each run once and are then frozen)
"""
import os, sys, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import bw_p4_common as P

EF = f"{P.P4}/envflag"
WIN = f"{EF}/win"
TAB = f"{EF}/tables"
FIG = f"{EF}/figures"
FLAG = f"{TAB}/quiet_flag_env.csv"
CAL0, CAL1 = "2010-01-10", "2011-12-31"
VAL0, VAL1 = "2012-01-01", "2013-12-31"
PAD = 3
CGRID = np.round(np.arange(0.30, 0.90 + 1e-9, 0.02), 2)
KGRID = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48]
NWIN = 575
AREA_KM = 30.0
SSETS = ("CN4", "all")
now = lambda: os.popen("TZ=America/Los_Angeles date").read().strip()


def load_days():
    """per (sset, series, date): coverage and the sorted per-window tremor scores m = min(coh, s3) of in-area windows."""
    fs = sorted(glob.glob(f"{WIN}/env_[0-9][0-9][0-9][0-9].parquet"))
    years = [int(os.path.basename(f)[4:8]) for f in fs]
    out = []
    for f in fs:
        A = pd.read_parquet(f, columns=["date", "series", "sset", "n_valid", "coh", "s3", "dfam_km"])
        A["use"] = A.n_valid >= 3
        A["m"] = np.where(A.use & (A.dfam_km <= AREA_KM), np.minimum(A.coh, A.s3), np.nan)
        for (ss, se, d), g in A.groupby(["sset", "series", "date"], sort=False):
            m = g.m.to_numpy(float); m = np.sort(m[np.isfinite(m)])
            out.append(dict(sset=ss, series=se, date=d, cov=float(g.use.sum()) / NWIN, m=m))
    return pd.DataFrame(out), years


def n_at(D, C):
    """N(d) at threshold C; NaN where coverage < 0.5."""
    n = np.array([len(m) - np.searchsorted(m, C, side="left") for m in D.m], float)
    return np.where(D["cov"].to_numpy() >= 0.5, n, np.nan)


def calendar():
    cal = pd.DataFrame(dict(date=[str(d.date()) for d in pd.date_range("2003-01-01", "2013-12-31")]))
    Q = P.pnsn_quiet("2003-01-01", "2013-12-31")[["date", "quiet", "n_tremor_day"]].rename(columns={"quiet": "pnsn_quiet"})
    cal = cal.merge(Q, on="date", how="left")
    cal["bostock_near"] = P.bostock_near(cal.date.to_numpy(), PAD)
    return cal


def flag(N, k, bnear):
    """N over the full calendar (NaN unknown) -> active, env_quiet, quiet_env3."""
    known = np.isfinite(N)
    act = known & (np.nan_to_num(N, nan=-1) >= k)
    near = np.convolve(act.astype(int), np.ones(2 * PAD + 1, int), mode="same") > 0
    eq = known & ~near
    return act, eq, eq & ~bnear


def sens_spec(pq, fq, m):
    pq = pq[m].astype(bool); fq = fq[m].astype(bool)
    sens = float((~fq[~pq]).mean()) if (~pq).any() else np.nan
    spec = float(fq[pq].mean()) if pq.any() else np.nan
    return sens, spec, int(pq.sum()), int((~pq).sum())


def auc(score, y):
    """Mann-Whitney AUC of score for y == True vs False (NaN scores dropped)."""
    s = np.asarray(score, float); y = np.asarray(y, bool); ok = np.isfinite(s)
    s, y = s[ok], y[ok]
    if y.all() or (~y).all():
        return np.nan
    from scipy.stats import rankdata
    r = rankdata(s); n1 = y.sum(); n0 = (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def max3(N):
    x = pd.Series(N).rolling(2 * PAD + 1, center=True, min_periods=1).max().to_numpy()
    return np.where(np.isfinite(N), x, np.nan)


def aucs(cal, N, lo, hi):
    m = ((cal.date >= lo) & (cal.date <= hi) & cal.pnsn_quiet.notna()).to_numpy()
    pq = cal.pnsn_quiet.to_numpy()[m].astype(bool); nt = cal.n_tremor_day.to_numpy()[m]
    n, n3 = N[m], max3(N)[m]
    big = (nt >= 50) | pq
    return dict(auc_N_pnsn_notquiet=auc(n, ~pq), auc_N_tremor_that_day=auc(n, nt > 0),
                auc_N_big_vs_quiet=auc(n[big], (nt >= 50)[big]), auc_max3N_pnsn_notquiet=auc(n3, ~pq))


def series_N(Dd, sset, series, cal, C):
    x = Dd[(Dd.sset == sset) & (Dd.series == series)].copy()
    x["N"] = n_at(x, C)
    return cal[["date"]].merge(x[["date", "N", "cov"]], on="date", how="left")


def calibrate(Dd, cal, sset, years):
    f = f"{TAB}/envflag_calibration_{sset}.json"
    if os.path.exists(f):
        return json.load(open(f))
    if not {2010, 2011} <= set(years):
        raise SystemExit("[envflag] calibration needs 2010 and 2011")
    m = ((cal.date >= CAL0) & (cal.date <= CAL1) & cal.pnsn_quiet.notna()).to_numpy()
    pq = cal.pnsn_quiet.to_numpy()
    bn = cal.bostock_near.to_numpy()
    rows = []
    for C in CGRID:
        N = series_N(Dd, sset, "data", cal, C).N.to_numpy(float)
        known = np.isfinite(N)
        for k in KGRID:
            act, eq, q3 = flag(N, k, bn)
            se, sp, npq, nnq = sens_spec(pq, q3, m)
            se_k, sp_k, _, _ = sens_spec(pq, q3, m & known)
            se_e, sp_e, _, _ = sens_spec(pq, eq, m)
            rows.append(dict(C=C, k=k, sens=se, spec=sp, J=se + sp - 1, sens_Nknown=se_k, spec_Nknown=sp_k,
                             sens_envonly=se_e, spec_envonly=sp_e, n_pnsn_quiet=npq, n_pnsn_notquiet=nnq,
                             n_flag_quiet=int(q3[m].sum())))
    R = pd.DataFrame(rows); R.to_csv(f"{TAB}/envflag_calibration_curve_{sset}.csv", index=False)
    b = R.sort_values(["J", "C", "k"], ascending=[False, False, False]).iloc[0]
    feas = R[(R.spec >= 0.9) & (R.sens >= 0.8)]
    s90 = R[R.spec >= 0.9].sort_values(["sens", "C", "k"], ascending=[False, False, False])
    N = series_N(Dd, sset, "data", cal, float(b.C)).N.to_numpy(float)
    out = dict(sset=sset, C=float(b.C), k=int(b.k), sens_cal=float(b.sens), spec_cal=float(b.spec), J=float(b.J),
               sens_cal_Nknown=float(b.sens_Nknown), spec_cal_Nknown=float(b.spec_Nknown),
               sens_cal_envonly=float(b.sens_envonly), spec_cal_envonly=float(b.spec_envonly),
               bars_feasible_on_grid=bool(len(feas)), n_feasible_points=int(len(feas)),
               best_sens_at_spec_ge_0p9=(dict(C=float(s90.iloc[0].C), k=int(s90.iloc[0].k), sens=float(s90.iloc[0].sens),
                                             spec=float(s90.iloc[0].spec)) if len(s90) else None),
               n_cal_days=int(m.sum()), n_cal_days_N_known=int((m & np.isfinite(N)).sum()),
               **{k_: v for k_, v in aucs(cal, N, CAL0, CAL1).items()},
               frozen_pt=now(), rule="(C,k) maximising sens+spec of quiet_env3 vs PNSN quiet on 2010-01-10..2011-12-31;"
                                      " ties -> larger C, then larger k")
    json.dump(out, open(f, "w"), indent=1)
    P.log(f"[envflag {sset}] calibrated C={b.C} k={b.k}: sens {b.sens:.3f} spec {b.spec:.3f} J {b.J:.3f}; "
          f"AUC(N, not quiet) {out['auc_N_pnsn_notquiet']:.3f}")
    return out


def sigma_p_days():
    f = f"{TAB}/sigmaP_days_2012_2013.csv"
    if os.path.exists(f):
        return pd.read_csv(f)
    g = f"{P.TAB}/qflag_sigmaP_days_2012_2013.csv"          # same function, written by the energy-flag validation
    if os.path.exists(g):
        S = pd.read_csv(g)
    else:
        import bw_p4_qflag as QF
        S = QF.sigma_p([2012, 2013])
    S.to_csv(f, index=False)
    return S


def validate(F, cal, sset, cj):
    f = f"{TAB}/envflag_validation_{sset}.json"
    if os.path.exists(f):
        return json.load(open(f))
    m = ((F.date >= VAL0) & (F.date <= VAL1) & F.pnsn_quiet.notna()).to_numpy()
    q = F[f"quiet_env3_{sset}"].to_numpy()
    N = F[f"N_{sset}"].to_numpy(float)
    se, sp, npq, nnq = sens_spec(F.pnsn_quiet.to_numpy(), q, m)
    se_k, sp_k, _, _ = sens_spec(F.pnsn_quiet.to_numpy(), q, m & np.isfinite(N))
    S = sigma_p_days()
    V = F[(F.date >= VAL0) & (F.date <= VAL1)][["date", "pnsn_quiet", f"quiet_env3_{sset}"]].merge(S, on="date", how="inner")
    V["month"] = V.date.str[:7]

    def rate(x):
        n = x.n_fd_data.sum()
        return (x.sp_data.sum() - x.sp_c1.sum()) / n if n else np.nan
    Vp = V[V.pnsn_quiet == True]; Vf = V[V[f"quiet_env3_{sset}"] == True]                     # noqa: E712
    r_p, r_f = rate(Vp), rate(Vf)
    months = []
    for mo, g in V.groupby("month"):
        gp, gf = g[g.pnsn_quiet == True], g[g[f"quiet_env3_{sset}"] == True]                  # noqa: E712
        months.append(dict(month=mo, n_days_pnsn=len(gp), n_days_env=len(gf), rate_pnsn=rate(gp), rate_env=rate(gf),
                           ratio=rate(gf) / rate(gp) if len(gp) and rate(gp) else np.nan))
    Mo = pd.DataFrame(months); Mo.to_csv(f"{TAB}/envflag_validation_monthly_{sset}.csv", index=False)
    el = Mo[(Mo.n_days_pnsn >= 5) & (Mo.n_days_env >= 5)]
    frac_ok = float(((el.ratio >= 0.9) & (el.ratio <= 1.1)).mean()) if len(el) else np.nan
    pooled = r_f / r_p if r_p else np.nan
    ok = dict(spec=sp >= 0.9, sens=se >= 0.8, sigmaP_pooled=bool(np.isfinite(pooled) and 0.9 <= pooled <= 1.1),
              sigmaP_months=bool(np.isfinite(frac_ok) and frac_ok >= 0.8))
    out = dict(sset=sset, C=cj["C"], k=cj["k"], sens=se, spec=sp, sens_Nknown=se_k, spec_Nknown=sp_k,
               n_pnsn_quiet=npq, n_pnsn_notquiet=nnq, n_flag_quiet=int(q[m].sum()),
               n_flag_quiet_but_pnsn_active=int((q[m] & ~F.pnsn_quiet.to_numpy()[m].astype(bool)).sum()),
               rate_pnsn=r_p, rate_env=r_f, pooled_ratio=pooled, n_months_eligible=int(len(el)),
               frac_months_within_10pct=frac_ok, checks=ok, pass_=bool(all(ok.values())),
               **aucs(F.rename(columns={}), N, VAL0, VAL1), written_pt=now())
    json.dump(out, open(f, "w"), indent=1, default=float)
    P.log(f"[envflag {sset}] validation 2012-13: sens {se:.3f} spec {sp:.3f} sigmaP ratio {pooled:.3f} months-ok {frac_ok}"
          f" -> {'PASS' if out['pass_'] else 'FAIL'}")
    return out


def sanity(F, years):
    """N on Bostock ETS window days vs days > 30 d from every window, per year; control beside it."""
    E = pd.read_csv(P.EPISODES)
    d = pd.to_datetime(F.date)
    inwin = np.zeros(len(F), bool)
    for r in E.itertuples():
        inwin |= ((d >= pd.Timestamp(r.start)) & (d <= pd.Timestamp(r.end))).to_numpy()
    far = ~P.bostock_near(F.date.to_numpy(), 30)
    rows = []
    for y in years:
        my = F.date.str[:4] == str(y)
        for lab, mm in (("ETS window", inwin), ("> 30 d from windows", far)):
            for ss in SSETS:
                x = F[my & mm]
                n = x[f"N_{ss}"]; c = x[f"Nc1_{ss}"]
                rows.append(dict(year=y, days=lab, sset=ss, n_days=int(n.notna().sum()), median_N=float(n.median()),
                                 mean_N=float(n.mean()), frac_active=float(x[f"active_{ss}"][n.notna()].mean())
                                 if n.notna().any() else np.nan, median_N_control=float(c.median()),
                                 mean_N_control=float(c.mean())))
    R = pd.DataFrame(rows); R.to_csv(f"{TAB}/sanity_ets_vs_far.csv", index=False)
    return R


def ets_map(ep0, ep1, tag, cj):
    """located tremor windows (main set, data) in one ETS vs PNSN tremor epicentres of the same days."""
    y = int(ep0[:4])
    f = f"{WIN}/env_{y}.parquet"
    if not os.path.exists(f):
        return None
    A = pd.read_parquet(f)
    A = A[(A.series == "data") & (A.date >= ep0) & (A.date <= ep1)]
    tr = pd.read_csv(P.C.PNSN, usecols=["time", "lat", "lon"]); tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    tr = tr[(tr.time >= ep0) & (tr.time < pd.Timestamp(ep1) + pd.Timedelta("1D"))]
    tr = tr[tr.lat.between(47.5, 49.6) & tr.lon.between(-125.2, -122.2)]
    tr["date"] = tr.time.dt.strftime("%Y-%m-%d")
    pm = tr.groupby("date")[["lat", "lon"]].median()
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from bw_p4env_detect import epi_km
    Fm = pd.read_csv("/home/jovyan/tremorferometry/lfe_catalog_repair/results/e5r_families_2012.csv")
    St = pd.read_csv("/home/jovyan/tremorferometry/svi_science/results/backward_scope/inventory_stations.csv").drop_duplicates("sta").set_index("sta")
    out = {}
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    for ax, ss in zip(axs, SSETS):
        C = cj[ss]["C"]
        w = A[(A.sset == ss) & (A.n_valid >= 3) & (A.coh >= C) & (A.s3 >= C) & (A.dfam_km <= AREA_KM)]
        wa = A[(A.sset == ss) & (A.n_valid >= 3) & (A.coh >= C) & (A.s3 >= C)]
        ax.scatter(tr.lon, tr.lat, s=2, c="0.6", alpha=0.4, label=f"PNSN tremor ({len(tr)})")
        sc = ax.scatter(wa.lon + np.random.default_rng(0).uniform(-.02, .02, len(wa)),
                        wa.lat + np.random.default_rng(1).uniform(-.02, .02, len(wa)), s=5,
                        c=pd.to_datetime(wa.date).dt.dayofyear, cmap="viridis", label=f"env-CC windows coh,s3 ≥ {C} ({len(wa)})")
        ax.scatter(Fm.lon, Fm.lat, s=8, marker="x", c="r", label="SVI families")
        stl = A[A.sset == ss].stas.iloc[0].split(",") if len(A[A.sset == ss]) else []
        ax.scatter([St.loc[s, "lon"] for s in stl], [St.loc[s, "lat"] for s in stl], marker="^", s=50, c="k", label="stations")
        ax.set_title(f"{tag}: {ss} set"); ax.set_xlim(-125.15, -122.25); ax.set_ylim(47.45, 49.6); ax.legend(fontsize=7, loc="lower left")
        ax.set_xlabel("lon")
        plt.colorbar(sc, ax=ax, label="day of year")
        dd = w.merge(pm, left_on="date", right_index=True, suffixes=("", "_p"))
        dist = epi_km(dd.lat.to_numpy(float), dd.lon.to_numpy(float), dd.lat_p.to_numpy(float), dd.lon_p.to_numpy(float))
        out[ss] = dict(n_windows_in_area=int(len(w)), n_windows_any_location=int(len(wa)),
                       median_km_to_day_pnsn_median=float(np.median(dist)) if len(dist) else np.nan,
                       p80_km=float(np.quantile(dist, .8)) if len(dist) else np.nan,
                       median_dlat_km=float(np.median((dd.lat - dd.lat_p) * 111.2)) if len(dd) else np.nan,
                       median_dlon_km=float(np.median((dd.lon - dd.lon_p) * 111.2 * np.cos(np.radians(48.5)))) if len(dd) else np.nan)
    axs[0].set_ylabel("lat")
    fig.tight_layout(); fig.savefig(f"{FIG}/ets_locations_{tag}.png", dpi=120); plt.close(fig)
    json.dump(dict(epoch=[ep0, ep1], **out), open(f"{TAB}/ets_locations_{tag}.json", "w"), indent=1)
    return out


def main():
    os.makedirs(TAB, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    Dd, years = load_days()
    cal = calendar()
    cj = {ss: calibrate(Dd, cal, ss, years) for ss in SSETS}
    F = cal.copy()
    bn = F.bostock_near.to_numpy()
    for ss in SSETS:
        C, k = cj[ss]["C"], cj[ss]["k"]
        x = series_N(Dd, ss, "data", cal, C); xc = series_N(Dd, ss, "c1", cal, C)
        N = x.N.to_numpy(float)
        act, eq, q3 = flag(N, k, bn)
        F[f"N_{ss}"] = N; F[f"cov_{ss}"] = x["cov"].to_numpy(); F[f"Nc1_{ss}"] = xc.N.to_numpy(float)
        F[f"active_{ss}"] = act; F[f"quiet_env_{ss}"] = eq; F[f"quiet_env3_{ss}"] = q3
        F[f"C_{ss}"] = C; F[f"k_{ss}"] = k
    F["year_computed"] = F.date.str[:4].astype(int).isin(years)
    F["quiet_bostock3"] = ~F.bostock_near
    F["quiet_bostock30"] = ~P.bostock_near(F.date.to_numpy(), 30)
    E = f"{P.TAB}/quiet_flag_cn4.csv"
    if os.path.exists(E):
        F = F.merge(pd.read_csv(E, usecols=["date", "Z", "quiet_opt3"]).rename(columns={"Z": "energy_Z", "quiet_opt3": "quiet_energy_opt3"}),
                    on="date", how="left")
    st = [pd.read_csv(f).rename(columns={"stations": "stations_cached", "n_sta": "n_sta_cached"})
          for f in sorted(glob.glob(f"{WIN}/env_*_stations.csv"))]
    if st:
        F = F.merge(pd.concat(st), on="date", how="left")
    F.to_csv(FLAG + ".tmp", index=False); os.replace(FLAG + ".tmp", FLAG)
    P.log(f"[envflag] flag written, years {years}; CN4 quiet_env3 days {int(F.quiet_env3_CN4.sum())}, "
          f"all {int(F.quiet_env3_all.sum())}")
    sanity(F, years)
    for ep0, ep1, tag in (("2010-08-16", "2010-09-08", "2010-08"), ("2012-08-30", "2012-09-25", "2012-09")):
        if not os.path.exists(f"{TAB}/ets_locations_{tag}.json") and int(ep0[:4]) in years:
            ets_map(ep0, ep1, tag, cj)
    ready = {2012, 2013} <= set(years) and all(os.path.exists(f"{P.BW}/{m}_DONE") for m in
                                                ("P2_2012", "P2_2013", "P5F_2012", "P5F_2013"))
    if ready:
        for ss in SSETS:
            validate(F, cal, ss, cj[ss])


if __name__ == "__main__":
    main()
