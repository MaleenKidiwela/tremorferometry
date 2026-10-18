#!/usr/bin/env python
"""SH_CHECKS (2026-10-01) CHECK 2: outside check of our 2025-26 Shelly-family rates against PNSN tremor.
Box chosen BEFORE looking at rates: circle of radius 30 km around the centroid of Shelly's 27 family locations
(results/families_south.csv, source Shelly2025 = his released Table S1 locations used by v3b); also 15 and 60 km.
  python scripts/sh_checks_pnsn.py  -> results/sh_checks/pnsn_*.csv + figure
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
WS = "/home/jovyan/tremorferometry/lfe_south"; RES = f"{WS}/results"; OUT = f"{RES}/sh_checks"
PN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
RADII = [30, 15, 60]          # km; 30 = pre-chosen primary
M0, M1 = "2018-01", "2026-05"
FIT1 = "2023-11"; JOIN0, JOIN1 = "2023-12", "2024-11"; EXT0 = "2024-12"; LASTPN = pd.Timestamp("2026-05-27 23:59")


def hav(la1, lo1, la2, lo2):
    R = 6371.0; p = np.radians
    a = np.sin(p(la2 - la1) / 2) ** 2 + np.cos(p(la1)) * np.cos(p(la2)) * np.sin(p(lo2 - lo1) / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def main():
    F = pd.read_csv(f"{RES}/families_south.csv"); F = F[F.source == "Shelly2025"]
    la0, lo0 = F.lat.mean(), F.lon.mean()
    months = pd.period_range(M0, M1, freq="M").astype(str)
    T = pd.DataFrame(index=months); T.index.name = "month"
    P = pd.read_csv(PN, parse_dates=["time"]); P = P[(P.time >= M0) & (P.time <= LASTPN)]
    P["d"] = hav(la0, lo0, P.lat.values, P.lon.values); P["month"] = P.time.dt.strftime("%Y-%m")
    for r in RADII:
        g = P[P.d <= r].groupby("month")
        T[f"pnsn_{r}km"] = g.size().reindex(months).fillna(0)
        T[f"pnsn_{r}km_nsta_med"] = g.num_stas.median().reindex(months)
    T["pnsn_days_in_month"] = [pd.Period(m).days_in_month for m in months]
    T.loc["2026-05", "pnsn_days_in_month"] = 27            # catalogue ends 2026-05-27
    from sh_v3_common import his_catalogue
    S = his_catalogue(); S["month"] = S.t.dt.strftime("%Y-%m")
    T["his_n"] = S.groupby("month").size().reindex(months)
    T.loc[T.index > "2024-11", "his_n"] = np.nan
    M = pd.read_csv(f"{RES}/sh_conf/monthly.csv"); M = M[M.fam.astype(str) == "all"].set_index("month")
    for c in ["days", "S1_rate", "sumP_same_rate", "sumP_any_rate", "cand_rate", "ctl_sumP_same_rate"]:
        T["ours_" + c] = M[c].reindex(months)
    fwd = T.index >= JOIN0
    for c in ["S1_rate", "sumP_same_rate", "sumP_any_rate"]:   # monthly count from daily rate (forward run is continuous)
        T.loc[~fwd, "ours_" + c] = np.nan
    T["ndays_m"] = [pd.Period(m).days_in_month for m in months]
    T["ours_S1_n"] = T.ours_S1_rate * T.ndays_m; T["ours_sumP_n"] = T.ours_sumP_same_rate * T.ndays_m
    T.to_csv(f"{OUT}/pnsn_monthly.csv")
    rows = []; lines = []
    for r in RADII:
        x = T[f"pnsn_{r}km"] * T.ndays_m / T.pnsn_days_in_month    # full-month equivalent
        T[f"pnsn_{r}km_fm"] = x
        fit = (T.index <= FIT1) & T.his_n.notna()
        xf, yf = x[fit].values, T.his_n[fit].values
        b, a = np.polyfit(xf, yf, 1); k0 = yf.sum() / xf.sum()
        r_fit = np.corrcoef(xf, yf)[0, 1]
        from scipy.stats import spearmanr
        rs_fit = spearmanr(xf, yf).correlation
        res_fit = yf - (a + b * xf); sd = res_fit.std(ddof=2)
        # join-year: his vs ours (scale ours to his units)
        j = (T.index >= JOIN0) & (T.index <= JOIN1)
        for ours in ["ours_sumP_n", "ours_S1_n"]:
            sc = T.his_n[j].sum() / T[ours][j].sum()
            y_ours = T[ours] * sc
            pred = a + b * x
            def blk(msk, lab, y):
                yy, pp, xx = y[msk].values, pred[msk].values, x[msk].values
                ok = np.isfinite(yy)
                yy, pp, xx = yy[ok], pp[ok], xx[ok]
                return dict(radius_km=r, ours=ours, period=lab, n_months=int(ok.sum()), obs_sum=yy.sum(), pred_sum=pp.sum(),
                            obs_over_pred=yy.sum() / pp.sum(), resid_mean=(yy - pp).mean(), resid_sd=(yy - pp).std(ddof=1) if ok.sum() > 1 else np.nan,
                            resid_mean_in_fit_sd=(yy - pp).mean() / sd, r_pearson=np.corrcoef(xx, yy)[0, 1] if ok.sum() > 2 else np.nan,
                            r_spearman=spearmanr(xx, yy).correlation if ok.sum() > 2 else np.nan, pnsn_sum=xx.sum())
            rows.append(blk(fit, "fit 2018-01..2023-11 (his)", T.his_n))
            rows.append(blk(j, "join 2023-12..2024-11 (his)", T.his_n))
            rows.append(blk(j, "join 2023-12..2024-11 (ours scaled)", y_ours))
            rows.append(blk(T.index >= EXT0, "ext 2024-12..2026-05 (ours scaled)", y_ours))
            rows.append(blk((T.index >= "2025-01") & (T.index <= "2025-12"), "2025 (ours scaled)", y_ours))
            rows.append(blk((T.index >= "2026-01") & (T.index <= "2026-05"), "2026-01..05 (ours scaled)", y_ours))
        lines.append(dict(radius_km=r, slope=b, intercept=a, ratio_his_over_pnsn=k0, r_fit=r_fit, rs_fit=rs_fit, resid_sd_fit=sd,
                          pnsn_2024=x[(T.index >= "2024-01") & (T.index <= "2024-12")].sum(),
                          pnsn_2025=x[(T.index >= "2025-01") & (T.index <= "2025-12")].sum(),
                          pnsn_2018_2023_per_yr=x[(T.index <= "2023-12")].sum() / 6,
                          pnsn_2026_jan_may=x[(T.index >= "2026-01") & (T.index <= "2026-05")].sum(),
                          pnsn_2024_jan_may=x[(T.index >= "2024-01") & (T.index <= "2024-05")].sum(),
                          pnsn_2025_jan_may=x[(T.index >= "2025-01") & (T.index <= "2025-05")].sum()))
    R = pd.DataFrame(rows); Lf = pd.DataFrame(lines)
    R.to_csv(f"{OUT}/pnsn_relation_tests.csv", index=False); Lf.to_csv(f"{OUT}/pnsn_fit.csv", index=False)
    T.to_csv(f"{OUT}/pnsn_monthly.csv")
    # year ratios
    yr = pd.DataFrame({"year": [m[:4] for m in T.index]}, index=T.index)
    Y = pd.concat([T, yr], axis=1)
    agg = {c: "sum" for c in ["his_n", "ours_S1_n", "ours_sumP_n"] + [f"pnsn_{r}km_fm" for r in RADII]}
    G = Y.groupby("year").agg(agg); G["n_months_ours"] = Y.groupby("year").ours_S1_n.count()
    G.to_csv(f"{OUT}/pnsn_yearly.csv")
    pd.set_option("display.width", 250)
    print(Lf.round(3).to_string()); print(R.round(3).to_string()); print(G.round(0).to_string())
    fig(T)


def fig(T):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    t = pd.PeriodIndex(T.index, freq="M").to_timestamp()
    f, ax = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    j = (T.index >= JOIN0) & (T.index <= JOIN1)
    sc = T.his_n[j].sum() / T.ours_sumP_n[j].sum(); sc1 = T.his_n[j].sum() / T.ours_S1_n[j].sum()
    ax[0].plot(t, T.pnsn_30km_fm, color="#444", lw=1.2, label="PNSN tremor windows ≤30 km (full-month equiv.)")
    ax[0].plot(t, T.pnsn_15km_fm, color="#999", lw=0.8, label="≤15 km"); ax[0].plot(t, T.pnsn_60km_fm / 3, color="#bbb", lw=0.8, ls="--", label="≤60 km ÷3")
    ax[0].set_ylabel("PNSN windows / month"); ax[0].legend(fontsize=8, loc="upper left")
    ax[1].plot(t, T.his_n, color="#1f77b4", lw=1.2, label="Shelly (2025) catalogue, events/month")
    ax[1].plot(t, T.ours_sumP_n * sc, color="#d62728", lw=1.2, label=f"ours ΣP(same) × {sc:.2f} (join-year scale)")
    ax[1].plot(t, T.ours_S1_n * sc1, color="#ff9896", lw=1.0, label=f"ours S1 9× count × {sc1:.2f}")
    ax[1].set_ylabel("events / month"); ax[1].legend(fontsize=8, loc="upper left")
    L = pd.read_csv(f"{OUT}/pnsn_fit.csv").set_index("radius_km").loc[30]
    pred = L.intercept + L.slope * T.pnsn_30km_fm
    ax[2].plot(t, T.his_n / pred, color="#1f77b4", lw=1, label="his / predicted-from-PNSN")
    ax[2].plot(t, T.ours_sumP_n * sc / pred, color="#d62728", lw=1, label="ours ΣP scaled / predicted-from-PNSN")
    ax[2].axhline(1, color="k", lw=0.6); ax[2].set_yscale("log"); ax[2].set_ylabel("obs / pred (30 km)"); ax[2].legend(fontsize=8)
    for a in ax:
        a.axvspan(pd.Timestamp("2023-12-01"), pd.Timestamp("2024-12-01"), color="0.9", zorder=0)
        a.axvline(pd.Timestamp("2025-09-17"), color="0.5", ls=":", lw=0.8); a.axvline(pd.Timestamp("2026-03-12"), color="0.5", ls=":", lw=0.8)
    ax[0].set_title("Shelly families vs PNSN tremor (grey band = join year; dotted = NC.KMR sensor swaps)", fontsize=10)
    f.tight_layout(); f.savefig(f"{OUT}/pnsn_vs_rates.png", dpi=130)
    f, ax = plt.subplots(1, 3, figsize=(13, 4))
    for a, r in zip(ax, RADII):
        x = T[f"pnsn_{r}km_fm"]; fit = (T.index <= FIT1) & T.his_n.notna(); j2 = (T.index >= JOIN0) & (T.index <= JOIN1); e = T.index >= EXT0
        a.scatter(x[fit], T.his_n[fit], s=12, color="#1f77b4", label="his 2018-2023-11")
        a.scatter(x[j2], T.his_n[j2], s=14, color="#2ca02c", label="his join yr")
        a.scatter(x[e], T.ours_sumP_n[e] * sc, s=14, color="#d62728", label="ours ΣP scaled 2024-12+")
        L = pd.read_csv(f"{OUT}/pnsn_fit.csv").set_index("radius_km").loc[r]; xx = np.linspace(0, x.max(), 10)
        a.plot(xx, L.intercept + L.slope * xx, color="k", lw=0.8); a.set_xlabel(f"PNSN windows/month ≤{r} km"); a.set_title(f"{r} km, r_fit={L.r_fit:.2f}", fontsize=9)
    ax[0].set_ylabel("LFE events / month"); ax[0].legend(fontsize=7)
    f.tight_layout(); f.savefig(f"{OUT}/pnsn_scatter.png", dpi=130)


if __name__ == "__main__":
    main()
