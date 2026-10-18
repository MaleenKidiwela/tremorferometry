#!/usr/bin/env python
"""SH_KMRFIX task 4: per-family daily rates v1 vs v2, 2018 -> 2026-09-28; join-year scaling as before (each family's
ours/his ratio over 2023-12-01..2024-11-30, all forward days); monthly figure with the KMR dates; v1 vs v2 for
2025-09-17..2026-09-28; step checks at 2025-07-25, 2025-08-26 and 2025-09-17.
Series: his = Shelly (2025) catalogue (to 2024-11-30); S1_v1 = results/sh_ext2026/catalogue_S1_9x; S1_v2 = catalogue_S1_9x_v2;
sumP_v1 = original model on the original candidates; sumP_v2a = original model up to 2025-07-24, v2 model after;
sumP_v2b = v2 model on all forward days.  Outputs results/sh_kmrfix/t4_*, figures/."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sh_v3_common import his_catalogue
RES = "/home/jovyan/tremorferometry/lfe_south/results"; OUT = f"{RES}/sh_kmrfix"; EXT = f"{RES}/sh_ext2026"
os.makedirs(f"{OUT}/figures", exist_ok=True)
FAMS = [str(i) for i in range(1, 28)]
J0, J1 = "2023-12-01", "2024-11-30"; V2A = "2025-07-25"; SW1 = "2025-09-17"
KD = [("2025-07-25", "KMR STS-2 drift starts"), ("2025-08-03", "KMR HHE dead"), ("2025-08-26", "KMR no data"),
      ("2025-09-17", "KMR swap: 16-station run"), ("2026-03-12", "KMR 2nd swap")]


def daily():
    days = pd.read_csv(f"{EXT}/days_processed.csv").day.astype(str).str[:10]
    grid = pd.MultiIndex.from_product([sorted(days), FAMS], names=["date", "fam"])
    out = pd.DataFrame(index=grid)
    for name, p in [("S1_v1", f"{EXT}/cat_S1_9x_internal.parquet"), ("S1_v2", f"{OUT}/catalogue_S1_9x_v2.parquet")]:
        C = pd.read_parquet(p); C["date"] = C.t.dt.strftime("%Y-%m-%d"); C["fam"] = C.fam.astype(str)
        out[name] = C.groupby(["date", "fam"]).size().reindex(grid, fill_value=0)
    A = pd.read_parquet(f"{RES}/sh_conf/candidates_conf.parquet", columns=["run", "series", "date", "fam", "P_same", "P_any"])
    A = A[(A.run == "fwd") & (A.series == "data")]; A["fam"] = A.fam.astype(str)
    out["sumP_v1"] = A.groupby(["date", "fam"]).P_same.sum().reindex(grid, fill_value=0)
    out["sumPany_v1"] = A.groupby(["date", "fam"]).P_any.sum().reindex(grid, fill_value=0)
    F = pd.read_parquet(f"{OUT}/conf_v2_rows.parquet"); F = F[(F.run == "fwd") & (F.series == "data")]
    F["pa"] = np.where(F.date < V2A, F.P_same_FULL, F.P_same_v2)
    out["sumP_v2a"] = F.groupby(["date", "fam"]).pa.sum().reindex(grid, fill_value=0)
    out["sumP_v2b"] = F.groupby(["date", "fam"]).P_same_v2.sum().reindex(grid, fill_value=0)
    out["sumPany_v2b"] = F.groupby(["date", "fam"]).P_any_v2.sum().reindex(grid, fill_value=0)
    out = out.reset_index()
    S = his_catalogue(); S["date"] = S.t.dt.strftime("%Y-%m-%d")
    S = S[S.date <= "2024-11-30"]
    H = S.groupby(["date", "fam"]).size().rename("his_n").reset_index()
    hd = pd.MultiIndex.from_product([pd.date_range("2018-01-01", "2024-11-30").strftime("%Y-%m-%d"), FAMS], names=["date", "fam"])
    H = H.set_index(["date", "fam"]).his_n.reindex(hd, fill_value=0).reset_index()
    out = H.merge(out, on=["date", "fam"], how="outer").sort_values(["date", "fam"]).reset_index(drop=True)
    return out


def main():
    D = daily()
    ser = ["S1_v1", "S1_v2", "sumP_v1", "sumP_v2a", "sumP_v2b"]
    jy = D[(D.date >= J0) & (D.date <= J1)]
    nd = jy[ser[0]].notna().groupby(jy.fam).sum()
    sc = pd.DataFrame({s: (jy.groupby("fam").his_n.sum() / 366.0) / (jy.groupby("fam")[s].sum() / nd) for s in ser})
    sc.to_csv(f"{OUT}/t4_join_scale_his_over_ours.csv")
    for s in ser: D[f"{s}_scaled"] = D[s] * D.fam.map(sc[s])
    D.to_csv(f"{OUT}/t4_daily_per_family.csv.gz", index=False, float_format="%.4f")
    # combined series: his to 2024-11-30, ours (scaled) after
    D["month"] = D.date.str[:7]
    Mo = D.groupby(["month", "fam"]).agg(his_n=("his_n", "sum"), his_days=("his_n", lambda v: v.notna().sum()),
                                         ours_days=("S1_v1", lambda v: v.notna().sum()), **{s: (f"{s}_scaled", "sum") for s in ser}).reset_index()
    for s in ser: Mo[f"{s}_rate"] = Mo[s] / Mo.ours_days.replace(0, np.nan)
    Mo["his_rate"] = Mo.his_n / Mo.his_days.replace(0, np.nan)
    Mo.to_csv(f"{OUT}/t4_monthly_per_family.csv", index=False, float_format="%.4f")
    T = Mo.groupby("month")[[f"{s}_rate" for s in ser] + ["his_rate"]].sum(min_count=1).reset_index()
    T.loc[T.month > "2024-11", "his_rate"] = np.nan
    T.to_csv(f"{OUT}/t4_monthly_total.csv", index=False, float_format="%.3f")
    # v1 vs v2, 2025-09-17 .. 2026-09-28 (unscaled counts; ratio identical scaled for S1 and v2a)
    P = D[D.date >= SW1].groupby("fam")[ser + ["sumPany_v1", "sumPany_v2b"]].sum()
    P.loc["all"] = P.sum()
    for a, b in [("S1_v2", "S1_v1"), ("sumP_v2a", "sumP_v1"), ("sumP_v2b", "sumP_v1"), ("sumPany_v2b", "sumPany_v1")]:
        P[f"{a}/{b}"] = P[a] / P[b]
    P["sumP_v2b_scaled/sumP_v1_scaled"] = [(P.loc[f, "sumP_v2b"] * sc.loc[f, "sumP_v2b"]) / (P.loc[f, "sumP_v1"] * sc.loc[f, "sumP_v1"]) if f != "all" else np.nan for f in P.index]
    P.to_csv(f"{OUT}/t4_v1_vs_v2_post_swap.csv", float_format="%.4f")
    print(P[[c for c in P.columns if "/" in c]].describe().round(3).to_string()); print(P.loc["all"].round(3).to_string())
    # step checks (network totals, unscaled = our own detector's counts)
    Tt = D[D.S1_v1.notna()].groupby("date")[ser].sum()
    Tt.index = pd.to_datetime(Tt.index)
    rows = []
    for d0 in ["2025-07-25", "2025-08-26", "2025-09-17"]:
        d0 = pd.Timestamp(d0)
        for w in (21, 45):
            for yr in (0, 1):
                dd = d0 - pd.DateOffset(years=yr)
                bef = Tt[(Tt.index >= dd - pd.Timedelta(days=w)) & (Tt.index < dd)]
                aft = Tt[(Tt.index >= dd) & (Tt.index < dd + pd.Timedelta(days=w))]
                if d0 == pd.Timestamp("2025-08-26") and w == 45 and yr == 0:
                    aft = aft[aft.index < pd.Timestamp("2025-09-17")]
                if d0 == pd.Timestamp("2025-07-25") and yr == 0:
                    aft = aft[aft.index < pd.Timestamp("2025-08-26")]
                for s in ["S1_v2", "sumP_v2b", "S1_v1", "sumP_v1"]:
                    b, a = bef[s], aft[s]
                    qb, qa = b[b <= b.median()], a[a <= a.median()]
                    rows.append(dict(date=str(d0.date()), year=dd.year, window_days=w, series=s, n_before=len(b), n_after=len(a),
                                     mean_before=b.mean(), mean_after=a.mean(), ratio_mean=a.mean() / b.mean(),
                                     median_before=b.median(), median_after=a.median(), ratio_median=a.median() / b.median(),
                                     ratio_quiet_half=qa.mean() / qb.mean()))
    ST = pd.DataFrame(rows); ST.to_csv(f"{OUT}/t4_step_checks.csv", index=False, float_format="%.3f")
    print(ST[ST.series.isin(["S1_v2", "sumP_v2b"])][["date", "year", "window_days", "series", "mean_before", "mean_after", "ratio_mean", "ratio_median", "ratio_quiet_half"]].round(2).to_string())
    # 2025-26 level vs 2024 same calendar months (v2 scaled, network total)
    figs(T, Mo, Tt)


def figs(T, Mo, Tt):
    BLUE, ORANGE, GREY, INK = "#2a78d6", "#eb6834", "#8c8c8c", "#222222"
    t = pd.to_datetime(T.month)
    fig, ax = plt.subplots(2, 1, figsize=(12, 8.5), sharex=False)
    a = ax[0]
    a.plot(t, T.his_rate, color=INK, lw=1.6, label="Shelly (2025) catalogue")
    a.plot(t, T.S1_v1_rate, color=GREY, lw=1.4, ls="--", label="ours S1 9x, v1 (with KMR after swap)")
    a.plot(t, T.S1_v2_rate, color=BLUE, lw=2, label="ours S1 9x, v2 (16 stations from 2025-09-17)")
    a.plot(t, T.sumP_v1_rate, color=GREY, lw=1.2, ls=":", label="ours ΣP(same), v1")
    a.plot(t, T.sumP_v2b_rate, color=ORANGE, lw=2, label="ours ΣP(same), v2 model")
    a.set_ylabel("events / day, all 27 families\n(ours scaled to his at 2023-12..2024-11)")
    a.set_title("Mendocino LFE rate, monthly: v1 vs v2 (KMR fix)", loc="left")
    b = ax[1]
    w = Tt[(Tt.index >= "2025-05-01") & (Tt.index <= "2026-01-31")].rolling(7, min_periods=4).mean()
    b.plot(w.index, w.S1_v1, color=GREY, lw=1.4, ls="--", label="S1 9x v1")
    b.plot(w.index, w.S1_v2, color=BLUE, lw=2, label="S1 9x v2")
    b.plot(w.index, w.sumP_v1, color=GREY, lw=1.2, ls=":", label="ΣP v1")
    b.plot(w.index, w.sumP_v2b, color=ORANGE, lw=2, label="ΣP v2")
    b.set_ylabel("events / day (7-day mean, unscaled)"); b.set_title("Around the KMR failure and swap (2025-05 .. 2026-01)", loc="left")
    for axx in ax:
        lo, hi = axx.get_ylim()
        for d, lab in KD:
            x = pd.Timestamp(d)
            if axx is ax[1] and x > pd.Timestamp("2026-01-31"): continue
            axx.axvline(x, color="#b0b0b0", lw=0.8, zorder=0)
            if axx is ax[1]: axx.text(x, hi * 0.98, " " + lab, rotation=90, va="top", fontsize=7, color="#555555")
        axx.grid(alpha=0.25, lw=0.5); axx.spines[["top", "right"]].set_visible(False)
        axx.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(f"{OUT}/figures/monthly_rate_total_v1_v2.png", dpi=130); plt.close(fig)
    # per family small multiples (monthly scaled rate)
    fig, axs = plt.subplots(7, 4, figsize=(14, 16), sharex=True); axs = axs.ravel()
    for k, f in enumerate(FAMS):
        g = Mo[Mo.fam == f].sort_values("month"); tt = pd.to_datetime(g.month); a = axs[k]
        a.plot(tt, g.his_rate, color=INK, lw=1); a.plot(tt, g.sumP_v1_rate, color=GREY, lw=1, ls="--"); a.plot(tt, g.sumP_v2b_rate, color=ORANGE, lw=1.4)
        a.plot(tt, g.S1_v2_rate, color=BLUE, lw=1.1)
        a.axvline(pd.Timestamp("2025-07-25"), color="#b0b0b0", lw=0.7); a.axvline(pd.Timestamp("2025-09-17"), color="#b0b0b0", lw=0.7)
        a.set_title(f"family {f}", fontsize=8, loc="left"); a.tick_params(labelsize=7); a.spines[["top", "right"]].set_visible(False)
    axs[27].axis("off")
    axs[27].legend(handles=[plt.Line2D([], [], color=INK), plt.Line2D([], [], color=GREY, ls="--"), plt.Line2D([], [], color=ORANGE, lw=1.4), plt.Line2D([], [], color=BLUE)],
                   labels=["Shelly (2025)", "ours ΣP v1", "ours ΣP v2", "ours S1 9x v2"], loc="center", frameon=False, fontsize=9)
    fig.suptitle("Per family, events/day per month (ours scaled at the join year); grey lines 2025-07-25, 2025-09-17", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{OUT}/figures/monthly_rate_per_family_v1_v2.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
