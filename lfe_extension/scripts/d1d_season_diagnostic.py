#!/usr/bin/env python
"""Track D, D1(d) — WHAT DRIVES THE WINTER / SUMMER QUIET-DAY RATIO (coordinator request 2026-09-17).

The floor λ is ~3.5x higher in DJF than in JJA (section 1c of RESULTS.md).  The bar is reported as NOT
met and nothing is tuned; this script asks WHY, from Track C's own outputs only (read-only).  Two
competing explanations:

  (A) a NOISE / MAD-adaptivity effect — winter storms raise the ambient level, and the detector's
      threshold is in units of each family-day's own MAD, so the rate only changes if the noise gets
      BURSTIER (heavier-tailed relative to its own MAD), not merely louder.  Signatures: higher daily
      e8, a SHALLOWER decay of detections with threshold in winter, a positive within-season slope of
      the rate against the daily noise term v, more HF-vetoed detections, more EQ flags.
  (B) a REAL inter-ETS activity effect — families genuinely fire more in winter.  Signature: the winter
      quiet-day detections register at a withheld station (tested separately by
      `d1d_season_stack_test.py`, the Track A step-3b test applied per season).

Outputs: results/trackD/season_diagnostic_{tail,noise,byfamily}.csv and figures/d1d_season_diagnostic.png"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_tiers  # noqa: E402

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TC = f"{EXT}/results/trackC"
TD = f"{EXT}/results/trackD"
FIG = f"{EXT}/figures"
SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
          9: "SON", 10: "SON", 11: "SON"}
GRID = [8.0, 8.25, 8.5, 8.75, 9.0, 9.5, 10.0, 11.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-prefix", default=f"{TD}/season_diagnostic")
    a = ap.parse_args()
    EPOCHS, TIERS, TSRC = d_tiers.load()
    FD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{TC}/family_days_*.parquet"))],
                   ignore_index=True).drop_duplicates(["fam", "date"])
    DT = pd.concat([pd.read_parquet(p, columns=["fam", "date", "net_mad", "mad_scale", "n_sites",
                                                "hf_pass", "eq_flag", "quiet_flag", "hfr"])
                    for p in sorted(glob.glob(f"{TC}/detections_*.parquet"))], ignore_index=True)
    for D in (FD, DT):
        D["epoch_id"] = d_tiers.epoch_of(D.date.to_numpy().astype(str), EPOCHS)
        D["season"] = pd.DatetimeIndex(D.date).month.map(SEASON)
    FD = FD[(FD.epoch_id != "") & FD.observable]
    Q = FD[FD.quiet_flag]
    DQ = DT[DT.quiet_flag & (DT.epoch_id != "")]

    # ---- (1) TAIL SHAPE: detections per quiet template-day vs threshold, by epoch and season
    rows = []
    for (eid, sea), g in Q.groupby(["epoch_id", "season"]):
        nd = len(g)
        d = DQ[(DQ.epoch_id == eid) & (DQ.season == sea)]
        for T in GRID:
            rows.append(dict(epoch_id=eid, season=sea, T=T, n_quiet_days=nd,
                             rate_all=float((d.net_mad >= T).sum()) / max(nd, 1),
                             rate_hfpass=float(((d.net_mad >= T) & d.hf_pass).sum()) / max(nd, 1),
                             rate_hfveto=float(((d.net_mad >= T) & ~d.hf_pass).sum()) / max(nd, 1),
                             rate_eq=float(((d.net_mad >= T) & d.eq_flag).sum()) / max(nd, 1)))
    TA = pd.DataFrame(rows)
    # log-linear decay slope over 8.25..10 (how fast the tail falls; a heavier tail = shallower)
    sl = []
    for (eid, sea), g in TA.groupby(["epoch_id", "season"]):
        g = g[(g["T"] >= 8.25) & (g["T"] <= 10.0) & (g.rate_hfpass > 0)]
        if len(g) >= 3:
            b = np.polyfit(g["T"], np.log10(g.rate_hfpass), 1)[0]
            sl.append(dict(epoch_id=eid, season=sea, tail_slope_dex_per_MAD=float(b),
                           n_quiet_days=int(g.n_quiet_days.iloc[0])))
    SL = pd.DataFrame(sl)
    TA = TA.merge(SL, on=["epoch_id", "season", "n_quiet_days"], how="left")
    TA.to_csv(f"{a.out_prefix}_tail.csv", index=False)

    # ---- (2) NOISE and MAD by season (the adaptivity question)
    MC = pd.read_parquet(f"{TD}/mc_family_day.parquet") if os.path.exists(f"{TD}/mc_family_day.parquet") else None
    N = []
    dn = [pd.read_parquet(p) for p in sorted(glob.glob(f"{TC}/daily_noise_*.parquet"))]
    if dn:
        DN = pd.concat(dn, ignore_index=True)
        DN = DN[DN.status.isin(["cached", "built"])].copy()
        DN["e8"] = DN[["e8_Z", "e8_N", "e8_E"]].mean(axis=1)
        DN["season"] = pd.DatetimeIndex(DN.date).month.map(SEASON)
        DN["epoch_id"] = d_tiers.epoch_of(DN.date.to_numpy().astype(str), EPOCHS)
        ref = DN.groupby("netsta").e8.median()
        DN["e8_rel"] = DN.e8 / DN.netsta.map(ref)
        N = DN.groupby(["epoch_id", "season"]).agg(
            e8_rel_median=("e8_rel", "median"), e8_rel_p90=("e8_rel", lambda v: v.quantile(.9)),
            n_station_days=("e8", "size")).reset_index()
    M = Q.groupby(["epoch_id", "season"]).agg(mad_scale_median=("mad_scale", "median"),
                                              n_sites_median=("n_sites", "median"),
                                              n_quiet_days=("fam", "size")).reset_index()
    if len(N):
        M = M.merge(N, on=["epoch_id", "season"], how="left")
    if MC is not None:
        mq = MC.merge(Q[["fam", "date"]].rename(columns={"fam": "family_id"}), on=["family_id", "date"])
        mq["season"] = pd.DatetimeIndex(mq.date).month.map(SEASON)
        M = M.merge(mq.groupby(["epoch_id", "season"]).agg(u_median=("u", "median"),
                                                           v_median=("v", "median"),
                                                           Mc_core_median=("Mc_T_core", "median")).reset_index(),
                    on=["epoch_id", "season"], how="left")
    M.to_csv(f"{a.out_prefix}_noise.csv", index=False)

    # ---- (3) within-season slope of the quiet-day rate against the daily noise term v
    slope_rows = []
    if MC is not None:
        q = Q[["fam", "date", "epoch_id", "season"]].rename(columns={"fam": "family_id"})
        q = q.merge(MC[["family_id", "date", "v", "u"]], on=["family_id", "date"], how="left")
        cnt = DQ[DQ.hf_pass].groupby(["fam", "date"]).apply(
            lambda g: pd.Series({f"n{T:g}": int((g.net_mad >= T).sum()) for T in (8.25, 8.75)}),
            include_groups=False).reset_index().rename(columns={"fam": "family_id"})
        q = q.merge(cnt, on=["family_id", "date"], how="left").fillna({"n8.25": 0, "n8.75": 0})
        for (eid, sea), g in q.groupby(["epoch_id", "season"]):
            if len(g) < 200 or g.v.std() == 0:
                continue
            for col, lab in (("n8.25", "T_ext"), ("n8.75", "T_core")):
                b = np.polyfit(g.v.fillna(0), g[col], 1)
                slope_rows.append(dict(epoch_id=eid, season=sea, tier=lab, n_family_days=len(g),
                                       rate=float(g[col].mean()), slope_per_dex_v=float(b[0]),
                                       corr=float(np.corrcoef(g.v.fillna(0), g[col])[0, 1])))
        q.to_csv(f"{a.out_prefix}_byfamily.csv", index=False)
    SLV = pd.DataFrame(slope_rows)
    if len(SLV):
        SLV.to_csv(f"{a.out_prefix}_vslope.csv", index=False)

    pd.set_option("display.width", 220)
    print(f"[tiers] {TSRC}")
    print("\n[tail shape — detections per quiet template-day vs threshold, HF-passing]")
    piv = TA.pivot_table(index=["epoch_id", "season"], columns="T", values="rate_hfpass")
    print(piv.round(4).to_string())
    print("\n[tail decay slope (dex per MAD, 8.25–10); more negative = lighter tail]")
    print(SL.round(3).to_string(index=False))
    print("\n[noise, MAD and coverage on quiet family-days]")
    print(M.round(4).to_string(index=False))
    if len(SLV):
        print("\n[within-season slope of the quiet-day count against the daily noise term v]")
        print(SLV.round(4).to_string(index=False))

    # ---- figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for (eid, sea), g in TA.groupby(["epoch_id", "season"]):
        if g.rate_hfpass.max() <= 0:
            continue
        ax[0].semilogy(g["T"], g.rate_hfpass.replace(0, np.nan), "o-", ms=3.5, label=f"{eid} {sea}")
    ax[0].set_xlabel("threshold (MAD)"); ax[0].set_ylabel("detections / quiet template-day")
    ax[0].set_title("(a) quiet-day tail by season", fontsize=10); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
    if "e8_rel_median" in M:
        o = ["DJF", "MAM", "JJA", "SON"]
        for eid, g in M.groupby("epoch_id"):
            g = g.set_index("season").reindex(o)
            ax[1].plot(o, g.e8_rel_median, "o-", label=f"{eid} noise e8 (rel)")
            ax[1].plot(o, g.mad_scale_median / M.mad_scale_median.median(), "s--",
                       label=f"{eid} MAD scale (rel)")
        ax[1].set_ylabel("relative to the record median")
        ax[1].set_title("(b) does the MAD track the noise?", fontsize=10)
        ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    if len(SLV):
        for tier, g in SLV.groupby("tier"):
            ax[2].scatter(g.rate, g.slope_per_dex_v, s=30, label=tier)
            for r in g.itertuples():
                ax[2].annotate(f"{r.epoch_id} {r.season}", (r.rate, r.slope_per_dex_v), fontsize=6,
                               xytext=(3, 3), textcoords="offset points")
        ax[2].axhline(0, color="k", lw=0.8)
        ax[2].set_xlabel("mean quiet-day detections / family-day")
        ax[2].set_ylabel("d(count) / d(log10 noise amplitude)")
        ax[2].set_title("(c) within-season noise dependence", fontsize=10)
        ax[2].legend(fontsize=7); ax[2].grid(alpha=.3)
    fig.suptitle("Track D — what drives the winter / summer quiet-day ratio", fontsize=12)
    fig.tight_layout(); fig.savefig(f"{FIG}/d1d_season_diagnostic.png", dpi=135); plt.close(fig)
    print(f"\n[write] {a.out_prefix}_{{tail,noise,byfamily,vslope}}.csv, {FIG}/d1d_season_diagnostic.png")


if __name__ == "__main__":
    main()
