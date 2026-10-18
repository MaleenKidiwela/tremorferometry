#!/usr/bin/env python
"""P4-0a step 3 - rank the Puget candidates by PREDICTED rho and admit sub-regions.

Predicted rho/rho_PGC comes from the north fit (results/p4/rho_calibration.json), applied per
station PER SUB-REGION: the distance term is the tremor-weighted mean HYPOCENTRAL distance to that
sub-region's tremor-active cells, which is the same construction the north fit used (mean
hypocentral distance to the 130 Bostock family locations).  The model is linear, so the mean of the
per-cell predictions equals the prediction at the mean distance - no extra assumption is made.

Admission bars (Merlin P4-0a, as relayed):
  predicted rho/rho_PGC >= 0.5
  tremor contrast >= +3 dB on >= 3 of the 4 tremor days
  far-tremor contrast <= +1 dB          (u_fartremor: 120-350 km tremor, none within 120 km)
  10-19 Hz contrast < 2-8 Hz contrast   (native rate, tremor is a low-frequency source)
  delivery >= 90 % of the planned probe days with coverage >= 0.9
  metadata parsed (response evaluated, channel dip read from StationXML)
The 2-8 Hz PSD bar is reported in BOTH readings (power <= +4.8 dB primary, amplitude <= +9.5 dB)
but is not a separate admission clause: the same noise measurement already enters predicted rho.

Sub-regions: R1 = 47.0-47.5 N, R2 = 48.0-48.3 N, R3 = 47.5-48.0 N (secondary, Merlin F8).
A sub-region PROCEEDS only if, after withholding ranks 3-4 by predicted rho, >= 4 admitted sites
(1-km merge) lie within 45 km hypocentral of >= 50 % of its tremor-active cells AND >= 2 withheld
sites exist.

Writes results/p4/shortlist_v1.csv, subregion_summary.csv, figures/u_shortlist_v1.png
"""
from __future__ import annotations
import os
import sys
import json
import argparse
import warnings
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402
import u_calibrate as UCAL                                             # noqa: E402

SUBREGIONS = {"R1": (47.0, 47.5), "R2": (48.0, 48.3), "R3": (47.5, 48.0)}
RHO_BAR = 0.5
CONTRAST_BAR_DB, CONTRAST_MIN_DAYS = 3.0, 3
FAR_BAR_DB = 1.0
DELIVERY_BAR = 0.90
COV_BAR = 0.9
PB_BACKFILL = ["B001", "B003", "B004", "B005", "B006", "B007", "B013", "B943", "B014",
               "B017", "B018", "B941"]


def db(x):
    x = np.asarray(x, float)
    return 10.0 * np.log10(np.where(x > 0, x, np.nan))


def main():
    ap = argparse.ArgumentParser()
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    CAL = json.load(open(os.path.join(U.P4, "rho_calibration.json")))
    beta = np.array(CAL["models"][CAL["primary"]]["beta"], float)
    cols = CAL["models"][CAL["primary"]]["cols"]
    print(f"[v1] calibration model {CAL['primary']} {cols} beta={np.round(beta,4).tolist()} "
          f"verdict={CAL['verdict']} LOO rho_s={CAL['models'][CAL['primary']]['spearman_loo']:+.3f}")

    TAG = os.environ.get("U_TAG", "")
    S = pd.read_csv(os.path.join(U.P4, os.environ.get("U_STATS_FILE", "station_day_stats.csv")))
    PD = pd.read_csv(os.path.join(U.P4, os.environ.get("U_DAYS_FILE", "probe_days.csv")))
    C = pd.read_csv(os.path.join(U.P4, "candidates.csv"), keep_default_na=False,
                    dtype={"loc": str})
    A = UCAL.aggregate(S, PD).merge(
        C.assign(netsta=C.net + "." + C.sta), on="netsta", how="left")
    INV = pd.read_csv(os.path.join(U.P4, "noise_inventory.csv"))
    A = A.merge(INV[["station", "psd_ratio_med", "psd_ratio_db", "pass_psd3_power",
                     "pass_psd3_amp", "eligible"]].rename(columns={"station": "netsta"}),
                on="netsta", how="left")

    # ---- delivery: delivered/usable over the days the station actually HAD AN EPOCH for.
    # (Dividing by all 8 would punish a station for days before it was installed, which is a
    # station-age statement, not a delivery one; the number of planned days is reported alongside.)
    Q = pd.read_csv(os.path.join(U.P4, os.environ.get("U_QUEUE_FILE", "probe_queue.csv")),
                    keep_default_na=False)
    plan = Q.groupby(Q.net + "." + Q.sta).size()
    ok = S[(S.status == "ok") & (S.coverage >= COV_BAR)].groupby("netsta").size()
    A["n_days_planned"] = A.netsta.map(plan).fillna(0).astype(int)
    A["n_days_ok"] = A.netsta.map(ok).fillna(0).astype(int)
    A["delivery"] = A.n_days_ok / A.n_days_planned.replace(0, np.nan)
    A["n_tremor_planned"] = A.netsta.map(
        Q[Q.date.isin(PD[PD.day_type == "tremor"].date)].groupby(Q.net + "." + Q.sta).size()
    ).fillna(0).astype(int)
    A["metadata_ok"] = (A.resp_mode != "none") & A.chan_az_dip.notna()

    # ---- far-tremor control
    FS = pd.read_csv(os.path.join(U.P4, "far_station_day_stats.csv"))
    FD = pd.read_csv(os.path.join(U.P4, "far_tremor_days.csv"))
    qd = list(PD[PD.day_type == "quiet"].date)
    qbase = S[(S.status == "ok") & S.date.isin(qd)].groupby("netsta").e6_3c.median()
    fmed = FS[FS.status == "ok"].groupby("netsta").e6_3c.median()
    A["far_contrast_db"] = db(A.netsta.map(fmed) / A.netsta.map(qbase))
    A["n_far_days"] = A.netsta.map(FS[FS.status == "ok"].groupby("netsta").size()).fillna(0).astype(int)

    # ---- 10-19 Hz control (native rate)
    hp = os.path.join(U.P4, os.environ.get("U_HF_FILE", "hf_station_day_stats.csv"))
    if os.path.exists(hp):
        H = pd.read_csv(hp)
        H = H[H.status == "ok"]
        piv = H.pivot_table(index="netsta", columns="day_type",
                            values=["e6_lf_3c", "e6_hf_3c"], aggfunc="median")
        hfc = pd.Series(db(piv[("e6_hf_3c", "tremor")].values
                           / piv[("e6_hf_3c", "quiet")].values), index=piv.index)
        lfc = pd.Series(db(piv[("e6_lf_3c", "tremor")].values
                           / piv[("e6_lf_3c", "quiet")].values), index=piv.index)
        A["hf_contrast_db"] = A.netsta.map(hfc)
        A["lf_contrast_native_db"] = A.netsta.map(lfc)
    else:
        A["hf_contrast_db"] = np.nan
        A["lf_contrast_native_db"] = np.nan
    A["hf_below_lf"] = A.hf_contrast_db < A.lf_contrast_native_db
    # ---- the contrast bar is specified in 2-8 Hz; the 20 Hz cache is band-passed 1-8 Hz and the
    # 1-2 Hz microseism dilutes it (north check: Spearman with rho 0.56 on the cache band vs 0.75
    # on the native 2-8 Hz band).  The per-tremor-day bar therefore uses the NATIVE 2-8 Hz
    # statistic; the cache-band contrast stays the calibration variable, because that is the band
    # the north fit was built in.
    if os.path.exists(hp):
        Hq = H[H.day_type == "quiet"].groupby("netsta").e6_lf_3c.median()
        nge, nt = {}, {}
        for ns, g in H[H.day_type == "tremor"].groupby("netsta"):
            base = Hq.get(ns, np.nan)
            v = db(g.e6_lf_3c.values / base) if np.isfinite(base) else np.array([])
            nge[ns] = int(np.sum(v >= CONTRAST_BAR_DB))
            nt[ns] = int(len(v))
        A["n_tremor_days_ge3db_native"] = A.netsta.map(nge).fillna(0).astype(int)
        A["n_tremor_days_native"] = A.netsta.map(nt).fillna(0).astype(int)
    else:
        A["n_tremor_days_ge3db_native"] = A.n_tremor_days_ge3db
        A["n_tremor_days_native"] = A.n_tremor
    # ---- the far bar's own noise floor: the quiet-day-to-quiet-day scatter of the SAME statistic.
    # Median max|deviation| from the quiet median is ~3.4 dB, so an absolute "<= +1 dB" clause is a
    # factor 3 inside the measurement noise.  Both readings are carried: the SPECIFIED absolute bar
    # (primary) and the INTENT bar (the far response must be well below the local one).
    qdev = {}
    for ns, g in S[(S.status == "ok") & S.date.isin(qd)].groupby("netsta"):
        if len(g) >= 3:
            qdev[ns] = float(np.max(np.abs(db(g.e6_3c.values / np.median(g.e6_3c.values)))))
    A["quiet_day_scatter_db"] = A.netsta.map(qdev)
    A["local_minus_far_db"] = A.contrast_3c_db - A.far_contrast_db

    # ---- predicted rho per sub-region
    CELLS = U.active_cells()
    out = []
    for R, (la0, la1) in SUBREGIONS.items():
        cc = CELLS[(CELLS.lat >= la0) & (CELLS.lat < la1 + (0.0001 if R != "R2" else 0.0001))]
        d = U.hyp_km(cc.lat.values[None, :], cc.lon.values[None, :], cc.depth_km.values[None, :],
                     A.lat.values[:, None], A.lon.values[:, None])
        w = cc.tremor_w.values[None, :]
        dtw = (d * w).sum(1) / w.sum()
        B = A.copy()
        B["sub_region"] = R
        B["n_cells_sub"] = len(cc)
        B["dist_field_km"] = dtw
        B["dist_nearest_sub_km"] = d.min(1)
        B["n_cells45"] = (d <= 45.0).sum(1)
        B["n_cells60"] = (d <= 60.0).sum(1)
        X = np.column_stack([np.ones(len(B))] + [B[c].values.astype(float) for c in cols])
        B["pred_rho"] = X @ beta
        out.append(B[B.dist_nearest_sub_km <= 60.0])
    V = pd.concat(out, ignore_index=True)

    V["pass_rho"] = V.pred_rho >= RHO_BAR
    # Merlin's original clause is ">= 4 of 6"; a station that was NOT re-probed still carries
    # only 4 tremor days and is judged on the ">= 3 of 4" reading it was measured under.
    V["contrast_need"] = np.where(V.n_tremor_days_native >= 6, 4, CONTRAST_MIN_DAYS)
    V["pass_contrast"] = V.n_tremor_days_ge3db_native >= V.contrast_need
    V["pass_contrast_cacheband"] = V.n_tremor_days_ge3db >= CONTRAST_MIN_DAYS
    V["pass_far"] = V.far_contrast_db <= FAR_BAR_DB                      # SPECIFIED (primary)
    V["pass_far_margin"] = V.local_minus_far_db >= 3.0                   # INTENT (noise-aware)
    V["pass_hf"] = V.hf_below_lf.fillna(False)
    V["pass_delivery"] = V.delivery >= DELIVERY_BAR
    V["enough_tremor_days"] = V.n_tremor_planned >= CONTRAST_MIN_DAYS
    V["admit"] = (V.pass_rho & V.pass_contrast & V.pass_far & V.pass_hf
                  & V.pass_delivery & V.metadata_ok)
    V["admit_margin"] = (V.pass_rho & V.pass_contrast & V.pass_far_margin & V.pass_hf
                         & V.pass_delivery & V.metadata_ok)
    V["admit_plus_psd"] = V.admit & V.pass_psd3_power.fillna(False)
    def _why(r):
        if not r.enough_tremor_days:
            return f"only {r.n_tremor_planned} tremor days in its epoch"
        for f, lab in [("metadata_ok", "metadata"), ("pass_delivery", "delivery"),
                       ("pass_rho", "predicted_rho"), ("pass_contrast", "tremor_contrast"),
                       ("pass_far", "far_tremor"), ("pass_hf", "10-19Hz")]:
            if not bool(getattr(r, f)):
                return lab
        return ""
    V["first_failing_bar"] = [_why(r) for r in V.itertuples()]
    V["is_pb_backfill"] = V.sta.isin(PB_BACKFILL)
    V = V.sort_values(["sub_region", "pred_rho"], ascending=[True, False]).reset_index(drop=True)

    # ---- site-level ranking, withholding, and the sub-region decision
    rows, summary = [], []
    ADMIT_COL = os.environ.get("U_ADMIT_COL", "admit")
    for R, g in V.groupby("sub_region"):
        la0, la1 = SUBREGIONS[R]
        cc = CELLS[(CELLS.lat >= la0) & (CELLS.lat < la1 + 0.0001)]
        adm = g[g[ADMIT_COL]].copy()
        # one entry per SITE: the best-predicted station of the site represents it
        site_rank = (adm.sort_values("pred_rho", ascending=False)
                     .drop_duplicates("site").reset_index(drop=True))
        site_rank["site_rank"] = np.arange(1, len(site_rank) + 1)
        withheld = set(site_rank[site_rank.site_rank.isin([3, 4])].site)
        summing = site_rank[~site_rank.site.isin(withheld)]
        if len(summing):
            d = U.hyp_km(cc.lat.values[:, None], cc.lon.values[:, None], cc.depth_km.values[:, None],
                         summing.lat.values[None, :], summing.lon.values[None, :])
            n45 = (d <= 45.0).sum(1)
            az = U.azimuth(cc.lat.values[:, None], cc.lon.values[:, None],
                           summing.lat.values[None, :], summing.lon.values[None, :])
            gaps = np.array([U.az_gap(az[i][d[i] <= 60.0]) for i in range(len(cc))])
            frac4 = float((n45 >= 4).mean())
            frac3 = float((n45 >= 3).mean())
        else:
            n45 = np.zeros(len(cc)); gaps = np.full(len(cc), 360.0); frac4 = frac3 = 0.0
        proceed = bool(len(summing) >= 4 and frac4 >= 0.5 and len(withheld) >= 2)
        summary.append(dict(sub_region=R, lat0=la0, lat1=la1, n_active_cells=len(cc),
                            n_probed=int(g.netsta.nunique()),
                            n_admitted_stations=int(adm.netsta.nunique()),
                            n_admitted_sites=int(len(site_rank)),
                            n_summing_sites=int(len(summing)), n_withheld_sites=len(withheld),
                            withheld=";".join(sorted(withheld)),
                            frac_cells_n45_ge4=frac4, frac_cells_n45_ge3=frac3,
                            n45_median=float(np.median(n45)),
                            gap_median_deg=float(np.median(gaps)),
                            verdict="PROCEED" if proceed else "STOP"))
        for r in g.itertuples():
            adm_r = bool(getattr(r, ADMIT_COL))
            role = ("withheld" if (adm_r and r.site in withheld)
                    else "admitted" if adm_r
                    else "backfill_only" if r.is_pb_backfill else "rejected")
            rows.append(dict(
                station=r.netsta, net=r.net, sta=r.sta, site=r.site, sub_region=R, role=role,
                sensor_type=r.sensor_type, band=r.band, fs=r.fs, chain=r.chain,
                pred_rho=r.pred_rho, dist_field_km=r.dist_field_km,
                dist_nearest_cell_km=r.dist_nearest_sub_km, n_cells45=r.n_cells45,
                noise_3c_db=r.noise_3c_db, noise_Z_db=r.noise_Z_db, noise_H_db=r.noise_H_db,
                psd_ratio_power_db=r.psd_ratio_db,
                pass_psd_power_4p8dB=bool(r.pass_psd3_power), pass_psd_amp_9p5dB=bool(r.pass_psd3_amp),
                contrast_3c_db=r.contrast_3c_db, n_tremor_days_ge3db=r.n_tremor_days_ge3db,
                far_contrast_db=r.far_contrast_db, n_far_days=r.n_far_days,
                hf_contrast_db=r.hf_contrast_db, lf_contrast_native_db=r.lf_contrast_native_db,
                delivery=r.delivery, n_days_ok=r.n_days_ok, coverage_med=r.coverage_med,
                metadata_ok=bool(r.metadata_ok), zdip_inverted=bool(r.zsign < 0),
                n_tremor_days_ge3db_native=r.n_tremor_days_ge3db_native,
                n_tremor_days_native=r.n_tremor_days_native,
                contrast_need=int(r.contrast_need),
                n_tremor_planned=r.n_tremor_planned, n_days_planned=r.n_days_planned,
                quiet_day_scatter_db=r.quiet_day_scatter_db,
                local_minus_far_db=r.local_minus_far_db,
                pass_rho=bool(r.pass_rho), pass_contrast=bool(r.pass_contrast),
                pass_contrast_cacheband=bool(r.pass_contrast_cacheband),
                pass_far=bool(r.pass_far), pass_far_margin=bool(r.pass_far_margin),
                pass_hf=bool(r.pass_hf), pass_delivery=bool(r.pass_delivery),
                first_failing_bar=r.first_failing_bar,
                admit=bool(r.admit), admit_margin=bool(r.admit_margin),
                admit_plus_psd=bool(r.admit_plus_psd), is_pb_backfill=bool(r.is_pb_backfill)))
    SL = pd.DataFrame(rows).sort_values(["sub_region", "role", "pred_rho"],
                                        ascending=[True, True, False])
    SL.to_csv(os.path.join(U.P4, f"shortlist_v1{TAG}{'' if ADMIT_COL == 'admit' else '_margin'}.csv"),
              index=False)
    SUM = pd.DataFrame(summary)
    SUM["admit_rule"] = ADMIT_COL
    SUM.to_csv(os.path.join(U.P4, f"subregion_summary{TAG}{'' if ADMIT_COL == 'admit' else '_margin'}.csv"),
               index=False)
    print(SUM.to_string(index=False))
    print(f"[v1] shortlist_v1.csv: {len(SL)} station x sub-region rows; "
          f"roles {SL.role.value_counts().to_dict()}")
    for R, g in SL.groupby("sub_region"):
        adm = g[g.role.isin(["admitted", "withheld"])]
        print(f"  {R}: " + ", ".join(f"{r.sta}({r.role[0]},{r.pred_rho:.2f})"
                                     for r in adm.sort_values('pred_rho', ascending=False).itertuples()))
    bars = ["pass_rho", "pass_contrast", "pass_contrast_cacheband", "pass_far",
            "pass_far_margin", "pass_hf", "pass_delivery", "enough_tremor_days"]
    uniq = V.drop_duplicates("netsta")
    print("[v1] per-bar pass counts over the " + str(len(uniq)) + " probed stations: "
          + str({b: int(uniq[b].sum()) for b in bars}))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    for i, R in enumerate(["R1", "R3", "R2"]):
        g = V[V.sub_region == R]
        adm = g[g.admit] if ADMIT_COL == "admit" else g[g.admit_margin]
        ax[i].scatter(g.dist_field_km, g.pred_rho, s=24, c="#bbbbbb", edgecolor="k", linewidth=0.3,
                      label=f"probed (n={len(g)})")
        ax[i].scatter(adm.dist_field_km, adm.pred_rho, s=52, c="#ff7f0e", edgecolor="k",
                      linewidth=0.5, label=f"admitted (n={len(adm)})")
        for r in adm.itertuples():
            ax[i].annotate(r.sta, (r.dist_field_km, r.pred_rho), fontsize=6,
                           xytext=(3, 3), textcoords="offset points")
        ax[i].axhline(RHO_BAR, color="k", ls="--", lw=1)
        la0, la1 = SUBREGIONS[R]
        v = SUM[SUM.sub_region == R].verdict.iloc[0]
        ax[i].set_title(f"{R}  {la0}-{la1} N  -> {v}")
        ax[i].set_xlabel("tremor-weighted mean hypocentral distance to the sub-region's cells (km)")
        ax[i].set_ylabel(r"predicted $\rho/\rho_{PGC}$")
        ax[i].legend(fontsize=8)
        ax[i].grid(alpha=0.25)
    fig.suptitle("P4-0a: Puget candidates ranked by predicted single-event registration", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(U.FIGS, f"u_shortlist_v1{TAG}.png"), dpi=150)
    print(f"[fig] -> {U.FIGS}/u_shortlist_v1.png")


if __name__ == "__main__":
    main()
