#!/usr/bin/env python
"""P4-0 step 5 - the inventory, the shortlist, the geometry it yields, and the figure.

Writes
  results/p4/station_probe_days.csv  one row per planned station-day: what was asked for, what
                                     arrived, and the per-day statistics.
  results/p4/noise_inventory.csv     one row per station: site, network, sensor, distance to the
                                     nearest tremor-active cell, 2-8 Hz velocity PSD ratio to
                                     CN.PGC on each of the 8 days and its median, the incoherent
                                     6-s tremor/quiet contrast per component, coverage.
  results/p4/shortlist_v0.csv        the stations passing BOTH clauses, ranked.
  results/p4/station_rank_all.csv    every station with the same columns and the pass flags, so a
                                     rejected station is never lost.
  results/p4/shortlist_cell_geometry.csv + shortlist_geometry.json
                                     the <=60 km (and <=45 km) SITE count per active cell that the
                                     shortlist yields, recomputed with the P3D geometry code, next
                                     to the all-candidate ceiling and the S15 baseline.
  figures/u_noise_inventory.png

SHORTLIST RULE (as specified, applied literally):
    psd_ratio_to_PGC (POWER, median over days) <= 3      AND
    tremor/quiet contrast (dB, 3-C) >= PGC's contrast on the SAME days.
Eligibility: >=2 tremor days and >=2 quiet days with coverage >= 0.8 and a usable response.
A second flag `pass_amp3` applies the same rule with the 3x read as an AMPLITUDE ratio (9.5 dB
rather than 4.8 dB) because "3x PGC" is ambiguous; both counts are reported.
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

MIN_COV = 0.8
MIN_DAYS_EACH = 2
S15 = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007",
       "B009", "B010", "B011", "B926"]


def db(x):
    x = np.asarray(x, float)
    return 10.0 * np.log10(np.where(x > 0, x, np.nan))


def load():
    S = pd.read_csv(os.path.join(U.P4, "station_day_stats.csv"), keep_default_na=False,
                    na_values=[""])
    Q = pd.read_csv(os.path.join(U.P4, "probe_queue.csv"), keep_default_na=False)
    PD = pd.read_csv(os.path.join(U.P4, "probe_days.csv"))
    C = pd.read_csv(os.path.join(U.P4, "candidates.csv"), keep_default_na=False, na_values=[""])
    M = pd.read_csv(os.path.join(U.P4, "download_manifest.csv"), keep_default_na=False)
    M = M.drop_duplicates(["net", "sta", "date"], keep="last")
    return S, Q, PD, C, M


def probe_days_table(S, Q, PD, M, out):
    T = Q.merge(PD[["date", "day_type", "n_tremor_50km", "n_tremor_75km"]], on="date", how="left")
    T = T.merge(M[["net", "sta", "date", "status", "bytes", "fs", "coverage_frac", "n_frag", "note"]]
                .rename(columns={"status": "download_status"}), on=["net", "sta", "date"], how="left")
    keep = ["netsta", "date", "status", "coverage", "n_seg", "n_win", "chain", "fs_native",
            "psd_Z", "psd_N", "psd_E", "psd_3c", "e6_Z", "e6_N", "e6_E", "e6_3c", "resp_mode", "err"]
    S2 = S[[c for c in keep if c in S.columns]].rename(columns={"status": "stat_status"})
    T["netsta"] = T.net + "." + T.sta
    T = T.merge(S2, on=["netsta", "date"], how="left")
    T = T.sort_values(["netsta", "day_type", "date"])
    T.to_csv(out, index=False)
    return T


def inventory(S, C, PD, out):
    ok = S[(S.status == "ok") & (S.coverage >= MIN_COV) & S.psd_3c.notna()].copy()
    pgc = ok[ok.netsta == U.PGC_KEY].set_index("date")
    if not len(pgc):
        raise RuntimeError("CN.PGC has no usable probe day - the reference is missing")
    dates = list(PD.sort_values(["day_type", "date"]).date)
    rows = []
    for ns, g in ok.groupby("netsta"):
        g = g.set_index("date")
        m = C[C.net + "." + C.sta == ns]
        m = m.iloc[0] if len(m) else None
        r = dict(station=ns, net=ns.split(".")[0], sta=ns.split(".")[1])
        if m is not None:
            r.update(site=m.site, sensor_type=m.sensor_type, sensor_desc=m.sensor_desc,
                     band=m.band, loc=m["loc"], fs=m.fs, chain=m.chain,
                     lat=m.lat, lon=m.lon, dist_nearest_cell_km=m.dmin_cell_km,
                     nearest_cell=m.nearest_cell, n_cells45=m.n_cells45, n_cells60=m.n_cells60,
                     tw_cells60=m.tw_cells60, in_roster183=m.in_roster183,
                     zdip_inverted=bool(m.zsign < 0), sensitivity=m.sensitivity_med)
        common = [d for d in dates if d in g.index and d in pgc.index]
        for d in dates:
            r[f"ratio_{d}"] = (float(g.psd_3c[d] / pgc.psd_3c[d])
                               if (d in g.index and d in pgc.index) else np.nan)
        rat = np.array([r[f"ratio_{d}"] for d in dates], float)
        r["psd_ratio_med"] = float(np.nanmedian(rat)) if np.isfinite(rat).any() else np.nan
        r["psd_ratio_db"] = float(db(r["psd_ratio_med"]))
        r["psd_3c_med"] = float(np.nanmedian(g.psd_3c.values))
        td = [d for d in common if PD.set_index("date").day_type[d] == "tremor"]
        qd = [d for d in common if PD.set_index("date").day_type[d] == "quiet"]
        r["n_days_tremor"], r["n_days_quiet"] = len(td), len(qd)
        r["coverage_med"] = float(np.nanmedian(g.coverage.values))
        r["resp_mode"] = g.resp_mode.mode().iloc[0] if len(g) else ""
        for lab in ["Z", "N", "E", "3c"]:
            c = f"e6_{lab}"
            if td and qd and c in g:
                r[f"contrast_{lab}_db"] = float(db(np.nanmedian(g.loc[td, c].values))
                                                - db(np.nanmedian(g.loc[qd, c].values)))
            else:
                r[f"contrast_{lab}_db"] = np.nan
            if td and qd and c in pgc:
                pt = [d for d in td if d in pgc.index]
                pq = [d for d in qd if d in pgc.index]
                r[f"pgc_contrast_{lab}_db"] = (float(db(np.nanmedian(pgc.loc[pt, c].values))
                                                     - db(np.nanmedian(pgc.loc[pq, c].values)))
                                               if pt and pq else np.nan)
            else:
                r[f"pgc_contrast_{lab}_db"] = np.nan
        r["contrast_margin_db"] = r["contrast_3c_db"] - r["pgc_contrast_3c_db"]
        r["eligible"] = bool(r["n_days_tremor"] >= MIN_DAYS_EACH
                             and r["n_days_quiet"] >= MIN_DAYS_EACH
                             and r["coverage_med"] >= MIN_COV
                             and r["resp_mode"] != "none")
        r["pass_psd3_power"] = bool(np.isfinite(r["psd_ratio_med"]) and r["psd_ratio_med"] <= 3.0)
        r["pass_psd3_amp"] = bool(np.isfinite(r["psd_ratio_med"]) and r["psd_ratio_med"] <= 9.0)
        r["pass_contrast"] = bool(np.isfinite(r["contrast_margin_db"]) and r["contrast_margin_db"] >= 0)
        r["pass_v0"] = bool(r["eligible"] and r["pass_psd3_power"] and r["pass_contrast"])
        r["pass_amp3"] = bool(r["eligible"] and r["pass_psd3_amp"] and r["pass_contrast"])
        rows.append(r)
    I = pd.DataFrame(rows)
    I = I.sort_values(["pass_v0", "contrast_margin_db"], ascending=[False, False])
    I.to_csv(out, index=False)
    return I


def geometry(I, log):
    """Recompute the <=60 km / <=45 km SITE count per active cell for a set of stations."""
    C = U.active_cells()
    ST = pd.read_csv(os.path.join(U.GEO, "stations_all_2019plus.csv"))
    def cells_for(keys, tag):
        D = ST[ST.key.isin(keys)][["key", "lat", "lon"]].copy()
        if not len(D):
            return pd.DataFrame(), {}
        sm = U.cluster_sites(D)
        D["site"] = [sm[k] for k in D.key]
        SI = D.groupby("site")[["lat", "lon"]].mean().reset_index()
        d = U.hyp_km(C.lat.values[:, None], C.lon.values[:, None], C.depth_km.values[:, None],
                     SI.lat.values[None, :], SI.lon.values[None, :])
        az = U.azimuth(C.lat.values[:, None], C.lon.values[:, None],
                       SI.lat.values[None, :], SI.lon.values[None, :])
        rows = []
        for i, c in enumerate(C.itertuples()):
            m60, m45 = d[i] <= 60.0, d[i] <= 45.0
            rows.append(dict(set=tag, cell_id=c.cell_id, lat=c.lat, lon=c.lon,
                             depth_km=c.depth_km, tremor_w=c.tremor_w,
                             n45=int(m45.sum()), n60=int(m60.sum()),
                             nearest_km=float(d[i].min()),
                             nearest_site=SI.site.values[int(d[i].argmin())],
                             gap_deg=U.az_gap(az[i][m60]),
                             sites60="+".join(sorted(SI.site.values[m60]))))
        G = pd.DataFrame(rows)
        s = dict(set=tag, n_sites=int(SI.site.nunique()), n_stations=int(len(D)),
                 cells=int(len(C)),
                 n60_median=float(G.n60.median()), n60_p10=float(G.n60.quantile(0.10)),
                 n60_p90=float(G.n60.quantile(0.90)),
                 n45_median=float(G.n45.median()), n45_p90=float(G.n45.quantile(0.90)),
                 frac_cells_n60_ge3=float((G.n60 >= 3).mean()),
                 frac_cells_n60_ge4=float((G.n60 >= 4).mean()),
                 frac_cells_n45_ge3=float((G.n45 >= 3).mean()),
                 gap_median=float(G.gap_deg.median()),
                 tw_weighted_n60=float((G.n60 * G.tremor_w).sum() / G.tremor_w.sum()))
        return G, s
    CAND = pd.read_csv(os.path.join(U.P4, "candidates.csv"))
    sets = {"shortlist_v0": set(I[I.pass_v0].station),
            "shortlist_amp3": set(I[I.pass_amp3].station),
            "eligible_probed": set(I[I.eligible].station),
            "all_probed": set(I.station),
            "all_candidates": set(CAND.key),
            "S15_baseline": {k for k in ST.key if k.split(".")[-1] in S15}}
    Gs, ss = [], []
    for tag, keys in sets.items():
        G, s = cells_for(keys, tag)
        if len(G):
            Gs.append(G)
            ss.append(s)
            log(f"  [geom] {tag:16s} {s['n_stations']:3d} stations / {s['n_sites']:3d} sites -> "
                f"n60 median {s['n60_median']:.0f} (p10 {s['n60_p10']:.0f}, p90 {s['n60_p90']:.0f}), "
                f"n45 median {s['n45_median']:.0f}, cells with >=3 sites <=60 km "
                f"{100*s['frac_cells_n60_ge3']:.0f}%, median gap {s['gap_median']:.0f} deg")
    return pd.concat(Gs, ignore_index=True), ss


def figure(I, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    E = I[I.eligible].copy()
    col = {"surface_broadband": "#1f77b4", "borehole": "#d62728",
           "surface_short_period": "#2ca02c", "reference": "#000000"}
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.6))
    for st, g in E.groupby("sensor_type"):
        ax[0].scatter(g.dist_nearest_cell_km, g.psd_ratio_db, s=26, alpha=0.85,
                      c=col.get(st, "#888888"), label=f"{st} (n={len(g)})",
                      edgecolor="k", linewidth=0.3)
    ax[0].axhline(10 * np.log10(3.0), color="k", ls="--", lw=1)
    ax[0].text(0.99, 10 * np.log10(3.0) + 0.6, "3x PGC power (4.8 dB)", transform=
               ax[0].get_yaxis_transform(), ha="right", fontsize=8)
    ax[0].axhline(10 * np.log10(9.0), color="k", ls=":", lw=1)
    ax[0].text(0.99, 10 * np.log10(9.0) + 0.6, "3x PGC amplitude (9.5 dB)", transform=
               ax[0].get_yaxis_transform(), ha="right", fontsize=8)
    ax[0].axhline(0, color="#888888", lw=0.8)
    ax[0].set_xlabel("hypocentral distance to the nearest tremor-active cell (km)")
    ax[0].set_ylabel("median 2-8 Hz velocity PSD relative to CN.PGC (dB)")
    ax[0].set_title("Site noise vs distance (P4 band 47.0-48.3 N)")
    ax[0].legend(fontsize=8, loc="upper left")
    ax[0].grid(alpha=0.25)
    sel = E[E.pass_v0]
    for st, g in E.groupby("sensor_type"):
        ax[1].scatter(g.psd_ratio_db, g.contrast_3c_db, s=26, alpha=0.85,
                      c=col.get(st, "#888888"), edgecolor="k", linewidth=0.3, label=st)
    ax[1].scatter(sel.psd_ratio_db, sel.contrast_3c_db, s=110, facecolor="none",
                  edgecolor="#ff7f0e", linewidth=1.6, label=f"shortlist_v0 (n={len(sel)})")
    pg = float(np.nanmedian(E.pgc_contrast_3c_db))
    ax[1].axhline(pg, color="k", ls="--", lw=1)
    ax[1].text(0.01, pg + 0.2, f"CN.PGC contrast {pg:.1f} dB", transform=
               ax[1].get_yaxis_transform(), fontsize=8)
    ax[1].axvline(10 * np.log10(3.0), color="k", ls="--", lw=1)
    ax[1].set_xlabel("median 2-8 Hz PSD relative to CN.PGC (dB)")
    ax[1].set_ylabel("incoherent 6-s energy: tremor days - quiet days (dB, 3-C)")
    ax[1].set_title("Registration proxy vs site noise")
    ax[1].legend(fontsize=8, loc="best")
    ax[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"[fig] -> {path}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    lines = []

    def log(m):
        print(m, flush=True)
        lines.append(m)

    S, Q, PD, C, M = load()
    T = probe_days_table(S, Q, PD, M, os.path.join(U.P4, "station_probe_days.csv"))
    log(f"[inv] station_probe_days.csv: {len(T)} planned station-days, "
        f"{int((T.get('stat_status') == 'ok').sum())} with statistics")
    I = inventory(S, C, PD, os.path.join(U.P4, "noise_inventory.csv"))
    log(f"[inv] noise_inventory.csv: {len(I)} stations, eligible {int(I.eligible.sum())}, "
        f"pass PSD<=3x(power) {int(I.pass_psd3_power.sum())}, pass contrast "
        f"{int(I.pass_contrast.sum())}, pass_v0 {int(I.pass_v0.sum())}, "
        f"pass_amp3 {int(I.pass_amp3.sum())}")
    I.to_csv(os.path.join(U.P4, "station_rank_all.csv"), index=False)
    SL = I[I.pass_v0].sort_values(["contrast_margin_db", "psd_ratio_med"],
                                  ascending=[False, True]).reset_index(drop=True)
    SL.insert(0, "rank", np.arange(1, len(SL) + 1))
    G, ss = geometry(I, log)
    G.to_csv(os.path.join(U.P4, "shortlist_cell_geometry.csv"), index=False)
    gsum = {s["set"]: s for s in ss}
    U.jdump(gsum, os.path.join(U.P4, "shortlist_geometry.json"))
    for k, v in gsum.get("shortlist_v0", {}).items():
        if k != "set":
            SL[f"geom_{k}"] = v
    SL.to_csv(os.path.join(U.P4, "shortlist_v0.csv"), index=False)
    log(f"[inv] shortlist_v0.csv: {len(SL)} stations / "
        f"{gsum.get('shortlist_v0', {}).get('n_sites', 0)} sites")
    figure(I, os.path.join(U.FIGS, "u_noise_inventory.png"))
    with open(os.path.join(U.P4, "u_inventory.log"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
