#!/usr/bin/env python
"""E2 figures (2026-09-10, resumed run).  Okabe-Ito colour-blind-safe palette throughout.

  e2_pgc2010_validation.png   recall/precision vs picker threshold with the shifted-time null, and
                              the pick-minus-predicted residual histograms that carry the signal
  e2_pgc2010_exampleday.png   one ETS day at PGC: 1-8 Hz envelope, picks, Lin/Bostock arrival times
  e2_assoc_rate_vs_tremor.png associated-event rate vs PNSN tremor, vi2010 and vi2017, + the
                              precision/recall trade-off against the station requirement
  e2_location_scatter.png     NonLinLoc epicentres vs Bostock family hypocentres, offsets, depths
  e2_offgt_ets_quiet.png      off-ground-truth stations: pick rate vs tremor + the degenerate-window
                              QC hole
Every panel states units; every claim-bearing panel carries its null or its reference.
"""
from __future__ import annotations

import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
FIG = "lfe_catalog_repair/figures"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
          purple="#CC79A7", yellow="#F0E442", sky="#56B4E9", black="#000000")
MC = {"cascadia": OI["blue"], "jcms": OI["orange"]}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "figure.dpi": 130,
                     "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False})


def savefig(fig, name):
    os.makedirs(FIG, exist_ok=True)
    p = f"{FIG}/{name}"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig] {p}")


# ---------------------------------------------------------------- 1. PGC validation
def fig_validation():
    v = pd.read_csv(f"{RES}/e2_validation_pgc2010.csv")
    ex = pd.read_csv(f"{RES}/e2_resid_excess.csv")
    h = pd.read_csv(f"{RES}/e2_resid_hist.csv")
    fig, ax = plt.subplots(2, 3, figsize=(10.5, 5.6))
    for j, gk in enumerate(["lin_hq", "bostock_near40"]):
        s = v[(v.phase == "S") & (v.tol_s == 3.0) & (v["gt"] == gk)]
        a = ax[j, 0]
        for m, g in s.groupby("model"):
            g = g.sort_values("threshold")
            a.plot(g.threshold, g.recall, "o-", color=MC[m], label=f"{m}")
            a.plot(g.threshold, g.recall_null, "--", color=MC[m], alpha=.6,
                   label=f"{m}, time-shifted null")
        a.set_yscale("log"); a.set_xlabel("LFEDetect S threshold")
        a.set_ylabel("recall (fraction of catalogue\nevents with a pick within 3 s)")
        a.set_title(f"({'ab'[j]}) recall vs {gk}, CN.PGC 2010-05..11")
        a.legend(fontsize=6)
        a = ax[j, 1]
        for m, g in s.groupby("model"):
            g = g.sort_values("threshold")
            a.plot(g.threshold, g.precision, "o-", color=MC[m], label=m)
            a.plot(g.threshold, g.precision_null, "--", color=MC[m], alpha=.6,
                   label=f"{m}, null")
        a.set_yscale("log"); a.set_xlabel("LFEDetect S threshold")
        a.set_ylabel("precision (fraction of picks\nwithin 3 s of a catalogue event)")
        a.set_title(f"({'cd'[j]}) precision vs {gk}")
        a.legend(fontsize=6)
        a = ax[j, 2]
        gk2 = "lin_hq" if gk == "lin_hq" else "bostock_stack"
        for m in ("cascadia", "jcms"):
            hh = h[(h.model == m) & (h.threshold == 0.2) & (h["gt"] == gk2)].sort_values("lag_s")
            if not len(hh):
                continue
            a.step(hh.lag_s, hh["count"], where="mid", color=MC[m], lw=.8, label=m)
            a.axhline(hh.bg_per_bin.iloc[0], color=MC[m], ls=":", lw=.8)
            pk = hh.peak_lag_s.iloc[0]
            a.axvline(pk, color=MC[m], ls="--", lw=.6)
            a.annotate(f"peak {pk:+.2f} s", (pk, hh["count"].max()), fontsize=6,
                       color=MC[m], ha="left", va="top",
                       xytext=(3, -2 - 10 * (m == "jcms")), textcoords="offset points")
        a.set_xlabel("pick time - predicted S arrival (s)")
        a.set_ylabel("pick count per 0.25 s bin")
        a.set_title(f"({'ef'[j]}) residual histogram, {gk2}, thr 0.2\n"
                    "(dotted = chance background)")
        a.legend(fontsize=6)
    fig.suptitle("E2: LFEDetect single-station validation at CN.PGC, 2010-05-01..2010-11-01 "
                 "(Lin 2024 high-quality subset and Bostock 2015)", fontsize=9)
    fig.tight_layout()
    savefig(fig, "e2_pgc2010_validation.png")


# ---------------------------------------------------------------- 2. example day
def fig_exampleday(day="2010-08-24", t_hour=None):
    import obspy
    import e2_validate_pgc2010 as V
    p = f"lfe_catalog_repair/traces/CN.PGC/{day[:4]}/{pd.Timestamp(day).dayofyear:03d}.mseed"
    if not os.path.exists(p):
        print(f"[skip] example day {p} missing")
        return
    st = obspy.read(p)
    st.detrend("demean"); st.filter("bandpass", freqmin=1, freqmax=8, zerophase=True)
    tr = [t for t in st if t.stats.channel.endswith(("N", "1", "E", "2"))]
    tr = (tr or list(st))[0]
    env = np.abs(tr.data.astype(float))
    t0 = tr.stats.starttime.timestamp
    tt = t0 + np.arange(len(env)) / tr.stats.sampling_rate
    pk = V.load_picks()
    bo, lin, C = V.gt_tables()
    off = pd.read_csv(f"{RES}/e2_pgc2010_offsets.csv").set_index("gt").peak_lag_s
    linq = lin[(lin.residual < 0.5) & (lin.N >= 4)]
    gts = {"Lin 2024 (residual<0.5, N>=4)": (linq.pred_1d.to_numpy() + off["lin_hq"], OI["green"]),
           "Bostock 2015 (+12.4 s conv.)": (bo.pred_stack.dropna().to_numpy()
                                            + off["bostock_stack"], OI["purple"])}
    # pick the hour with the most Lin events that day
    d0 = pd.Timestamp(day).timestamp()
    if t_hour is None:
        g = gts["Lin 2024 (residual<0.5, N>=4)"][0]
        g = g[(g >= d0) & (g < d0 + 86400)]
        t_hour = int(np.median(np.floor((g - d0) / 3600))) if len(g) else 12
    fig, ax = plt.subplots(2, 1, figsize=(10.5, 5.0))
    for k, (span, lab) in enumerate([(3600.0, "1 hour"), (600.0, "10 minutes")]):
        a = ax[k]
        ta = d0 + t_hour * 3600 + (0 if k == 0 else 1200)
        tb = ta + span
        m = (tt >= ta) & (tt < tb)
        a.plot((tt[m] - ta), env[m] / max(env[m].max(), 1e-9), color="0.55", lw=.4,
               label=f"{tr.stats.channel} 1-8 Hz |amplitude| (normalised)")
        for mo in ("cascadia", "jcms"):
            s = pk[(pk.model == mo) & (pk.phase == "S") & (pk.thr_first <= 0.3)
                   & (pk.prob >= 0.3)]
            s = s[(s.t >= ta) & (s.t < tb)]
            a.vlines(s.t - ta, 0, s.prob, color=MC[mo], lw=1.0,
                     label=f"{mo} S pick (height = probability)")
        for lab2, (g, col) in gts.items():
            gg = g[(g >= ta) & (g < tb)]
            a.vlines(gg - ta, -0.18, -0.02, color=col, lw=1.0, label=f"{lab2} predicted S")
        a.set_xlim(0, span); a.set_ylim(-0.2, 1.05)
        a.set_xlabel(f"seconds after {pd.Timestamp(ta, unit='s')} UTC")
        a.set_ylabel("pick probability /\nnormalised amplitude")
        a.set_title(f"({'ab'[k]}) CN.PGC {day}, {lab}: picks (thr 0.3) vs catalogue arrival times")
        if k == 0:
            a.legend(fontsize=6, ncol=3, loc="upper right")
    fig.tight_layout()
    savefig(fig, "e2_pgc2010_exampleday.png")


# ---------------------------------------------------------------- 3. association rate vs tremor
def fig_assoc():
    sc = pd.read_csv(f"{RES}/e2_assoc_scores.csv")
    sc["n"] = sc.tag.str.extract(r"_n(\d)$").astype(int)
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 5.6))
    for j, (setname, thr) in enumerate([("vi2010", 0.3), ("vi2017", 0.3)]):
        a = ax[j, 0]
        f = f"{RES}/e2_assoc_daily_{setname}_cascadia_t{thr}_n4.csv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, index_col=0, parse_dates=True)
        a.bar(d.index, d.events, color=OI["blue"], width=.8,
              label="associated events / day (cascadia, thr 0.3, >=4 stations)")
        a2 = a.twinx(); a2.grid(False)
        a2.plot(d.index, d.tremor, color=OI["red"], lw=1.2,
                label="PNSN tremor detections / day in the box")
        a.set_ylabel("associated LFE events / day")
        a2.set_ylabel("PNSN tremor detections / day", color=OI["red"])
        row = sc[(sc["set"] == setname) & (sc.model == "cascadia") & (sc.threshold == thr)
                 & (sc.n == 4)]
        a.set_title(f"({'ab'[j]}) {setname}: event rate vs tremor")
        h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
        a.legend(h1 + h2, l1 + l2, fontsize=6, loc="upper left")
        for lab in a.get_xticklabels():
            lab.set_rotation(30); lab.set_ha("right")
        a = ax[j, 1]
        s = sc[(sc["set"] == setname) & (sc.tol_s == 2.0) & (sc["gt"] == "bostock")]
        if not len(s) and os.path.exists(f"{RES}/e2_assoc_rates.csv"):
            # vi2017 has no ground truth (the Lin catalogue ends 2017-02-21): show instead how the
            # tremor contrast responds to the station requirement -- the only reference available.
            rr = pd.read_csv(f"{RES}/e2_assoc_rates.csv")
            rr = rr[rr.tag.str.startswith(setname) & rr.tag.str.contains("_n")]
            rr["model"] = np.where(rr.tag.str.contains("cascadia"), "cascadia", "jcms")
            rr["thr"] = rr.tag.str.extract(r"_t([\d.]+)_").astype(float)
            rr["n"] = rr.tag.str.extract(r"_n(\d)$").astype(int)
            for (m, t2), g in rr.groupby(["model", "thr"]):
                g = g.sort_values("n")
                a.plot(g.n, g.ets_quiet, "o-", color=MC[m],
                       alpha=0.35 + 0.65 * (t2 == 0.3), label=f"{m}, pick thr {t2}")
            a.set_xlabel("minimum number of stations per event")
            a.set_ylabel("ETS / quiet event-rate ratio")
            a.set_xticks([3, 4, 5])
            a.legend(fontsize=5, ncol=2)
            a.set_title(f"(d) {setname}: tremor contrast vs the station requirement\n"
                        "(no LFE catalogue overlaps this window)")
            continue
        if len(s):
            for m, g in s.groupby("model"):
                for thr2, gg in g.groupby("threshold"):
                    gg = gg.sort_values("n")
                    a.plot(gg.recall, gg.precision, "o-", color=MC[m],
                           alpha=0.35 + 0.3 * (thr2 == 0.3) + 0.35 * (thr2 == 0.5),
                           label=f"{m}, pick thr {thr2}")
                    for _, r in gg.iterrows():
                        a.annotate(f"N>={int(r.n)}", (r.recall, r.precision), fontsize=5,
                                   xytext=(2, 2), textcoords="offset points")
            a.set_xscale("log"); a.set_xlabel("recall vs Bostock 2015 (within 2 s)")
            a.set_ylabel("precision vs Bostock 2015")
            a.legend(fontsize=5, ncol=2)
        a.set_title(f"({'cd'[j]}) {setname}: precision-recall vs the station requirement\n"
                    "(Bostock times shifted by the measured +12.4 s convention offset)")
    fig.suptitle("E2: PyOcto association of LFEDetect picks, southern Vancouver Island",
                 fontsize=9)
    fig.tight_layout()
    savefig(fig, "e2_assoc_rate_vs_tremor.png")


# ---------------------------------------------------------------- 4. locations
def fig_locations():
    f = sorted(glob.glob(f"{RES}/e2_loc_compare_*.csv"))
    if not f:
        print("[skip] no e2_loc_compare_*.csv yet")
        return
    lc = pd.read_csv(f[0])
    tag = os.path.basename(f[0])[len("e2_loc_compare_"):-4]
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    fam = bo.groupby(["lat", "lon", "dep"]).size().reset_index(name="n")
    stas = pd.read_csv(f"{RES}/e2_stations.csv")
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.9))
    a = ax[0]
    hb = a.hexbin(lc.lon, lc.lat, gridsize=45, cmap="Greys", bins="log", mincnt=1)
    plt.colorbar(hb, ax=a, label="located events per cell (log)")
    a.plot(fam.lon, fam.lat, "^", ms=3.5, color=OI["purple"], label="Bostock family hypocentres")
    st = stas[stas.sta.isin(["PGC", "SNB", "VGZ", "LZB", "PFB", "YOUB", "NLLB", "B011",
                             "B926", "LRIV"])]
    a.plot(st.lon, st.lat, "*", ms=7, color=OI["orange"], label="pilot stations")
    a.set_xlabel("longitude (deg)"); a.set_ylabel("latitude (deg)")
    a.set_title(f"(a) NonLinLoc epicentres\n{tag}")
    a.legend(fontsize=6, loc="lower left")
    a = ax[1]
    rng = np.random.default_rng(1)
    R_KM = 111.195; CL = np.cos(np.radians(48.6))
    rlat = rng.uniform(lc.lat.min(), lc.lat.max(), 4000)
    rlon = rng.uniform(lc.lon.min(), lc.lon.max(), 4000)
    Dn = np.array([np.hypot((fam.lon.values - rlon[k]) * R_KM * CL,
                            (fam.lat.values - rlat[k]) * R_KM).min() for k in range(len(rlat))])
    a.hist(Dn, bins=np.arange(0, 121, 4), histtype="step", color="0.4", ls="--", density=True,
           label=f"NULL: uniform random point to nearest family\nmedian {np.median(Dn):.1f} km")
    for col, lab, c in (("dh_nearest_family", "to nearest Bostock family", OI["purple"]),
                        ("dh_lin", "to the time-matched Lin event (<=2 s)", OI["green"]),
                        ("dh_bostock", "to the time-matched Bostock event (<=2 s)", OI["blue"])):
        v = lc[col].dropna()
        if len(v) < 5:
            continue
        a.hist(v, bins=np.arange(0, 121, 4), histtype="step", color=c, density=True,
               label=f"{lab}\nmedian {np.median(v):.1f} km (n={len(v)})")
    a.set_xlabel("horizontal offset (km)"); a.set_ylabel("probability density (1/km)")
    a.set_title("(b) epicentre offsets"); a.legend(fontsize=6)
    a = ax[2]
    a.hist(lc.depth, bins=np.arange(0, 81, 2), histtype="stepfilled", color=OI["sky"],
           alpha=.7, label="NonLinLoc depth")
    a.hist(lc.slab2_km.dropna(), bins=np.arange(0, 81, 2), histtype="step", color=OI["red"],
           label="Slab2 interface under the epicentre")
    a.axvspan(25, 45, color=OI["yellow"], alpha=.25, label="25-45 km interface band")
    a.set_xlabel("depth (km)"); a.set_ylabel("events per 2 km bin")
    a.set_title(f"(c) depth: median {lc.depth.median():.1f} km, formal ERZ median "
                f"{lc.err_z_km.median():.1f} km\n(S-only picks: depth is not resolved)")
    a.legend(fontsize=6)
    fig.suptitle("E2: NonLinLoc locations of the PyOcto-associated LFE events (Crosson 1-D, "
                 "S-weighted)", fontsize=9)
    fig.tight_layout()
    savefig(fig, "e2_location_scatter.png")


# ---------------------------------------------------------------- 5. off ground truth
def fig_offgt():
    if not os.path.exists(f"{RES}/e2_offgt_summary.csv"):
        print("[skip] no e2_offgt_summary.csv yet")
        return
    su = pd.read_csv(f"{RES}/e2_offgt_summary.csv")
    da = pd.read_csv(f"{RES}/e2_offgt_daily.csv", parse_dates=["date"])
    stas = [s for s in ("WISH", "B943", "YACT") if (su.station == s).any()]
    ncol = len(stas)
    fig, ax = plt.subplots(2, max(ncol, 1), figsize=(3.6 * max(ncol, 1), 5.4), squeeze=False)
    for j, sta in enumerate(stas):
        a = ax[0, j]
        d = da[(da.station == sta) & (da.model == "cascadia") & (da.threshold == 0.3)]
        a.plot(d.date, d.picks, lw=.7, color=OI["blue"], label="S picks / day (cascadia, thr 0.3)")
        a2 = a.twinx(); a2.grid(False)
        a2.plot(d.date, d.tremor, lw=.9, color=OI["red"], alpha=.8,
                label="PNSN tremor / day within 40 km")
        r = su[(su.station == sta) & (su.model == "cascadia") & (su.threshold == 0.3)]
        txt = (f"ETS/quiet {r.ets_quiet_ratio.iloc[0]:.2f}\nrho {r.spearman_rho.iloc[0]:+.2f}"
               if len(r) else "")
        a.set_title(f"({'abc'[j]}) {sta}: pick rate vs tremor\n{txt}")
        a.set_ylabel("S picks / day"); a2.set_ylabel("tremor / day", color=OI["red"])
        for lab in a.get_xticklabels():
            lab.set_rotation(30); lab.set_ha("right")
        if j == 0:
            h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
            a.legend(h1 + h2, l1 + l2, fontsize=6, loc="upper left")
        a = ax[1, j]
        if os.path.exists(f"{RES}/e2_offgt_degen.csv"):
            dg = pd.read_csv(f"{RES}/e2_offgt_degen.csv")
            g = dg[(dg.station == sta) & (dg.model == "cascadia")].sort_values("threshold")
            if len(g):
                w = 0.35
                x = np.arange(len(g))
                a.bar(x - w / 2, 100 * g.pick_gap, w, color=OI["blue"],
                      label="picks whose 60 s window is gap-affected")
                a.bar(x + w / 2, 100 * g.win_gap, w, color="0.6",
                      label="all 60 s windows of the day (baseline)")
                a.set_xticks(x); a.set_xticklabels([f"thr {t}" for t in g.threshold])
                a.set_ylabel("percent of picks / windows")
                a.set_title(f"({'def'[j]}) {sta}: degenerate windows -- picks are DEPLETED\n"
                            f"(0 % of picks in a fully empty window; day coverage "
                            f"{g.day_coverage_mean.iloc[0]:.4f})")
                if j == 0:
                    a.legend(fontsize=6)
        else:
            a.set_axis_off()
    fig.suptitle("E2: off-ground-truth stations, 2017-2018 (no LFE catalogue exists there)",
                 fontsize=9)
    fig.tight_layout()
    savefig(fig, "e2_offgt_ets_quiet.png")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for name, fn in (("validation", fig_validation), ("exampleday", fig_exampleday),
                     ("assoc", fig_assoc), ("locations", fig_locations), ("offgt", fig_offgt)):
        if which in ("all", name):
            try:
                fn()
            except Exception as e:
                import traceback
                print(f"[fail] {name}: {e}")
                traceback.print_exc()
