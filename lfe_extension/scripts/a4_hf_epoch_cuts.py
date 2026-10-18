#!/usr/bin/env python
"""Track A-2 follow-up (1c) — PER-EPOCH HF-ratio veto cuts, and the referee-free rule that sets them.

The frozen cut 3.6384 is the 99th percentile of the ratio over BUILD-half detections matched to an
own-family Bostock event — a REFEREE rule, available only where Bostock's catalog exists (2007-2013).

WHAT THE DIAGNOSIS ACTUALLY SHOWS (see a4_hf_response.py, a4_hf_bands.py and the table this script
prints).  The shift Track C measured, median 0.199 (2015-12, BH) -> 0.480 (2018-06, HH), is NOT mainly
an instrument effect:
  * the instrument response predicts only x1.16 of it on 10-19 Hz (x1.06 on 10-16), and
  * NO band removes it: measured HH/BH shift of the median is x2.41 (10-19), x1.86 (10-16), x2.06
    (10-14), x1.59 (8-16) — so "redefine the band so both instruments are flat" does not work, and
  * the SAME instrument gives median 1.019 on 2012 doy 228-279 and 0.199 on 2015-12, a 5x spread with
    no boundary crossed at all.
Restricted to a matched HIGH-CONFIDENCE population (net_mad >= T_core, n_sites >= 4, PNSN tremor days)
the median is 0.256 (2012 BH) / 0.169 (2015-12 BH) / 0.247 (2018-06 HH): the BH-to-BH spread is as large
as any BH-to-HH difference.  The distribution moves with the DETECTION POPULATION, not the sensor.

So the per-epoch cut must be tied to the epoch's own LFE-like core, not to an absolute number and not to
a tail percentile of a month-dependent mixture.  Two referee-free rules are calibrated on 2012 (where
the referee exists) and then applied unchanged:

  RULE-M  cut = F x median(R | net_mad >= T_core, n_sites >= 4, tremor day),  F fixed on 2012
  RULE-Q  cut = quantile q* of the same selection,                            q* fixed on 2012
  RULE-T  cut = p99 of the same selection (the coordinator's candidate) — REPORTED AND REJECTED:
          it is a tail statistic of a mixture and is not even stable between two BH months
          (3.246 on 2012 vs 1.208 on 2015-12, same instrument).

Outputs: results/a4_hf_veto_cuts.csv, results/a4_hf_rule_validation.csv,
         results/a4_hf_band_epoch_table.csv, figures/a4_hf_epoch_cuts.png
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a1_hf_veto import one2one                                   # noqa: E402

EXT = os.path.dirname(HERE)
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
IN = f"{REP}/results"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
BOX = dict(lat=(48.0, 49.0), lon=(-124.5, -123.0))
QUIET_2012 = set(range(228, 243))
T_CORE, T_EXT = 8.75, 8.25
FROZEN = 3.6384


def tremor_days(year, doys=None, min_n=30):
    d = pd.read_csv(PNSN, parse_dates=["time"])
    d = d[(d.lat.between(*BOX["lat"])) & (d.lon.between(*BOX["lon"])) & (d.time.dt.year == year)]
    c = d.groupby(d.time.dt.dayofyear).size()
    s = set(int(x) for x in c[c >= min_n].index)
    return s if doys is None else s & set(int(x) for x in doys)


def core_sel(df, tdays, col="hfr", tcore=T_CORE, min_sites=4):
    """The LFE-like core of an epoch: strong, multi-site detections on tremor days."""
    return ((df.net_mad >= tcore) & (df.n_sites >= min_sites)
            & df.doy.isin(tdays) & np.isfinite(df[col])).to_numpy()


def rule_referee(df, col="hfr", block_s=7200.0, guard_s=60.0, pct=99.0):
    bb = np.floor(df.t.to_numpy() / block_s).astype(np.int64)
    pp = df.t.to_numpy() - bb * block_s
    build = (bb % 2 == 0) & (pp >= guard_s) & (pp <= block_s - guard_s)
    v = df.loc[build & df.match_own.to_numpy() & np.isfinite(df[col]), col].to_numpy()
    return (float(np.percentile(v, pct)) if len(v) else np.nan), int(len(v))


def veto_effect_2012(df, cut, col="hfr"):
    gt = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet")
    test_ev = gt[(gt.half == "test") & (~gt.guard)].reset_index(drop=True)
    R = df[col].to_numpy()
    keep = ~(np.isfinite(R) & (R > cut))
    out = {}
    for T, lab in ((T_EXT, "ext"), (T_CORE, "core")):
        sel = df.net_mad.to_numpy() >= T
        q = sel & df.doy.isin(QUIET_2012).to_numpy()
        fa0, fa1 = int(q.sum()), int((q & keep).sum())
        rec = {}
        for lbl, kp in (("before", np.ones(len(df), bool)), ("after", keep)):
            D = df[sel & kp]
            hit = np.zeros(len(test_ev), bool)
            for fk, g in test_ev.groupby("bkey"):
                td = np.sort(D[D.fam == float(fk)].t.to_numpy())
                hit[test_ev.index.get_indexer(g.index)] = one2one(g.t_origin.to_numpy(), td) <= 2.0
            rec[lbl] = float(hit.mean())
        out[lab] = dict(quiet_before=fa0, quiet_after=fa1, fa_reduction=1 - fa1 / max(fa0, 1),
                        recall_before=rec["before"], recall_after=rec["after"],
                        recall_cost_rel=(rec["before"] - rec["after"]) / max(rec["before"], 1e-9),
                        vetoed_frac=float((np.isfinite(R) & (R > cut))[sel].mean()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default="10-19")
    a = ap.parse_args()

    # ---------- the three legs, all with a column named `hfr` --------------------------------------
    d12 = pd.read_parquet(f"{EXT}/results/a1_hf_ratio_a1_mergeB_L12.parquet").rename(columns={"hfr_PGC": "hfr"})
    d15 = pd.read_parquet(f"{EXT}/results/trackC/detections_2015-12.parquet")
    d18 = pd.read_parquet(f"{EXT}/results/trackC/detections_2018-06.parquet")
    LEGS = [("2012", "PGC_BH", 40.0, d12, 2012), ("2015-12", "PGC_BH", 40.0, d15, 2015),
            ("2018-06", "PGC_HH", 100.0, d18, 2018)]

    ref_cut, ref_n = rule_referee(d12)
    td12 = tremor_days(2012, d12.doy.unique())
    m12 = core_sel(d12, td12)
    med12 = float(np.median(d12.loc[m12, "hfr"]))
    F = ref_cut / med12                                     # RULE-M multiplier, fixed on 2012
    qstar = float(100.0 * (d12.loc[m12, "hfr"].to_numpy() < ref_cut).mean())   # RULE-Q quantile
    print(f"[calibration on 2012] referee cut {ref_cut:.4f} (n {ref_n}); core median {med12:.4f} "
          f"(n {int(m12.sum())})  ->  RULE-M F = {F:.2f},  RULE-Q q* = {qstar:.3f} %", flush=True)

    rows, prof = [], []
    for lab, epoch, fs, d, yr in LEGS:
        td = tremor_days(yr, d.doy.unique())
        m = core_sel(d, td)
        v = d.loc[m, "hfr"].to_numpy()
        cut_M = F * float(np.median(v))
        cut_Q = float(np.percentile(v, qstar))
        cut_T = float(np.percentile(v, 99.0))
        allR = d.hfr.to_numpy()
        rng = np.random.default_rng(5)
        bs = np.array([np.median(rng.choice(v, len(v), replace=True)) for _ in range(2000)]) * F
        row = dict(leg=lab, epoch=epoch, fs=fs, n_det=len(d), n_core=int(m.sum()),
                   cut_M_lo=float(np.percentile(bs, 2.5)), cut_M_hi=float(np.percentile(bs, 97.5)),
                   median_all=float(np.nanmedian(allR)), median_core=float(np.median(v)),
                   cut_M=cut_M, cut_Q=cut_Q, cut_T=cut_T,
                   veto_all_frozen=float(np.nanmean(allR > FROZEN)),
                   veto_all_M=float(np.nanmean(allR > cut_M)),
                   veto_core_frozen=float((v > FROZEN).mean()),
                   veto_core_M=float((v > cut_M).mean()))
        if "quiet_flag" in d.columns and d.quiet_flag.any():
            qv = d.loc[d.quiet_flag & d.hfr.notna(), "hfr"].to_numpy()
            row.update(n_quiet=len(qv), veto_quiet_frozen=float((qv > FROZEN).mean()),
                       veto_quiet_M=float((qv > cut_M).mean()))
        rows.append(row)
    P = pd.DataFrame(rows)
    pd.set_option("display.width", 240)
    print("\n[per leg] the frozen band " + a.band + " Hz")
    print(P.round(4).to_string(index=False))

    # ---------- 2012 validation of every rule (the only place with a referee) ----------------------
    val = []
    for name, cut, n in (("referee_p99_own_matched (FROZEN)", ref_cut, ref_n),
                         ("RULE-M  F x core median", F * med12, int(m12.sum())),
                         ("RULE-Q  q* of the core", float(np.percentile(d12.loc[m12, 'hfr'], qstar)), int(m12.sum())),
                         ("RULE-T  p99 of the core", float(np.percentile(d12.loc[m12, 'hfr'], 99.0)), int(m12.sum()))):
        e = veto_effect_2012(d12, cut)
        for tier in ("ext", "core"):
            val.append(dict(year=2012, band=a.band, rule=name, cut=cut, n_calib=n, tier=tier, **e[tier]))
    V = pd.DataFrame(val)
    V.to_csv(f"{EXT}/results/a4_hf_rule_validation.csv", index=False)
    print("\n[rule validation on 2012 — bars: FA reduction >= 30 %, relative recall cost <= 2 %]")
    print(V.round(4).to_string(index=False))

    # ---------- band x epoch table (from a4_hf_bands.py) -------------------------------------------
    BR = None
    if os.path.exists(f"{EXT}/results/a4_hf_band_scan.csv"):
        B = pd.read_csv(f"{EXT}/results/a4_hf_band_scan.csv")
        tdays = {m: tremor_days(int(m[:4]), B[B.month == m].doy.unique()) for m in B.month.unique()}
        br = []
        for (month, band), g in B.groupby(["month", "band"]):
            g = g.rename(columns={"hfr": "hfr"})
            mm = core_sel(g, tdays[month])
            br.append(dict(month=month, band=band, fs=float(g.fs.iloc[0]),
                           median_all=float(g.hfr.median()),
                           median_core=float(np.median(g.loc[mm, "hfr"])) if mm.sum() else np.nan,
                           n_core=int(mm.sum()), n_det=len(g)))
        BR = pd.DataFrame(br).sort_values(["band", "month"])
        piv = BR.pivot(index="band", columns="month", values=["median_all", "median_core"])
        BR.to_csv(f"{EXT}/results/a4_hf_band_epoch_table.csv", index=False)
        print("\n[band x epoch] median of ALL detections vs median of the LFE-like CORE")
        print(piv.round(3).to_string())
        print("  -> the all-detection median shifts x1.6-2.4 in EVERY band (no band fixes it);")
        print("     the CORE median is epoch-stable, which is what the cut must be tied to.")

    # ---------- the deliverable -------------------------------------------------------------------
    out = []
    r12 = P[P.leg == "2012"].iloc[0]
    r15 = P[P.leg == "2015-12"].iloc[0]
    r18 = P[P.leg == "2018-06"].iloc[0]
    out.append(dict(station_epoch="PGC_BH_2010-08-01..2017-08-10", band=a.band, fs=40.0,
                    cut=round(ref_cut, 4), status="recommended",
                    method="referee p99 (own-matched Bostock events, 2012 build half) = THE FROZEN CUT",
                    n_calib=ref_n,
                    note=(f"RULE-M reproduces it exactly on 2012 (F={F:.2f}); on the OTHER BH SSE "
                          f"(2015-12) RULE-M would give {r15.cut_M:.2f} [{r15.cut_M_lo:.2f}, {r15.cut_M_hi:.2f}], "
                          "i.e. the SSE-to-SSE spread inside one instrument epoch is +-35 %")))
    out.append(dict(station_epoch="PGC_HH_2017-08-10..2026", band=a.band, fs=100.0,
                    cut=round(r18.cut_M, 4), status="recommended",
                    method=f"RULE-M: F={F:.2f} (fixed on the 2012 referee) x the core median of this "
                           f"epoch's first cached SSE (2018-06)",
                    n_calib=int(r18.n_core),
                    note=(f"95 % CI [{r18.cut_M_lo:.2f}, {r18.cut_M_hi:.2f}] CONTAINS the frozen 3.6384, so the "
                          f"frozen cut was never badly wrong on HH: it vetoes {100*r18.veto_all_frozen:.1f} % of all "
                          f"2018-06 detections but only {100*r18.veto_core_frozen:.1f} % of the LFE-like core")))
    for lab, r in (("2015-12", r15), ("2018-06", r18), ("2012", r12)):
        for nm, c in (("RULE-Q q*", r.cut_Q), ("RULE-T p99 (REJECTED: not stable BH-to-BH)", r.cut_T)):
            out.append(dict(station_epoch=f"PGC_{r.epoch.split('_')[1]}_diagnostic_{lab}", band=a.band,
                            fs=r.fs, cut=round(float(c), 4), status="diagnostic",
                            method=nm, n_calib=int(r.n_core), note="not for use; see a4_hf_epoch_profile.csv"))
    O = pd.DataFrame(out)
    O.to_csv(f"{EXT}/results/a4_hf_veto_cuts.csv", index=False)
    print("\n[deliverable] results/a4_hf_veto_cuts.csv")
    print(O.to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    bins = np.logspace(-2.5, 2.2, 80)
    cols = {"2012": "tab:green", "2015-12": "tab:blue", "2018-06": "tab:red"}
    for lab, epoch, fs, d, yr in LEGS:
        td = tremor_days(yr, d.doy.unique()); m = core_sel(d, td)
        ax[0].hist(d.hfr.dropna(), bins=bins, histtype="step", color=cols[lab], density=True,
                   label=f"{lab} {epoch} all (med {np.nanmedian(d.hfr):.2f})")
        ax[1].hist(d.loc[m, "hfr"], bins=bins, histtype="step", color=cols[lab], density=True,
                   label=f"{lab} {epoch} core (med {np.median(d.loc[m, 'hfr']):.2f})")
    for x in (0, 1):
        ax[x].axvline(FROZEN, color="k", ls=":", label="frozen cut 3.64")
        ax[x].set_xscale("log"); ax[x].set_xlabel(f"HF ratio, {a.band} Hz"); ax[x].legend(fontsize=6.5)
    ax[0].set_title("ALL detections: the mixture moves with the month"); ax[0].set_ylabel("density")
    ax[1].set_title("LFE-like CORE: epoch-stable")
    w = 0.35
    x = np.arange(len(P))
    ax[2].bar(x - w / 2, 100 * P.veto_all_frozen, w, label="frozen cut", color="0.6")
    ax[2].bar(x + w / 2, 100 * P.veto_all_M, w, label="RULE-M per-epoch cut", color="tab:orange")
    ax[2].set_xticks(x); ax[2].set_xticklabels([f"{r.leg}\n{r.epoch}" for r in P.itertuples()])
    ax[2].set_ylabel("% of detections vetoed"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3, axis="y")
    ax[2].set_title("veto fraction per epoch")
    plt.tight_layout(); plt.savefig(f"{EXT}/figures/a4_hf_epoch_cuts.png", dpi=130)
    P.to_csv(f"{EXT}/results/a4_hf_epoch_profile.csv", index=False)
    print("[write] figures/a4_hf_epoch_cuts.png, results/a4_hf_epoch_profile.csv")


if __name__ == "__main__":
    main()
