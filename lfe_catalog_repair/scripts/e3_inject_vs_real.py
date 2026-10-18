#!/usr/bin/env python
"""E3 step 3, interpretation — put the INJECTED recall(SNR) and the REAL recall(magnitude) on one
axis, and work out how far the false-alarm rate F must fall for the N>=4 chance rate to become
negligible.

Three pieces.

1. The REAL curve, recomputed here rather than quoted: PGC S picks vs Bostock events with
   epi_PGC < 40 km, in six equal-count magnitude bins, scored by the chance-free +-1 s PEAK EXCESS
   around the family's own measured time offset (E3 step 5), at pick threshold 0.1 -- the same
   estimator and threshold as the audit's 0.056/0.098/0.142/0.179/0.241/0.281.

2. THE CONVERSION, stated explicitly.  D3 measured that the MEDIAN Lin LFE at PGC has amplitude
   s/n = 0.078 of the ambient 1-8 Hz noise, from rms(signal)/rms(noise) = sqrt(1 + N (s/n)^2) on a
   130-event single-cell stack at L = 6 s.  That is an RMS-over-6-s measure, which is exactly the
   `snr_rms` recorded for every injection.  Bostock's LFE magnitudes are moment magnitudes from
   template amplitudes, so amplitude scales as the moment: s/n(M) = 0.078 * 10^(1.5 (M - M_ref)),
   with M_ref = the median magnitude of the Bostock < 40 km population at PGC.  Two caveats carried
   into the plot: (a) the 0.078 anchor is a median LIN event, taken here as the median BOSTOCK
   event of the same population; (b) D3's mixed-source stack gave 0.032, a lower bound, so the
   anchor is uncertain by a factor ~2.4 -- shown as a horizontal error band, not hidden.

3. F-SCALING.  The step-1 null rate is measured at three pick thresholds, i.e. at three values of
   the per-station false-alarm rate F.  Fitting log(null events/day) = a + p log(network S
   picks/day) gives the chance exponent p; the real rate is (observed - null); and the factor by
   which F must fall for null <= 10 % of real follows as (null/(0.1*real))^(1/p).

Outputs results/e3_inject_vs_real.csv, results/e3_F_scaling.csv
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_resid_hist as RH  # noqa: E402

SN_ANCHOR = 0.078          # D3: median Lin LFE at PGC, rms(6 s)/rms(noise)
SN_ANCHOR_LO = 0.032       # D3 mixed-source stack, a lower bound


def epoch(s):
    return pd.to_datetime(s).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9


def peak_excess(gt, pk, half=1.0, maxlag=30.0, binw=0.25, bg_excl=5.0):
    if len(gt) < 10 or len(pk) == 0:
        return np.nan, np.nan
    d = RH.pairs_within(np.sort(gt), np.sort(pk), maxlag)
    e = np.arange(-maxlag, maxlag + binw / 2, binw)
    h, _ = np.histogram(d, bins=e)
    c = 0.5 * (e[1:] + e[:-1])
    bg = float(np.median(h[np.abs(c) > bg_excl]))
    m = np.abs(c) <= half
    exc = float(h[m].sum() - bg * m.sum())
    return exc, exc / np.sqrt(max(h[m].sum(), 1))


def main():
    # ---------------- 1. the real recall-vs-magnitude curve ----------------
    pk = pd.read_parquet(f"{RES}/e2_picks_pgc2010_cascadia.parquet")
    pk = pk[(pk.phase == "S") & (pk.sta == "PGC")]
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    acc = fo[fo.accept.astype(bool)]
    med = float(acc.offset_s.median())
    offmap = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    bo = bo[bo.epi_PGC < 40.0].copy()
    bo["t_pred"] = (epoch(bo.time) + bo.bkey.astype(int).map(offmap).fillna(med).to_numpy()
                    + bo.ttS_PGC.to_numpy())
    M_ref = float(bo.mag.median())
    print(f"[real] {len(bo)} Bostock events with epi_PGC < 40 km; magnitude "
          f"{bo.mag.min():.2f}..{bo.mag.max():.2f}, median M_ref = {M_ref:.3f}")
    q = np.quantile(bo.mag, np.linspace(0, 1, 7))
    q[0] -= 1e-6; q[-1] += 1e-6
    rows = []
    for thr in (0.1, 0.3):
        # a pick exists at threshold t iff thr_first <= t <= prob (e2_pick docstring)
        s = np.sort(epoch(pk.time[(pk.thr_first <= thr) & (pk.prob >= thr)]))
        for i in range(6):
            g = bo[(bo.mag > q[i]) & (bo.mag <= q[i + 1])]
            exc, sig = peak_excess(g.t_pred.to_numpy(), s)
            Mm = float(g.mag.median())
            rows.append(dict(kind="real", threshold=thr, bin=i + 1, m_lo=q[i], m_hi=q[i + 1],
                             mag_median=Mm, n_gt=len(g), excess=exc, sigma=sig,
                             recall=exc / len(g) if len(g) else np.nan,
                             snr_rms=SN_ANCHOR * 10 ** (1.5 * (Mm - M_ref)),
                             snr_rms_lo=SN_ANCHOR_LO * 10 ** (1.5 * (Mm - M_ref))))
    real = pd.DataFrame(rows)
    print("\n[real curve, recomputed with the per-family Bostock offsets]")
    print(real[["threshold", "bin", "m_lo", "m_hi", "mag_median", "n_gt", "excess", "sigma",
                "recall", "snr_rms"]].round(4).to_string(index=False))

    # ---------------- the injected curve ----------------
    inj = pd.read_csv(f"{RES}/e3_inject_recall.csv")
    ic = inj[(inj.tol == 1.0)].groupby(["station", "model", "threshold", "level"]).agg(
        n=("n_inj", "sum") if "n_inj" in inj.columns else ("recall", "size"),
        snr_rms=("snr_rms", "median"), recall=("recall", "mean"),
        baseline=("baseline", "mean"), bias_s=("bias_s", "mean"),
        scatter_s=("scatter_s", "mean")).reset_index()
    ic["kind"] = "injected"
    print("\n[injected curve, tol 1 s, PGC cascadia]")
    print(ic[(ic.station == "PGC") & (ic.model == "cascadia") & ic.threshold.isin([0.1, 0.3])]
          [["threshold", "level", "snr_rms", "recall", "baseline", "bias_s", "scatter_s"]]
          .round(4).to_string(index=False))
    both = pd.concat([real, ic], ignore_index=True)
    both.to_csv(f"{RES}/e3_inject_vs_real.csv", index=False)

    # ---------------- 3. F scaling ----------------
    nr = pd.read_csv(f"{RES}/e3_nullevents_rates.csv")
    nr = nr[(nr.set == "vi2010") & (nr.model == "cascadia")]
    fa = pd.read_csv(f"{RES}/e3_inject_falsealarm.csv")
    fa = fa[(fa.run == "base") & (fa.model == "cascadia") & (fa.station == "PGC")]
    Fq = fa.groupby("threshold").n_S.mean()
    # network S picks/day at each threshold, from the pick table itself
    pv = pd.read_parquet(f"{RES}/e2_picks_vi2010_cascadia.parquet")
    pv = pv[pv.phase == "S"]
    nday = pv.day.nunique()
    frows = []
    for thr in (0.2, 0.3, 0.5):
        npk = int(((pv.thr_first <= thr) & (pv.prob >= thr)).sum())
        for n in (3, 4, 5):
            r = nr[(nr.threshold == thr) & (nr.n_min == n)]
            if not len(r):
                continue
            r = r.iloc[0]
            frows.append(dict(threshold=thr, n_min=n, net_picks_per_day=npk / nday,
                              pgc_quietday_picks=float(Fq.get(thr, np.nan)),
                              obs_per_day=r.obs_per_day, null_per_day=r.null_all_per_day,
                              real_per_day=r.obs_per_day - r.null_all_per_day,
                              ratio=r.ratio_all))
    F = pd.DataFrame(frows)
    for n, g in F.groupby("n_min"):
        p, a = np.polyfit(np.log(g.net_picks_per_day), np.log(g.null_per_day), 1)
        need = (g.null_per_day / (0.1 * g.real_per_day)).clip(lower=1.0)
        F.loc[g.index, "chance_exponent_p"] = p
        F.loc[g.index, "null_reduction_needed"] = need
        F.loc[g.index, "F_reduction_needed"] = need ** (1.0 / p)
    F.to_csv(f"{RES}/e3_F_scaling.csv", index=False)
    print("\n[F scaling: how far the per-station false-alarm rate must fall for the chance rate "
          "to be <= 10 % of the real rate]")
    print(F.round(3).to_string(index=False))
    print(f"\n[write] {RES}/e3_inject_vs_real.csv, {RES}/e3_F_scaling.csv")


def _run():
    main()
    burst_split()


def burst_split():
    """E3 step 3, follow-up: the injected curve sits FAR BELOW the real recall-vs-magnitude curve
    at matched s/n.  Neither of the audit's two branches covers that direction, and the obvious
    candidate explanation is that the injections are ISOLATED single LFEs in quiet noise whereas
    real LFEs arrive in bursts, so LFEDetect's 60 s window usually contains several of them and the
    operative signal-to-noise is the burst, not the event.  This tests it with no new run: split
    the real Bostock population by how many OTHER Bostock events (any family, epi_PGC < 40 km) lie
    within +-30 s, and recompute the same +-1 s peak-excess recall in matched magnitude bins.
    """
    pk = pd.read_parquet(f"{RES}/e2_picks_pgc2010_cascadia.parquet")
    pk = pk[(pk.phase == "S") & (pk.sta == "PGC")]
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    acc = fo[fo.accept.astype(bool)]
    med = float(acc.offset_s.median())
    omap = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    bo = bo[bo.epi_PGC < 40.0].copy()
    bo["t_pred"] = (epoch(bo.time) + bo.bkey.astype(int).map(omap).fillna(med).to_numpy()
                    + bo.ttS_PGC.to_numpy())
    t = np.sort(bo.t_pred.to_numpy())
    nb = (np.searchsorted(t, t + 30.0) - np.searchsorted(t, t - 30.0) - 1)
    bo = bo.assign(t_pred=t, n_within_30s=nb)
    print(f"\n[burst] neighbours within +-30 s: median {np.median(nb):.0f}, "
          f"{float((nb <= 1).mean()):.2f} isolated (<=1), {float((nb >= 10).mean()):.2f} dense "
          f"(>=10)")
    M_ref = float(bo.mag.median())
    q = np.quantile(bo.mag, np.linspace(0, 1, 4)); q[0] -= 1e-6; q[-1] += 1e-6
    rows = []
    for thr in (0.1,):
        s = np.sort(epoch(pk.time[(pk.thr_first <= thr) & (pk.prob >= thr)]))
        for lab, sel in (("isolated (<=1 in +-30 s)", bo.n_within_30s <= 1),
                         ("moderate (2-9)", bo.n_within_30s.between(2, 9)),
                         ("dense (>=10)", bo.n_within_30s >= 10)):
            for i in range(3):
                g = bo[sel & (bo.mag > q[i]) & (bo.mag <= q[i + 1])]
                if len(g) < 50:
                    continue
                exc, sig = peak_excess(g.t_pred.to_numpy(), s)
                Mm = float(g.mag.median())
                rows.append(dict(group=lab, mag_bin=i + 1, mag_median=Mm, n_gt=len(g),
                                 n_within_30s_med=float(g.n_within_30s.median()),
                                 excess=exc, sigma=sig, recall=exc / len(g),
                                 snr_rms=SN_ANCHOR * 10 ** (1.5 * (Mm - M_ref))))
    B = pd.DataFrame(rows)
    B.to_csv(f"{RES}/e3_real_recall_by_burst.csv", index=False)
    print("[burst] real +-1 s peak-excess recall at PGC thr 0.1, by local Bostock density")
    print(B.round(4).to_string(index=False))
    print(f"[write] {RES}/e3_real_recall_by_burst.csv")


if __name__ == "__main__":
    _run()
