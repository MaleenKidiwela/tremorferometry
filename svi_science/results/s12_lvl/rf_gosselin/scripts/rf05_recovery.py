#!/usr/bin/env python
"""rf05_recovery: (a) recovery curve of the ETS-related change — sliding windows of RF quantities vs days from the
ETS midpoint (centres -60..+120 d every 3 d, widths 7 and 15 d), anomaly = window mean - mean of all RFs of the
analysed period, with a scrambled-control 95 % band (RF dates shuffled among events, 5000 draws);
(b) local slip-front timing: per station and episode (2010+, PNSN tremor), the time tremor first comes within 30 km of
the station and the median time of tremor within 30 km, relative to the GPS midpoint.
Usage: python -P rf05_recovery.py --tag repro [--stations ...] [--ets SVI_ALBH]
Outputs: tables/<tag>/recovery_<STA>.csv, figures/<tag>/recovery_<STA>.png, tables/slipfront_<ets>.csv
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
W = f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin"
Q = [("cc_sv_clean", "ΔCC SV"), ("cc_sh_clean", "ΔCC SH"), ("dvs_top", "ΔdVs top (km/s)"), ("dvs_bot", "ΔdVs bottom (km/s)")]


def curve(rel, x, centres, w):
    ok = np.isfinite(x)
    mu = x[ok].mean()
    out = np.full(len(centres), np.nan); n = np.zeros(len(centres), int)
    for i, c in enumerate(centres):
        m = ok & (rel >= c - w / 2) & (rel < c + w / 2)
        n[i] = m.sum()
        if n[i] >= 3:
            out[i] = x[m].mean() - mu
    return out, n


def run(tag, key, nboot=5000):
    f = f"{W}/tables/{tag}/{key}_events.csv"
    if not os.path.exists(f):
        return
    M = pd.read_csv(f)
    rel = M.rel_days.values
    centres = np.arange(-60, 121, 3)
    rng = np.random.default_rng(7)
    rows = {"centre": centres}
    fig, axs = plt.subplots(2, 2, figsize=(12, 7), sharex=True)
    for ax, (q, lab) in zip(axs.flat, Q):
        x = M[q].values.astype(float)
        for w, col in ((15, "k"), (7, "tab:orange")):
            c, n = curve(rel, x, centres, w)
            rows[f"{q}_w{w}"] = c; rows[f"{q}_w{w}_n"] = n
            if w == 15:
                ok = np.isfinite(x); xv = x[ok]; rv = rel[ok]
                B = np.empty((nboot, len(centres)))
                for b in range(nboot):
                    B[b] = curve(rv, rng.permutation(xv), centres, w)[0]
                lo, hi = np.nanpercentile(B, 2.5, 0), np.nanpercentile(B, 97.5, 0)
                rows[f"{q}_w15_lo95"] = lo; rows[f"{q}_w15_hi95"] = hi
                ax.fill_between(centres, lo, hi, color="0.85", label="scrambled control 95 %")
            ax.plot(centres, c, "-o", ms=2.5, color=col, lw=1.5 if w == 15 else 0.8, label=f"{w}-day window")
        ax.axvline(0, color="0.5", ls=":"); ax.axhline(0, color="0.5", lw=0.5)
        ax.set_title(lab, fontsize=9)
    axs[0, 0].legend(fontsize=7)
    for ax in axs[1]:
        ax.set_xlabel("window centre, days from ETS midpoint")
    fig.suptitle(f"{key} {tag}: recovery curve (anomaly vs mean of all RFs)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{W}/figures/{tag}/recovery_{key}.png", dpi=110); plt.close(fig)
    pd.DataFrame(rows).to_csv(f"{W}/tables/{tag}/recovery_{key}.csv", index=False)


def slipfront(ets, stations):
    T = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    T["t"] = pd.to_datetime(T.time).astype("datetime64[ns]")
    E = pd.read_csv(f"{W}/tables/ets_mid_{ets}.csv"); E["mid"] = pd.to_datetime(E["mid"])
    S = pd.read_csv(f"{ROOT}/data/rf_teleseis/stations.csv")
    rows = []
    for r in S.itertuples():
        if f"{r.net}.{r.sta}" not in stations:
            continue
        d = np.hypot((T.lon - r.lon) * 111.19 * np.cos(np.deg2rad(r.lat)), (T.lat - r.lat) * 111.19)
        Tn = T[d < 30]
        for e in E.itertuples():
            if e.mid < pd.Timestamp("2009-06-01"):
                continue
            m = Tn[(Tn.t >= e.mid - pd.Timedelta(days=40)) & (Tn.t <= e.mid + pd.Timedelta(days=40))]
            if len(m) < 10:
                rows.append(dict(station=f"{r.net}.{r.sta}", episode=e.episode, mid=e.mid.date(), n_trem30=len(m)))
                continue
            rr = (m.t - e.mid) / pd.Timedelta(days=1)
            rows.append(dict(station=f"{r.net}.{r.sta}", episode=e.episode, mid=e.mid.date(), n_trem30=len(m),
                             first_d=round(float(np.percentile(rr, 2)), 1), median_d=round(float(np.median(rr)), 1),
                             last_d=round(float(np.percentile(rr, 98)), 1)))
    F = pd.DataFrame(rows); F.to_csv(f"{W}/tables/slipfront_{ets}.csv", index=False)
    print(F.groupby("station")[["first_d", "median_d", "last_d", "n_trem30"]].median().to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="repro"); ap.add_argument("--ets", default="SVI_ALBH")
    ap.add_argument("--stations", default="CN.PGC,CN.SNB,CN.VGZ,CN.LZB,CN.PFB")
    ap.add_argument("--slipfront", action="store_true")
    a = ap.parse_args()
    st = a.stations.split(",")
    if a.slipfront:
        slipfront(a.ets, st)
    else:
        for k in st:
            run(a.tag, k)
