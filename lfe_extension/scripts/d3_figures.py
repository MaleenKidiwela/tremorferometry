#!/usr/bin/env python
"""Track D, D3 -- figures: one space-time plot per SSE (along-strike position vs time, with the
leading-edge track, the first-burst main front and the secondary-front windows) and one summary figure
across every SSE."""
from __future__ import annotations
import argparse, glob, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TD = f"{EXT}/results/trackD"
TH = f"{TD}/thomas"
FIG = f"{EXT}/figures"


def per_sse(tag, title, out):
    fcsv = f"{TH}/e5_thomas_families_{tag}.csv"
    ccsv = f"{TH}/cat_{tag}.csv"
    if not (os.path.exists(fcsv) and os.path.exists(ccsv)):
        return False
    F = pd.read_csv(fcsv); C = pd.read_csv(ccsv)
    xmap = dict(zip(F.family_id.astype(int), F.x_strike))
    C["x"] = C.family_id.astype(int).map(xmap)
    C["time"] = pd.to_datetime(C.time)
    C = C[C.x.notna()]
    if len(C) < 20:
        return False
    t0 = C.time.min()
    C["day"] = (C.time - t0).dt.total_seconds() / 86400.0
    fig, ax = plt.subplots(figsize=(11, 5.2))
    rng = np.random.default_rng(3)
    ax.scatter(C.day, C.x + rng.normal(0, 0.7, len(C)), s=5, c="0.35", alpha=0.45, lw=0,
               label=f"events (n={len(C)})")
    p = f"{TH}/e5_thomas_leadtrack_{tag}_mcall.csv"
    if os.path.exists(p):
        L = pd.read_csv(p)
        if len(L):
            d = (pd.to_datetime(L.t_mid, unit="s") - t0).dt.total_seconds() / 86400.0
            ax.plot(d, L.x_lead, "-", color="tab:blue", lw=2, label="6-h leading edge (10th pct)")
            ax.plot(d, L.x_mean, "--", color="tab:cyan", lw=1.2, label="6-h centroid")
    p = f"{TH}/e5_thomas_firstburst_{tag}_mcall.csv"
    if os.path.exists(p):
        B = pd.read_csv(p)
        if len(B) and "t_fb" in B:
            d = (pd.to_datetime(B.t_fb, unit="s") - t0).dt.total_seconds() / 86400.0
            ax.scatter(d, B.x, s=42, facecolor="none", edgecolor="tab:red", lw=1.3,
                       label=f"family first burst (n={len(B)})")
    p = f"{TH}/e5_thomas_secondary_{tag}_mcall.csv"
    if os.path.exists(p):
        S = pd.read_csv(p)
        # the Theil-Sen line of each accepted secondary-front window, anchored on its own mean position
        for i, r in enumerate(S.itertuples()):
            da = (pd.Timestamp(r.t_start, unit="s") - t0).total_seconds() / 86400.0
            db = (pd.Timestamp(r.t_end, unit="s") - t0).total_seconds() / 86400.0
            dm = (pd.Timestamp(r.t_mid, unit="s") - t0).total_seconds() / 86400.0
            ax.plot([da, db], [r.x_mean + r.speed_kmh * 24 * (da - dm),
                               r.x_mean + r.speed_kmh * 24 * (db - dm)],
                    color="tab:orange", lw=2.4, alpha=0.9,
                    label="secondary front (Theil-Sen)" if i == 0 else None)
    ax.set_xlabel(f"days from {t0:%Y-%m-%d %H:%M} UTC")
    ax.set_ylabel("along-strike position (km, + = NW)")
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(out, dpi=135); plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="core")
    ap.add_argument("--window-set", default="main")
    a = ap.parse_args()
    os.makedirs(FIG, exist_ok=True)
    S = pd.read_csv(f"{TD}/thomas_per_sse.csv")
    S = S[S.window_set.isin([a.window_set, "referee"])]
    made = []
    for r in S.itertuples():
        if r.tier != f"cat_{a.tier}" and r.tier != "bostock":
            continue
        ttl = (f"{r.sse_id} ({r.t0} .. {r.t1}), {r.tier}"
               + (f" -- first-burst front {r.front_firstburst_kmday:+.2f} km/day"
                  if hasattr(r, "front_firstburst_kmday") and np.isfinite(getattr(r, "front_firstburst_kmday", np.nan)) else ""))
        suf = "" if a.window_set == "main" else f"_{a.window_set}"
        out = f"{FIG}/d3_sse_{r.sse_id}_{a.tier}{suf}.png"
        if per_sse(r.tag, ttl, out):
            made.append(out)
    # ---- summary
    M = S[(S.tier == f"cat_{a.tier}") & S.sse_id.ne("REF2012")
          & (S.window_set == a.window_set)].copy()
    M = M[M.n_events.notna()].sort_values("t0")
    suf2 = "" if a.window_set == "main" else f"_{a.window_set}"
    if len(M):
        fig, ax = plt.subplots(4, 1, figsize=(12, 11), sharex=True)
        x = np.arange(len(M))
        ax[0].bar(x, M.n_events, color="0.45"); ax[0].set_ylabel("events in window")
        ax[0].set_title(f"Thomas harness across every SSE window, cat_{a.tier} "
                        f"({a.window_set} window set)", fontsize=12)
        v = M.front_firstburst_kmday.to_numpy(float)
        lo = M.get("front_firstburst_ci_lo", pd.Series(np.nan, index=M.index)).to_numpy(float)
        hi = M.get("front_firstburst_ci_hi", pd.Series(np.nan, index=M.index)).to_numpy(float)
        ax[1].errorbar(x, v, yerr=[np.clip(v - lo, 0, None), np.clip(hi - v, 0, None)], fmt="o",
                       color="tab:red", ms=5, capsize=3)
        ax[1].axhline(-4.695, color="k", ls="--", lw=1, label="Bostock 2012 (-4.70 km/day)")
        ax[1].axhspan(-4.695 * 1.2, -4.695 * 0.8, color="k", alpha=0.08, label="+-20 %")
        ax[1].set_ylabel("main-front speed (km/day)"); ax[1].legend(fontsize=8)
        ax[2].bar(x - 0.2, M.n_secondary, width=0.4, color="tab:orange", label="observed")
        ax[2].bar(x + 0.2, M.n_secondary_null_circ, width=0.4, color="0.6", label="circular null (mean)")
        if "n_secondary_null_dayshuf" in M:
            ax[2].plot(x, M.n_secondary_null_dayshuf, "kv", ms=5, label="day-shuffle null (mean)")
        ax[2].set_ylabel("secondary fronts"); ax[2].legend(fontsize=8)
        ax[3].bar(x, M.reactivation_interval_d_p50, color="tab:green")
        ax[3].axhline(0.704, color="k", ls="--", lw=1, label="Bostock 2012 median (0.70 d)")
        ax[3].set_ylabel("re-activation interval\nmedian (d)"); ax[3].legend(fontsize=8)
        ax[3].set_xticks(x); ax[3].set_xticklabels(M.sse_id, rotation=70, fontsize=7)
        for b in ax:
            b.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(f"{FIG}/d3_thomas_summary_{a.tier}{suf2}.png", dpi=135); plt.close(fig)
        made.append(f"{FIG}/d3_thomas_summary_{a.tier}{suf2}.png")
    print(f"[write] {len(made)} figures:\n  " + "\n  ".join(made))


if __name__ == "__main__":
    main()
