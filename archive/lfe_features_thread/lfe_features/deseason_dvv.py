"""Per-family deseasoning of the rolling coda dv/v.

For each family's daily dv/v, fit [const + linear trend + annual + semi-annual harmonics] by
cc-weighted least squares and SUBTRACT only the harmonic (seasonal) terms -> residual keeps the
mean/trend/transients but removes the annual cycle. Reveals any coherent non-seasonal signal.

Usage: python lfe_features/deseason_dvv.py --csv data/daily_dvv_B011p70_2to4_roll.csv --sta B011
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
P = 365.25


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True); ap.add_argument("--sta", required=True)
    a = ap.parse_args()
    d = pd.read_csv(a.csv); d["date"] = pd.to_datetime(d.date)
    d["pct"] = d.dvv * (100.0 if d.dvv.abs().median() < 0.1 else 1.0)
    w = d.cc_max.values if "cc_max" in d.columns else np.ones(len(d))
    t0 = d.date.min()
    out = []
    seas_amp = []
    for fam, g in d.groupby("patch"):
        g = g.sort_values("date")
        t = (g.date - t0).dt.days.values.astype(float)
        y = g.pct.values
        wj = (g.cc_max.values if "cc_max" in g else np.ones(len(g)))
        if len(g) < 120 or (t.max() - t.min()) < 365:
            g = g.assign(resid=np.nan, seasonal=np.nan); out.append(g); continue
        X = np.column_stack([np.ones_like(t), t,
                             np.cos(2*np.pi*t/P), np.sin(2*np.pi*t/P),
                             np.cos(4*np.pi*t/P), np.sin(4*np.pi*t/P)])
        W = np.sqrt(wj)
        beta, *_ = np.linalg.lstsq(X * W[:, None], y * W, rcond=None)
        seasonal = X[:, 2:] @ beta[2:]          # harmonic terms only
        resid = y - seasonal
        a1 = np.hypot(beta[2], beta[3]); a2 = np.hypot(beta[4], beta[5])
        seas_amp.append(a1)
        g = g.assign(resid=resid, seasonal=seasonal)
        out.append(g)
    D = pd.concat(out)
    D.to_csv(a.csv.replace(".csv", "_deseason.csv"), index=False)

    # plot residual: faint per-family + cross-family median (>=10 families/day)
    fig, ax = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    # (top) raw with fitted seasonal removed shown as median comparison
    medraw = D.groupby("date").filter(lambda x: len(x) >= 10).groupby("date").pct.median()
    medres = D.dropna(subset=["resid"]).groupby("date").filter(lambda x: len(x) >= 10).groupby("date").resid.median()
    ax[0].plot(medraw.index, medraw.values, color="#888", lw=0.8, label="raw dv/v median")
    ax[0].set_ylabel("dv/v (%)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25); ax[0].axhline(0, color="k", lw=.5, ls="--")
    ax[0].set_title(f"{a.sta} dv/v (2-4s, 52 families) — raw (top) vs deseasoned residual (bottom); "
                    f"median annual amp removed = {np.median(seas_amp):.3f}%", fontsize=9)
    for fam, g in D.dropna(subset=["resid"]).groupby("patch"):
        ax[1].plot(g.date, g.resid, lw=0.4, alpha=0.10, color="#888")
    ax[1].plot(medres.index, medres.values, color="#b30000", lw=1.3, label="cross-family residual median")
    ax[1].axhline(0, color="k", lw=.5, ls="--")
    lo, hi = medres.quantile(.01), medres.quantile(.99); pad = 0.3*(hi-lo) + 0.02
    ax[1].set_ylim(lo-pad, hi+pad)
    ax[1].set_ylabel("residual dv/v (%)"); ax[1].set_xlabel("date"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(f"figures/dvv_{a.sta}p70_2to4_deseason.png", dpi=130)
    print(f"median annual seasonal amplitude removed: {np.median(seas_amp):.3f} %")
    print(f"residual cross-family median: std {medres.std():.3f}%  range {medres.min():.3f}..{medres.max():.3f}%")
    print(f"saved figures/dvv_{a.sta}p70_2to4_deseason.png + {a.csv.replace('.csv','_deseason.csv')}")


if __name__ == "__main__":
    main()
