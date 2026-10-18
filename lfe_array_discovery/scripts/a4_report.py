#!/usr/bin/env python
"""A4 report: real vs null separation per subarray, every criterion a column, no cuts.

Reads <dir>/<SUB>_<cond>/loc_table.csv (cond real / nullT / nullR), writes
<dir>/all_locations.csv, <dir>/separation.csv, <dir>/curves.csv, <dir>/separation.png, <dir>/REPORT.md.

AUC = P(real location scores better than a null location) (ties 1/2); a missing value (no
seed / no detections / not located) counts as the worst value.  Curves: number of locations whose
value is at least as good as x, real and null, on a grid of the pooled values."""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
TRUTH = os.path.join(HERE, "../results/ab_royer/truth/royer2014_templates_digitized.csv")

METRICS = [  # column, direction ('hi' better / 'lo' better), meaning
    ("hyp_dist_km", "lo", "Hyp2000 epicentre of the final stack to the seed location (km)"),
    ("hyp_d_truth_km", "lo", "Hyp2000 epicentre to the NEAREST Royer template (km)"),
    ("frac_S_1s", "hi", "fraction of stack stations whose horizontal envelope peak is within 1 s of predicted S"),
    ("frac_P_1s", "hi", "same for the vertical peak before S vs predicted P"),
    ("dS_med_abs", "lo", "median |observed - predicted S| (s)"),
    ("moveout_corr", "hi", "corr(observed S, predicted S) across stations"),
    ("n_picks_S", "hi", "unambiguous (R>=1.5) S picks on the Hilbert-rotated stack"),
    ("n_picks_P", "hi", "unambiguous (R>=1.5) P picks"),
    ("snr_med", "hi", "median station stack SNR (envelope peak / median, 30 s stack)"),
    ("seed_excess", "hi", "seed coherence: mean top-5 similarity minus median candidate similarity"),
    ("det_per_day_fin", "hi", "final matched-filter detections per day (8 x MAD)"),
    ("nmad_med_fin", "hi", "median nmad of the final detections"),
    ("cc_med_fin", "hi", "median network cc of the final detections"),
    ("day_max_nmad_med_fin", "hi", "median over days of the day's maximum nmad"),
]


def auc(r, n, d):
    r = np.asarray(r, float)
    n = np.asarray(n, float)
    worst = -np.inf if d == "hi" else np.inf
    r = np.where(np.isfinite(r), r, worst)
    n = np.where(np.isfinite(n), n, worst)
    if d == "lo":
        r, n = -r, -n
    if not len(r) or not len(n):
        return np.nan
    gt = (r[:, None] > n[None]).mean()
    eq = (r[:, None] == n[None]).mean()
    return float(gt + 0.5 * eq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--title", default="A4 location-seeded positive control")
    a = ap.parse_args()
    tabs = []
    for p in sorted(glob.glob(f"{a.dir}/*_*/loc_table.csv")):
        tabs.append(pd.read_csv(p))
    if not tabs:
        raise SystemExit("no loc_table.csv")
    T = pd.concat(tabs, ignore_index=True)
    tr = pd.read_csv(TRUTH)
    tr = tr[tr.locator == "Hyp2000"]
    if "hyp_lat" in T:
        d = np.full(len(T), np.nan)
        for i, r in enumerate(T.itertuples()):
            if np.isfinite(r.hyp_lat):
                t = tr[tr.subarray == r.subarray]
                dd = np.hypot(111.19 * (t.lat - r.hyp_lat),
                              111.19 * np.cos(np.radians(r.hyp_lat)) * (t.lon - r.hyp_lon))
                d[i] = dd.min()
        T["hyp_d_truth_km"] = d
    T["grp"] = T.subarray + ("/" + T.seed_mode if "seed_mode" in T else "")
    T.to_csv(f"{a.dir}/all_locations.csv", index=False)
    rows, curves = [], []
    for sub, S in T.groupby("grp"):
        R = S[S.cond == "real"]
        for nc in ("nullT", "nullR"):
            N = S[S.cond == nc]
            if not len(N):
                continue
            for m, d, _ in METRICS:
                if m not in S:
                    continue
                rows.append(dict(subarray=sub, null=nc, metric=m, better=d, n_real=len(R), n_null=len(N),
                                 real_median=float(R[m].median()), null_median=float(N[m].median()),
                                 real_finite=int(R[m].notna().sum()), null_finite=int(N[m].notna().sum()),
                                 auc=auc(R[m], N[m], d)))
        for m, d, _ in METRICS:
            if m not in S:
                continue
            v = S[m].dropna().values
            if not v.size:
                continue
            grid = np.unique(np.quantile(v, np.linspace(0, 1, 21)))
            for x in grid:
                row = dict(subarray=sub, metric=m, better=d, x=float(x))
                for c in ("real", "nullT", "nullR"):
                    vv = S[S.cond == c][m]
                    row[f"n_{c}"] = int(((vv >= x) if d == "hi" else (vv <= x)).sum())
                    row[f"of_{c}"] = int((S.cond == c).sum())
                curves.append(row)
    Sep = pd.DataFrame(rows)
    Cv = pd.DataFrame(curves)
    Sep.to_csv(f"{a.dir}/separation.csv", index=False)
    Cv.to_csv(f"{a.dir}/curves.csv", index=False)

    # figure: cumulative curves for 4 key metrics per subarray
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        keys = ["hyp_dist_km", "frac_S_1s", "snr_med", "seed_excess"]
        subs = sorted(T.grp.unique())
        fig, ax = plt.subplots(len(subs), len(keys), figsize=(4 * len(keys), 3.2 * len(subs)), squeeze=False)
        col = {"real": "#1f5fbf", "nullT": "#c0392b", "nullR": "#e59f3a"}
        for i, sub in enumerate(subs):
            for j, m in enumerate(keys):
                A_ = ax[i, j]
                if m not in T:
                    continue
                d = dict((x[0], x[1]) for x in METRICS)[m]
                for c in ("real", "nullT", "nullR"):
                    v = T[(T.grp == sub) & (T.cond == c)][m].values.astype(float)
                    if not v.size:
                        continue
                    n = v.size
                    fv = np.sort(v[np.isfinite(v)])
                    if d == "hi":
                        fv = fv[::-1]
                    A_.step(fv, np.arange(1, fv.size + 1) / n, where="post", color=col[c],
                            label=f"{c} (n={n}, finite {fv.size})")
                A_.set_title(f"{sub}: {m} ({'lower' if d == 'lo' else 'higher'} = better)", fontsize=9)
                A_.set_ylabel("fraction of locations at least this good", fontsize=8)
                if m == "hyp_dist_km":
                    A_.set_xlim(0, 30)
                A_.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(f"{a.dir}/separation.png", dpi=110)
    except Exception as e:
        print("figure failed:", e)

    def md(df):
        h = "| " + " | ".join(df.columns) + " |\n|" + "---|" * len(df.columns) + "\n"
        return h + "\n".join("| " + " | ".join(f"{v:.3g}" if isinstance(v, float) else str(v) for v in r) + " |"
                             for r in df.itertuples(index=False))

    lines = [f"# {a.title}\n",
             f"Written by `scripts/a4_report.py` from `{a.dir}`.  Method, choices and caveats: "
             "`results/a4_seed/A4_LOG.md` and the docstring of `scripts/a4_seed.py`.\n",
             "Conditions: **real** = Royer & Bostock (2014) Fig 4 Hyp2000 template locations, real data; "
             "**nullT** = same locations, every station reading a different (day, in-day offset) -- the "
             "no-family null; **nullR** = the template cloud rotated + shifted (a location-specificity "
             "control; its points are a median 2-5 km from some real template, so it is NOT family-free).\n",
             "AUC = P(a real location scores better than a null location); 0.5 = no separation.  Missing "
             "values (no seed, no detections, not located) count as the worst value.  Nothing is cut.\n",
             "## Counts per subarray and condition\n"]
    g = T.groupby(["grp", "cond"])
    cnt = g.agg(n=("lid", "size"),
                seeded=("seed_score", lambda x: int(x.notna().sum())) if "seed_score" in T else ("lid", "size"),
                located=("hyp_dist_km", lambda x: int(x.notna().sum())) if "hyp_dist_km" in T else ("lid", "size"),
                within_3km=("hyp_dist_km", lambda x: int((x <= 3).sum())) if "hyp_dist_km" in T else ("lid", "size"),
                within_5km=("hyp_dist_km", lambda x: int((x <= 5).sum())) if "hyp_dist_km" in T else ("lid", "size"),
                frac_S_1s_med=("frac_S_1s", "median"), snr_med=("snr_med", "median"),
                det_per_day_med=("det_per_day_fin", "median")).reset_index()
    lines.append(md(cnt) + "\n")
    lines.append("## Real vs null, criterion by criterion (AUC)\n")
    if len(Sep):
        P = Sep.pivot_table(index=["subarray", "metric", "better"], columns="null",
                            values=["auc", "real_median", "null_median"]).reset_index()
        P.columns = ["_".join([c for c in col if c]) for col in P.columns]
        lines.append(md(P) + "\n")
    lines.append("## Metric definitions\n")
    lines += [f"- `{m}` ({'lower' if d == 'lo' else 'higher'} = better): {t}" for m, d, t in METRICS]
    lines.append("\nFigure: `separation.png` (cumulative fraction of locations at least as good as x). "
                 "Curve table: `curves.csv`.  Per-location table: `all_locations.csv`.\n")
    open(f"{a.dir}/REPORT.md", "w").write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
