#!/usr/bin/env python
"""2-8 Hz eye-review sheets from an ab_restack28 stacks npz, in the layout of
results/ab_royer/NW/hypo_qc28/review/sheet_*.png (2026-09-26).

6 families per sheet, sorted by number of detections (descending); per family two panels (N, Z), one
trace per pool station, each trace scaled to its own max; stations sorted by distance from the
reference station = the station with the highest stack SNR (envelope peak / median, ab_grow.envelope),
nearest at the bottom; time axis -12..32 s (the LONG_PRE / LONG_LEN window).  Writes refs.csv.

  python scripts/ab_review_sheets.py --stacks .../stacks_loc.npz --det .../det_it14.parquet --out .../review
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import ab_common as A  # noqa: E402
import ab_grow as G  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stacks", required=True)
    ap.add_argument("--det", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    z = np.load(a.stacks)
    ST = {}
    for k in z.files:
        g, s = k.split("_", 1)
        ST.setdefault(int(g), {})[s] = z[k]
    n = pd.read_parquet(a.det).groupby("grp").size()
    fams = sorted(ST, key=lambda g: -int(n.get(g, 0)))
    refs = []
    for sh in range(0, len(fams), a.per):
        chunk = fams[sh:sh + a.per]
        fig, axs = plt.subplots(len(chunk), 2, figsize=(15, 3.6 * len(chunk)), squeeze=False)
        for r, g in enumerate(chunk):
            T = {s: x for s, x in ST[g].items() if np.isfinite(x).all() and np.abs(x).max() > 0}
            snr = {}
            for s, x in T.items():
                env = G.envelope(x)
                med = float(np.median(env))
                snr[s] = float(env.max() / med) if med > 0 else 0.0
            ref = max(snr, key=snr.get)
            refs.append(dict(grp=g, sheet=sh // a.per, ref=ref))
            la, lo = A.COORD[ref][:2]
            order = sorted(T, key=lambda s: -A.epi_km(la, lo, *A.COORD[s][:2]))   # nearest at bottom
            for c, (ci, cn) in enumerate([(1, "N"), (0, "Z")]):
                ax = axs[r, c]
                for i, s in enumerate(order):
                    y = T[s][ci]
                    tt = np.arange(y.size) / A.FS - G.LONG_PRE
                    mx = np.abs(y).max()
                    if mx > 0:
                        ax.plot(tt, y / mx * 0.45 + (len(order) - 1 - i), "k", lw=0.5)
                ax.set_yticks(range(len(order)))
                ax.set_yticklabels(order[::-1], fontsize=5)
                ax.set_xlim(-G.LONG_PRE, G.LONG_LEN - G.LONG_PRE)
                ax.set_ylim(-1, len(order))
                ax.grid(axis="x", alpha=0.3)
                ax.set_title(f"fam {g} ({int(n.get(g, 0))} det) {cn} 2-8 Hz, ref {ref}", fontsize=7)
        plt.tight_layout()
        plt.savefig(f"{a.out}/sheet_{sh // a.per:02d}.png", dpi=70)
        plt.close()
    pd.DataFrame(refs).to_csv(f"{a.out}/refs.csv", index=False)
    A.log(f"[sheets] {len(fams)} families -> {a.out}")


if __name__ == "__main__":
    main()
