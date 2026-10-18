#!/usr/bin/env python
"""Gate figures: one PASS, one FAIL, one UNTESTABLE example.

  python scripts/g_figures.py --examples pos_bost_N200:posb200_f1 alias_N50:alias_f10 ...
Each figure: (a) the Gate-1 withheld-station stack with the S and pre-arrival windows, (b) the
Gate-1 null distributions (random-time and moveout-shuffled) with the observed value, (c) the Gate-0
per-site stack cc against the per-station cross-family chance level, (d) the Gate-2 map.
"""
from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402
import g_gates as GG                                             # noqa: E402
import g_locate as L                                             # noqa: E402


def one_figure(ctx, prop, tag, out):
    rng = np.random.default_rng(ctx.cfg.seed)
    r0 = GG.gate0(ctx, prop, rng)
    r1 = GG.gate1(ctx, prop, rng)
    r2 = GG.gate2(ctx, prop, rng)
    best = max(r1["per_station"], key=lambda x: x.get("excess_db", -99) or -99)
    sta = best["sta"]
    _d, ttS, _ = ctx.pred_S(sta, prop.lat, prop.lon, prop.dep)
    st, n = ctx.store.stack(ctx.netsta(sta), prop.times + ttS)
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.4))
    # (a) stack
    t = (np.arange(G.NW) - G.IDX_S) / G.FS
    if st is not None:
        for i, c in enumerate("ZNE"):
            ax[0].plot(t, st[i] / np.abs(st).max() * 0.8 + (1 - i), lw=.7, label=c)
        ax[0].axvspan(-1, 5, color="tab:orange", alpha=.2)
        ax[0].axvspan(-12, -6, color="tab:grey", alpha=.2)
        ax[0].axvline(0, color="k", lw=.8)
    ax[0].set_xlim(-20, 20); ax[0].set_xlabel("t - predicted S (s)")
    ax[0].set_title(f"Gate 1 stack at withheld {sta}\nN={best.get('n_used')} d={best.get('dist_km',0):.0f} km "
                    f"excess {best.get('excess_db', np.nan):+.2f} dB")
    ax[0].legend(fontsize=7, loc="lower right")
    # (b) nulls
    nm, ns = best.get("null_mean_db", 0.0), best.get("null_sd_db", 1.0)
    x = np.linspace(nm - 4 * ns, max(best.get("sdb_db", 0), nm + 4 * ns) + 1, 200)
    ax[1].plot(x, np.exp(-(x - nm) ** 2 / (2 * ns ** 2)), color="tab:grey",
               label=f"random-time null ({best.get('n_null')})")
    ax[1].axvline(best.get("null_p99_db", np.nan), color="tab:grey", ls="--", label="null p99")
    ax[1].axvline(nm + (best.get("moveout_p95_excess_db") or 0), color="tab:red", ls=":",
                  label="moveout-shuffled p95")
    ax[1].axvline(best.get("sdb_db", np.nan), color="tab:blue", lw=2, label="observed")
    ax[1].set_xlabel("S/pre-arrival ratio (dB)"); ax[1].set_yticks([])
    ax[1].legend(fontsize=7); ax[1].set_title(f"Gate 1 nulls -> {r1['verdict']}")
    # (c) gate 0
    P = pd.DataFrame(r0.get("cc_per_station", []))
    if len(P):
        col = ["tab:red" if b else "tab:grey" for b in P.registers]
        ax[2].bar(range(len(P)), P.cc, color=col)
        ax[2].bar([], [], color="tab:red", label="site registers (S excess >= 2 dB)")
        ax[2].bar([], [], color="tab:grey", label="does not register")
        ax[2].plot(range(len(P)), P.cc_p95, "k_", ms=14, label="cross-family p95")
        ax[2].set_xticks(range(len(P))); ax[2].set_xticklabels(P.sta, rotation=90, fontsize=7)
    ax[2].set_ylabel("stack cc to nearest accepted family"); ax[2].legend(fontsize=7)
    ax[2].set_title(f"Gate 0 -> {r0['verdict']}  (d={r0['min_dist_km']:.2f} km to f{r0['nearest_bkey']}, "
                    f"HF {r0.get('hf_median', np.nan):.2f})")
    # (d) gate 2
    sll = ctx.site_ll
    ax[3].scatter([sll[s][1] for s in ctx.cfg.sites], [sll[s][0] for s in ctx.cfg.sites],
                  marker="^", c="g", s=60)
    for s in ctx.cfg.sites:
        ax[3].annotate(s, (sll[s][1], sll[s][0]), fontsize=6)
    ax[3].scatter([prop.lon], [prop.lat], marker="*", s=180, c="k", label="proposal")
    if "interface" in r2:
        ax[3].scatter([r2["interface"]["lon"]], [r2["interface"]["lat"]], s=60, c="tab:blue",
                      label="Gate 2 (interface)")
        ax[3].scatter([r2["free_depth"]["lon"]], [r2["free_depth"]["lat"]], s=40, c="tab:orange",
                      label="free depth")
    ax[3].legend(fontsize=7)
    ax[3].set_title(f"Gate 2 -> {r2['verdict']}\n{r2.get('reason', '')}"[:90])
    ax[3].set_xlabel("lon"); ax[3].set_ylabel("lat")
    fig.suptitle(f"{prop.pid} [{tag}]  N={len(prop.times)}  ->  overall "
                 f"{GG.evaluate(ctx, prop, gates=())['verdict'] if False else ''}"
                 f"{'FAIL' if 'FAIL' in (r0['verdict'], r1['verdict'], r2['verdict']) else ('UNTESTABLE' if 'UNTESTABLE' in (r0['verdict'], r1['verdict'], r2['verdict']) else 'PASS')}",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=115)
    plt.close(fig)
    print("->", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--examples", nargs="+", required=True, help="set:pid[:tag]")
    a = ap.parse_args()
    ctx = GG.GateContext(GG.GateConfig())
    for ex in a.examples:
        parts = ex.split(":")
        st, pid = parts[0], parts[1]
        tag = parts[2] if len(parts) > 2 else st
        props = {p.pid: p for p in GG.load_proposals(f"{G.GATES}/proposals_{st}.parquet")}
        one_figure(ctx, props[pid], tag, f"{G.FIGS}/g_example_{tag}_{pid}.png")


if __name__ == "__main__":
    main()
