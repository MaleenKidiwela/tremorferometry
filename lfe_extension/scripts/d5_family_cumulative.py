#!/usr/bin/env python
"""Track D v1.1 regeneration of Track E's per-family cumulative figure (coordinator request
2026-09-17).  Verbatim copy of `scripts/e_family_cumulative.py` — same construction, same palette —
with FOUR things parameterised so it can be pointed at catalog v1.1 without touching Track E's outputs:
`--catalog` (the v1.1 catalog), `--realness` (Track E's re-scored p_real; optional — without it the
high-purity curve is omitted and the figure says so), `--family-table`, `--outdir` / `--grid` / `--daily`.

Original docstring — per-family cumulative LFE count, 2011-01-01 .. 2026-09-15.

One panel per Bostock family (130, incl. the 14 template-less as empty panels):
  * our cat_core cumulative count 2014->2026                      solid,  blue 450   #2a78d6
  * our cat_ext  cumulative count                                 solid,  blue 300   #6da7ec
  * the high-purity subset p_real >= 0.8                          solid,  orange     #eb6834
  * cat_core rescaled by the per-(family, epoch) recall           dotted, blue 600   #184f95
    (MODEL-BASED: count / recall from results/trackD/completeness_family_epoch.csv, drawn
     continuing Bostock's 2013 total so the 2013->2014 join is rate-comparable)
  * Bostock's own 2011-2013 events, at HIS completeness           dashed, grey       #52514e
SSE windows shaded; the network-epoch boundaries 2017-03 / 2019-01 / 2021-04 as thin verticals.

Colours: one blue ORDINAL ramp for the three readings of our own detector (same measure, different
threshold / correction), one categorical hue for the high-purity subset, neutral grey for the
external catalog.  Validated with the dataviz palette validator (categorical pair PASS, ordinal
ramp PASS, light mode).

  results/trackE/family_cumulative_daily.parquet
  figures/family_cumulative/<family_id>.png          (130)
  figures/e_family_cumulative_grid.png               (12 x 11)
"""
from __future__ import annotations
import os, json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
RES = f"{EXT}/results"
OUT = f"{RES}/trackE"
FIG = f"{EXT}/figures"
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
BOST = f"{REP}/results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"

C_EXT, C_CORE, C_CORR, C_HP, C_BOST = "#6da7ec", "#2a78d6", "#184f95", "#eb6834", "#52514e"
C_SSE, C_EPOCH = "#f0efec", "#c3c2b7"
T0, T1 = "2011-01-01", "2026-09-15"
RECALL_FLOOR = 0.05   # never divide by a smaller recall (see the rescaling note)
EPOCH_LINES = ["2017-03-01", "2019-01-01", "2021-04-01"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored", default="")
    ap.add_argument("--catalog", default=f"{RES}/trackD/lfe_extension_catalog_v1_1.parquet")
    ap.add_argument("--realness", default=f"{OUT}/realness_v1_1.parquet")
    ap.add_argument("--family-table", default=f"{RES}/trackD/family_table_v1_1.csv")
    ap.add_argument("--outdir", default=f"{FIG}/family_cumulative_v1_1")
    ap.add_argument("--grid", default=f"{FIG}/d5_family_cumulative_grid_v1_1.png")
    ap.add_argument("--daily", default=f"{RES}/trackD/family_cumulative_daily_v1_1.parquet")
    ap.add_argument("--meta", default=f"{RES}/trackD/family_cumulative_meta_v1_1.json")
    ap.add_argument("--label", default="v1.1")
    ap.add_argument("--hp", type=float, default=0.8)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    if a.scored:
        S = pd.read_parquet(a.scored, columns=["family_id", "date", "epoch_id", "cat_core", "cat_ext",
                                               "p_real"])
        HAVE_HP = True
    else:
        S = pd.read_parquet(a.catalog, columns=["family_id", "t", "date", "epoch_id",
                                                "cat_core", "cat_ext"])
        HAVE_HP = os.path.exists(a.realness)
        if HAVE_HP:
            # join Track E's re-scored p_real on the event key the catalog and the scorer share
            R = pd.read_parquet(a.realness, columns=["family_id", "t", "p_real"])
            n0 = len(S)
            S = S.merge(R, on=["family_id", "t"], how="left")
            n_hit = int(S.p_real.notna().sum())
            print(f"[hp] realness joined: {n_hit}/{n0} events matched "
                  f"({100 * n_hit / max(n0, 1):.1f} %) from {a.realness}", flush=True)
            if n_hit < 0.5 * n0:
                print("[hp] WARNING: fewer than half the events carry a score — "
                      "the high-purity curve is dropped", flush=True)
                HAVE_HP = False
        else:
            print(f"[hp] {a.realness} not present -> the high-purity curve is OMITTED", flush=True)
        if not HAVE_HP:
            S["p_real"] = np.nan
    S["date"] = pd.to_datetime(S.date)
    FT = pd.read_csv(a.family_table)
    CM = pd.read_csv(f"{RES}/trackD/completeness_family_epoch.csv")
    CM = CM[CM.tier == "T_core"][["epoch_id", "family_id", "recall"]]
    rec = {(int(r.family_id), str(r.epoch_id)): float(r.recall) for r in CM.itertuples()}
    rec_epoch_med = CM.groupby("epoch_id").recall.median().to_dict()
    SSE = pd.read_csv(f"{RES}/trackD/sse_windows.csv")

    # ---- Bostock 2011-2013, mapped to families by his own family hypocentre (exact, 130/130) ----
    B = pd.read_csv(BOST)
    B["time"] = pd.to_datetime(dict(year=B.year, month=B.month, day=B.day, hour=B.hour,
                                    minute=B["min"], second=B.sec.astype(int)), errors="coerce")
    key = {(round(r.lat_bostock, 4), round(r.lon_bostock, 4), round(r.depth_bostock, 4)):
           int(r.family_id) for r in FT.itertuples()}
    B["family_id"] = [key.get((round(la, 4), round(lo, 4), round(dp, 4)), -1)
                      for la, lo, dp in zip(B.lat, B.lon, B.dep)]
    nB = int((B.family_id < 0).sum())
    B = B[(B.family_id >= 0) & B.time.between(T0, "2013-12-31 23:59:59")]
    B["date"] = B.time.dt.normalize()

    days = pd.date_range(T0, T1, freq="D")
    fams = sorted(FT.family_id.astype(int))
    di = {d: i for i, d in enumerate(days)}

    def daily(D, mask=None):
        d = D[mask] if mask is not None else D
        g = d.groupby(["family_id", "date"]).size()
        M = np.zeros((len(fams), len(days)))
        if len(g):
            fi = {f: i for i, f in enumerate(fams)}
            for (f, dt), n in g.items():
                if f in fi and dt in di:
                    M[fi[f], di[dt]] = n
        return M

    D_core = daily(S, S.cat_core.to_numpy())
    D_ext = daily(S, S.cat_ext.to_numpy())
    D_hp = (daily(S, (S.cat_core & (S.p_real >= a.hp)).to_numpy()) if HAVE_HP
            else np.zeros_like(D_core))
    D_bos = daily(B)
    # completeness-corrected core: count / recall(family, epoch of that day)
    ep_of_day = np.empty(len(days), object)
    for i, d in enumerate(days):
        ym = d.strftime("%Y-%m")
        ep_of_day[i] = ("E0" if ym <= "2017-02" else "E1" if ym <= "2018-12"
                        else "E2" if ym <= "2021-03" else "E3")
    D_corr = np.zeros_like(D_core)
    fam_med = CM.groupby("family_id").recall.median().to_dict()
    for i, f in enumerate(fams):
        r = np.array([rec.get((f, ep_of_day[j]), np.nan) for j in range(len(days))])
        bad = ~np.isfinite(r) | (r < 0.01)
        if bad.any():
            fb = fam_med.get(f, np.nan)
            if not np.isfinite(fb) or fb < 0.01:
                fb = np.nan
            r = np.where(bad, fb if np.isfinite(fb) else
                         np.array([rec_epoch_med.get(ep_of_day[j], np.nan) for j in range(len(days))]), r)
        # the rescaling is floored at recall 0.05: below that the division is a 20x+ extrapolation
        # of a model, not a measurement (a few 3-5-template families sit there).
        r = np.clip(r, RECALL_FLOOR, None)
        with np.errstate(invalid="ignore", divide="ignore"):
            D_corr[i] = np.where(np.isfinite(r), D_core[i] / r, np.nan)
    D_corr = np.nan_to_num(D_corr, nan=0.0)

    K = dict(cum_core=np.cumsum(D_core, 1), cum_ext=np.cumsum(D_ext, 1),
             cum_hp=np.cumsum(D_hp, 1), cum_core_corrected=np.cumsum(D_corr, 1),
             cum_bostock=np.cumsum(D_bos, 1))
    # Bostock's catalog stops at 2013-12-31 -> the line stops there
    bos_end = di[pd.Timestamp("2013-12-31")]
    bos_tot = K["cum_bostock"][:, bos_end].copy()
    K["cum_bostock"][:, bos_end + 1:] = np.nan
    # our curves start at 2014-01-01
    ours_start = di[pd.Timestamp("2014-01-01")]
    for k in ("cum_core", "cum_ext", "cum_hp", "cum_core_corrected"):
        K[k][:, :ours_start] = np.nan
    K["cum_core_corrected_joined"] = K["cum_core_corrected"] + bos_tot[:, None]

    rows = []
    for i, f in enumerate(fams):
        rows.append(pd.DataFrame(dict(family=f, date=days,
                                      **{k: v[i] for k, v in K.items()})))
    T = pd.concat(rows, ignore_index=True)
    T.to_parquet(a.daily, index=False)

    FTi = FT.set_index("family_id")
    tl = set(FT[FT.template_less.astype(str).str.lower() == "true"].family_id.astype(int))

    def panel(ax, i, f, small):
        for r in SSE.itertuples():
            ax.axvspan(pd.Timestamp(r.t0), pd.Timestamp(r.t1), color=C_SSE, lw=0, zorder=0)
        for e in EPOCH_LINES:
            ax.axvline(pd.Timestamp(e), color=C_EPOCH, lw=0.7, zorder=1)
        if f in tl:
            ax.text(0.5, 0.5, "template-less\n(no detector)", ha="center", va="center",
                    fontsize=6 if small else 11, color="#8a8a85", transform=ax.transAxes)
        ax.plot(days, K["cum_bostock"][i], color=C_BOST, lw=1.2, ls=(0, (5, 2)), zorder=3)
        if K["cum_ext"][i][-1] > 0:            # a family with no detector has no curve of ours
            ax.plot(days, K["cum_ext"][i], color=C_EXT, lw=1.2, zorder=4)
            ax.plot(days, K["cum_core"][i], color=C_CORE, lw=1.4, zorder=5)
            if K["cum_core"][i][-1] > 0:
                ax.plot(days, K["cum_core_corrected_joined"][i], color=C_CORR, lw=1.1,
                        ls=(0, (1, 1.6)), zorder=6)
            if HAVE_HP:
                ax.plot(days, K["cum_hp"][i], color=C_HP, lw=1.4, zorder=7)
        ax.set_xlim(pd.Timestamp(T0), pd.Timestamp(T1))
        ax.margins(y=0.05)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=5 if small else 9, length=2, colors="#52514e")
        ax.grid(axis="y", color="#e8e7e3", lw=0.5, zorder=0)
        ax.set_axisbelow(True)

    LEG = [Line2D([], [], color=C_CORE, lw=1.6, label="ours, cat_core"),
           Line2D([], [], color=C_EXT, lw=1.6, label="ours, cat_ext")] + \
          ([Line2D([], [], color=C_HP, lw=1.6, label=f"high purity (p_real >= {a.hp:g})")]
           if HAVE_HP else []) + \
          [
           Line2D([], [], color=C_CORR, lw=1.4, ls=(0, (1, 1.6)),
                  label="cat_core / recall (model-based rescaling, joined to Bostock 2013)"),
           Line2D([], [], color=C_BOST, lw=1.4, ls=(0, (5, 2)),
                  label="Bostock (2015) 2011-2013, at HIS completeness")]

    for i, f in enumerate(fams):
        fig, ax = plt.subplots(figsize=(8.2, 4.2))
        panel(ax, i, f, small=False)
        m = FTi.loc[f]
        rr = " ".join(f"{e}:{rec.get((f, e), float('nan')):.2f}" for e in ("E0", "E1", "E2", "E3"))
        ax.set_title(f"family {f} — {m.lat_bostock:.3f}, {m.lon_bostock:.3f}, "
                     f"{m.depth_bostock:.1f} km — {int(m.n_template_stations)} template stations\n"
                     f"per-epoch T_core recall   {rr}", fontsize=10, color="#0b0b0b", loc="left")
        ax.set_ylabel("cumulative LFE count", fontsize=9, color="#52514e")
        ax.legend(handles=LEG, fontsize=6.5, frameon=False, loc="upper left")
        fig.tight_layout()
        fig.savefig(f"{a.outdir}/{f}.png", dpi=130)
        plt.close(fig)

    nr, nc = 12, 11
    fig, axs = plt.subplots(nr, nc, figsize=(24, 26), sharex=True)
    for j, ax in enumerate(axs.ravel()):
        if j >= len(fams):
            ax.axis("off"); continue
        f = fams[j]
        panel(ax, j, f, small=True)
        m = FTi.loc[f]
        ax.set_title(f"F{f}  {m.depth_bostock:.0f} km  n{int(m.n_template_stations)}",
                     fontsize=6, pad=2, color="#0b0b0b")
        if j % nc:
            pass
    fig.suptitle(f"Per-family cumulative LFE count, 2011-2026 — lfe_extension catalog {a.label}"
                 + (" + Track E realness score" if HAVE_HP else " (no realness score available)")
                 + "; Bostock (2015) 2011-2013 for reference (his completeness)",
                 fontsize=13, y=0.995)
    fig.legend(handles=LEG, fontsize=10, frameon=False, ncol=3, loc="upper center",
               bbox_to_anchor=(0.5, 0.982))
    fig.tight_layout(rect=(0, 0, 1, 0.972))
    fig.savefig(a.grid, dpi=110)
    plt.close(fig)

    meta = dict(label=a.label, catalog=a.catalog, realness=(a.realness if HAVE_HP else None),
                hp_curve_present=bool(HAVE_HP),
                n_families=len(fams), template_less=len(tl), bostock_unmapped=nB,
                bostock_events_2011_2013=int(len(B)),
                totals={k: float(np.nansum(v[:, -1])) for k, v in K.items()},
                bostock_total_2013=float(bos_tot.sum()))
    json.dump(meta, open(a.meta, "w"), indent=1)
    print("[cumulative]", json.dumps(meta, indent=1), flush=True)


if __name__ == "__main__":
    main()
