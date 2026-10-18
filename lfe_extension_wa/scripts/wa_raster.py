#!/usr/bin/env python
"""Per-family time raster of headline-selected detections (core network, prod run), 2007-2026. 2026-10-06.
Rule: NW/ArmB P_own_WAVE_PNET >= 0.75 (all days); SW P_own_WAVE_PNET >= 0.5 on tremor days; OlyA P_nbr3_WAVE >= 0.75.
Colour = day class (tremor <=30 km / quiet = covered by PNSN catalogue & no tremor <=60 km / unknown = no catalogue / near)."""
import glob, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.dates as mdates
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/scripts")
import wa_prod_common as P, wa_prod_post as W
CUT = float(sys.argv[1]) if len(sys.argv) > 1 else None   # uniform cut-off, all day classes (None = headline rules)
TAG = "" if CUT is None else f"_P{CUT:g}"
parts = []
for p in sorted(glob.glob(f"{P.OUT}/conf/P/prod/P_*.parquet")):
    D = pd.read_parquet(p, columns=["date", "t", "net", "series", "fid", "reg", "role", "P_own_WAVE_PNET", "P_nbr3_WAVE"])
    D = D[(D.series == "data") & (D.net == "core") & (~D.role.astype(str).str.startswith("sealed"))]
    if not len(D): continue
    D = W.day_class(D)
    if CUT is None:
        sel = ((D.reg.isin(["NW", "ArmB"])) & (D.P_own_WAVE_PNET >= 0.75)) | ((D.reg == "SW") & (D.P_own_WAVE_PNET >= 0.5) & D.tremor30.astype(bool)) \
              | ((D.reg == "OlyA") & (D.P_nbr3_WAVE >= 0.75))
    else:
        sel = ((D.reg.isin(["NW", "ArmB", "SW"])) & (D.P_own_WAVE_PNET >= CUT)) | ((D.reg == "OlyA") & (D.P_nbr3_WAVE >= CUT))
    parts.append(D[sel][["date", "t", "fid", "reg", "tremor30", "tremor60"]])
S = pd.concat(parts, ignore_index=True)
cov = (S.date >= "2009-08-06") & (S.date <= "2026-05-27")
S["cls"] = np.where(S.tremor30.astype(bool), "tremor", np.where(S.tremor60.astype(bool), "near", np.where(cov, "quiet", "unknown")))
S["time"] = pd.to_datetime(S.t.astype(float), unit="s")
F = pd.read_csv(f"{P.OUT}/setup/families.csv")
S.to_parquet(f"{P.OUT}/figures/raster_selected{TAG}.parquet", index=False)
col = {"tremor": "#1f77b4", "near": "#9467bd", "quiet": "#d62728", "unknown": "#7f7f7f"}
names = {k: v for k, v in {"NW": "NW (Royer), PNET P ≥ 0.75", "SW": "SW (Royer), PNET P ≥ 0.5, tremor days", "ArmB": "Arm B (Chestler), PNET P ≥ 0.75",
         "OlyA": "Olympic Arm A, nbr3 P ≥ 0.75 (group-level)"}.items()}
if CUT is not None:
    names = {"NW": f"NW (Royer), PNET P ≥ {CUT:g}, all days", "SW": f"SW (Royer), PNET P ≥ {CUT:g}, all days",
             "ArmB": f"Arm B (Chestler), PNET P ≥ {CUT:g}, all days", "OlyA": f"Olympic Arm A, nbr3 P ≥ {CUT:g} (group-level)"}
for reg in ["NW", "SW", "ArmB", "OlyA"]:
    f = F[F.reg == reg].sort_values("lat_use", ascending=False).reset_index(drop=True)
    ypos = {fid: i for i, fid in enumerate(f.fid)}
    x = S[S.reg == reg]
    h = max(4, 0.11 * len(f) + 1.5)
    fig, ax = plt.subplots(figsize=(15, h))
    for c in ["tremor", "near", "unknown", "quiet"]:
        g = x[x.cls == c]
        ax.scatter(g.time, g.fid.map(ypos), s=2, c=col[c], lw=0, alpha=0.6, label=f"{c} ({len(g):,})", rasterized=True)
    ax.axvspan(pd.Timestamp("2007-01-01"), pd.Timestamp("2011-12-31"), color="0.9", zorder=-1)
    ax.text(pd.Timestamp("2007-03-01"), -1.5, "array era", fontsize=8, color="0.4")
    ax.set_ylim(len(f) - 0.5 + 1, -2.5); ax.set_xlim(pd.Timestamp("2007-01-01"), pd.Timestamp("2026-10-01"))
    step = 1 if len(f) <= 45 else 5
    ax.set_yticks(range(0, len(f), step)); ax.set_yticklabels([str(int(v)) for v in f.fid[::step]], fontsize=6)
    ax.set_ylabel("family id (sorted north → south)"); ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    nfam = x.fid.nunique(); n10 = (x[x.date >= "2012-01-01"].groupby("fid").size() >= 10).sum()
    ax.set_title(f"{names[reg]} — {len(x):,} detections, {nfam} of {len(f)} families; {n10} with ≥ 10 in 2012–2026 (core network)", loc="left")
    ax.legend(loc="lower left", markerscale=6, fontsize=8, ncol=4, title="day class")
    fig.tight_layout(); out = f"{P.OUT}/figures/raster_{reg}{TAG}.png"; fig.savefig(out, dpi=120); plt.close(fig); print(out, len(x))
