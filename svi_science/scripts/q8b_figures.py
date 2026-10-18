"""Q8b figures: PNSN space-time examples with RTRs, tides summary."""
import json, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd, q8_rtr_lib as L
OUT = f"{L.OUT}/tables"; FIG = f"{L.OUT}/figures"
C_DATA, C_FWD, C_GREY, C_ROY, C_B = "#2a78d6", "#eda100", "#b5b4ae", "#e87ba4", "#4a3aa7"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
S = pd.read_parquet(f"{OUT}/streaks_all_abs_pnsn_d2.parquet"); S = S[S["mode"] == "data"]
B = pd.read_parquet(f"{OUT}/streaks_all_abs_d2.parquet"); B = B[(B.cat == "bostock") & (B["mode"] == "data") & (B.cls == "RTR")]
Ro = pd.read_csv(f"{OUT}/royer2015_rtr_list.csv", parse_dates=["start_utc", "end_utc"])
LD = {"pnsn_svi": qd.pnsn("SVI"), "pnsn_nw": qd.pnsn("NW"), "pnsn_sw": qd.pnsn("SW")}
def pt(ts):
    return pd.DatetimeIndex(ts).tz_localize("UTC").tz_convert("America/Los_Angeles").tz_localize(None)
def panel(ax, cat, ep, title, lim=None, bostock=False):
    ev, E = LD[cat]; g = ev[ev.ep == ep]
    ax.scatter(pt(g.time), g.x, s=2, c=C_GREY, lw=0, rasterized=True)
    q = S[(S.cat == cat) & (S.ep == ep) & S.cls.isin(["RTR", "forward"])]
    for s in q.itertuples():
        ax.plot(pt([s.start, s.end]), [s.x_start, s.x_end], color=C_DATA if s.cls == "RTR" else C_FWD, lw=2.5, solid_capstyle="round")
    if bostock:
        for b in B.itertuples():
            ax.axvspan(*pt([b.start, b.end]), color=C_B, alpha=0.18, lw=0)
        for r in Ro.itertuples():
            ax.axvline(pt([r.start_utc])[0], color=C_ROY, lw=1, ls="--")
    if lim:
        ax.set_xlim(*pt([pd.Timestamp(lim[0]), pd.Timestamp(lim[1])]))
    ax.set_title(title, loc="left", fontsize=9); ax.set_ylabel("along strike (km, + = NW)"); ax.set_xlabel("time (PT)")
    ax.tick_params(axis="x", labelrotation=30)
fig, axs = plt.subplots(3, 2, figsize=(13, 12))
panel(axs[0, 0], "pnsn_svi", "2012-08", "PNSN SVI 2012-09: blue RTR, yellow forward;\nviolet = Bostock-LFE RTRs, pink dashes = Royer",
      ("2012-09-08", "2012-09-16"), bostock=True)
panel(axs[0, 1], "pnsn_svi", "2011-08", "PNSN SVI 2011-08 (same marks)", ("2011-08-20", "2011-09-03"), bostock=True)
panel(axs[1, 0], "pnsn_svi", "2023-04", "PNSN SVI 2023-04")
panel(axs[1, 1], "pnsn_svi", "2020-11", "PNSN SVI 2020-11")
panel(axs[2, 0], "pnsn_nw", "2018-05", "PNSN WA-NW 2018-05")
panel(axs[2, 1], "pnsn_sw", "2021-09", "PNSN WA-SW 2021-09 (NW episode window)")
fig.tight_layout(); fig.savefig(f"{FIG}/q8b_pnsn_examples.png", dpi=125); plt.close(fig)

Q = json.load(open(f"{L.OUT}/results_q8b.json")); Q8 = json.load(open(f"{L.OUT}/results_questions_d2.json"))
d = Q["d_tides_svi"]; d8 = Q8["d_tides"]
rows = [("LFE RTRs\n2003–13", d8["ours_on_bostock_2003_2013"]["frac_all_positive"], d8["ours_on_bostock_2003_2013"]["expected_all_positive"], d8["ours_on_bostock_2003_2013"]["n"]),
        ("LFE RTRs\n2014–26", d8["ours_svi_2014_2026"]["frac_all_positive"], d8["ours_svi_2014_2026"]["expected_all_positive"], d8["ours_svi_2014_2026"]["n"]),
        ("tremor RTRs\n2010–13", d["RTR_2010-13"]["frac_all"], d["RTR_2010-13"]["exp_all"], d["RTR_2010-13"]["n"]),
        ("tremor RTRs\n2014–26", d["RTR_2014-26"]["frac_all"], d["RTR_2014-26"]["exp_all"], d["RTR_2014-26"]["n"]),
        ("tremor forward\nsurges", d["forward"]["frac_all"], d["forward"]["exp_all"], d["forward"]["n"]),
        ("tremor non-migr.\nsurges", d["surge_only"]["frac_all"], d["surge_only"]["exp_all"], d["surge_only"]["n"])]
fig, ax = plt.subplots(figsize=(10, 4)); xs = np.arange(len(rows))
ax.bar(xs - 0.18, [r[1] for r in rows], 0.34, color=C_DATA, label="observed")
ax.bar(xs + 0.18, [r[2] for r in rows], 0.34, color=C_GREY, label="expected (random start, same length)")
for i, r in enumerate(rows):
    ax.text(i - 0.18, r[1] + 0.02, f"n={r[3]}", ha="center", fontsize=7)
ax.axhline(0.72, color=C_ROY, ls="--", lw=1); ax.text(len(rows) - 0.5, 0.74, "Royer 72 %", fontsize=7, ha="right")
ax.set_xticks(xs); ax.set_xticklabels([r[0] for r in rows], fontsize=8); ax.set_ylim(0, 1.1)
ax.set_ylabel("fraction wholly in positive up-dip shear"); ax.set_title("SVI: bursts during positive tidal up-dip shear (LFE catalogue vs PNSN tremor)", loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig(f"{FIG}/q8b_tides.png", dpi=130); plt.close(fig)
