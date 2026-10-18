"""Q8c figures: example episodes (all candidates; P-weighted draw 0; scrambled), tides summary."""
import json, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd, q8_rtr_lib as L
OUT = f"{L.OUT}/tables"; FIG = f"{L.OUT}/figures"
C_DATA, C_FWD, C_GREY, C_SCR, C_B = "#2a78d6", "#eda100", "#b5b4ae", "#eb6834", "#4a3aa7"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
S = pd.read_parquet(f"{OUT}/streaks_all_abs_q8c_d2.parquet")
B = pd.read_parquet(f"{OUT}/streaks_all_abs_d2.parquet"); B = B[(B.cat == "bostock") & (B["mode"] == "data") & (B.cls == "RTR")]
def pt(ts):
    return pd.DatetimeIndex(ts).tz_localize("UTC").tz_convert("America/Los_Angeles").tz_localize(None)
ev, E = qd.allcand_svi("raw", "data"); evb, Eb = qd.back_2010("raw", "data")
def panel(ax, evx, cat, var_pw, ep, title, lim=None, bost=False):
    g = evx[evx.ep == ep]
    ax.scatter(pt(g.time), g.x, s=1, c=C_GREY, lw=0, alpha=0.5, rasterized=True)
    for mode, var, col, lw in [("pw", var_pw, C_DATA, 2.5), ("data", var_pw.replace("raw", "raw"), None, 0)]:
        pass
    q = S[(S.cat == cat) & (S["var"] == var_pw) & (S["mode"] == ("pw" if cat == "ac_svi" else "data")) & (S.ep == ep) & S.cls.isin(["RTR", "forward"])]
    if cat == "ac_svi":
        q = q[q.rep == 0]
    for s in q.itertuples():
        ax.plot(pt([s.start, s.end]), [s.x_start, s.x_end], color=C_DATA if s.cls == "RTR" else C_FWD, lw=2.5, solid_capstyle="round")
    qs = S[(S.cat == cat) & (S["var"] == var_pw.replace("data", "c1")) & (S["mode"] == ("pw" if cat == "ac_svi" else "data")) & (S.ep == ep) & (S.cls == "RTR")]
    if cat == "ac_svi":
        qs = qs[qs.rep == 0]
    for s in qs.itertuples():
        ax.plot(pt([s.start, s.end]), [s.x_start, s.x_end], color=C_SCR, lw=2, ls=(0, (2, 1)))
    if bost:
        for b in B[B.ep == ep].itertuples():
            ax.axvspan(*pt([b.start, b.end]), color=C_B, alpha=0.2, lw=0)
    if lim:
        ax.set_xlim(*pt([pd.Timestamp(lim[0]), pd.Timestamp(lim[1])]))
    ax.set_ylim(-60, 50); ax.set_title(title, loc="left", fontsize=9); ax.set_ylabel("along strike (km, + = NW)"); ax.set_xlabel("time (PT)")
    ax.tick_params(axis="x", labelrotation=30)
fig, axs = plt.subplots(3, 2, figsize=(13, 12))
panel(axs[0, 0], ev, "ac_svi", "raw|data", "2018-06", "SVI 2018-06, all candidates; P-weighted draw: blue RTR, yellow forward;\nred dashed = RTR in scrambled-control candidates")
panel(axs[0, 1], ev, "ac_svi", "raw|data", "2023-04", "SVI 2023-04 (same marks)")
panel(axs[1, 0], ev, "ac_svi", "raw|data", "2024-08", "SVI 2024-08")
panel(axs[1, 1], ev, "ac_svi", "raw|data", "2015-12", "SVI 2015-12")
panel(axs[2, 0], evb, "ac_back", "once|data", "2012-08", "2012 forward-detector candidates (counted once): blue RTR;\nviolet = Bostock-LFE RTRs", ("2012-09-06", "2012-09-21"), bost=True)
panel(axs[2, 1], evb, "ac_back", "once|data", "2010-08", "2010 forward-detector candidates (counted once)", ("2010-08-20", "2010-09-04"), bost=True)
fig.tight_layout(); fig.savefig(f"{FIG}/q8c_examples.png", dpi=120); plt.close(fig)
Q = json.load(open(f"{L.OUT}/results_q8c.json"))["questions"]; P = json.load(open(f"{L.OUT}/results_q8c_pw.json"))["tides_svi_pw"]
rows = [("SVI 14–26\nP-weighted", P["data_RTR_RTR"]["frac_all"], P["data_RTR_RTR"]["exp_all"], P["scrambled_RTR_RTR"]["frac_all"], P["scrambled_RTR_RTR"]["exp_all"]),
        ("SVI 14–26\nraw", Q["ac_svi|raw"]["tides"]["RTR"]["frac_all"], Q["ac_svi|raw"]["tides"]["RTR"]["exp_all"], Q["ac_svi|raw"]["tides"]["scrambled_RTR"]["frac_all"], Q["ac_svi|raw"]["tides"]["scrambled_RTR"]["exp_all"]),
        ("SVI 14–26\ncounted once", Q["ac_svi|once"]["tides"]["RTR"]["frac_all"], Q["ac_svi|once"]["tides"]["RTR"]["exp_all"], Q["ac_svi|once"]["tides"]["scrambled_RTR"]["frac_all"], Q["ac_svi|once"]["tides"]["scrambled_RTR"]["exp_all"]),
        ("2010–13 fwd-det.\nraw", Q["ac_back|raw"]["tides"]["RTR"]["frac_all"], Q["ac_back|raw"]["tides"]["RTR"]["exp_all"], Q["ac_back|raw"]["tides"]["scrambled_RTR"]["frac_all"], Q["ac_back|raw"]["tides"]["scrambled_RTR"]["exp_all"]),
        ("2010–13 fwd-det.\ncounted once", Q["ac_back|once"]["tides"]["RTR"]["frac_all"], Q["ac_back|once"]["tides"]["RTR"]["exp_all"], Q["ac_back|once"]["tides"]["scrambled_RTR"]["frac_all"], Q["ac_back|once"]["tides"]["scrambled_RTR"]["exp_all"])]
fig, ax = plt.subplots(figsize=(11, 4)); xs = np.arange(len(rows)); w = 0.2
ax.bar(xs - 1.5 * w, [r[1] for r in rows], w, color=C_DATA, label="data RTRs, observed")
ax.bar(xs - 0.5 * w, [r[2] for r in rows], w, color="#9ec5f4", label="data RTRs, expected")
ax.bar(xs + 0.5 * w, [r[3] for r in rows], w, color=C_SCR, label="scrambled-control RTRs, observed")
ax.bar(xs + 1.5 * w, [r[4] for r in rows], w, color="#f3b49a", label="scrambled-control RTRs, expected")
ax.set_xticks(xs); ax.set_xticklabels([r[0] for r in rows], fontsize=8); ax.set_ylim(0, 1.05)
ax.set_ylabel("fraction wholly in positive up-dip shear"); ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper right")
ax.set_title("All-candidate RTRs: tidal timing, data vs scrambled-control candidates", loc="left")
fig.tight_layout(); fig.savefig(f"{FIG}/q8c_tides.png", dpi=130); plt.close(fig)
