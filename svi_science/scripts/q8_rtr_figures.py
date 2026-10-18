"""Q8 RTR figures: example space-time plots with RTRs marked + summary panels."""
import json
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd
import q8_rtr_lib as L

OUT = f"{L.OUT}/tables"
FIG = f"{L.OUT}/figures"
C_DATA, C_C1, C_C2, C_ROY, C_FWD, C_GREY = "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#eda100", "#b5b4ae"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def spacetime(ax, ev, E, ep, S, title, royer=None, t_lim=None, sst=None, pt=True):
    r = E.set_index("ep").loc[ep]
    g = ev[ev.ep == ep]
    tz = pd.Timedelta(hours=-7) if pt else pd.Timedelta(0)
    tt = g.time + tz
    ax.scatter(tt, g.x, s=1.5 if len(g) > 5000 else 5, c=C_GREY, lw=0, rasterized=True)
    q = S[(S.ep == ep)]
    for s in q.itertuples():
        col = C_DATA if s.cls == "RTR" else C_FWD
        if s.cls not in ("RTR", "forward"):
            continue
        ax.plot([s.start + tz, s.end + tz], [s.x_start, s.x_end], color=col, lw=2.5, solid_capstyle="round")
    if royer is not None:
        for rr in royer.itertuples():
            ax.axvspan(rr.start_utc + tz, rr.end_utc + tz, color=C_ROY, alpha=0.25, lw=0)
    if t_lim is not None:
        ax.set_xlim(t_lim[0] + tz, t_lim[1] + tz)
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_ylabel("along strike (km, + = NW)")
    ax.tick_params(axis="x", labelrotation=30)


def main():
    st = json.load(open(f"{L.OUT}/results_step1_d2.json"))
    SST = st["S_star_frozen"]
    S = pd.read_parquet(f"{OUT}/streaks_all_abs_d2.parquet")
    Ro = pd.read_csv(f"{OUT}/royer2015_rtr_list.csv", parse_dates=["start_utc", "end_utc"])
    evB, EB = qd.bostock(); evS, ES = qd.ours_svi("head"); evN, EN = qd.ours_nw("head")
    sel = lambda c, v: S[(S.cat == c) & (S["var"] == v) & (S["mode"] == "data") & (S.S >= SST) & S.cls.isin(["RTR", "forward"])]
    Bs, Ss, Ns = sel("bostock", "all"), sel("svi", "head"), sel("nw", "head")

    # ---- example episodes
    fig, axs = plt.subplots(3, 2, figsize=(13, 12))
    spacetime(axs[0, 0], evB, EB, "2003-02", Bs, "Bostock 2003 (Mar 1–5 PT): blue = our RTRs, yellow = forward-moving surges, pink = Royer RTRs",
              Ro[Ro.year == 2003], (pd.Timestamp("2003-03-01"), pd.Timestamp("2003-03-06")))
    spacetime(axs[0, 1], evB, EB, "2012-08", Bs, "Bostock 2012 (Sep 8–13 PT)", Ro[Ro.year == 2012],
              (pd.Timestamp("2012-09-08"), pd.Timestamp("2012-09-14")))
    spacetime(axs[1, 0], evB, EB, "2005-09", Bs, "Bostock 2005 whole episode", Ro[Ro.year == 2005])
    spacetime(axs[1, 1], evS, ES, "2024-08", Ss, "Ours SVI 2024-08 (P ≥ 0.75, counted once): 3 RTRs")
    spacetime(axs[2, 0], evS, ES, "2020-11", Ss, "Ours SVI 2020-11: 2 RTRs")
    nw_best = Ns[Ns.cls == "RTR"].ep.value_counts().index[0] if (Ns.cls == "RTR").any() else "2025-09"
    spacetime(axs[2, 1], evN, EN, nw_best, Ns, f"Ours WA-NW {nw_best}")
    for a in axs.flat:
        a.set_xlabel("time (PT)")
    fig.tight_layout(); fig.savefig(f"{FIG}/q8_examples.png", dpi=130); plt.close(fig)

    # per-episode panels for Bostock Royer episodes
    eps = ["2003-02", "2004-07", "2005-09", "2008-05", "2010-08", "2011-08", "2012-08"]
    fig, axs = plt.subplots(len(eps), 1, figsize=(13, 3.0 * len(eps)))
    yr = {"2003-02": 2003, "2004-07": 2004, "2005-09": 2005, "2008-05": 2008, "2010-08": 2010, "2011-08": 2011, "2012-08": 2012}
    for a, ep in zip(axs, eps):
        q = Ro[Ro.year == yr[ep]]
        lim = (q.start_utc.min() - pd.Timedelta(days=2), q.end_utc.max() + pd.Timedelta(days=2))
        spacetime(a, evB, EB, ep, Bs, f"Bostock {ep}: our RTRs (blue), forward surges (yellow), Royer RTRs (pink)", q, lim)
    fig.tight_layout(); fig.savefig(f"{FIG}/q8_bostock_royer_episodes.png", dpi=110); plt.close(fig)

    # ---- summary panels
    CV = pd.read_csv(f"{OUT}/threshold_curves_d2.csv")
    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    for a, (c, v, m, lab) in zip(axs[0], [("bostock", "all", "data", "Bostock 2003–13 (10 episodes)"),
                                         ("svi", "head", "data", "Ours SVI 2014–26 (10)"), ("nw", "head", "data", "Ours WA-NW 2015–26 (14)")]):
        q = CV[(CV.cat == c) & (CV["var"] == v) & (CV["mode"] == m)]
        a.plot(q.S, q.data, color=C_DATA, lw=2, label="data")
        a.plot(q.S, q.c1, color=C_C1, lw=2, label="time-shuffled control")
        a.plot(q.S, q.c2, color=C_C2, lw=2, label="position-permuted control")
        
        a.set_yscale("symlog", linthresh=1); a.set_xlabel("surge score S"); a.set_ylabel("RTRs per run with score ≥ S"); a.set_title(lab, loc="left")
        a.set_xlim(2, 100)
    axs[0, 0].legend(frameon=False, fontsize=8)
    rv = pd.read_csv(f"{OUT}/royer_recall_vs_S_d2.csv")
    a = axs[1, 0]
    a.plot(rv.S, rv.recall, "-o", color=C_DATA, lw=2, ms=4, label="Royer RTRs matched (±6 h)")
    a.plot(rv.S, rv.chance, "-o", color=C_GREY, lw=2, ms=4, label="chance: our list shifted ±1–2 days")
    a.axhline(0.7, color="k", lw=1, ls=":")
    a.set_xlabel("score threshold"); a.set_ylabel("fraction of Royer's 64"); a.set_ylim(0, 1); a.legend(frameon=False, fontsize=8)
    a.set_title("Step 1: match to Royer 2015", loc="left")
    R = pd.read_csv(f"{OUT}/rtr_list_frozen_d2.csv")
    cols = {"bostock": C_DATA, "svi": C_C1, "nw": C_C2}
    names = {"bostock": "Bostock 03–13", "svi": "ours SVI", "nw": "ours WA-NW"}
    a = axs[1, 1]
    for k in cols:
        q = R[R.cat == k]
        if len(q):
            a.hist(q.v_ev.abs(), bins=np.arange(0, 21, 1), histtype="step", lw=2, color=cols[k], density=True, label=f"{names[k]} (n={len(q)})")
    rq = pd.read_csv(f"{OUT}/royer2015_rtr_list.csv")
    a.hist(rq.v_kmh, bins=np.arange(0, 21, 1), histtype="step", lw=2, ls="--", color=C_ROY, density=True, label="Royer list (n=64)")
    a.set_xlabel("RTR speed (km/h)"); a.set_ylabel("density"); a.legend(frameon=False, fontsize=8); a.set_title("Speeds", loc="left")
    a = axs[1, 2]
    for k in cols:
        q = R[R.cat == k]
        if len(q):
            a.hist(q.days_after_onset, bins=np.arange(-3, 35, 2), histtype="step", lw=2, color=cols[k], density=True, label=names[k])
    a.set_xlabel("days after episode onset"); a.set_ylabel("density"); a.legend(frameon=False, fontsize=8); a.set_title("Timing", loc="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/q8_summary.png", dpi=130); plt.close(fig)

    # tides panel
    try:
        Q = json.load(open(f"{L.OUT}/results_questions_d2.json"))
        d = Q["d_tides"]
        labs = [("royer_list_2003_2012", "Royer's list"), ("ours_on_bostock_2003_2013", "our RTRs 03–13"),
                ("ours_svi_2014_2026", "our RTRs 14–26"), ("bostock_forward_streaks", "forward surges 03–13")]
        labs = [z for z in labs if z[0] in d and d[z[0]].get("n", 0) > 0]
        fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
        xs = np.arange(len(labs))
        for a, key, ek, ttl in [(axs[0], "frac_all_positive", "expected_all_positive", "whole RTR during positive up-dip shear"),
                                (axs[1], "frac_mean_positive", "expected_mean_positive", "mean up-dip shear during RTR > 0")]:
            a.bar(xs - 0.18, [d[k][key] for k, _ in labs], 0.34, color=C_DATA, label="observed")
            a.bar(xs + 0.18, [d[k][ek] for k, _ in labs], 0.34, color=C_GREY, label="expected (random start, same length)")
            for i, (k, _) in enumerate(labs):
                a.text(i - 0.18, d[k][key] + 0.02, f"n={d[k]['n']}", ha="center", fontsize=7)
            a.set_xticks(xs); a.set_xticklabels([l for _, l in labs], fontsize=8); a.set_ylim(0, 1.12); a.set_title(ttl, loc="left")
        axs[0].axhline(0.72, color=C_ROY, ls="--", lw=1); axs[0].text(len(labs) - 0.5, 0.74, "Royer 72 %", fontsize=7, ha="right")
        axs[0].legend(frameon=False, fontsize=8)
        fig.tight_layout(); fig.savefig(f"{FIG}/q8_tides.png", dpi=130); plt.close(fig)
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
