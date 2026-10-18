"""C figures: rate vs time since front (per depth class), inter-event times, burst sizes. Reads results/c_decay/."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import special

R = "/home/jovyan/tremorferometry/svi_science/results/c_decay"
FIG = "/home/jovyan/tremorferometry/svi_science/figures"
CL = {"shallow": "#86b6ef", "mid": "#2a78d6", "deep": "#0d366b"}
OBS, SUR, UNI = "#2a78d6", "#8a8984", "#eb6834"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e5e1", "grid.linewidth": 0.6, "axes.edgecolor": "#52514e"})

E = pd.read_parquet(f"{R}/events_lag.parquet")
CT = pd.read_csv(f"{R}/cells_duration.csv")
res = json.load(open(f"{R}/result.json"))
regions = [("SVI", "Southern Vancouver Island (10 big ETS)"), ("NW", "Washington-NW (14 big ETS)")]

# ---------------- figure 1: rate vs time since front
fig, ax = plt.subplots(2, 3, figsize=(13, 7.2))
bins = np.arange(-5, 30.01, 0.5); ctr = bins[:-1] + 0.25
for i, (reg, title) in enumerate(regions):
    eall = E[(E.region == reg) & (E.run == "headline")]
    e = eall.dropna(subset=["lag"])
    neps = e.ep.nunique()
    a = ax[i, 0]
    for cl in ["shallow", "mid", "deep"]:
        x = e[e.dcl == cl].lag
        h, _ = np.histogram(x, bins)
        a.semilogy(ctr, np.where(h > 0, h / 0.5 / neps, np.nan), color=CL[cl], lw=2, label=f"{cl} (n={len(x)})")
    a.axvline(0, color="#52514e", lw=0.8, ls="--")
    a.set_xlabel("days since tremor front reached the patch's along-strike position")
    a.set_ylabel("LFEs per day per episode (stack)")
    a.set_title(f"{title}\nstacked rate by depth class", loc="left")
    a.legend(frameon=False, fontsize=8)
    # single patch-episode shape: daily counts aligned on each cell's peak PT day, / peak count, median over cells
    a = ax[i, 1]
    c = CT[(CT.region == reg) & (CT.run == "headline")].dropna(subset=["D80"])
    dd = np.arange(-4, 16)
    for cl in ["shallow", "mid", "deep"]:
        Ms = []
        for _, r in c[c.dcl == cl].iterrows():
            g = eall[(eall.ep == r.ep) & (eall.patch == r.patch)]
            day = pd.to_datetime(g.tpt).dt.floor("D")
            vc = day.value_counts(); pk = vc.idxmax()
            Ms.append([vc.get(pk + pd.Timedelta(days=int(k)), 0) / vc.max() for k in dd])
        if not Ms:
            continue
        Ms = np.array(Ms)
        a.plot(dd, np.median(Ms, 0), "o-", color=CL[cl], lw=2, ms=4, label=f"{cl} ({len(Ms)} patch-episodes)")
        a.fill_between(dd, np.quantile(Ms, 0.25, 0), np.quantile(Ms, 0.75, 0), color=CL[cl], alpha=0.12, lw=0)
    a.set_xlabel("days from the patch's busiest day (PT)")
    a.set_ylabel("daily LFEs / busiest-day count")
    a.set_title("one patch-episode at a time\n(median, IQR band)", loc="left")
    a.legend(frameon=False, fontsize=8)
    # duration vs down-dip
    a = ax[i, 2]
    c = c.dropna(subset=["D80_n30"])
    for cl in ["shallow", "mid", "deep"]:
        cc = c[c.dcl == cl]
        a.scatter(cc.dd, cc.D80_n30, s=22, color=CL[cl], edgecolor="white", lw=0.8, label=cl)
    d = res["duration"][f"{reg}_headline"]
    a.set_xlabel("patch down-dip position (km, + = deeper)")
    a.set_ylabel("active duration D80, days")
    a.set_title(f"duration vs depth\nwithin-episode ρ = {d['rho_D80n30_dd']:+.2f}, p = {d['p']:.2f}", loc="left")
    a.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/c_rate_vs_front.png", dpi=150); plt.close(fig)

# ---------------- figure 2: inter-event times
Z = np.load(f"{R}/sizes_iets.npz")
fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
for a, key, title in [(ax[0], "SVI_headline", "SVI headline"), (ax[1], "NW_headline", "Washington-NW headline"),
                      (ax[2], "SVI_scrambled", "SVI scrambled control (P ≥ 0.5)")]:
    nb = np.logspace(-4, 1.5, 45)
    for k, col, lab in [("obs", OBS, "observed"), ("day", SUR, "shuffled within each day"), ("uni", UNI, "uniform (Poisson)")]:
        v = Z[f"{key}_iet_{k}"]
        if len(v) == 0:
            continue
        h, _ = np.histogram(v, nb, density=True)
        a.loglog(np.sqrt(nb[1:] * nb[:-1]), np.where(h > 0, h, np.nan), color=col, lw=2 if k == "obs" else 1.5, label=lab)
    q = res["q2"].get(key, {})
    if "day" in q:
        a.set_title(f"{title}\nmedian CV {q['day']['obs_median_cv']:.2f} vs day-shuffle {q['day']['sur_median_cv']:.2f} "
                    f"(Poisson 1)", loc="left")
    a.set_xlabel("inter-event time / patch-episode mean")
    a.set_ylabel("probability density")
    a.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/c_iet.png", dpi=150); plt.close(fig)

# ---------------- figure 3: burst sizes
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
for a, key, title in [(ax[0], "SVI_headline", "SVI"), (ax[1], "NW_headline", "Washington-NW")]:
    for k, col, lab in [("obs", OBS, "observed"), ("day", SUR, "shuffled within each day")]:
        s = np.sort(Z[f"{key}_size_{k}"])
        if len(s) == 0:
            continue
        u = np.arange(2, s.max() + 1)
        ccdf = np.array([(s >= x).mean() for x in u])
        a.loglog(u, ccdf, "o-" if k == "obs" else "-", color=col, ms=4, lw=1.5, label=f"{lab} (n={len(s)})")
    f = res["burst_sizes"][key]["obs"]
    s = Z[f"{key}_size_obs"]; u = np.arange(2, s.max() + 1)
    pl = special.zeta(f["alpha"], u) / special.zeta(f["alpha"], 2)
    ge = f["q"] ** (u - 2)
    a.loglog(u, pl, ":", color="#52514e", lw=1.2, label=f"power-law fit α={f['alpha']:.2f}")
    a.loglog(u, ge, "--", color="#52514e", lw=1.2, label="exponential fit")
    a.set_ylim(1e-4 if len(s) > 1000 else 1e-3, 1.2)
    a.set_xlabel("burst size (events linked by gaps < 0.1 × mean)"); a.set_ylabel("fraction of bursts ≥ size")
    a.set_title(f"{title}: Vuong z = {f['z']:+.1f} → {f['verdict']}", loc="left")
    a.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/c_burst_sizes.png", dpi=150); plt.close(fig)
print("ok")
