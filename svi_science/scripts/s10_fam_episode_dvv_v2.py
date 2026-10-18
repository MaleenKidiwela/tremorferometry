"""Family per-episode dv/v, v2 (user 2026-10-11): one reference per station over ALL time (no sensor epochs),
2-4 s after picked S only, Z component, P>=0.5 weighted episode stacks; points with CC <= 0.6 dropped.
dv/v = -eps, S-anchored stretch (current sampled at tS+(t-tS)(1+eps) vs reference). Plot 2: episode stacks as
wiggles, y = episode date, x = time after predicted S."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"; FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
STAS = ["CN.PGC", "PB.B011", "PB.B009", "PB.B010", "PB.B926", "CN.VGZ", "CN.SNB"]
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet"); rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet"); ser = pd.read_csv(f"{PIL}/tables/series.csv")
EPS = np.linspace(-0.01, 0.01, 801)
def stretch(cur, ref, tau, tS):
    m = (tau >= tS + 2) & (tau <= tS + 4); t = tau[m]; r = ref[m] - ref[m].mean(); best = (-2, 0)
    for e in EPS:
        x = np.interp(tS + (t - tS) * (1 + e), tau, cur); x = x - x.mean()
        c = float(x @ r / (np.linalg.norm(x) * np.linalg.norm(r) + 1e-12))
        if c > best[0]: best = (c, e)
    return -best[1] * 100, best[0]
res, stacks = [], {}
for ns in STAS:
    z = np.load(f"{PIL}/win/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}; tau = z["tau"]
    g = rows.merge(rs[rs.netsta == ns][["row_id"]], on="row_id"); g = g[g.row_id.isin(pos)].copy()
    if len(g) == 0: continue
    g["date"] = pd.to_datetime(g.date); g = g.reset_index(drop=True)
    W = z["W"][[pos[r] for r in g.row_id], 0, :].astype(np.float64); W[~np.isfinite(W)] = 0
    ok = np.abs(W).sum(1) > 0; g = g[ok].reset_index(drop=True); W = W[ok]
    tS = float(ser[(ser.fam == FAM) & (ser.netsta == ns)].tS_pick.median())
    ref = (W * g.P.to_numpy()[:, None]).sum(0) / g.P.sum()
    days = np.sort(g.date.unique()); ep = np.r_[0, np.cumsum(np.diff(days).astype("timedelta64[D]").astype(int) > 20)]
    g["ep"] = g.date.map(dict(zip(days, ep)))
    for e, gg in g.groupby("ep"):
        w = gg.P.to_numpy()
        if w.sum() < 10: continue
        st = (W[gg.index] * w[:, None]).sum(0) / w.sum()
        dvv, cc = stretch(st, ref, tau, tS)
        res.append(dict(station=ns, date=gg.date.min(), n=len(gg), sumP=w.sum(), dvv=dvv, cc=cc))
        stacks[(ns, gg.date.min())] = st
    stacks[(ns, "tau")] = tau; stacks[(ns, "tS")] = tS
T = pd.DataFrame(res); T.to_csv(f"{PIL}/tables/fam{FAM}_episode_dvv_v2.csv", index=False)
K = T[T.cc > 0.6]
fig, axs = plt.subplots(len(STAS) + 1, 1, figsize=(13, 2.0 * (len(STAS) + 1)), sharex=True)
for ax, ns in zip(axs, STAS):
    h = K[K.station == ns]
    sc = ax.scatter(h.date, h.dvv, c=h.cc, cmap="viridis", vmin=0.6, vmax=1, s=45, edgecolor="k", lw=0.4)
    ax.axhline(0, color="grey", lw=0.6); ax.set_ylim(-0.6, 0.6); ax.set_ylabel(f"{ns.split('.')[1]}\ndv/v (%)", fontsize=9)
    ax.text(0.005, 0.85, f"{len(h)} episodes with CC > 0.6 (of {(T.station == ns).sum()})", transform=ax.transAxes, fontsize=8)
p = K.groupby(K.date.dt.to_period("M")).apply(lambda g: pd.Series(dict(date=g.date.min(), dvv=np.average(g.dvv, weights=g.cc ** 2),
        cc=g.cc.mean(), nst=g.station.nunique())), include_groups=False)
ax = axs[-1]; ax.scatter(p.date, p.dvv, c=p.cc, cmap="viridis", vmin=0.6, vmax=1, s=20 + 15 * p.nst, edgecolor="k", lw=0.5)
ax.axhline(0, color="grey", lw=0.6); ax.set_ylim(-0.4, 0.4); ax.set_ylabel("all stations\ndv/v (%)", fontsize=9)
ax.text(0.005, 0.85, "CC²-weighted mean over stations; size = number of stations", transform=ax.transAxes, fontsize=8)
fig.colorbar(sc, ax=axs, shrink=0.4, pad=0.01, label="CC with all-time reference")
axs[0].set_title(f"Family {FAM}: dv/v per episode, 2–4 s after picked S (Z), one all-time reference per station, CC > 0.6", fontsize=10)
fig.savefig(f"{PIL}/figures/s10_fam{FAM}_episode_dvv_v2.png", dpi=115, bbox_inches="tight")
# plot 2: wiggles, y = episode date
SHOW = ["CN.PGC", "PB.B011"]
fig, axs = plt.subplots(1, len(SHOW), figsize=(7 * len(SHOW), 11), sharey=True)
for ax, ns in zip(axs, SHOW):
    tau = stacks[(ns, "tau")]; tS = stacks[(ns, "tS")]; m = (tau >= -1) & (tau <= 5)
    keys = sorted([k for k in stacks if k[0] == ns and not isinstance(k[1], str)], key=lambda k: k[1])
    for k in keys:
        s = stacks[k][m]; s = s / (np.abs(s).max() + 1e-12); y = pd.Timestamp(k[1]).year + pd.Timestamp(k[1]).dayofyear / 365.25
        ax.plot(tau[m], y - 0.45 * s, "k", lw=0.6); ax.fill_between(tau[m], y, y - 0.45 * s, where=s > 0, color="k", lw=0)
    ax.axvspan(tS + 2, tS + 4, color="orange", alpha=0.2); ax.axvline(tS, color="orange", lw=1)
    ax.set_xlim(-1, 5); ax.set_xlabel("time after predicted S (s); orange = picked S and 2–4 s window")
    ax.set_title(f"family {FAM} at {ns}: one stack per episode (Z, 2–8 Hz)", fontsize=10)
axs[0].set_ylabel("episode date (year)"); axs[0].invert_yaxis()
fig.tight_layout(); fig.savefig(f"{PIL}/figures/s10_fam{FAM}_episode_wiggles.png", dpi=115)
print(T.groupby("station").agg(n=("dvv", "size"), kept=("cc", lambda c: (c > 0.6).sum()), cc_med=("cc", "median")).to_string())
