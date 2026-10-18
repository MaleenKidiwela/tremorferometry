"""Family per-episode dv/v at PGC only (user 2026-10-11): one all-time PGC reference (no sensor epochs), 2-4 s after
picked S, Z, P>=0.5-weighted episode stacks, points with CC <= 0.6 dropped. dv/v = -eps (S-anchored stretch).
Plot 2: episode stacks as wiggles, y = episode date, x = time after predicted S."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"; FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10; NS = "CN.PGC"
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet"); rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet"); ser = pd.read_csv(f"{PIL}/tables/series.csv")
z = np.load(f"{PIL}/win/{NS}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}; tau = z["tau"]
g = rows.merge(rs[rs.netsta == NS][["row_id"]], on="row_id"); g = g[g.row_id.isin(pos)].reset_index(drop=True)
g["date"] = pd.to_datetime(g.date)
W = z["W"][[pos[r] for r in g.row_id], 0, :].astype(np.float64); W[~np.isfinite(W)] = 0
ok = np.abs(W).sum(1) > 0; g = g[ok].reset_index(drop=True); W = W[ok]
tS = float(ser[(ser.fam == FAM) & (ser.netsta == NS)].tS_pick.median())
ref = (W * g.P.to_numpy()[:, None]).sum(0) / g.P.sum()
EPS = np.linspace(-0.01, 0.01, 801); m24 = (tau >= tS + 2) & (tau <= tS + 4); t24 = tau[m24]; r = ref[m24] - ref[m24].mean()
def stretch(cur):
    best = (-2, 0.0)
    for e in EPS:
        x = np.interp(tS + (t24 - tS) * (1 + e), tau, cur); x = x - x.mean()
        c = float(x @ r / (np.linalg.norm(x) * np.linalg.norm(r) + 1e-12))
        if c > best[0]: best = (c, e)
    return -best[1] * 100, best[0]
def align(cur, maxlag=0.25):
    """shift an episode stack so its S wave (tS-1..tS+2 s) lines up with the reference (cc, sub-sample); same idea as the pilot realign."""
    mS = (tau >= tS - 1) & (tau <= tS + 2); a = ref[mS] - ref[mS].mean(); best = (-2, 0.0)
    for dt in np.arange(-maxlag, maxlag + 1e-9, 0.002):
        x = np.interp(tau[mS] + dt, tau, cur); x = x - x.mean(); c = float(x @ a / (np.linalg.norm(x) * np.linalg.norm(a) + 1e-12))
        if c > best[0]: best = (c, dt)
    return np.interp(tau + best[1], tau, cur), best[1]
days = np.sort(g.date.unique()); ep = np.r_[0, np.cumsum(np.diff(days).astype("timedelta64[D]").astype(int) > 20)]
g["ep"] = g.date.map(dict(zip(days, ep)))
res, stacks = [], []
for e, gg in g.groupby("ep"):
    w = gg.P.to_numpy()
    if w.sum() < 10: continue
    st0 = (W[gg.index] * w[:, None]).sum(0) / w.sum(); st, sh = align(st0); dvv, cc = stretch(st); dvv0, cc0 = stretch(st0)
    res.append(dict(date=gg.date.min(), n=len(gg), sumP=round(w.sum(), 1), S_shift_ms=round(sh * 1000, 1), dvv=dvv, cc=cc, dvv_unaligned=dvv0)); stacks.append((gg.date.min(), st))
T = pd.DataFrame(res); T.to_csv(f"{PIL}/tables/fam{FAM}_PGC_episode_dvv.csv", index=False)
K = T[T.cc > 0.6]
fig, ax = plt.subplots(figsize=(13, 4))
sc = ax.scatter(K.date, K.dvv, c=K.cc, cmap="viridis", vmin=0.6, vmax=1, s=60, edgecolor="k", lw=0.5)
ax.axhline(0, color="grey", lw=0.6); ax.set_ylabel("dv/v (%)"); fig.colorbar(sc, ax=ax, label="CC with all-time PGC reference")
ax.set_title(f"Family {FAM} at PGC: dv/v per episode, 2–4 s after picked S (Z); one all-time reference; S-aligned; CC > 0.6 ({len(K)} of {len(T)} episodes)", fontsize=10)
fig.tight_layout(); fig.savefig(f"{PIL}/figures/s10_fam{FAM}_PGC_episode_dvv.png", dpi=120)
fig, ax = plt.subplots(figsize=(8, 12)); mw = (tau >= -1) & (tau <= 5)
for d0, s in stacks:
    s = s[mw] / (np.abs(s[mw]).max() + 1e-12); y = d0.year + (d0.dayofyear - 1) / 365.25
    ax.plot(tau[mw], y - 0.45 * s, "k", lw=0.6); ax.fill_between(tau[mw], y, y - 0.45 * s, where=s > 0, color="k", lw=0)
ax.axvspan(tS + 2, tS + 4, color="orange", alpha=0.2); ax.axvline(tS, color="orange", lw=1)
ax.set_xlim(-1, 5); ax.invert_yaxis(); ax.set_ylabel("episode date (year)"); ax.set_xlabel("time after predicted S (s); orange = picked S and 2–4 s window")
ax.set_title(f"Family {FAM} at PGC: one stack per episode (Z, 2–8 Hz)", fontsize=10)
fig.tight_layout(); fig.savefig(f"{PIL}/figures/s10_fam{FAM}_PGC_episode_wiggles.png", dpi=120)
print(T.round(3).to_string(index=False))
