"""Are family-10 P and S arrival times stable over 2003-2026? Per-episode P-weighted stacks (P>=0.5, unaligned windows,
time relative to predicted S from the catalogue). Lag of P (Z, window around the P arrival) and S (horizontals N+E,
window around S) vs the station-epoch grand stack, by cross-correlation with parabolic sub-sample refinement.
S-P change = S lag - P lag (independent of origin time). User request 2026-10-11. Plot: no error bars, colour = CC."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"; FAM = 10
STAS = ["CN.PGC", "PB.B011"]
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet"); rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet")
def lag_cc(a, b, fs, maxlag=0.5):
    a = a - a.mean(); b = b - b.mean(); L = int(maxlag * fs); best = (-2, 0)
    cs = []
    for k in range(-L, L + 1):
        x = a[max(0, k):len(a) + min(0, k)]; y = b[max(0, -k):len(b) + min(0, -k)]
        c = float(x @ y / (np.linalg.norm(x) * np.linalg.norm(y) + 1e-12)); cs.append(c)
    cs = np.array(cs); i = int(np.argmax(cs)); d = 0.0
    if 0 < i < len(cs) - 1:
        den = cs[i - 1] - 2 * cs[i] + cs[i + 1]; d = 0.5 * (cs[i - 1] - cs[i + 1]) / den if den != 0 else 0.0
    return (i - L + d) / fs, cs[i]
out = []
for ns in STAS:
    z = np.load(f"{PIL}/win/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}; tau = z["tau"]; fs = 1 / (tau[1] - tau[0])
    g = rows.merge(rs[rs.netsta == ns][["row_id", "epoch"]], on="row_id"); g = g[g.row_id.isin(pos)].copy()
    g["date"] = pd.to_datetime(g.date); W = z["W"][[pos[r] for r in g.row_id]].astype(np.float64); W[~np.isfinite(W)] = 0
    days = np.sort(g.date.unique()); ep = np.r_[0, np.cumsum(np.diff(days).astype("timedelta64[D]").astype(int) > 20)]
    g["ep"] = g.date.map(dict(zip(days, ep))); g = g.reset_index(drop=True)
    for epoch, ge in g.groupby("epoch"):
        w = ge.P.to_numpy(); Gz = (W[ge.index, 0] * w[:, None]).sum(0) / w.sum()
        Gh = [(W[ge.index, c] * w[:, None]).sum(0) / w.sum() for c in (1, 2)]
        # P window: around the P peak on Z in [-9, -3] s; S window: [-0.5, 1.5] s on horizontals
        mP = (tau > -9) & (tau < -3); iP = np.flatnonzero(mP)[np.argmax(np.abs(Gz[mP]))]; pw = (tau > tau[iP] - 1.0) & (tau < tau[iP] + 1.5)
        sw = (tau > -0.5) & (tau < 1.5)
        for e, gg in ge.groupby("ep"):
            if gg.P.sum() < 15: continue
            ww = gg.P.to_numpy(); ez = (W[gg.index, 0] * ww[:, None]).sum(0) / ww.sum()
            lp, cp = lag_cc(ez[pw], Gz[pw], fs)
            ls, cs_ = [], []
            for c in (1, 2):
                eh = (W[gg.index, c] * ww[:, None]).sum(0) / ww.sum(); l, cc = lag_cc(eh[sw], Gh[c - 1][sw], fs); ls.append(l); cs_.append(cc)
            out.append(dict(station=ns, epoch=epoch, date=gg.date.min(), n=len(gg), P_lag_ms=lp * 1000, P_cc=cp,
                            S_lag_ms=np.mean(ls) * 1000, S_cc=np.mean(cs_), SmP_ms=(np.mean(ls) - lp) * 1000, P_time=tau[iP]))
T = pd.DataFrame(out); T.to_csv(f"{PIL}/tables/fam10_PS_timing_by_episode.csv", index=False)
fig, axs = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
for ns, mk in zip(STAS, ["s", "o"]):
    h = T[T.station == ns]
    for ax, col, ccol, lab in [(axs[0], "P_lag_ms", "P_cc", "P arrival shift (ms, Z)"), (axs[1], "S_lag_ms", "S_cc", "S arrival shift (ms, horizontals)"),
                               (axs[2], "SmP_ms", "P_cc", "S−P change (ms)")]:
        sc = ax.scatter(h.date, h[col], c=h[ccol], cmap="viridis", vmin=0.3, vmax=1, marker=mk, s=40, edgecolor="k", lw=0.3, label=ns)
        ax.axhline(0, color="grey", lw=0.6); ax.set_ylabel(lab)
axs[0].legend(fontsize=8); fig.colorbar(sc, ax=axs, shrink=0.5, label="CC with grand stack")
axs[0].set_title("Family 10: per-episode arrival-time shifts vs grand stack (each station-epoch its own reference); 25 ms = 1 sample at 40 Hz", fontsize=10)
out_png = f"{PIL}/figures/s10_fam10_PS_timing.png"; fig.savefig(out_png, dpi=120); print(out_png)
for ns in STAS:
    h = T[T.station == ns]; print(ns, f"episodes {len(h)}; P shift sd {h.P_lag_ms.std():.1f} ms (range {h.P_lag_ms.min():.0f}..{h.P_lag_ms.max():.0f}); "
          f"S shift sd {h.S_lag_ms.std():.1f} ms (range {h.S_lag_ms.min():.0f}..{h.S_lag_ms.max():.0f}); S-P change sd {h.SmP_ms.std():.1f} ms "
          f"(range {h.SmP_ms.min():.0f}..{h.SmP_ms.max():.0f}); median P cc {h.P_cc.median():.2f}, S cc {h.S_cc.median():.2f}")
