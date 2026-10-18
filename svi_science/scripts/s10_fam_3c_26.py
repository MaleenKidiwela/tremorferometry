"""Family 10 at PGC: daily + per-episode dv/v with all 3 components and a 2-6 s window (user 2026-10-12).
Same pipeline as the pilot/SVD runs: P>=0.5 rows, per-event 3-C realignment on the S window (s10_measure.realign),
P-weighted day stacks per component, SVD-Wiener per component per epoch (s10_svd_measure.svdw), reference =
realigned epoch grand stack. Joint 3-C S-anchored stretch: one eps maximising the cc of the concatenated
(Z,N,E) windows tS+2..tS+6 s. Compared with Z-only 2-4 s through the same code. dv/v = -eps. No error bars."""
import sys, numpy as np, pandas as pd, matplotlib
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import s10_common as S, s10_measure as M, s10_svd_measure as V
W1 = float(sys.argv[1]) if len(sys.argv) > 1 else 6.0; TAG = f"3c2{int(W1)}"
matplotlib.use("Agg"); import matplotlib.pyplot as plt
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"; FAM = 10; NS = "CN.PGC"
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet"); rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet"); rs = rs[rs.netsta == NS][["row_id", "epoch"]]
ser = pd.read_csv(f"{PIL}/tables/series.csv"); ser = ser[(ser.fam == FAM) & (ser.netsta == NS)].set_index("epoch")
z = np.load(f"{PIL}/win/{NS}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}; tau = z["tau"]
EPS = np.linspace(-0.01, 0.01, 801)
def stretch(cur, ref, tS, w0, w1, comps):
    m = (tau >= tS + w0) & (tau <= tS + w1); t = tau[m]
    R = np.concatenate([(ref[c][m] - ref[c][m].mean()) / (ref[c][m].std() + 1e-12) for c in comps])
    best = (-2, 0.0)
    for e in EPS:
        X = np.concatenate([(lambda x: (x - x.mean()) / (x.std() + 1e-12))(np.interp(tS + (t - tS) * (1 + e), tau, cur[c])) for c in comps])
        cc = float(X @ R / (np.linalg.norm(X) * np.linalg.norm(R) + 1e-12))
        if cc > best[0]: best = (cc, e)
    return -best[1] * 100, best[0]
out = []
for ep_id, rr in rs.groupby("epoch"):
    g = rows[rows.row_id.isin(set(rr.row_id)) & rows.row_id.isin(pos)].reset_index(drop=True)
    if len(g) < 30: continue
    g["date"] = pd.to_datetime(g.date); W = z["W"][[pos[r] for r in g.row_id]].astype(np.float64)
    tS = float(ser.loc[ep_id, "tS_pick"]); w = g.P.to_numpy()
    G = np.array([np.nansum(np.nan_to_num(W[:, c]) * w[:, None], 0) / w.sum() for c in range(3)])
    Wa = M.realign(W, G, tS)[0]
    Ga = np.array([np.nansum(np.nan_to_num(Wa[:, c]) * w[:, None], 0) / w.sum() for c in range(3)])
    days = sorted(g.date.unique()); D = np.zeros((len(days), 3, len(tau))); keep = []
    for k, d in enumerate(days):
        gg = g.index[g.date == d]; ww = w[gg]
        if ww.sum() < 3: continue
        D[k] = [np.nansum(np.nan_to_num(Wa[gg, c]) * ww[:, None], 0) / ww.sum() for c in range(3)]; keep.append(k)
    D = D[keep]; days = [days[k] for k in keep]
    noise = np.flatnonzero(tau < tS - 8)
    Df = np.stack([V.svdw(D[:, c], noise) for c in range(3)], 1)
    for k, d in enumerate(days):
        rec = dict(epoch=ep_id, date=d)
        for lab, X in (("raw", D[k]), ("svd", Df[k])):
            rec[f"{lab}_3c"], rec[f"{lab}_cc3c26"] = stretch(X, Ga, tS, 2, W1, (0, 1, 2))
            rec[f"{lab}_z24"], rec[f"{lab}_ccz24"] = stretch(X, Ga, tS, 2, 4, (0,))
        out.append(rec)
T = pd.DataFrame(out).sort_values("date"); T.to_csv(f"{PIL}/tables/fam10_PGC_{TAG}_daily.csv", index=False)
dd = np.sort(T.date.unique()); epn = np.r_[0, np.cumsum(np.diff(dd).astype("timedelta64[D]").astype(int) > 20)]; T["ep"] = T.date.map(dict(zip(dd, epn)))
def per_ep(col, cc):
    K = T[(T[cc] > 0.8) & (T[col].abs() < 0.99)]
    return K, K.groupby("ep").apply(lambda q: pd.Series(dict(date=q.date.min(), dvv=np.average(q[col], weights=q[cc]), cc=q[cc].mean(), nd=len(q))), include_groups=False)
fig, axs = plt.subplots(2, 2, figsize=(16, 8), sharex=True)
for j, lab in enumerate(["raw", "svd"]):
    for i, (col, cc, name) in enumerate([(f"{lab}_z24", f"{lab}_ccz24", "Z only, 2–4 s"), (f"{lab}_3c", f"{lab}_cc3c26", f"3 components, 2–{W1:g} s")]):
        K, E = per_ep(col, cc); ax = axs[i, j]
        ax.scatter(K.date, K[col], c=K[cc], cmap="viridis", vmin=0.8, vmax=1, s=10, alpha=0.6, lw=0)
        ax.plot(E.date, E.dvv, "-", color="k", lw=1.2, zorder=3)
        sc = ax.scatter(E.date, E.dvv, c=E.cc, cmap="viridis", vmin=0.8, vmax=1, s=25 + 10 * E.nd, edgecolor="k", lw=0.6, zorder=4)
        ax.axhline(0, color="grey", lw=0.6, ls=":"); ax.set_ylim(-0.5, 0.5)
        ax.set_title(f"{'raw' if lab == 'raw' else 'SVD'}, {name}\nday sd {K[col].std():.2f} %, episode sd {E.dvv.std():.3f} %, {len(K)} days CC>0.8", fontsize=9)
        ax.set_ylabel("dv/v (%)")
fig.colorbar(sc, ax=axs, shrink=0.5, label="CC with epoch reference")
fig.suptitle("Family 10 at PGC: small dots = days (CC > 0.8), big dots + line = per-episode CC-weighted mean", fontsize=11)
outp = f"{PIL}/figures/s10_fam10_PGC_{TAG}.png"; fig.savefig(outp, dpi=115, bbox_inches="tight"); print(outp)
for lab in ["raw", "svd"]:
    for col, cc in [(f"{lab}_z24", f"{lab}_ccz24"), (f"{lab}_3c", f"{lab}_cc3c26")]:
        K, E = per_ep(col, cc); within = K.groupby("ep")[col].std().median()
        print(f"{col:10s}: days CC>0.8 {len(K):3d}, median CC {K[cc].median():.3f}, day sd {K[col].std():.3f} %, within-episode sd {within:.3f} %, episode-mean sd {E.dvv.std():.3f} %")
