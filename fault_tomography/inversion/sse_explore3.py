#!/usr/bin/env python
"""Exploration round 3: MECHANISM of the station-specific ETS-locked response, and the remaining predictors.
D1 shrinkage-release test (Merlin mechanism e): per family, ETS contrast vs (PRE level - family all-time mean). SVD-Wiener
   shrinks low-SNR rows toward the all-time reference; higher ETS SNR releases it -> contrast ~ -(PRE - mean).
D2 detection-count test (mechanism b): per family, contrast vs log(n_det ON / n_det PRE) from the long_window npz n_det.
H  station response vs geometry: distance to footprint, up-dip/down-dip of the slip patch, distance to tremor, instrument.
I  amplitude scaling: per-family contrast vs episode footprint slip (pooled Spearman, shift null).
J  local vs regional timing: per station, r(station family-mean dv/v, LOCAL slip-in-window at nearest footprint cells)
   vs r(.., SEGMENT-mean slip); lag of max |xcorr|.
L  tremor vs slip as predictor: per station, r(dv/v 30-d change, local tremor count in trailing 30 d) vs r(.., local slip).
F  3-5 s coda as primary observable: pooled N+C episode contrast of the family 3-5 s product with both nulls.
Outputs: res_catalog_sse/explore3_*.csv, figures/sse_explore3.png; all printed.
"""
import numpy as np, pandas as pd, os, glob
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(41)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
clat, clon, dep = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"]; gi = pd.Series(np.arange(G), index=days); inc = zc["slip_thrust_inc30"].astype(float) * 1e3
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok & ep.seg.isin(["N", "C", "S"])].reset_index(drop=True); zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
E0, E1 = int(gi[pd.Timestamp("2011-04-01")]), int(gi[pd.Timestamp("2026-04-01")]); span = E1 - E0; wrap = lambda i: E0 + (i - E0) % span
pairs = pd.read_csv(f"{RD}/pairs.csv"); sll = pairs.groupby("tag")[["sta_lat", "sta_lon"]].first()
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam", "cell", "date", "dvv", "cc"]); cpos = {c: i for i, c in enumerate(cells)}
fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos); fd["v"] = -fd.dvv
tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"]); tr["t"] = pd.to_datetime(tr.time, errors="coerce"); tr = tr.dropna(subset=["t"]); tr["d"] = tr.t.dt.floor("D"); TXY = km(tr.lat.values, tr.lon.values)
win = np.asarray((days >= "2011-01-01") & (days <= "2026-07-20"))
def deseason(y): s = pd.Series(y, index=days); return (s - s.groupby(s.index.month).transform("mean")).values
def d30(y): o = np.full_like(y, np.nan); o[30:] = y[30:] - y[:-30]; return o
def corr(a, b): m = np.isfinite(a) & np.isfinite(b) & win; return np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 500 else np.nan
# episode geometry
EP = []
for _, r in okep.iterrows():
    fp = zf[f"fp_{r.key}"]; EP.append(dict(key=r.key, seg=r.seg, ion=int(gi[pd.Timestamp(r.fp_onset)]), fp=fp, cxy=XY[fp].mean(0), size=r.fp_slip_sum_mm, updip_dir=None))
# per-family contrasts + covariates, all stations with >=4 episodes within 150 km
recs = []; sta_series = {}
for tag, g in fd.groupby("tag"):
    if tag not in sll.index: continue
    sxy = km(np.array([sll.sta_lat[tag]]), np.array([sll.sta_lon[tag]]))[0]
    eps = [e for e in EP if np.sqrt(((XY[e["fp"]] - sxy) ** 2).sum(1)).min() <= 150]
    if len(eps) < 4: continue
    piv = g.pivot_table(index="date", columns="fam", values="v").reindex(days); allmean = piv.mean(0); pivd = piv - piv.groupby(piv.index.month).transform("mean")
    ccp = g.pivot_table(index="date", columns="fam", values="cc").reindex(days)
    sta_series[tag] = dict(y=pivd.mean(axis=1).values, xy=sxy, eps=eps, pivd=pivd.values.astype(np.float32), piv=piv.values.astype(np.float32), allmean=allmean.values.astype(np.float32))
    for fam in piv.columns:
        y = pivd[fam].values.astype(float); yr = piv[fam].values.astype(float); c = ccp[fam].values.astype(float)
        for e in eps:
            i = e["ion"]; on = slice(max(i + 15, 0), i + 61); pre = slice(max(i - 90, 0), i - 29)
            if np.isfinite(y[on]).mean() < .6 or np.isfinite(y[pre]).mean() < .6: continue
            recs.append(dict(tag=tag, fam=fam, key=e["key"], contrast=np.nanmean(y[on]) - np.nanmean(y[pre]), pre_rel=np.nanmean(yr[pre]) - allmean[fam],
                             dcc=np.nanmean(c[on]) - np.nanmean(c[pre]), size=e["size"], dist=float(np.sqrt(((XY[e["fp"]] - sxy) ** 2).sum(1)).min()),
                             updip=float(sxy[0] - e["cxy"][0]) < 0))   # station west (up-dip) of the footprint centre
R = pd.DataFrame(recs); R.to_csv(f"{SD}/explore3_family_episode_contrasts.csv", index=False)
print(f"family-episode contrasts: {len(R):,} rows, {R.tag.nunique()} stations, {R.key.nunique()} episodes", flush=True)
from scipy.stats import spearmanr
# D1 shrinkage release
rs = spearmanr(R.contrast, R.pre_rel).correlation
# mean-reversion BASELINE: the same statistic with episode times rigidly shifted (contrast = ON-PRE always anticorrelates with PRE level for pure noise)
base = []
for k in range(60):
    sh = rng.randint(90, span - 90); cc_, pr_ = [], []
    for tag, s_ in sta_series.items():
        pivd, piv, allmean = s_["pivd"], s_["piv"], s_["allmean"]
        for e in s_["eps"]:
            i = wrap(e["ion"] + sh); on = slice(max(i + 15, 0), i + 61); pre = slice(max(i - 90, 0), i - 29)
            a = np.nanmean(pivd[on], 0); b = np.nanmean(pivd[pre], 0); br = np.nanmean(piv[pre], 0) - allmean
            m = np.isfinite(a) & np.isfinite(b); cc_.append(a[m] - b[m]); pr_.append(br[m])
    base.append(spearmanr(np.concatenate(cc_), np.concatenate(pr_)).correlation)
print(f"D1 SHRINKAGE-RELEASE: Spearman(contrast, PRE level - all-time mean) at TRUE episodes = {rs:+.3f}; at SHIFTED times (mean-reversion baseline) {np.mean(base):+.3f} +- {np.std(base):.3f} -> excess {rs - np.mean(base):+.3f} ({(rs - np.mean(base)) / np.std(base):+.1f} sd); mechanism predicts a significant NEGATIVE excess")
st_ = R.groupby("tag").apply(lambda d: spearmanr(d.contrast, d.pre_rel).correlation if len(d) > 30 else np.nan); print(f"   per station: median {st_.median():+.2f}, fraction < -0.3: {(st_ < -0.3).mean():.2f}, > +0.3: {(st_ > 0.3).mean():.2f}")
# cc change
rc = spearmanr(R.contrast, R.dcc).correlation; print(f"   Spearman(contrast, cc change ON-PRE) = {rc:+.3f}")
# D2 n_det from long_window npz for the strongest-responding stations
tags_nd = ["ANGUS", "B004p90f40", "RADR", "B033p90f40", "B935p90f40", "LOPEZ", "B024p90f40", "BFSB", "B928p90f40", "BTB", "B011p90f40", "PGChhz", "B018p90f40", "STOR"]
nd_rows = []
for tag in tags_nd:
    fn = f"data/long_window_daily_{tag}_Z.npz"
    if not os.path.exists(fn) or tag not in sta_series: continue
    z = np.load(fn, allow_pickle=True); nd = pd.DataFrame(dict(fam=np.asarray(z["patches"]).astype(str), date=pd.to_datetime(np.asarray(z["dates"]).astype(str)), n=np.asarray(z["n_det"], float)))
    piv = nd.pivot_table(index="date", columns="fam", values="n", aggfunc="sum").reindex(days)
    sub = R[R.tag == tag]
    for _, r in sub.iterrows():
        if r.fam not in piv.columns: continue
        i = [e["ion"] for e in sta_series[tag]["eps"] if e["key"] == r.key][0]; y = piv[r.fam].values
        on, pre = np.nansum(y[max(i + 15, 0):i + 61]), np.nansum(y[max(i - 90, 0):i - 29])
        nd_rows.append(dict(tag=tag, fam=r.fam, key=r.key, contrast=r.contrast, lr=np.log((on + 1) / (pre * 46 / 61 + 1))))
N_ = pd.DataFrame(nd_rows)
if len(N_):
    rn = spearmanr(N_.contrast, N_.lr).correlation; per = N_.groupby("tag").apply(lambda d: pd.Series(dict(n=len(d), r=spearmanr(d.contrast, d.lr).correlation, med_lr=d.lr.median(), mean_c=d.contrast.mean())))
    print(f"D2 DETECTION COUNT: Spearman(contrast, log n_det ON/PRE) pooled {rn:+.3f} over {len(N_):,} family-episodes; median log-ratio {N_.lr.median():.2f} (= x{np.exp(N_.lr.median()):.1f} detections during ON)")
    print(per.round(3).to_string())
# H geometry + instrument
S = R.groupby("tag").agg(c=("contrast", "mean"), fpos=("contrast", lambda v: (v > 0).mean()), n=("contrast", "size"), dist=("dist", "median"), updip=("updip", "mean")).reset_index()
S["borehole"] = S.tag.str.contains("p90f40"); S["lat"] = S.tag.map(sll.sta_lat); S["lon"] = S.tag.map(sll.sta_lon)
S["dist_tremor"] = [np.sqrt(((TXY - km(np.array([sll.sta_lat[t]]), np.array([sll.sta_lon[t]]))[0]) ** 2).sum(1)).min() for t in S.tag]
print(f"\nH GEOMETRY ({len(S)} stations): Spearman(c, dist to footprint) {spearmanr(S.c, S.dist).correlation:+.2f} | (c, lat) {spearmanr(S.c, S.lat).correlation:+.2f} | (|c|, dist to nearest tremor) {spearmanr(S.c.abs(), S.dist_tremor).correlation:+.2f} | up-dip stations mean c {S[S.updip > 0.5].c.mean():+.3f} (n={(S.updip > 0.5).sum()}) vs down-dip {S[S.updip <= 0.5].c.mean():+.3f} | borehole {S[S.borehole].c.mean():+.3f} vs broadband {S[~S.borehole].c.mean():+.3f}; |c| borehole {S[S.borehole].c.abs().mean():.3f} vs broadband {S[~S.borehole].c.abs().mean():.3f}")
# I amplitude scaling
ri = spearmanr(R.contrast.abs(), R["size"]).correlation; rs2 = R.groupby(["tag"]).apply(lambda d: spearmanr(d.contrast, d["size"]).correlation if len(d) > 30 else np.nan)
print(f"I AMPLITUDE: Spearman(|contrast|, episode slip) pooled {ri:+.3f}; per-station Spearman(contrast, slip) median {rs2.median():+.2f} (sign-aware; a response scaling with slip would show |median| >> 0)")
# J local vs regional timing, L tremor vs slip
jl = []
for tag, s in sta_series.items():
    y = d30(s["y"]); sxy = s["xy"]; d = np.sqrt(((XY - sxy) ** 2).sum(1)); near = np.isfinite(dep) & (dep >= 30) & (dep <= 45) & (d <= 60)
    if near.sum() < 3: continue
    local = np.nanmean(inc[:, near], 1); seg = (clat >= sll.sta_lat[tag] - 1.2) & (clat <= sll.sta_lat[tag] + 1.2) & np.isfinite(dep) & (dep >= 30) & (dep <= 45); regional = np.nanmean(inc[:, seg], 1)
    tl = np.sqrt(((TXY - sxy) ** 2).sum(1)) <= 60; tc = pd.Series(1.0, index=tr.d[tl]).groupby(level=0).sum().reindex(days, fill_value=0).rolling(30, min_periods=1).sum().values
    lags = np.arange(-90, 91, 5); xl = [corr(y, np.roll(local, -L)) for L in lags]; xr = [corr(y, np.roll(regional, -L)) for L in lags]
    if not np.isfinite(xl).any(): continue
    nul = [abs(corr(y, np.roll(local, rng.randint(365, G - 365)))) for _ in range(100)]
    jl.append(dict(tag=tag, r_local=corr(y, local), r_regional=corr(y, regional), r_tremor=corr(y, tc), lag_local=int(lags[int(np.nanargmax(np.abs(xl)))]), rmax_local=xl[int(np.nanargmax(np.abs(xl)))],
                   null95=np.nanpercentile(nul, 95), r_slip_vs_tremor=corr(local, tc)))
JL = pd.DataFrame(jl); JL.to_csv(f"{SD}/explore3_station_timing.csv", index=False)
sig = JL[JL.rmax_local.abs() > JL.null95]
print(f"\nJ/L LOCAL vs REGIONAL vs TREMOR ({len(JL)} stations): |r_local| > |r_regional| at {100*(JL.r_local.abs() > JL.r_regional.abs()).mean():.0f}%; stations above their local null95: {len(sig)} ({100*len(sig)/len(JL):.0f}%); "
      f"lag of max |xcorr| among those: median {sig.lag_local.median():+.0f} d, IQR {sig.lag_local.quantile(.25):+.0f}..{sig.lag_local.quantile(.75):+.0f}; sign of r_local among them: {100*(sig.r_local > 0).mean():.0f}% positive")
print(f"   tremor as predictor: |r_tremor| > |r_local slip| at {100*(JL.r_tremor.abs() > JL.r_local.abs()).mean():.0f}% of stations; median |r_tremor| {JL.r_tremor.abs().median():.3f} vs |r_local| {JL.r_local.abs().median():.3f}; corr(local slip, local tremor) median {JL.r_slip_vs_tremor.median():.2f}")
# F 3-5 s as primary observable (pooled N+C footprint contrast, both nulls)
rows = []
for f in glob.glob("data/daily_dvv_*_Z_3to5.csv"):
    tag = os.path.basename(f).split("daily_dvv_")[1].split("_Z_3to5")[0]
    x = pd.read_csv(f, usecols=["patch", "date", "dvv", "cc_max"]); x = x[x.cc_max > 0.7]; x["tag"] = tag; x["date"] = pd.to_datetime(x.date); x["v"] = -x.dvv * 100; rows.append(x.rename(columns={"patch": "fam"})[["tag", "fam", "date", "v"]])
f35 = pd.concat(rows).merge(fd[["tag", "fam", "ci"]].drop_duplicates(), on=["tag", "fam"]); c35 = f35.groupby(["ci", "date"]).v.mean().reset_index(); ok35 = c35.date.isin(days)
T35 = np.full((G, M), np.nan); T35[gi[c35.date[ok35]].values, c35.ci[ok35].values] = c35.v[ok35].values
T24 = np.full((G, M), np.nan); c24 = fd.groupby(["ci", "date"]).v.mean().reset_index(); ok24 = c24.date.isin(days); T24[gi[c24.date[ok24]].values, c24.ci[ok24].values] = c24.v[ok24].values
month = days.month.values
def deseason_grid(X, eps):
    mask = np.zeros(G, bool)
    for e in eps: mask[max(e["ion"] - 30, 0):e["ion"] + 121] = True
    Xd = X.copy()
    for m in range(1, 13): sel = (month == m) & ~mask; Xd[month == m] -= np.nanmean(X[sel], 0)[None]
    return Xd
PRIM = [e for e in EP if e["seg"] in ("N", "C")]
def fp_contrast(X, e, shift=0):
    i = wrap(e["ion"] + shift) if shift else e["ion"]; on = X[max(i + 15, 0):i + 61][:, e["fp"]]; pre = X[max(i - 90, 0):i - 29][:, e["fp"]]
    a, b = np.nanmean(on, 0), np.nanmean(pre, 0); m = np.isfinite(a) & np.isfinite(b); return np.nanmean(a[m] - b[m]) if m.sum() >= 3 else np.nan
for name, X in (("family 2-4 s", T24), ("family 3-5 s", T35)):
    Xd = deseason_grid(X, PRIM); s = np.array([fp_contrast(Xd, e) for e in PRIM]); K, N, A = int((s < 0).sum()), int(np.isfinite(s).sum()), np.nanmean(s)
    nA = [np.nanmean([fp_contrast(Xd, e, rng.randint(90, span - 90)) for e in PRIM]) for _ in range(500)]
    nB = [np.nanmean([fp_contrast(Xd, e, int(rng.choice([-1, 1]) * rng.randint(1, 15) * 365 + rng.randint(-15, 16))) for e in PRIM]) for _ in range(500)]
    print(f"F {name} pooled N+C footprint contrast: mean {A:+.3f}% (K={K}/{N} decrease) | Null A p {np.mean(np.abs(nA) >= abs(A)):.3f} (sd {np.std(nA):.3f}) | Null B p {np.mean(np.abs(nB) >= abs(A)):.3f}")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(17, 4.5))
ax[0].scatter(R.pre_rel, R.contrast, s=2, alpha=.15, color="#1d4ed8"); ax[0].set_xlim(-1, 1); ax[0].set_ylim(-1, 1); ax[0].axhline(0, color="k", lw=.5); ax[0].axvline(0, color="k", lw=.5)
ax[0].set_xlabel("PRE level − family all-time mean (%)"); ax[0].set_ylabel("ETS contrast ON−PRE (%)"); ax[0].set_title(f"D1 shrinkage-release test: Spearman {rs:+.2f}", fontsize=9)
if len(N_): ax[1].scatter(N_.lr, N_.contrast, s=4, alpha=.3, color="#7c3aed"); ax[1].axhline(0, color="k", lw=.5); ax[1].set_xlabel("log(n_det ON / PRE)"); ax[1].set_ylabel("contrast (%)"); ax[1].set_title(f"D2 detection count: Spearman {rn:+.2f}", fontsize=9)
ax[2].scatter(JL.r_regional, JL.r_local, s=12, color="#059669"); ax[2].plot([-.3, .3], [-.3, .3], "k--", lw=.5); ax[2].axhline(0, color="k", lw=.5); ax[2].axvline(0, color="k", lw=.5); ax[2].set_xlabel("r(station dv/v, segment slip)"); ax[2].set_ylabel("r(station dv/v, local slip <=60 km)"); ax[2].set_title("J local vs regional slip as predictor", fontsize=9)
fig.suptitle("Exploration round 3: mechanism tests for the station-specific ETS-locked response", fontsize=10); fig.savefig("figures/sse_explore3.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_explore3.png")
