#!/usr/bin/env python
"""Does the interface dbeta/beta change in the same 125-km patch when the geodetic model says slow slip occurred?
Merlin-prescribed test (2026-09-02, bounded ruling; notes/2026-09-02_Notes.md):
  predictor  p_q(t) = thrust-positive 30-d slip increment  -[cum(t) - cum(t-30)]  (mm; stored sign has ETS NEGATIVE),
             projected on the 10 DEEP 125-km patches (transmission_T.npz patterns), finite-cell renormalised.
  response   y_q(t) = <pat_q, MODEL(t)> / <pat_q, pat_q> over finite cells (frame valid if >= 60 % pattern weight
             finite), 2011-01-01..2026-07-20, DESEASONALISED by monthly climatology (raw kept as sensitivity).
  episodes   local maxima of p_q >= 3 x RMS(negative side of p_q), >= 90-d apart, t* in [2011-04-01, 2026-04-01];
             dropped if PRE/ON overlap +-35 d of an era_table boundary of a station giving >= 10 % of the patch rows.
  contrast   s = mean y_q over ON=[t*-10, t*+35] - mean over PRE=[t*-90, t*-30]  (>= 60 % valid frames each).
  statistic  P1 = pooled count K of s < 0 (velocity decrease), TWO-SIDED vs a rigid circular-shift null of the
             episode times (one delta per along-strike group: south q10-12 / central q13-15 / north q16-19; 4000
             draws); P2 = sum (2K_q - N_q)^2 / N_q (sign-agnostic); binomial and episode-scramble nulls as checks;
             leave-one-patch-out; per-patch Sidak.
  controls   same contrast on cc and rows/cell-day (source-side / roster), SITES of footprint stations (near-station
             response, separate line), day-of-year concentration of episodes (seasonal phase lock), raw vs deseasoned.
  readout    DETECTION p<0.01 + robustness; BOUND p>=0.05 -> 95 % upper bound on |contrast| (estimator units);
             else AMBIGUOUS.  Strain invariants dropped (linear in the same time functions).
Outputs: res_catalog_sse/episode_test.npz, episode_table.csv, figures/sse_episode_composite.png; prints the readout.
"""
import numpy as np, pandas as pd, os
RD = "fault_tomography/inversion/res_catalog"; SD = "fault_tomography/inversion/res_catalog_sse"
rng = np.random.RandomState(2026)
NDRAW = 4000
T0, T1 = pd.Timestamp("2011-01-01"), pd.Timestamp("2026-07-20")
E0, E1 = pd.Timestamp("2011-04-01"), pd.Timestamp("2026-04-01")
PRE, ON, POST = (-90, -30), (-10, 35), (60, 120)

# ---------- load ----------
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True)
mdays = pd.to_datetime(zp["days"]); MODEL = np.asarray(zp["MODEL"], float); cells = zp["cell"]; M = MODEL.shape[1]
SITES = np.asarray(zp["SITES"], float); sta_tags = list(zp["sta_tags"])
zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True)
PAT = np.asarray(zt["patterns"], float); band = zt["patch_band"]; plat = np.asarray(zt["patch_lat"], float)
zs = np.load(f"{SD}/slip_at_cells.npz", allow_pickle=True)
assert list(zs["cell"]) == list(cells)
grid = pd.date_range(T0, T1, freq="D"); G = len(grid)
cum = pd.DataFrame(zs["slip_cum"].T, index=pd.to_datetime(zs["t_cum"])).groupby(level=0).mean()
cum = cum.reindex(pd.date_range(cum.index.min(), cum.index.max(), freq="D")).interpolate(limit_area="inside")
inc = -(cum - cum.shift(30)) * 1000.0                                        # thrust-positive 30-d increment, mm
P_cell = inc.reindex(grid).values                                            # (G, M)
Y_cell = pd.DataFrame(MODEL, index=mdays).reindex(grid).values               # (G, M), NaN on missing days

# ---------- gate (Merlin step 1): whole-mesh northern zone vs monthly tremor ----------
cum_all = pd.DataFrame(zs["slip_cum_all"].T, index=pd.to_datetime(zs["t_cum"])).groupby(level=0).mean()
zone = (zs["patch_lat"] >= 47.3) & (zs["patch_lat"] <= 48.8) & (-zs["patch_h"] >= 28) & (-zs["patch_h"] <= 45)
inc_all = -(cum_all - cum_all.shift(30)).iloc[:, zone] * 1000.0
zmax = inc_all.resample("ME").max().mean(axis=1)
trem = pd.read_csv(f"{RD}/ets_north.csv", index_col=0, parse_dates=True).tremor
j = zmax.index.intersection(trem.index)
from scipy.stats import spearmanr
rho = spearmanr(trem.loc[j], zmax.loc[j]).correlation
neg_rms = np.sqrt(np.mean(np.minimum(inc_all.values, 0) ** 2))
top5 = trem.loc[j].nlargest(5)
gate_ok = rho >= 0.25 and all(zmax.loc[t] >= 3 * neg_rms for t in top5.index)
print(f"GATE: Spearman(monthly tremor, zone max thrust-positive 30-d increment) = {rho:+.2f} (bar >= +0.25); "
      f"negative-side RMS {neg_rms:.1f} mm; five largest tremor months: "
      + ", ".join(f"{t.strftime('%Y-%m')}: {zmax.loc[t]:.1f} mm" for t in top5.index) + f" -> {'PASS' if gate_ok else 'FAIL'}")
if not gate_ok:
    raise SystemExit("gate failed: predictor sign/definition problem — stop (Merlin step 1 failure signature)")

# ---------- patch functionals ----------
deep = [q for q in np.argsort(plat) if band[q] == "deep"]                    # south -> north
groups = {"south": [], "central": [], "north": []}
for k, q in enumerate(deep):
    groups["south" if k < 3 else "central" if k < 6 else "north"].append(q)
def functional(X, q, minw=0.6):
    w = PAT[q]; fin = np.isfinite(X); num = np.nansum(np.where(fin, X, 0) * w[None, :], 1)
    den = (fin * (w ** 2)[None, :]).sum(1); wfin = (fin * w[None, :]).sum(1) / w.sum()
    out = num / np.where(den > 0, den, np.nan); out[wfin < minw] = np.nan
    return out
Y = {q: functional(Y_cell, q) for q in deep}
P = {q: functional(P_cell, q) for q in deep}
def deseason(y):
    s = pd.Series(y, index=grid); clim = s.groupby(s.index.month).transform("mean")
    return (s - clim).values
Yd = {q: deseason(Y[q]) for q in deep}

# ---------- episodes ----------
def episodes(p):
    neg = np.sqrt(np.nanmean(np.minimum(np.nan_to_num(p), 0) ** 2)); thr = 3 * neg
    pp = np.nan_to_num(p, nan=-np.inf)
    cand = [i for i in range(1, G - 1) if pp[i] >= thr and pp[i] >= pp[i - 1] and pp[i] >= pp[i + 1]
            and E0 <= grid[i] <= E1]
    cand = sorted(cand, key=lambda i: -pp[i]); keep = []
    for i in cand:
        if all(abs(i - k) >= 90 for k in keep): keep.append(i)
    return sorted(keep), neg
# era exclusion: stations giving >= 10 % of a patch's rows
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "cell", "date", "cc"])
cpos = {c: i for i, c in enumerate(cells)}; fd["ci"] = fd.cell.map(cpos)
et = pd.read_csv(f"{RD}/era_table.csv")
BND = {r.tag: [pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan"] for _, r in et.iterrows()}
rows_ct = fd.groupby(["tag", "ci"]).size().reset_index(name="n")
def patch_boundaries(q):
    w = PAT[q][rows_ct.ci.values]; share = (rows_ct.n.values * w); tot = share.sum()
    st = pd.Series(share).groupby(rows_ct.tag.values).sum() / max(tot, 1e-9)
    return [b for t in st[st >= 0.10].index for b in BND.get(t, [])], list(st[st >= 0.10].index)
# per-cell-day cc and rows (controls) and SITES footprint
cd = fd.groupby(["ci", "date"]).agg(cc=("cc", "mean"), n=("cc", "size")).reset_index()
CC_cell = np.full((G, M), np.nan); N_cell = np.full((G, M), np.nan)
gi = pd.Series(np.arange(G), index=grid); ok = cd.date.isin(grid)
CC_cell[gi[cd.date[ok]].values, cd.ci[ok].values] = cd.cc[ok].values
N_cell[gi[cd.date[ok]].values, cd.ci[ok].values] = cd.n[ok].values
CC = {q: functional(CC_cell, q) for q in deep}; NR = {q: functional(N_cell, q) for q in deep}
S_full = pd.DataFrame(SITES, index=mdays, columns=sta_tags).reindex(grid)
def sites_series(q):
    _, stas = patch_boundaries(q); stas = [s for s in stas if s in S_full.columns]
    return S_full[stas].mean(axis=1).values if stas else np.full(G, np.nan), stas

# ---------- window means via cumsums ----------
def cs(x):
    v = np.nan_to_num(x); f = np.isfinite(x).astype(float)
    return np.concatenate([[0], np.cumsum(v)]), np.concatenate([[0], np.cumsum(f)])
def wmean(c, i0, i1, minfrac=0.6):
    """mean over grid index [i0, i1] inclusive, None if < minfrac valid; i0/i1 arrays"""
    i0 = np.clip(i0, 0, G - 1); i1 = np.clip(i1, 0, G - 1)
    n = c[1][i1 + 1] - c[1][i0]; s = c[0][i1 + 1] - c[0][i0]; L = i1 - i0 + 1
    return np.where(n >= minfrac * L, s / np.maximum(n, 1), np.nan)
def contrasts(c, idx, win_on=ON, win_pre=PRE):
    idx = np.asarray(idx)
    return wmean(c, idx + win_on[0], idx + win_on[1]) - wmean(c, idx + win_pre[0], idx + win_pre[1])

EP, tab = {}, []
for q in deep:
    idx, neg = episodes(P[q]); bnds, stas = patch_boundaries(q)
    keep = []
    for i in idx:
        t = grid[i]
        near = any((t + pd.Timedelta(days=PRE[0]) - pd.Timedelta(days=35) <= b <= t + pd.Timedelta(days=ON[1]) + pd.Timedelta(days=35)) for b in bnds)
        if not near: keep.append(i)
    EP[q] = np.array(keep, int)
    tab.append(dict(q=q, lat=plat[q], n_cand=len(idx), n_kept=len(keep), neg_rms_mm=neg, thr_mm=3 * neg,
                    top_stations=",".join(stas[:6]), n_bnd=len(bnds)))
ept = pd.DataFrame(tab); print("\nEPISODES per deep patch:\n" + ept.round(2).to_string(index=False))

# ---------- observed statistics ----------
def stats_for(series_by_q, ep_by_q):
    K = 0; N = 0; per = {}
    for q in deep:
        idx = ep_by_q[q]
        if len(idx) == 0: per[q] = (0, 0, np.nan, np.array([])); continue
        s = contrasts(cs(series_by_q[q]), idx); s = s[np.isfinite(s)]
        per[q] = (int((s < 0).sum()), len(s), float(np.median(s)) if len(s) else np.nan, s)
        K += per[q][0]; N += per[q][1]
    return K, N, per
def P2_of(per):
    return sum((2 * k - n) ** 2 / n for k, n, _, _ in per.values() if n >= 6)
K_obs, N_obs, per_obs = stats_for(Yd, EP); P2_obs = P2_of(per_obs)
K_raw, N_raw, per_raw = stats_for(Y, EP)

# ---------- nulls ----------
span = (E1 - E0).days; e0 = int(gi[E0])
def shift_eps(ep_by_q, deltas):
    out = {}
    for g, qs in groups.items():
        for q in qs:
            out[q] = e0 + ((ep_by_q[q] - e0 + deltas[g]) % span) if len(ep_by_q[q]) else ep_by_q[q]
    return out
def null_shift(series_by_q, ep_by_q, ndraw=NDRAW):
    Ks, P2s = np.zeros(ndraw), np.zeros(ndraw)
    for d in range(ndraw):
        deltas = {g: rng.randint(90, span - 90) for g in groups}
        K, N, per = stats_for(series_by_q, shift_eps(ep_by_q, deltas)); Ks[d] = K - N / 2; P2s[d] = P2_of(per)
    return Ks, P2s
def null_scramble(series_by_q, ep_by_q, ndraw=NDRAW):
    Ks = np.zeros(ndraw)
    lo, hi = int(gi[E0]), int(gi[E1])
    for d in range(ndraw):
        ep2 = {}
        for q in deep:
            n = len(ep_by_q[q]); pts = []
            while len(pts) < n:
                c = rng.randint(lo, hi)
                if all(abs(c - p) >= 90 for p in pts): pts.append(c)
            ep2[q] = np.array(sorted(pts), int)
        K, N, _ = stats_for(series_by_q, ep2); Ks[d] = K - N / 2
    return Ks
def two_sided_p(null_dev, obs_dev):
    return float((np.abs(null_dev) >= abs(obs_dev) - 1e-12).mean())
Ks, P2s = null_shift(Yd, EP)
p1 = two_sided_p(Ks, K_obs - N_obs / 2); p2 = float((P2s >= P2_obs - 1e-12).mean())
from scipy.stats import binomtest
p_binom = binomtest(K_obs, N_obs, 0.5).pvalue if N_obs else np.nan
Ks_scr = null_scramble(Yd, EP); p1_scr = two_sided_p(Ks_scr, K_obs - N_obs / 2)
Ks_raw, _ = null_shift(Y, EP, 1000); p1_raw = two_sided_p(Ks_raw, K_raw - N_raw / 2)
# per-patch Sidak, leave-one-patch-out
perp = {}
for q in deep:
    k, n, a, s = per_obs[q]
    if n < 6: perp[q] = np.nan; continue
    dev = np.array([stats_for({q: Yd[q]}, {q: shift_eps({q: EP[q]}, {g: rng.randint(90, span - 90) for g in groups})[q]})[0] - n / 2
                    for _ in range(1000)])
    perp[q] = two_sided_p(dev, k - n / 2)
sidak = {q: (1 - (1 - p) ** 10 if np.isfinite(p) else np.nan) for q, p in perp.items()}
loo = {}
for q in deep:
    EPl = {k: (v if k != q else np.array([], int)) for k, v in EP.items()}
    Kl, Nl, _ = stats_for(Yd, EPl); Ksl, _ = null_shift(Yd, EPl, 1000); loo[q] = two_sided_p(Ksl, Kl - Nl / 2)
# group consistency
gsign = {g: sum(2 * per_obs[q][0] - per_obs[q][1] for q in qs) for g, qs in groups.items()}
# controls
def control_p(series_by_q):
    K, N, per = stats_for(series_by_q, EP); Kn, _ = null_shift(series_by_q, EP, 1000)
    return K, N, two_sided_p(Kn, K - N / 2)
ctl_cc = control_p(CC); ctl_n = control_p(NR)
SI = {}; foot = {}
for q in deep: SI[q], foot[q] = sites_series(q)
ctl_sites = control_p(SI)
doy = {q: (lambda a: float(np.abs(np.mean(np.exp(1j * 2 * np.pi * a / 365.25)))) if len(a) else np.nan)(grid[EP[q]].dayofyear.values) for q in deep}
# amplitude bound (pooled precision-weighted contrast and its shift-null 95 % interval)
A_obs = np.concatenate([per_obs[q][3] for q in deep if per_obs[q][1]])
amp_null = []
for d in range(1000):
    K, N, per = stats_for(Yd, shift_eps(EP, {g: rng.randint(90, span - 90) for g in groups}))
    amp_null.append(np.nanmean(np.concatenate([per[q][3] for q in deep if per[q][1]])))
amp_null = np.array(amp_null); amp_obs = float(np.nanmean(A_obs))
ub95 = float(np.nanpercentile(np.abs(amp_null - amp_null.mean()), 95) + abs(amp_obs))

# ---------- readout ----------
print(f"\nP1 pooled sign count: K = {K_obs}/{N_obs} contrasts with velocity DECREASE (deseasonalised) | shift-null p = {p1:.4f} "
      f"| scramble-null p = {p1_scr:.4f} | binomial p = {p_binom:.4f} | raw (undeseasonalised) K = {K_raw}/{N_raw}, p = {p1_raw:.4f}")
print(f"P2 sign-agnostic: {P2_obs:.1f}, shift-null p = {p2:.4f}")
print("group sums of (2K-N): " + ", ".join(f"{g} {v:+d}" for g, v in gsign.items()))
print("per patch: " + " | ".join(f"{plat[q]:.1f}N K={per_obs[q][0]}/{per_obs[q][1]} med={per_obs[q][2]:+.2f}% p={perp[q]:.2f} sidak={sidak[q]:.2f} doy={doy[q]:.2f}" for q in deep))
print("leave-one-patch-out pooled p: " + ", ".join(f"{plat[q]:.1f}N {loo[q]:.3f}" for q in deep))
print(f"CONTROLS (same contrast at true episodes vs shift null): cc K={ctl_cc[0]}/{ctl_cc[1]} p={ctl_cc[2]:.3f} | rows/cell-day K={ctl_n[0]}/{ctl_n[1]} p={ctl_n[2]:.3f} "
      f"| SITES (near-station, footprint stations) K={ctl_sites[0]}/{ctl_sites[1]} p={ctl_sites[2]:.3f}")
print(f"AMPLITUDE: pooled mean contrast {amp_obs:+.3f} % (estimator units, uncalibrated); shift-null sd {amp_null.std():.3f} %; 95% upper bound on |contrast| {ub95:.3f} %")
verdict = "AMBIGUOUS"
robust = (p1 < 0.01 and sum(np.sign(v) == np.sign(K_obs - N_obs / 2) for v in gsign.values() if v != 0) >= 2
          and all(v < 0.05 for v in loo.values()) and ctl_cc[2] >= 0.05 and ctl_n[2] >= 0.05)
if robust: verdict = "DETECTION"
elif p1 >= 0.05 and abs(amp_obs) <= np.nanpercentile(np.abs(amp_null - amp_null.mean()), 95): verdict = "BOUND"
print(f"\nVERDICT (pre-registered): {verdict}")
if verdict == "DETECTION":
    Kp, Np, perpost = stats_for(Yd, EP) if False else (None, None, None)
    post = np.concatenate([contrasts(cs(Yd[q]), EP[q], POST, PRE) for q in deep if len(EP[q])])
    print(f"  POST-PRE mean {np.nanmean(post):+.3f}% vs ON-PRE {amp_obs:+.3f}% -> {'STEP' if abs(np.nanmean(post)) > 0.5 * abs(amp_obs) else 'PULSE'}")
    lags = list(range(-90, 91, 15)); lagK = []
    for L in lags:
        K, N, _ = stats_for(Yd, {q: EP[q] + L for q in deep}); lagK.append(K - N / 2)
    print("  lag scan (K-N/2): " + " ".join(f"{L:+d}d:{v:+.0f}" for L, v in zip(lags, lagK)))

# ---------- composite figure ----------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 5, figsize=(17, 6), sharex=True)
lag = np.arange(-150, 181)
for k, q in enumerate(deep):
    a = ax[k // 5, k % 5]; idx = EP[q]
    if len(idx):
        Ym = np.array([[Yd[q][i + l] if 0 <= i + l < G else np.nan for l in lag] for i in idx])
        Pm = np.array([[P[q][i + l] if 0 <= i + l < G else np.nan for l in lag] for i in idx])
        a.plot(lag, np.nanmean(Ym, 0), color="#1d4ed8", lw=1.6, label="δβ/β composite (deseason.)")
        a.fill_between(lag, np.nanmean(Ym, 0) - np.nanstd(Ym, 0) / np.sqrt(len(idx)), np.nanmean(Ym, 0) + np.nanstd(Ym, 0) / np.sqrt(len(idx)), color="#bfdbfe", alpha=.7)
        a2 = a.twinx(); a2.plot(lag, np.nanmean(Pm, 0), color="#f59e0b", lw=1.2); a2.set_ylim(0, None); a2.tick_params(labelsize=7)
        a.axvspan(PRE[0], PRE[1], color="#94a3b8", alpha=.15); a.axvspan(ON[0], ON[1], color="#f59e0b", alpha=.12)
    a.axhline(0, color="k", lw=.5); a.set_title(f"deep {plat[q]:.1f}N  N={per_obs[q][1]}  K={per_obs[q][0]}  p={perp[q]:.2f}", fontsize=9)
    a.set_ylim(-0.6, 0.6); a.tick_params(labelsize=7)
    if k % 5 == 0: a.set_ylabel("δβ/β (%)", fontsize=8)
    if k // 5 == 1: a.set_xlabel("days from model slip peak", fontsize=8)
fig.suptitle(f"Episode composites: interface δβ/β (2-4 s joint inversion, deep 125-km patches, deseasonalised) around modeled slow-slip peaks (orange = thrust 30-d slip increment)  |  pooled K={K_obs}/{N_obs}, shift-null p={p1:.3f} -> {verdict}", fontsize=9)
fig.savefig("figures/sse_episode_composite.png", dpi=110, bbox_inches="tight")
rows = []
for q in deep:
    for i, s in zip(EP[q], contrasts(cs(Yd[q]), EP[q]) if len(EP[q]) else []):
        rows.append(dict(patch_lat=plat[q], t_peak=str(grid[i].date()), slip_inc_mm=P[q][i], contrast_pct=s))
pd.DataFrame(rows).to_csv(f"{SD}/episode_table.csv", index=False)
np.savez_compressed(f"{SD}/episode_test.npz", K=K_obs, N=N_obs, p1=p1, p1_scr=p1_scr, p_binom=p_binom, P2=P2_obs, p2=p2,
                    K_raw=K_raw, p1_raw=p1_raw, null_K=Ks, per_patch_p=np.array([perp[q] for q in deep]), loo=np.array([loo[q] for q in deep]),
                    ctl=np.array([ctl_cc[2], ctl_n[2], ctl_sites[2]]), amp_obs=amp_obs, amp_null_sd=amp_null.std(), ub95=ub95, verdict=verdict,
                    deep=np.array(deep), plat=plat[deep], gate_rho=rho)
print("wrote figures/sse_episode_composite.png, res_catalog_sse/episode_table.csv, episode_test.npz")
