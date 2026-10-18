#!/usr/bin/env python
"""Repeating (overlapping) footprints through time (user 2026-09-10: "in patches that overlap, look how the response may be
changing over time"). Sites = footprint clusters from sse_patch_pattern_explore.py (patch_pattern_features.csv: site id per patch).
For every site with >= 3 post-2011 substantive patches (Mw >= MWMIN): the sequence of responses (DUR-PRE patch median, POST-PRE,
60-d DUR-PRE) vs date; per-site Spearman trend with time and with cycle index, p from the per-patch shift-null draws; consecutive-
event differences; sign consistency. Aggregate: distribution of per-site trend rho vs null; pooled within-site trend (responses
centred per site, regressed on years since the site's first event), by terrane.
Outputs: res_catalog_sse/site_response_time.csv (per site), site_response_pairs.csv (consecutive pairs), figures/sse_site_response_time.png
Env: MWMIN (6.0).
"""
import os, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy.stats import spearmanr, theilslopes
SD = "fault_tomography/inversion/res_catalog_sse"; MWMIN = float(os.environ.get("MWMIN", "6.0")); rng = np.random.RandomState(3)
F = pd.read_csv(f"{SD}/patch_pattern_features.csv", parse_dates=["start"]); zn = np.load(f"{SD}/patch_prepost_null.npz", allow_pickle=True); pids = [int(x) for x in zn["pids"]]; pos = {p: i for i, p in enumerate(pids)}
PN = "interface δβ/β (corrected)"; NUL = {k: zn[f"{PN}|{k}"] for k in ("d_median_dur", "d_median_post")}; ND = NUL["d_median_dur"].shape[0]
TERR = ["Klamath/Franciscan (<43.2)", "Siletzia core (43.2-46.8)", "Olympic/Crescent (46.8-48.6)", "Vancouver Is. (>48.6)"]
Fs = F[F.Mw >= MWMIN].copy(); Fs["t_yr"] = (Fs.start - pd.Timestamp("2011-01-01")).dt.days / 365.25
sites = [s for s, n in Fs.site.value_counts().items() if n >= 3]; rows = []; pairs = []; centred = []; centred_null = []
for s in sites:
    e = Fs[Fs.site == s].sort_values("start"); ii = [pos[p] for p in e.id]; t = e.t_yr.values; y = e.d_median_dur.values; yp = e.d_median_post.values; ok = np.isfinite(y)
    def sp(a, b):
        m = np.isfinite(a) & np.isfinite(b); return spearmanr(a[m], b[m]).correlation if m.sum() >= 3 else np.nan
    rho = sp(t[ok], y[ok]) if ok.sum() >= 3 else np.nan; nr = np.array([sp(t[ok], NUL["d_median_dur"][d, ii][ok]) for d in range(ND)]); nr = nr[np.isfinite(nr)]
    slope = theilslopes(y[ok], t[ok])[0] if ok.sum() >= 3 else np.nan; nsl = np.array([theilslopes(NUL["d_median_dur"][d, ii][ok], t[ok])[0] for d in range(ND) if np.isfinite(NUL["d_median_dur"][d, ii][ok]).all() and ok.sum() >= 3])
    rows.append(dict(site=s, n=len(e), lat=e.lat.mean(), terrane=e.terrane.iloc[0], seg=e.seg.iloc[0], first=e.start.iloc[0].date(), last=e.start.iloc[-1].date(), mean_dur=np.nanmean(y), sd_dur=np.nanstd(y), frac_pos=np.nanmean(y > 0), mean_post=np.nanmean(yp),
                     trend_rho=rho, trend_p=float(np.mean(np.abs(nr) >= abs(rho))) if len(nr) > 20 else np.nan, theil_slope_pct_per_yr=slope, slope_null_sd=nsl.std() if len(nsl) > 20 else np.nan, slope_p=float(np.mean(np.abs(nsl) >= abs(slope))) if len(nsl) > 20 else np.nan,
                     first_resp=y[0], last_resp=y[-1], responses=";".join(f"{d.date()}:{v:+.2f}" for d, v in zip(e.start, y)), n_pre_overlap=int(e.pre_overlap.sum())))
    for (i1, r1), (i2, r2) in zip(e.iloc[:-1].iterrows(), e.iloc[1:].iterrows()):
        pairs.append(dict(site=s, terrane=r1.terrane, id_prev=r1.id, id_next=r2.id, dt_days=(r2.start - r1.start).days, resp_prev=r1.d_median_dur, resp_next=r2.d_median_dur, d_resp=r2.d_median_dur - r1.d_median_dur, Mw_prev=r1.Mw, Mw_next=r2.Mw))
    centred.append(np.column_stack([t - t[0], y - np.nanmean(y), [TERR.index(e.terrane.iloc[0])] * len(e)])); centred_null.append(np.stack([NUL["d_median_dur"][d, ii] - np.nanmean(NUL["d_median_dur"][d, ii]) for d in range(ND)]))
S = pd.DataFrame(rows).sort_values("lat"); S.to_csv(f"{SD}/site_response_time.csv", index=False); PPairs = pd.DataFrame(pairs); PPairs.to_csv(f"{SD}/site_response_pairs.csv", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 200)
print(f"Mw>={MWMIN:g}: {len(sites)} sites with >=3 events ({len(Fs)} patches, {Fs.site.isin(sites).sum()} on these sites)")
print(S[["site", "n", "lat", "terrane", "first", "last", "mean_dur", "sd_dur", "frac_pos", "trend_rho", "trend_p", "theil_slope_pct_per_yr", "slope_p", "n_pre_overlap"]].round(3).to_string(index=False))
print(f"\nsites with trend p<0.05: {(S.trend_p < 0.05).sum()} of {len(S)} (expected ~{0.05 * len(S):.1f}); with slope p<0.05: {(S.slope_p < 0.05).sum()}")
# pooled within-site trend (centred), overall and per terrane, with null
Cc = np.vstack(centred); Cn = np.concatenate(centred_null, axis=1)          # (ND, total)
def pooled(vals, tt): ok = np.isfinite(vals); return theilslopes(vals[ok], tt[ok])[0] if ok.sum() > 5 else np.nan
tt = Cc[:, 0]; obs = pooled(Cc[:, 1], tt); nul = np.array([pooled(Cn[d], tt) for d in range(ND)]); print(f"\npooled within-site trend (responses centred per site): {obs:+.3f} %/yr (null sd {np.nanstd(nul):.3f}, p {np.nanmean(np.abs(nul) >= abs(obs)):.2f})")
for k, terr in enumerate(TERR):
    m = Cc[:, 2] == k
    if m.sum() < 8: continue
    o = pooled(Cc[m, 1], tt[m]); n_ = np.array([pooled(Cn[d][m], tt[m]) for d in range(ND)]); print(f"   {terr:30s} n={m.sum():3d} slope {o:+.3f} %/yr (p {np.nanmean(np.abs(n_) >= abs(o)):.2f})")
# consecutive pairs: does the change between consecutive responses depend on the time gap or on the previous response?
ok = np.isfinite(PPairs.d_resp)
print(f"\nconsecutive pairs at a site: {ok.sum()} | corr(d_resp, dt) {spearmanr(PPairs.dt_days[ok], PPairs.d_resp[ok]).correlation:+.2f} | corr(resp_next, resp_prev) {spearmanr(PPairs.resp_prev[ok], PPairs.resp_next[ok]).correlation:+.2f} (negative ~ -0.5 is what pure noise + centring would give)")
# ---------- figure ----------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
COL = dict(zip(TERR, ["#ea580c", "#059669", "#1d4ed8", "#7c3aed"])); big = S[S.n >= 4].sort_values("lat", ascending=False); n = len(big); nc = 4; nr_ = int(np.ceil(n / nc)) + 1
fig = plt.figure(figsize=(18, 3.3 * nr_), constrained_layout=True); gs = fig.add_gridspec(nr_, nc)
for k, (_, r) in enumerate(big.iterrows()):
    a = fig.add_subplot(gs[k // nc, k % nc]); e = Fs[Fs.site == r.site].sort_values("start"); c = COL[r.terrane]
    a.errorbar(e.start, e.d_median_dur, yerr=1.96 * e.d_median_dur_nullsd, fmt="o-", color=c, ms=5, lw=1, ecolor=c, alpha=.85, label="DUR − PRE")
    a.plot(e.start, e.d_median_post, "s--", mfc="none", color=c, ms=5, lw=.8, alpha=.7, label="POST − PRE")
    for _, q in e[e.pre_overlap].iterrows(): a.plot(q.start, q.d_median_dur, "x", color="k", ms=9, mew=1.5)
    a.axhline(0, color="k", lw=.5); a.set_title(f"site {r.site} · {r.lat:.1f}N · {r.terrane.split(' (')[0]} · n {r.n}\ntrend rho {r.trend_rho:+.2f} (p {r.trend_p:.2f}), slope {r.theil_slope_pct_per_yr:+.2f} %/yr", fontsize=8); a.tick_params(labelsize=7); a.set_ylim(-2.6, 2.6)
    if k == 0: a.legend(fontsize=7)
a = fig.add_subplot(gs[-1, :2]); a.scatter(S.lat, S.trend_rho, s=25 * S.n, c=[COL[t] for t in S.terrane], edgecolor="k"); a.axhline(0, color="k", lw=.5); a.set_xlabel("site latitude"); a.set_ylabel("per-site trend rho (response vs time)"); a.set_title("per-site trend of the response with time (size = n events; black x in panels = PRE overlapped by another patch)", fontsize=9)
a = fig.add_subplot(gs[-1, 2:]); a.scatter(PPairs.dt_days / 365.25, PPairs.d_resp, c=[COL[t] for t in PPairs.terrane], s=18, alpha=.8); a.axhline(0, color="k", lw=.5); a.set_xlabel("years between consecutive events at a site"); a.set_ylabel("change of response, next − previous (%)"); a.set_title("consecutive events at the same site: response change vs recurrence time", fontsize=9)
fig.suptitle(f"Repeating footprints (sites) through time: dv/v response per event (Mw ≥ {MWMIN:g}; bars = per-patch null 95 %)", fontsize=11)
fig.savefig("figures/sse_site_response_time.png", dpi=130); print("wrote figures/sse_site_response_time.png, site_response_time.csv, site_response_pairs.csv")
