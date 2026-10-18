#!/usr/bin/env python
"""rf07_samemonth (written Wed 2026-10-14 ~2:55 pm PT; PLAN Amendment 2 + acf_pgc PLAN §10 C2, written before running).
Re-scores the existing per-RF CC_SV tables (tables/<tag>/<STA>_events.csv, column cc_sv_clean) with
 (1) a SAME-MONTH scrambled control: real ETS midpoints moved to other years, same calendar date +-U{-15..15} d,
     rejected if < 60 d from any real midpoint, outside the station's RF span, or < 30 d from a date already drawn;
     RF values stay on their own dates. 20,000 draws for bins / H, first 5,000 also for the recovery curve.
 (2) a seasonal regression CC = a + annual + semiannual, fitted on RFs > 60 d from any real midpoint (and on all RFs);
     season-predicted Delta and season-removed Delta.
 (3) pre-ETS window H = mean[-21,0) - mean[-60,-25) d (lead test C2).
Delta = mean(post [0,W)) - mean(pre [-W,0)), nearest-midpoint assignment, as rf04. Recovery as rf05 (15-d window).
Usage: python -P rf07_samemonth.py --runs full:CN.PGC,repro:CN.PGC,... [--ndraw 20000] [--ncurve 5000] [--out part1]
Out: tables/samemonth/<out>_{bins,curve,season,lead}.csv ; figures/samemonth/recovery_<tag>_<STA>.png
"""
import argparse, json, os, time
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
W = f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin"
BINS = [15, 30, 45, 60, 90, 120, 150]
CENTRES = np.arange(-60, 121, 3)
YR = 365.25
ETS_OF = {"repro": "SVI_ALBH", "full": "SVI_ALBH", "oos2019": "SVI_ALBH", "ext_svi": "SVI_ALBH", "ext_wa": "WA_SEAT"}


def days(t):
    return (np.asarray(t).astype("datetime64[ns]") - np.datetime64("1990-01-01", "ns")) / np.timedelta64(1, "D")


def nearest_rel(td, md):
    """td (n,), md (..., k) -> rel (..., n) days from the nearest midpoint."""
    if md.ndim == 1:
        r = td[:, None] - md[None, :]
        return r[np.arange(len(td)), np.argmin(np.abs(r), 1)]
    r = td[None, :, None] - md[:, None, :]          # (nb, n, k)
    j = np.argmin(np.abs(np.where(np.isfinite(r), r, np.inf)), 2)
    return np.take_along_axis(r, j[..., None], 2)[..., 0]


def bin_delta(x, rel, Wd):
    pre = (rel >= -Wd) & (rel < 0); post = (rel >= 0) & (rel < Wd)
    npre, npost = pre.sum(-1), post.sum(-1)
    d = (post * x).sum(-1) / np.maximum(npost, 1) - (pre * x).sum(-1) / np.maximum(npre, 1)
    return np.where((npre >= 3) & (npost >= 3), d, np.nan)


def win_mean_diff(x, rel, a0, a1, b0, b1):
    A = (rel >= a0) & (rel < a1); B = (rel >= b0) & (rel < b1)
    na, nb = A.sum(-1), B.sum(-1)
    d = (A * x).sum(-1) / np.maximum(na, 1) - (B * x).sum(-1) / np.maximum(nb, 1)
    return np.where((na >= 3) & (nb >= 3), d, np.nan)


def curve(x, rel, mu, w=15):
    out = []
    for c in CENTRES:
        m = (rel >= c - w / 2) & (rel < c + w / 2)
        n = m.sum(-1)
        out.append(np.where(n >= 3, (m * x).sum(-1) / np.maximum(n, 1) - mu, np.nan))
    return np.stack(out, -1)


def draw_mids(rng, real_all_d, scr_src_d, years_d, tmin, tmax):
    """One same-month draw. real_all_d: all real midpoints (days); scr_src_d: midpoints to move; years_d: Jan-1 (days) of
    years with RFs. Returns array (len(scr_src_d),) with NaN for dropped episodes."""
    out = np.full(len(scr_src_d), np.nan); chosen = []
    for i, m in enumerate(scr_src_d):
        y0 = np.datetime64("1990-01-01") + np.timedelta64(int(m), "D")
        doy = int(m - days(y0.astype("datetime64[Y]").astype("datetime64[D]"))[()])
        for _ in range(200):
            c = years_d[rng.integers(len(years_d))] + doy + rng.integers(-15, 16)
            if c < tmin or c > tmax:
                continue
            if np.min(np.abs(real_all_d - c)) < 60:
                continue
            if chosen and np.min(np.abs(np.array(chosen) - c)) < 30:
                continue
            out[i] = c; chosen.append(c); break
    return out


def season_design(td):
    w = 2 * np.pi * td / YR
    return np.c_[np.ones_like(td), np.cos(w), np.sin(w), np.cos(2 * w), np.sin(2 * w)]


def pct_low(B, o):
    B = B[np.isfinite(B)]
    return float(np.mean(B <= o)) if len(B) and np.isfinite(o) else np.nan


def pct_high(B, o):
    B = B[np.isfinite(B)]
    return float(np.mean(B >= o)) if len(B) and np.isfinite(o) else np.nan


def rsd(B):
    B = B[np.isfinite(B)]; m = np.median(B)
    return 1.4826 * np.median(np.abs(B - m))


def run_one(tag, key, ndraw, ncurve, seed):
    t0 = time.time()
    R = json.load(open(f"{W}/tables/{tag}/{key}_result.json"))
    ets = ETS_OF[tag]
    M = pd.read_csv(f"{W}/tables/{tag}/{key}_events.csv")
    M["t"] = pd.to_datetime(M["t"]).astype("datetime64[ns]")
    E = pd.read_csv(f"{W}/tables/ets_mid_{ets}.csv")
    allmid = pd.to_datetime(E["mid"]).values.astype("datetime64[ns]")
    runmid = np.sort(allmid[(allmid >= np.datetime64(R["t0"])) & (allmid < np.datetime64(R["t1"]))])
    td_all = days(M.t.values); tmin, tmax = td_all.min(), td_all.max()
    x_all = M.cc_sv_clean.values.astype(float)
    ok = np.isfinite(x_all); x = x_all[ok]; td = td_all[ok]
    real_all_d = np.sort(days(allmid)); run_d = days(runmid)
    rel = nearest_rel(td, run_d)
    # obs Delta per bin (must match rf04)
    obs = {Wd: float(bin_delta(x, rel, Wd)) for Wd in BINS}
    stored = {r["W"]: r for r in R["quantities"]["cc_sv_clean"]}
    chk = max(abs(obs[Wd] - stored[Wd]["d"]) for Wd in BINS if np.isfinite(stored[Wd]["d"]))
    H = float(win_mean_diff(x, rel, -21, 0, -60, -25))
    mu = x.mean()
    obs_curve = curve(x, rel, mu)
    # seasonal regression
    rel_allreal = nearest_rel(td, real_all_d)
    quiet = np.abs(rel_allreal) > 60
    G = season_design(td)
    cq = np.linalg.lstsq(G[quiet], x[quiet], rcond=None)[0]
    ca = np.linalg.lstsq(G, x, rcond=None)[0]
    rss1 = np.sum((x[quiet] - G[quiet] @ cq) ** 2); rss0 = np.sum((x[quiet] - x[quiet].mean()) ** 2)
    nq = quiet.sum(); F = ((rss0 - rss1) / 4) / (rss1 / (nq - 5)); pF = float(stats.f.sf(F, 4, nq - 5))
    model = G @ cq; resid = x - model
    seas_pred = {Wd: float(bin_delta(model, rel, Wd)) for Wd in (15, 30, 150)}
    seas_pred_all = {Wd: float(bin_delta(G @ ca, rel, Wd)) for Wd in (15, 30, 150)}
    obs_res = {Wd: float(bin_delta(resid, rel, Wd)) for Wd in (15, 30, 150)}
    H_res = float(win_mean_diff(resid, rel, -21, 0, -60, -25))
    # monthly means (descriptive)
    mon = pd.Series(x, index=M.t.values[ok]).groupby(pd.DatetimeIndex(M.t.values[ok]).month).mean()
    # draws
    rng = np.random.default_rng(seed)
    src = run_d[(run_d >= tmin) & (run_d <= tmax)]
    years = np.unique(M.t.values.astype("datetime64[Y]")).astype("datetime64[D]")
    years_d = days(years)
    SM = {Wd: np.empty(ndraw) for Wd in BINS}; SMres = {Wd: np.empty(ndraw) for Wd in (15, 30, 150)}
    SH = np.empty(ndraw); SHres = np.empty(ndraw); SC = np.full((ncurve, len(CENTRES)), np.nan)
    ndrop = 0; batch = 250
    for s in range(0, ndraw, batch):
        nb = min(batch, ndraw - s)
        MD = np.stack([draw_mids(rng, real_all_d, src, years_d, tmin, tmax) for _ in range(nb)])
        ndrop += int(np.isnan(MD).sum())
        r = nearest_rel(td, MD)                                    # (nb, n)
        for Wd in BINS:
            SM[Wd][s:s + nb] = bin_delta(x[None, :], r, Wd)
        for Wd in (15, 30, 150):
            SMres[Wd][s:s + nb] = bin_delta(resid[None, :], r, Wd)
        SH[s:s + nb] = win_mean_diff(x[None, :], r, -21, 0, -60, -25)
        SHres[s:s + nb] = win_mean_diff(resid[None, :], r, -21, 0, -60, -25)
        if s < ncurve:
            k = min(nb, ncurve - s)
            SC[s:s + k] = curve(x[None, :], r[:k], mu)
    # old (value-shuffle) control for H, residual Delta (bins come from the stored rf04 results)
    rng2 = np.random.default_rng(seed + 1); nold = 5000
    OH = np.empty(nold); OR15 = np.empty(nold)
    for i in range(0, nold, 500):
        P = np.argsort(rng2.random((500, len(x))), 1)
        OH[i:i + 500] = win_mean_diff(x[P], rel[None, :], -21, 0, -60, -25)
        OR15[i:i + 500] = bin_delta(resid[P], rel[None, :], 15)
    rows = []
    for Wd in BINS:
        B = SM[Wd]
        rows.append(dict(tag=tag, station=key, quantity="dCC_SV", W=Wd, obs=obs[Wd], stored_rf04=stored[Wd]["d"],
                         n_pre=stored[Wd]["n_pre"], n_post=stored[Wd]["n_post"],
                         sm_median=float(np.nanmedian(B)), sm_rsd=rsd(B), sm_lo95=float(np.nanpercentile(B, 2.5)),
                         sm_hi95=float(np.nanpercentile(B, 97.5)), sm_pct_low=pct_low(B, obs[Wd]),
                         sm_pct_high=pct_high(B, obs[Wd]),
                         old_pct_low=stored[Wd]["pct_below"], old_pct_high=stored[Wd]["pct_above"],
                         welch_p=stored[Wd]["p"]))
    for Wd in (15, 30, 150):
        B = SMres[Wd]
        rows.append(dict(tag=tag, station=key, quantity="dCC_SV_season_removed", W=Wd, obs=obs_res[Wd],
                         sm_median=float(np.nanmedian(B)), sm_rsd=rsd(B), sm_lo95=float(np.nanpercentile(B, 2.5)),
                         sm_hi95=float(np.nanpercentile(B, 97.5)), sm_pct_low=pct_low(B, obs_res[Wd]),
                         sm_pct_high=pct_high(B, obs_res[Wd]),
                         old_pct_low=pct_low(OR15, obs_res[Wd]) if Wd == 15 else np.nan,
                         old_pct_high=pct_high(OR15, obs_res[Wd]) if Wd == 15 else np.nan))
    bins = pd.DataFrame(rows)
    lead = pd.DataFrame([dict(tag=tag, station=key, H=H, sm_median=float(np.nanmedian(SH)), sm_rsd=rsd(SH),
                              sm_pct_high=pct_high(SH, H), old_pct_high=pct_high(OH, H),
                              H_season_removed=H_res, sm_pct_high_res=pct_high(SHres, H_res),
                              n_m21_0=int(((rel >= -21) & (rel < 0)).sum()), n_m60_m25=int(((rel >= -60) & (rel < -25)).sum()))])
    amp = lambda c1, c2: float(np.hypot(c1, c2))
    seas = pd.DataFrame([dict(tag=tag, station=key, n_rf=len(x), n_quiet=int(nq), mean_cc=float(mu),
                              annual_amp=amp(*cq[1:3]), semiannual_amp=amp(*cq[3:5]),
                              annual_peak_doy=float((np.degrees(np.arctan2(cq[2], cq[1])) % 360) / 360 * YR),
                              F=float(F), F_p=pF, annual_amp_allfit=amp(*ca[1:3]), semiannual_amp_allfit=amp(*ca[3:5]),
                              season_pred_d15=seas_pred[15], season_pred_d30=seas_pred[30], season_pred_d150=seas_pred[150],
                              season_pred_d15_allfit=seas_pred_all[15], season_pred_d150_allfit=seas_pred_all[150],
                              obs_d15=obs[15], obs_d150=obs[150],
                              **{f"month_{m:02d}": float(mon.get(m, np.nan) - mu) for m in range(1, 13)})])
    lo, hi = np.nanpercentile(SC, 2.5, 0), np.nanpercentile(SC, 97.5, 0)
    fo = f"{W}/tables/{tag}/recovery_{key}.csv"
    oldc = pd.read_csv(fo) if os.path.exists(fo) else pd.DataFrame({c: np.full(len(CENTRES), np.nan) for c in
                                                                    ["cc_sv_clean_w15_lo95", "cc_sv_clean_w15_hi95", "cc_sv_clean_w15"]})
    cur = pd.DataFrame(dict(tag=tag, station=key, centre=CENTRES, anomaly=obs_curve, sm_median=np.nanmedian(SC, 0),
                            sm_lo95=lo, sm_hi95=hi, old_lo95=oldc["cc_sv_clean_w15_lo95"].values,
                            old_hi95=oldc["cc_sv_clean_w15_hi95"].values, rf05_anomaly=oldc["cc_sv_clean_w15"].values))
    # figure: recovery curve with both bands (CC series; no error bars)
    os.makedirs(f"{W}/figures/samemonth", exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw=dict(width_ratios=[2, 1]))
    a = ax[0]
    a.fill_between(CENTRES, cur.old_lo95, cur.old_hi95, color="0.88", label="old control 95 % (values shuffled)")
    a.plot(CENTRES, cur.sm_lo95, color="tab:purple", lw=1, label="same-month control 95 %")
    a.plot(CENTRES, cur.sm_hi95, color="tab:purple", lw=1)
    a.plot(CENTRES, cur.sm_median, color="tab:purple", lw=1, ls=":", label="same-month control median")
    a.plot(CENTRES, cur.anomaly, "k-o", ms=3, lw=1.5, label="observed (15-d window)")
    a.axvline(0, color="0.5", ls=":"); a.axhline(0, color="0.5", lw=0.5)
    a.set_xlabel("window centre, days from ETS midpoint"); a.set_ylabel("CC_SV anomaly vs mean of all RFs")
    a.legend(fontsize=7); a.set_title(f"{key} {tag}: recovery curve, same-month scrambled control", fontsize=9)
    b = ax[1]
    mm = np.arange(1, 13); mv = [seas[f"month_{m:02d}"].iloc[0] for m in mm]
    b.bar(mm, mv, color="0.6", label="monthly mean − overall")
    tt = np.linspace(0, YR, 200)
    b.plot(0.5 + tt / YR * 12, season_design(tt) @ cq - cq[0] - 0 * mu, color="tab:red", label="annual+semiannual fit (quiet RFs)")
    mids_m = pd.DatetimeIndex(runmid).month
    for m in mids_m:
        b.axvline(m, color="tab:blue", alpha=0.15, lw=3)
    b.set_xlabel("month (blue = ETS midpoints)"); b.set_title(f"season: annual {seas.annual_amp.iloc[0]:.3f}, semi {seas.semiannual_amp.iloc[0]:.3f}, F p {pF:.3g}", fontsize=8)
    b.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{W}/figures/samemonth/recovery_{tag}_{key}.png", dpi=110); plt.close(fig)
    print(f"{tag} {key}: n {len(x)} ets {len(run_d)} (scrambled {len(src)}), check |obs-rf04| {chk:.2e}, dropped "
          f"{ndrop / ndraw:.2f}/draw, d15 {obs[15]:+.3f} sm pct {pct_low(SM[15], obs[15]):.4f} (median {np.nanmedian(SM[15]):+.3f}),"
          f" H {H:+.3f} pct {pct_high(SH, H):.4f}, {time.time() - t0:.0f} s", flush=True)
    return bins, cur, seas, lead


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True); ap.add_argument("--ndraw", type=int, default=20000)
    ap.add_argument("--ncurve", type=int, default=5000); ap.add_argument("--out", default="part")
    ap.add_argument("--seed", type=int, default=20261014)
    a = ap.parse_args()
    od = f"{W}/tables/samemonth"; os.makedirs(od, exist_ok=True)
    acc = [[], [], [], []]
    for i, rk in enumerate(a.runs.split(",")):
        tag, key = rk.split(":")
        for j, df in enumerate(run_one(tag, key, a.ndraw, a.ncurve, a.seed + i)):
            acc[j].append(df)
        for j, nm in enumerate(["bins", "curve", "season", "lead"]):
            pd.concat(acc[j]).to_csv(f"{od}/{a.out}_{nm}.csv", index=False)


if __name__ == "__main__":
    main()
