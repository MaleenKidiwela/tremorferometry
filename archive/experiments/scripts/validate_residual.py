#!/usr/bin/env python
"""Decisive validation of the noise-corrected residual dv/v (Merlin step 4 + 3).
Residual = fwd_common_mode - beta*rev_common_mode (beta = regression, deseasoned).
ETS-epoch PERMUTATION test: does the residual move at real tremor-episode onsets beyond chance?
  statistic = mean residual co-onset (0..30 d) - pre (-60..-20 d), averaged over episodes & stations.
  null = 1000 draws of equal-count FAKE onsets avoiding +-60 d of real episodes.
beta-guard: re-fit beta EXCLUDING episode days; confirm the excursion survives (<30% change).
"""
import numpy as np, pandas as pd
rng = np.random.RandomState(7)
STATIONS = [("B926", 48.82, -124.13), ("B011", 48.65, -123.45)]

def cm(csv, fams):
    d = pd.read_csv(csv); d = d[d.patch.isin(fams) & (d.cc_max >= 0.6)]; d["date"] = pd.to_datetime(d.date)
    m = d.groupby("date").dvv.median(); n = d.groupby("date").patch.nunique()
    return (m[n >= 3].rolling(15, center=True, min_periods=5).median()) * 100

def deseason(s):
    s = s.dropna(); t = (s.index - s.index[0]).days.values.astype(float); yr = 365.25
    X = np.column_stack([np.ones_like(t), t, np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr),
                         np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
    b, *_ = np.linalg.lstsq(X, s.values, rcond=None)
    return pd.Series(s.values - X @ b, index=s.index)

# tremor-episode onsets near a station (>=120-day separation)
cat = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
cat["t"] = pd.to_datetime(cat.time, errors="coerce"); cat = cat.dropna(subset=["t"])
def onsets(la, lo):
    c = cat[cat.lat.between(la-0.7, la+0.7) & cat.lon.between(lo-0.7, lo+0.7)]
    daily = c.groupby(c.t.dt.floor("D")).size()
    idx = pd.date_range(daily.index.min(), daily.index.max())
    rate = daily.reindex(idx, fill_value=0).rolling(15, center=True, min_periods=5).mean()
    thr = rate.quantile(0.80)
    hot = rate > thr
    ons = []
    for d in rate.index[hot.values]:
        if not ons or (d - ons[-1]).days >= 120:
            ons.append(d)
    return pd.DatetimeIndex(ons), rate

def excursion(resid, ons):
    vals = []
    for o in ons:
        pre = resid[(resid.index >= o - pd.Timedelta("60D")) & (resid.index < o - pd.Timedelta("20D"))]
        co = resid[(resid.index >= o) & (resid.index < o + pd.Timedelta("30D"))]
        if len(pre) >= 10 and len(co) >= 10:
            vals.append(co.mean() - pre.mean())
    return np.array(vals)

# build residuals + real excursions
resids, real_ons, all_days = {}, {}, {}
for S, la, lo in STATIONS:
    fams = set(pd.read_csv(f"data/{S.lower()}_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    f = deseason(cm(f"data/daily_dvv_{S}p90f40_Z_2to4.csv", fams))
    r = deseason(cm(f"data/daily_dvv_{S}p90f40rev_Z_2to4.csv", fams))
    j = pd.concat([f.rename("f"), r.rename("r")], axis=1).dropna()
    beta = np.polyfit(j.r, j.f, 1)[0]
    resids[S] = (j.f - beta * j.r)
    real_ons[S], _ = onsets(la, lo); all_days[S] = resids[S].index
    print(f"{S}: beta={beta:.2f}, residual std {resids[S].std():.3f}%, {len(real_ons[S])} episode onsets")

# observed statistic (mean excursion across both stations)
obs = np.concatenate([excursion(resids[S], real_ons[S]) for S, _, _ in STATIONS])
obs_stat = obs.mean(); frac_pos = (obs > 0).mean()
print(f"\nOBSERVED co-onset excursion: mean {obs_stat:+.4f}% over {len(obs)} episodes, {frac_pos*100:.0f}% positive")

# permutation null: fake onsets avoiding +-60d of real
NPERM = 1000
null = np.empty(NPERM)
for p in range(NPERM):
    fake_ex = []
    for S, _, _ in STATIONS:
        days = all_days[S]; k = len(real_ons[S]); ro = real_ons[S]
        cand = days[(days > days[0] + pd.Timedelta("60D")) & (days < days[-1] - pd.Timedelta("30D"))]
        ok = [d for d in cand if (np.abs((ro - d).days) >= 60).all()]
        if len(ok) < k: continue
        pick = pd.DatetimeIndex(rng.choice(ok, k, replace=False))
        fake_ex.append(excursion(resids[S], pick))
    null[p] = np.concatenate(fake_ex).mean() if fake_ex else np.nan
null = null[~np.isnan(null)]
p95 = np.quantile(null, 0.95); p975 = np.quantile(null, 0.975); p_val = (null >= obs_stat).mean()
print(f"NULL (n={len(null)}): mean {null.mean():+.4f}%, 95th {p95:+.4f}%, 97.5th {p975:+.4f}%")
print(f"-> p(null >= observed) = {p_val:.3f}   {'SIGNIFICANT (outside 95% null)' if obs_stat>p95 else 'inside null band (not significant)'}")

# beta-guard: re-fit beta excluding episode days
print("\n=== beta-guard (beta fit excluding +-30d of episodes) ===")
guard = []
for S, _, _ in STATIONS:
    fams = set(pd.read_csv(f"data/{S.lower()}_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    f = deseason(cm(f"data/daily_dvv_{S}p90f40_Z_2to4.csv", fams)); r = deseason(cm(f"data/daily_dvv_{S}p90f40rev_Z_2to4.csv", fams))
    j = pd.concat([f.rename("f"), r.rename("r")], axis=1).dropna()
    epi = np.zeros(len(j), bool)
    for o in real_ons[S]:
        epi |= (j.index >= o - pd.Timedelta("30D")) & (j.index < o + pd.Timedelta("40D"))
    beta_ne = np.polyfit(j.r[~epi], j.f[~epi], 1)[0]
    res_ne = j.f - beta_ne * j.r
    guard.append(excursion(res_ne, real_ons[S]))
gstat = np.concatenate(guard).mean()
print(f"excursion with episode-excluded beta: {gstat:+.4f}% (was {obs_stat:+.4f}%, change {100*abs(gstat-obs_stat)/abs(obs_stat):.0f}%)")
print("VERDICT:", "residual signal SURVIVES guards + permutation" if (obs_stat>p95 and abs(gstat-obs_stat)/abs(obs_stat)<0.3)
      else "NOT decisively validated -> product is a bound (per Merlin step 4/5 failure branch)")
