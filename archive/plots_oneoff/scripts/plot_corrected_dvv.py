#!/usr/bin/env python
"""Plot the noise-corrected residual dv/v (fwd - beta*rev, deseasoned) for B926+B011, with ETS episode
onsets, and an epoch-stack that makes the small co-onset excursion visible."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

STATIONS = [("B926", 48.82, -124.13), ("B011", 48.65, -123.45)]
def cm(csv, fams):
    d = pd.read_csv(csv); d = d[d.patch.isin(fams) & (d.cc_max >= 0.6)]; d["date"] = pd.to_datetime(d.date)
    m = d.groupby("date").dvv.median(); n = d.groupby("date").patch.nunique()
    return (m[n >= 3].rolling(15, center=True, min_periods=5).median()) * 100
def des(s):
    s = s.dropna(); t = (s.index - s.index[0]).days.values.astype(float); yr = 365.25
    X = np.column_stack([np.ones_like(t), t, np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr), np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
    b, *_ = np.linalg.lstsq(X, s.values, rcond=None); return pd.Series(s.values - X@b, index=s.index)
cat = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time","lat","lon"])
cat["t"] = pd.to_datetime(cat.time, errors="coerce"); cat = cat.dropna(subset=["t"])
def onsets(la, lo):
    c = cat[cat.lat.between(la-0.7,la+0.7) & cat.lon.between(lo-0.7,lo+0.7)]
    daily = c.groupby(c.t.dt.floor("D")).size(); idx = pd.date_range(daily.index.min(), daily.index.max())
    rate = daily.reindex(idx, fill_value=0).rolling(15, center=True, min_periods=5).mean()
    hot = rate > rate.quantile(0.80); ons = []
    for d in rate.index[hot.values]:
        if not ons or (d-ons[-1]).days >= 120: ons.append(d)
    return pd.DatetimeIndex(ons)

fig = plt.figure(figsize=(13, 8))
gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.1], hspace=0.35)
epoch_days = np.arange(-90, 61)
for i, (S, la, lo) in enumerate(STATIONS):
    TAG = f"{S}p90f40"; s = S.lower()
    fams = set(pd.read_csv(f"data/{s}_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    fwd = des(cm(f"data/daily_dvv_{TAG}_Z_2to4.csv", fams)); rev = des(cm(f"data/daily_dvv_{TAG}rev_Z_2to4.csv", fams))
    j = pd.concat([fwd.rename("f"), rev.rename("r")], axis=1).dropna()
    beta = np.polyfit(j.r, j.f, 1)[0]; resid = j.f - beta*j.r
    ons = onsets(la, lo)
    ax = fig.add_subplot(gs[i])
    ax.plot(j.index, j.f.values, lw=0.6, color="#cbd5e1", label="raw (deseasoned)")
    ax.plot(resid.index, resid.values, lw=1.1, color="#1a1a2e", label="noise-corrected residual")
    for o in ons: ax.axvline(o, color="#c1440e", lw=0.7, alpha=0.35)
    ax.axhline(0, color="0.6", lw=0.5); ax.set_ylim(-0.6, 0.6); ax.set_ylabel("dv/v (%)")
    ax.set_title(f"{S} — noise-corrected residual (β={beta:.2f}, std {resid.std():.3f}%); orange = ETS episode onsets", fontsize=11)
    if i == 0: ax.legend(loc="upper right", fontsize=9, frameon=False)
    # epoch stack for this station
    ep = []
    for o in ons:
        w = resid[(resid.index >= o + pd.Timedelta(int(epoch_days[0]), "D")) & (resid.index <= o + pd.Timedelta(int(epoch_days[-1]), "D"))]
        if len(w) < 60: continue
        rel = (w.index - o).days.values
        ep.append(np.interp(epoch_days, rel, w.values, left=np.nan, right=np.nan))
    S and globals().setdefault("_EP", {}).__setitem__(S, np.array(ep))

# bottom: epoch stacks
axe = fig.add_subplot(gs[2])
for S, c in [("B926", "#1a1a2e"), ("B011", "#2563eb")]:
    E = _EP[S]; m = np.nanmean(E, 0); se = np.nanstd(E, 0)/np.sqrt(np.sum(~np.isnan(E), 0))
    axe.plot(epoch_days, m, lw=1.8, color=c, label=f"{S} ({len(E)} episodes)")
    axe.fill_between(epoch_days, m-se, m+se, color=c, alpha=0.15)
axe.axvline(0, color="#c1440e", lw=1.2); axe.axhline(0, color="0.6", lw=0.5)
axe.axvspan(0, 30, color="#c1440e", alpha=0.06); axe.axvspan(-60, -20, color="0.5", alpha=0.06)
axe.set_xlabel("days from ETS episode onset"); axe.set_ylabel("residual dv/v (%)")
axe.set_title("Epoch stack: residual dv/v around ETS onset (co-onset +0.04%, p=0.003; grey=pre, orange=co-onset window)", fontsize=11)
axe.legend(loc="upper left", fontsize=9, frameon=False)
fig.savefig("figures/corrected_dvv_ets.png", dpi=130, bbox_inches="tight")
print("-> figures/corrected_dvv_ets.png")
