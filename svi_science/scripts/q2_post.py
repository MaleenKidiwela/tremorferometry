"""Q2 post-checks (NOT pre-registered; labelled post-hoc in RESULTS.md):
0. fix: thinned-Bostock columns were merged by episode name only and leaked onto NW-E rows with the same name -> blank them.
1. NaN-safe recomputation of T9 (SVI), T10 ratio, T11 tremor, duration-vs-speed in SVI-O (an SVI-O episode has no source with >= 30).
2. Floor-15 (W_dur15) and patch (W_dur_patch) versions of T2-T5.
3. Random thinning of Bostock leaves no source with >= 30 events -> redo random thinning with the 15-event floor and patches.
"""
import json, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q2_band as Q   # NB: re-opens q2.log (main log was copied to q2_main.log first)

OUT = Q.OUT
EP = pd.read_csv(f"{OUT}/episodes.csv")
for c in [c for c in EP.columns if c.endswith("_thinM") or c.endswith("_thinR")]:
    EP.loc[EP.era != "SVI-B 2003-13", c] = np.nan
res = {}

# 3. random thinning with floor 15 and patches
EV = Q.load_svi_b(); target = int(EP[EP.era == "SVI-O 2014-26"].n.median()); rows = []
for ep, g in EV.groupby("ep"):
    g = g.reset_index(drop=True)
    for k in range(50):
        gi = g.iloc[np.sort(Q.rng.choice(len(g), min(target, len(g)), replace=False))].reset_index(drop=True)
        r, _ = Q.episode_measure(gi, full=False, nsub_d=50, nsub_s=10)
        rows.append(dict(ep=ep, draw=k, v=r["v"], W_dur15=r["W_dur15"], W_dur_patch=r["W_dur_patch"], n_units15=r["n_units15"]))
    gm = g.sort_values("mag", ascending=False).iloc[:target].sort_values("time").reset_index(drop=True)
    r, _ = Q.episode_measure(gm, full=False, nsub_d=200, nsub_s=10)
    rows.append(dict(ep=ep, draw=-1, v=r["v"], W_dur15=r["W_dur15"], W_dur_patch=r["W_dur_patch"], n_units15=r["n_units15"]))
TH = pd.DataFrame(rows); TH.to_csv(f"{OUT}/thinning_post.csv", index=False)
tr = TH[TH.draw >= 0].groupby("ep")[["W_dur15", "W_dur_patch", "n_units15"]].median().add_suffix("_thinR")
tm = TH[TH.draw == -1].set_index("ep")[["W_dur15", "W_dur_patch"]].add_suffix("_thinM")
b = EP.era == "SVI-B 2003-13"
for c in tr.columns: EP.loc[b, c] = EP.loc[b, "ep"].map(tr[c]).values
for c in tm.columns: EP.loc[b, c] = EP.loc[b, "ep"].map(tm[c]).values
bb = EP[b & EP.single_front]
res["thin_ratio_post"] = {f"{c}_{s}": float(np.nanmedian(bb[f"{c}_{s}"] / bb[c])) for c in ("W_dur15", "W_dur_patch") for s in ("thinR", "thinM")}
res["thin_ratio_post"]["median_units15_after_random_thinning"] = float(bb.n_units15_thinR.median())

# 2. column versions of T2-T5
res["columns"] = {}
for col in ("W_dur15", "W_dur_patch", "W_dur_arr", "W_front"):
    d = {}
    for reg, g in EP[EP.single_front].groupby("region"):
        d[f"trend_{reg}"] = Q.trend(g.year.values, g[col].values)
        d[f"median_{reg}"] = float(g[col].median())
        d[f"CV_{reg}"] = float(g[col].std() / g[col].mean())
    sf = EP[EP.single_front]
    d["eras_SVI"] = Q.mw(sf[sf.era == "SVI-B 2003-13"][col], sf[sf.era == "SVI-O 2014-26"][col])
    d["eras_NW"] = Q.mw(sf[sf.era == "NW-E 2007-14"][col], sf[sf.era == "NW-O 2015-26"][col])
    d["regions"] = Q.mw(sf[sf.region == "SVI"][col], sf[sf.region == "NW"][col])
    if col in ("W_dur15", "W_dur_patch"):
        d["eras_SVI_thinR"] = Q.mw(sf[sf.era == "SVI-B 2003-13"][f"{col}_thinR"], sf[sf.era == "SVI-O 2014-26"][col])
        d["eras_SVI_thinM"] = Q.mw(sf[sf.era == "SVI-B 2003-13"][f"{col}_thinM"], sf[sf.era == "SVI-O 2014-26"][col])
    d["by_era_median"] = sf.groupby("era")[col].median().to_dict()
    res["columns"][col] = d

# 1. NaN-safe
sf = EP[EP.single_front & EP.W_dur.notna()]
res["T9_SVI_nan_safe"] = dict(SE=sf[(sf.region == "SVI") & (sf.v < 0)].W_dur.tolist(), NW=sf[(sf.region == "SVI") & (sf.v > 0)].W_dur.tolist(),
                              mw=Q.mw(sf[(sf.region == "SVI") & (sf.v < 0)].W_dur, sf[(sf.region == "SVI") & (sf.v > 0)].W_dur))
for reg in ("SVI", "NW"):
    g = EP[EP.single_front & (EP.region == reg)]
    res[f"T9_{reg}_W_dur15"] = Q.mw(g[g.v < 0].W_dur15, g[g.v > 0].W_dur15)
S = pd.read_csv(f"{OUT}/s13.csv")
res["T10_ratio_nan_safe"] = dict(median_ratio_W_dur=float(np.nanmedian(S.W_dur_s13 / S.W_dur)), n_valid=int((S.W_dur_s13 / S.W_dur).notna().sum()),
                                 ratios=(S.W_dur_s13 / S.W_dur).round(2).tolist())
m = EP.single_front & EP.W_dur.notna() & EP.W_dur_trem.notna() & (EP.p_speed_trem < 0.01) & (EP.v_trem.abs() >= 2)
r = stats.spearmanr(EP[m].W_dur, EP[m].W_dur_trem)
res["T11_nan_safe"] = dict(n=int(m.sum()), rho=float(r[0]), p=float(r[1]), median_ratio=float(np.median(EP[m].W_dur_trem / EP[m].W_dur)),
                           median_trem_by_era=EP[m].groupby("era").W_dur_trem.median().to_dict(),
                           median_lfe_by_era=EP[m].groupby("era").W_dur.median().to_dict(),
                           rho_v=float(stats.spearmanr(EP[m].v, EP[m].v_trem)[0]),
                           trem_eras_SVI=Q.mw(EP[m & (EP.era == "SVI-B 2003-13")].W_dur_trem, EP[m & (EP.era == "SVI-O 2014-26")].W_dur_trem),
                           trem_eras_NW=Q.mw(EP[m & (EP.era == "NW-E 2007-14")].W_dur_trem, EP[m & (EP.era == "NW-O 2015-26")].W_dur_trem))
g = EP[EP.single_front & (EP.era == "SVI-O 2014-26") & EP.medD80.notna()]
r = stats.spearmanr(g.v.abs(), g.medD80); res["SVI_O_D80_vs_v"] = dict(n=len(g), rho=float(r[0]), p=float(r[1]))
g = EP[EP.single_front & (EP.region == "SVI")]
r = stats.spearmanr(g.v.abs(), g.medD80_n15); res["SVI_all_D80n15_vs_v"] = dict(n=len(g), rho=float(r[0]), p=float(r[1]))
# Cochran with floor-15 is not possible (no sigma); report CV per era
res["CV_by_era"] = {e: dict(W_dur=float(x.W_dur.std() / x.W_dur.mean()), W_dur15=float(x.W_dur15.std() / x.W_dur15.mean()))
                    for e, x in EP[EP.single_front].groupby("era")}
EP.to_csv(f"{OUT}/episodes.csv", index=False)
json.dump(res, open(f"{OUT}/result_post.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
