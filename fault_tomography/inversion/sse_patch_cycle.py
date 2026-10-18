#!/usr/bin/env python
"""PATCH-CYCLE analysis (user hint 2026-09-04): each episode's slip patch is a fixed spatial unit; its dv/v is followed
through the WHOLE record and folded on the patch's OWN slip cycle.
For each model-episode footprint (patch): patch-mean slip-in-window series -> the patch's own slip events (peaks >= 3x
negative-side RMS, >= 90 d apart, 2010-2026); patch-mean dv/v series (interface deseasoned, family 2-4 s, family 3-5 s;
also the difference patch-minus-ring). Folding: (a) time since the patch's last slip peak (0..450 d, 30-d bins);
(b) cycle phase between consecutive peaks (10 bins). Aggregate over patches (each patch weighted equally, distinct
patches merged when >50 % cell overlap), null = rigid time shift of the dv/v series relative to the slip history (300
draws). Statistics: phase-curve amplitude (max-min), the slope of dv/v vs time-since-slip (loading trend), and the drop
across the slip (bin 0 minus last bin). Per-patch versions too.  Outputs: res_catalog_sse/patch_cycle.csv,
figures/sse_patch_cycle.png
"""
import numpy as np, pandas as pd
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(2026)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
clat, clon, dep = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"]; gi = pd.Series(np.arange(G), index=days); inc = zc["slip_thrust_inc30"].astype(float) * 1e3
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
win = np.asarray((days >= "2010-06-01") & (days <= "2026-07-20")); month = days.month.values
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok]; zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
mask = np.zeros(G, bool)
for i in okep.ipeak: mask[max(int(i) - 45, 0):int(i) + 91] = True
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); md = pd.to_datetime(zp["days"]); Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=md).reindex(days).values
F24 = np.load(f"{SD}/product_family24.npz", allow_pickle=True)["MODEL"].astype(float); F35 = np.load(f"{SD}/product_family35.npz", allow_pickle=True)["MODEL"].astype(float)
def deseason(y):
    o = y.copy()
    for m in range(1, 13): s = (month == m) & ~mask & np.isfinite(y); o[month == m] -= np.nanmean(y[s]) if s.any() else 0
    return o
def rmean(X, cs): v = X[:, cs]; o = np.nanmean(v, 1); o[np.isfinite(v).mean(1) < 0.5] = np.nan; return o
# distinct patches (merge footprints with >50 % overlap, keep the larger)
P = []
for _, r in okep.sort_values("n_fp", ascending=False).iterrows():
    fp = zf[f"fp_{r.key}"]
    if any(len(np.intersect1d(fp, p["fp"])) > 0.5 * min(len(fp), len(p["fp"])) for p in P): continue
    d = np.sqrt(((XY[:, None, :] - XY[fp][None, :, :]) ** 2).sum(-1)).min(1); ring = np.where(np.isfinite(dep) & (dep >= 30) & (dep <= 45) & (d >= 100) & (d <= 300))[0]
    P.append(dict(key=r.key, seg=r.seg, fp=fp, ring=ring))
print(f"distinct patches: {len(P)} (from {len(okep)} episodes)", flush=True)
TSB = np.arange(0, 451, 30); NB = len(TSB) - 1; PHB = np.linspace(0, 1, 11)
def fold(y, peaks):
    """(a) mean dv/v by time-since-last-peak bin, (b) by phase bin between consecutive peaks"""
    ts = np.full(G, np.nan); ph = np.full(G, np.nan)
    for a, b in zip(peaks[:-1], peaks[1:]):
        idx = np.arange(a, b); ts[idx] = idx - a; ph[idx] = (idx - a) / (b - a)
    fa = np.array([np.nanmean(y[win & (ts >= TSB[k]) & (ts < TSB[k + 1])]) for k in range(NB)]); fb = np.array([np.nanmean(y[win & (ph >= PHB[k]) & (ph < PHB[k + 1])]) for k in range(10)])
    return fa - np.nanmean(fa), fb - np.nanmean(fb)
PRODS = {"interface": deseason, "family 2-4 s": None, "family 3-5 s": None}
res = {}; rows = []
for pname, X in (("interface", Y), ("family 2-4 s", F24), ("family 3-5 s", F35)):
    FA, FB, FAr, FBr, nulA, nulB, nulD = [], [], [], [], [], [], []
    for p in P:
        s = rmean(inc, p["fp"]); neg = np.sqrt(np.nanmean(np.minimum(np.nan_to_num(s), 0) ** 2)); thr = 3 * neg; v = np.nan_to_num(s, nan=-np.inf)
        cand = sorted([i for i in range(1, G - 1) if v[i] >= thr and v[i] >= v[i - 1] and v[i] >= v[i + 1] and win[i]], key=lambda i: -v[i]); pk = []
        for i in cand:
            if all(abs(i - k) >= 90 for k in pk): pk.append(i)
        pk = sorted(pk)
        if len(pk) < 4: continue
        y = deseason(rmean(X, p["fp"])); yr = deseason(rmean(X, p["ring"])) if len(p["ring"]) >= 10 else np.full(G, np.nan)
        fa, fb = fold(y, pk); far, fbr = fold(y - yr, pk); FA.append(fa); FB.append(fb); FAr.append(far); FBr.append(fbr)
        na, nb, nd = [], [], []
        for _ in range(100):
            sh = rng.randint(365, G - 365); ys = np.roll(y, sh); a_, b_ = fold(ys, pk); na.append(a_); nb.append(b_); nd.append(fold(np.roll(y - yr, sh), pk)[0])
        nulA.append(np.array(na)); nulB.append(np.array(nb)); nulD.append(np.array(nd))
        amp = np.nanmax(fa) - np.nanmin(fa); slope = np.polyfit(TSB[:-1][np.isfinite(fa)], fa[np.isfinite(fa)], 1)[0] * 100 if np.isfinite(fa).sum() > 3 else np.nan
        namp = [np.nanmax(a_) - np.nanmin(a_) for a_ in na]; rows.append(dict(product=pname, patch=p["key"], seg=p["seg"], n_cycles=len(pk) - 1, recurrence_d=float(np.median(np.diff(pk))), amp=amp, p_amp=float(np.mean(np.array(namp) >= amp)), slope_per100d=slope, drop_at_slip=fa[0] - fa[-1] if np.isfinite(fa[0]) and np.isfinite(fa[-1]) else np.nan))
    FA, FB, FAr = np.array(FA), np.array(FB), np.array(FAr); NA, NB_, ND = np.array(nulA), np.array(nulB), np.array(nulD)      # (patches, draws, bins)
    ca, cb, car = np.nanmean(FA, 0), np.nanmean(FB, 0), np.nanmean(FAr, 0); nca, ncb, ncar = np.nanmean(NA, 0), np.nanmean(NB_, 0), np.nanmean(ND, 0)   # (draws, bins) aggregated null
    amp = np.nanmax(ca) - np.nanmin(ca); namp = np.nanmax(nca, 1) - np.nanmin(nca, 1); sl = np.polyfit(TSB[:-1], ca, 1)[0] * 100; nsl = [np.polyfit(TSB[:-1], r_, 1)[0] * 100 for r_ in nca]
    ampr = np.nanmax(car) - np.nanmin(car); nampr = np.nanmax(ncar, 1) - np.nanmin(ncar, 1)
    res[pname] = dict(ca=ca, cb=cb, car=car, lo=np.nanpercentile(nca, 2.5, 0), hi=np.nanpercentile(nca, 97.5, 0), lob=np.nanpercentile(ncb, 2.5, 0), hib=np.nanpercentile(ncb, 97.5, 0), n=len(FA))
    print(f"{pname:13s}: {len(FA)} patches | time-since-slip curve: amplitude {amp:.3f}% (null {np.mean(namp):.3f} +- {np.std(namp):.3f}, p {np.mean(namp >= amp):.3f}); slope {sl:+.4f} %/100 d (null sd {np.std(nsl):.4f}, p {np.mean(np.abs(nsl) >= abs(sl)):.3f}); bin0-minus-last {ca[0]-ca[-1]:+.3f} | patch-minus-ring amplitude {ampr:.3f} (p {np.mean(nampr >= ampr):.3f})", flush=True)
    print("   curve (0-30 ... 420-450 d):", np.round(ca, 3)); print("   phase curve (0..1):", np.round(cb, 3))
T = pd.DataFrame(rows); T.to_csv(f"{SD}/patch_cycle.csv", index=False)
print(f"\nper-patch amplitude p<0.05: " + ", ".join(f"{pn}: {(T[T['product']==pn].p_amp<0.05).sum()}/{(T['product']==pn).sum()}" for pn in res))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(17, 4.2))
for k, (pname, r) in enumerate(res.items()):
    a = ax[k]; a.fill_between(TSB[:-1] + 15, r["lo"], r["hi"], color="#e5e7eb", label="null 95%"); a.plot(TSB[:-1] + 15, r["ca"], "o-", color="#1d4ed8", label="patch dv/v"); a.plot(TSB[:-1] + 15, r["car"], "s--", color="#059669", label="patch − ring")
    a.axhline(0, color="k", lw=.5); a.set_xlabel("days since the patch's last slip peak"); a.set_ylabel("dv/v rel. cycle mean (%)"); a.set_title(f"{pname}: {r['n']} patches folded on their own slip cycles", fontsize=9); a.legend(fontsize=7)
fig.suptitle("Patch-cycle folding: dv/v of each slip patch through its own inter-slip interval", fontsize=10); fig.savefig("figures/sse_patch_cycle.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_patch_cycle.png, patch_cycle.csv")
