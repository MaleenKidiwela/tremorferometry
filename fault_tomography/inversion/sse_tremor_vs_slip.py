#!/usr/bin/env python
"""PNSN tremor vs Gualandi slip (user 2026-09-10). Tremor (catalogs/pnsn_tremor_cascadia_full.csv, 2010-2026) is assigned to the
nearest model mesh triangle (<= 40 km) -> daily counts per triangle -> 30-d trailing counts (same frame as slip_thrust_inc30).
(1) Latitude-time: 0.25-deg bands, tremor 30-d count vs area-weighted slip in window; per-band Spearman over time and the lag of
    maximum cross-correlation (slip vs tremor, -60..+60 d).
(2) Per slip patch (unbridged catalog, substantive, 2010+): tremor detections on the footprint per day in PRE(-120..-31) / DUR /
    POST(+31..+120); tremor-confirmed = DUR rate >= 3x PRE rate and >= 20 detections; rates by segment / Mw class / broad pulse;
    the MIRROR (back-slip) catalog as control.
(3) Per tremor episode (tremor_episodes.csv): fraction whose footprint has modeled slip > 5 mm/30 d at the peak, by segment.
(4) Scaling: DUR tremor count vs Mw (log-log) by segment.  (5) Timing: tremor onset on the footprint vs slip start (trailing frame:
    slip start should LAG tremor onset by 0-30 d).  (6) Position: tremor centroid during DUR vs footprint centroid (along-strike /
    across-strike km).
Outputs: res_catalog_sse/tremor_vs_slip_patches.csv, tremor_vs_slip_bands.csv, tremor_vs_slip_summary.txt, figures/sse_tremor_vs_slip.png
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, theilslopes
SD = "fault_tomography/inversion/res_catalog_sse"; REPO = "/home/jovyan/cascadia-sse-strain"
zl = np.load(f"{SD}/slip_patches_5mm.npz", allow_pickle=True); L = zl["label"].T; inc = zl["inc30_mm"].T; days = pd.to_datetime(zl["days"]); G = len(days); gi = pd.Series(np.arange(G), index=days)
plat, plon, area = zl["patch_lat"], zl["patch_lon"], zl["patch_area_km2"]; NP = len(plat)
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); TXY = km(plat, plon)
tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"]); tr["date"] = pd.to_datetime(tr.time).dt.normalize(); tr = tr[(tr.date >= days[0]) & (tr.date <= days[-1])]
d_, j_ = cKDTree(TXY).query(km(tr.lat.values, tr.lon.values)); ok = d_ <= 40; tr = tr[ok]; tri = j_[ok]; ti = gi[tr.date].values
TR = np.zeros((NP, G)); np.add.at(TR, (tri, ti), 1); TR30 = np.cumsum(TR, 1); TR30[:, 30:] -= TR30[:, :-30].copy()      # detections within (t-30, t]
t0 = int(gi[pd.Timestamp("2010-01-01")]); out = []
out.append(f"tremor: {len(tr)} detections assigned to mesh triangles (<=40 km) 2010-2026 of {ok.size} in catalog window; {(~ok).sum()} dropped (>40 km from mesh)")
# ---------- (1) latitude bands ----------
bands = np.arange(40.0, 50.0, 0.25); rows_b = []; TRB = []; SLB = []
for b in bands:
    m = (plat >= b) & (plat < b + 0.25)
    if m.sum() == 0: TRB.append(np.full(G, np.nan)); SLB.append(np.full(G, np.nan)); continue
    tb = TR30[m].sum(0); sb = np.nansum(inc[m] * area[m, None], 0) / area[m].sum(); TRB.append(tb); SLB.append(sb)
    x, y = tb[t0:], sb[t0:]; okk = np.isfinite(x) & np.isfinite(y); rho = spearmanr(x[okk], y[okk]).correlation
    lags = np.arange(-60, 61, 5); cc = [np.corrcoef(x[okk][max(0, -lg):len(x[okk]) - max(0, lg)], y[okk][max(0, lg):len(y[okk]) - max(0, -lg)])[0, 1] for lg in lags]   # positive lag: slip lags tremor
    rows_b.append(dict(band=b, n_tri=int(m.sum()), tremor_total=int(TR[m, t0:].sum()), slip_pos_days=int((sb[t0:] > 5).sum()), spearman=rho, best_lag_d=int(lags[int(np.nanargmax(cc))]), best_cc=float(np.nanmax(cc)), cc_lag0=float(cc[list(lags).index(0)])))
B = pd.DataFrame(rows_b); B.to_csv(f"{SD}/tremor_vs_slip_bands.csv", index=False); TRB = np.array(TRB); SLB = np.array(SLB)
out.append(f"bands: tremor-vs-slip Spearman median {B.spearman.median():.2f} (north >=46: {B[B.band>=46].spearman.median():.2f}, central 44-46: {B[(B.band>=44)&(B.band<46)].spearman.median():.2f}, south <44: {B[B.band<44].spearman.median():.2f}); best lag median {B.best_lag_d.median():+.0f} d (slip lags tremor if positive), IQR {B.best_lag_d.quantile(.25):+.0f}..{B.best_lag_d.quantile(.75):+.0f}")
# ---------- (2) per slip patch ----------
cat = pd.read_csv(f"{SD}/slip_patches_5mm.csv", parse_dates=["start", "end"]); seg = lambda la: "north" if la >= 46 else ("central" if la >= 44 else "south")
def patch_rows(cat, Lab, tag):
    R = []
    for _, r in cat[(cat.start >= "2010-06-01") & ~cat.minor].iterrows():
        ts, te = int(gi[r.start]), int(gi[r.end]); fp = np.unique(np.where(Lab[:, ts:te + 1] == r.id)[0])
        w = lambda a, b: TR[fp][:, max(a, 0):min(b, G - 1) + 1]
        pre, dur, post = w(ts - 120, ts - 31), w(ts, te), w(te + 31, te + 120); rp, rd, rpo = pre.sum() / max(pre.shape[1], 1), dur.sum() / max(dur.shape[1], 1), post.sum() / max(post.shape[1], 1)
        # tremor onset on footprint: first day in [ts-60, te] with 7-d footprint count >= max(10, 3x PRE 7-d mean)
        c7 = np.convolve(TR[fp].sum(0), np.ones(7), "same"); thr = max(10, 3 * 7 * rp); cand = np.where(c7[max(ts - 60, 0):te + 1] >= thr)[0]; onset = int(cand[0] + max(ts - 60, 0)) if len(cand) else -1
        # tremor centroid during DUR vs footprint centroid
        cnt = TR[fp][:, ts:te + 1].sum(1); cen_t = (TXY[fp] * cnt[:, None]).sum(0) / cnt.sum() if cnt.sum() > 0 else np.full(2, np.nan); cen_f = (TXY[fp] * area[fp, None]).sum(0) / area[fp].sum()
        R.append(dict(catalog=tag, id=r.id, start=r.start.date(), end=r.end.date(), seg=seg(r.lat_centroid), Mw=r.Mw, broad_pulse=bool(r.broad_pulse), n_tri=len(fp), dur_d=te - ts + 1, pre_rate=rp, dur_rate=rd, post_rate=rpo, dur_total=int(dur.sum()),
                      ratio_dur_pre=rd / rp if rp > 0 else (np.inf if rd > 0 else np.nan), confirmed=bool(rd >= 3 * rp and dur.sum() >= 20), tremor_onset_minus_slip_start_d=(onset - ts) if onset >= 0 else np.nan,
                      dx_km=cen_t[0] - cen_f[0], dy_km=cen_t[1] - cen_f[1]))
    return pd.DataFrame(R)
P = patch_rows(cat, L, "thrust"); zm = np.load(f"{SD}/slip_patches_5mm_mirror.npz", allow_pickle=True); Pm = patch_rows(pd.read_csv(f"{SD}/slip_patches_5mm_mirror.csv", parse_dates=["start", "end"]), zm["label"].T, "mirror")
ALL = pd.concat([P, Pm]); ALL.to_csv(f"{SD}/tremor_vs_slip_patches.csv", index=False)
def rate(df): return f"{df.confirmed.mean() * 100:.0f} % of {len(df)}"
out.append("tremor-confirmed slip patches (DUR tremor rate >= 3x PRE and >= 20 detections), substantive, 2010+:")
for mw, lab in ((6.2, "Mw>=6.2"), (6.0, "Mw>=6.0"), (0, "all substantive")):
    q = P[(P.Mw >= mw) & ~P.broad_pulse]; qm = Pm[(Pm.Mw >= mw) & ~Pm.broad_pulse]
    out.append(f"  {lab:16s}: thrust {rate(q)} [north {rate(q[q.seg=='north'])}, central {rate(q[q.seg=='central'])}, south {rate(q[q.seg=='south'])}] | MIRROR control {rate(qm)} [north {rate(qm[qm.seg=='north'])}, south {rate(qm[qm.seg=='south'])}]")
out.append(f"  broad pulses (thrust, any Mw): {rate(P[P.broad_pulse])}")
q = P[(P.Mw >= 6.0) & ~P.broad_pulse & P.confirmed]; out.append(f"timing (confirmed, Mw>=6.0, n {q.tremor_onset_minus_slip_start_d.notna().sum()}): tremor onset minus slip start median {q.tremor_onset_minus_slip_start_d.median():+.0f} d, IQR {q.tremor_onset_minus_slip_start_d.quantile(.25):+.0f}..{q.tremor_onset_minus_slip_start_d.quantile(.75):+.0f} (expect -30..0: slip start is the trailing-window threshold crossing)")
out.append(f"position (confirmed, Mw>=6.0): tremor centroid minus footprint centroid: E-W median {q.dx_km.median():+.0f} km (IQR {q.dx_km.quantile(.25):+.0f}..{q.dx_km.quantile(.75):+.0f}), N-S median {q.dy_km.median():+.0f} km (IQR {q.dy_km.quantile(.25):+.0f}..{q.dy_km.quantile(.75):+.0f}); positive E-W = tremor DOWN-DIP (east) of the slip footprint")
for s in ("north", "central", "south"):
    qq = q[q.seg == s]
    if len(qq) >= 5: out.append(f"   {s}: E-W offset median {qq.dx_km.median():+.0f} km, N-S {qq.dy_km.median():+.0f} km, timing {qq.tremor_onset_minus_slip_start_d.median():+.0f} d (n {len(qq)})")
# ---------- (3) tremor episodes -> modeled slip ----------
E = pd.read_csv(f"{SD}/tremor_episodes.csv"); out.append(f"tremor episodes (band-defined, {len(E)}): modeled slip in footprint at peak > 5 mm/30 d in {(E.model_slip_mm > 5).mean() * 100:.0f} % [" + ", ".join(f"{s}: {(E[E.seg == s].model_slip_mm > 5).mean() * 100:.0f} % of {len(E[E.seg == s])}" for s in E.seg.unique()) + f"]; median modeled slip at tremor peaks {E.model_slip_mm.median():.1f} mm/30 d")
# ---------- (4) scaling ----------
q = P[(P.Mw >= 6.0) & ~P.broad_pulse & (P.dur_total > 0)]; sl = theilslopes(np.log10(q.dur_total), q.Mw); out.append(f"scaling: log10(tremor count during slip) vs Mw slope {sl[0]:.2f} per Mw unit (n {len(q)}; Spearman {spearmanr(q.Mw, q.dur_total).correlation:+.2f})" + "; per segment Spearman " + ", ".join(f"{s} {spearmanr(q[q.seg == s].Mw, q[q.seg == s].dur_total).correlation:+.2f}" for s in ("north", "central", "south") if (q.seg == s).sum() >= 5))
open(f"{SD}/tremor_vs_slip_summary.txt", "w").write("\n".join(out)); print("\n".join(out))
# ---------- figure ----------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig = plt.figure(figsize=(20, 14), constrained_layout=True); gs = fig.add_gridspec(3, 3, height_ratios=[1.3, 1.3, 1])
a = fig.add_subplot(gs[0, :]); im = a.pcolormesh(days[t0:], bands + .125, np.log10(TRB[:, t0:] + 1), cmap="magma_r", vmin=0, vmax=3, shading="nearest"); plt.colorbar(im, ax=a, pad=.01, label="log10(1 + tremor detections in 30 d)")
sub = cat[(cat.start >= "2010-06-01") & ~cat.minor & ~cat.broad_pulse & (cat.Mw >= 6.0)]
for _, r in sub.iterrows(): a.plot([r.start, r.end], [r.lat_min, r.lat_min], color="#38bdf8", lw=1); a.plot([r.start, r.end], [r.lat_max, r.lat_max], color="#38bdf8", lw=1); a.plot([r.start, r.start], [r.lat_min, r.lat_max], color="#38bdf8", lw=1); a.plot([r.end, r.end], [r.lat_min, r.lat_max], color="#38bdf8", lw=1)
a.set_ylabel("latitude"); a.set_title("PNSN tremor (30-d trailing count per 0.25° band) with modeled slip patches (boxes: Mw ≥ 6.0, unbridged catalog)", fontsize=10)
a = fig.add_subplot(gs[1, :]); im = a.pcolormesh(days[t0:], bands + .125, SLB[:, t0:], cmap="RdBu_r", vmin=-15, vmax=15, shading="nearest"); plt.colorbar(im, ax=a, pad=.01, label="modeled thrust slip in 30 d (mm, band mean)")
a.set_ylabel("latitude"); a.set_title("Gualandi model: thrust slip within the trailing 30 d (band area-weighted mean); note the symmetric ± oscillation in the south", fontsize=10)
a = fig.add_subplot(gs[2, 0]); a.plot(B.spearman, B.band + .125, "o-", color="#1d4ed8", label="Spearman(tremor, slip) over time"); a.plot(B.cc_lag0, B.band + .125, "s--", color="#64748b", ms=4, label="Pearson at lag 0"); a.axvline(0, color="k", lw=.5); a.set_xlabel("correlation"); a.set_ylabel("band latitude"); a.legend(fontsize=8); a.set_title("per-band tremor–slip correlation 2010-26", fontsize=9)
a2 = a.twiny(); a2.plot(B.best_lag_d, B.band + .125, "^", color="#ea580c", ms=5, label="lag of max CC (d; + = slip lags)"); a2.set_xlabel("best lag (days, + = slip lags tremor)", color="#ea580c"); a2.set_xlim(-60, 60)
a = fig.add_subplot(gs[2, 1]); COL = dict(north="#1d4ed8", central="#059669", south="#ea580c")
for s in COL:
    qq = P[(P.Mw >= 6.0) & ~P.broad_pulse & (P.seg == s)]; a.plot(qq.Mw, qq.dur_total + 1, "o", color=COL[s], ms=5, alpha=.7, label=f"{s} (confirmed {qq.confirmed.mean() * 100:.0f} %)")
    qm_ = Pm[(Pm.Mw >= 6.0) & ~Pm.broad_pulse & (Pm.seg == s)]; a.plot(qm_.Mw, qm_.dur_total + 1, "x", color=COL[s], ms=4, alpha=.5)
a.set_yscale("log"); a.set_xlabel("patch Mw"); a.set_ylabel("tremor detections on footprint during slip (+1)"); a.set_title("tremor count vs Mw (o thrust patches, x mirror/back-slip control)", fontsize=9); a.legend(fontsize=8)
a = fig.add_subplot(gs[2, 2]); q = P[(P.Mw >= 6.0) & ~P.broad_pulse & P.confirmed]
a.hist(q.tremor_onset_minus_slip_start_d.dropna(), bins=np.arange(-62, 63, 7), color="#1d4ed8", alpha=.7); a.axvline(0, color="k", lw=.8); a.axvspan(-30, 0, color="#e5e7eb", alpha=.6); a.set_xlabel("tremor onset − slip start (days)"); a.set_ylabel("confirmed patches"); a.set_title("timing: tremor onset relative to the slip start (grey = expected −30..0 from the trailing window)", fontsize=9)
fig.suptitle("PNSN tremor vs Gualandi (2025) slip, 2010-2026", fontsize=12); fig.savefig("figures/sse_tremor_vs_slip.png", dpi=120); print("wrote figures/sse_tremor_vs_slip.png, tremor_vs_slip_{patches,bands}.csv, tremor_vs_slip_summary.txt")
