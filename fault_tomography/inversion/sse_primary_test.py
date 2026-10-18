#!/usr/bin/env python
"""Merlin steps 3-5 (notes/2026-09-03_merlin_sse_plan_ruling.md): does interface dbeta/beta change in the footprint of a
modeled (tremor-admitted) slow-slip episode?  Pre-registered primaries P1/P2/P3 with Nulls A and B and four controls.

Grid: model days (2007-01-02..2026-08-18, daily). dbeta/beta = res_catalog/inversion_4d_daily.npz MODEL reindexed.
Episodes: res_catalog_sse/episodes.csv (ok = tremor-admitted, glitch-free, footprint >= 5 cells, footprint duration >= 7 d),
  footprints + per-cell LOCAL onsets from episodes_footprints.npz. PRIMARY set = segments N + C; S = secondary; SC excluded.
Deseason (episode-masked): per cell monthly climatology from days NOT within [-30, +90] d of any ok episode peak in the
  cell's segment; a segment whose months retain < 40 % of days is flagged inseparable.
Response: footprint-mean deseasoned dbeta/beta as a function of lag from LOCAL onset (each footprint cell aligned on its
  own onset, then averaged; PRIMARY) and from the footprint onset (SECONDARY).
P1 = mean over ON [+15, +60] minus PRE [-90, -30] (trailing lag already included), >= 60 % valid frames each; pooled sign
  count K (s < 0 = velocity decrease) and mean amplitude.  P2 = P1 minus the same contrast on the RING (cells >= 100 km
  from every footprint cell, <= 300 km, deep 30-45 km, patch_dist <= 15 km, own segment quiet within +-45 d; >= 10 cells
  spanning >= 100 km).  P3 = pooled OLS of footprint-mean 30-d dbeta/beta change on the footprint-mean 30-d change of the
  6 strain components (trailing frame, 10 km) + footprint-mean slip within (t-30,t], at lags -90..+150 step 15 d; R^2.
Null A = rigid circular shift of the segment's episode times (one delta per segment, uniform [90 d, span-90 d], wrap in
  2011-04..2026-04); Null B = integer-year shift (1..14 yr, random sign) + uniform +-15 d jitter (keeps seasonal phase);
  4000 draws each; two-sided p on K - N/2 and on the mean.  Controls at the true catalog vs Null A (1000): same contrast on
  footprint-mean cc, rows per cell-day, SITES of footprint stations; day-of-year concentration R of the episodes.
Also: leave-one-segment-out, episode-onset alignment, raw (undeseasoned), strict catalog, gateD0.3 product, and the
  exploratory lag composite with a max-statistic FWER over lags.  Step-3 stress acceptance (footprint-median shear change
  at 0.5/10/20 km, local orientation, raw-amplitude windows) is printed first.
Outputs: res_catalog_sse/primary_test.npz, primary_episode_table.csv, figures/sse_primary_composite.png
"""
import numpy as np, pandas as pd, sys, os
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src"); from sse_strain import strain as st
rng = np.random.RandomState(11); NDRAW = int(os.environ.get("NDRAW", "4000")); NCTL = 1000
PRODUCT = os.environ.get("DVV_PRODUCT", f"{RD}/inversion_4d_daily.npz")
PRE, ON = (-90, -30), (15, 60); LAGS = np.arange(-150, 201)
SEGLAT = {"S": (40.0, 42.5), "SC": (42.5, 45.0), "C": (45.0, 47.3), "N": (47.3, 49.5)}
E0, E1 = pd.Timestamp("2011-04-01"), pd.Timestamp("2026-04-01")

# ---------- data ----------
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True)
days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells); clat, clon = zc["cell_lat"], zc["cell_lon"]
dep = zc["iface_depth_km"]; pdist = zc["patch_dist_km"]; gi = pd.Series(np.arange(G), index=days)
zp = np.load(PRODUCT, allow_pickle=True); assert list(zp["cell"]) == list(cells)
mdays = pd.to_datetime(zp["days"]); Y = pd.DataFrame(np.asarray(zp["MODEL"], float), index=mdays).reindex(days).values
SITES = pd.DataFrame(np.asarray(zp["SITES"], float), index=mdays, columns=list(zp["sta_tags"])).reindex(days)
ep = pd.read_csv(f"{SD}/episodes.csv"); zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
ep["key"] = ep.seg + "_" + ep.peak
segcell = np.full(M, "", object)
for s, (a, b) in SEGLAT.items(): segcell[(clat >= a) & (clat < b)] = s
i0, i1 = int(gi[E0]), int(gi[E1]); span = i1 - i0
def wrap(i): return i0 + (i - i0) % span
# strain (trailing, 10 km) 30-d changes per cell and slip within window
zg = np.load(f"{SD}/strain_gf_offsets.npz", allow_pickle=True); E10 = zg["E"][1].astype(float)
A_trail = np.nan_to_num(zc["A_trail"].astype(float)); A_raw = np.nan_to_num(zc["A_raw"].astype(float))
dAt30 = np.zeros_like(A_trail); dAt30[30:] = A_trail[30:] - A_trail[:-30]
inc30 = zc["slip_thrust_inc30"].astype(float)                     # (G, M) m within (t-30, t]
# family-day controls: cc and rows per cell-day; station rows per cell
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "cell", "date", "cc"]); cpos = {c: i for i, c in enumerate(cells)}
fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos)
cd = fd.groupby(["ci", "date"]).agg(cc=("cc", "mean"), n=("cc", "size")).reset_index(); ok = cd.date.isin(days)
CC = np.full((G, M), np.nan); NR = np.full((G, M), np.nan)
CC[gi[cd.date[ok]].values, cd.ci[ok].values] = cd.cc[ok].values; NR[gi[cd.date[ok]].values, cd.ci[ok].values] = cd.n[ok].values
sta_rows = fd.groupby(["ci", "tag"]).size().reset_index(name="n"); del fd

# ---------- episode-masked deseason ----------
okep = ep[ep.ok]; peaks = {s: okep[okep.seg == s].ipeak.values for s in SEGLAT}
mask_days = {s: np.zeros(G, bool) for s in SEGLAT}
for s, pk in peaks.items():
    for i in pk: mask_days[s][max(i - 30, 0):i + 91] = True
Yd = Y.copy(); month = days.month.values
for s in SEGLAT:
    cs = np.where(segcell == s)[0]
    if not len(cs): continue
    free = ~mask_days[s]
    ret = [free[month == m].mean() for m in range(1, 13)]
    if min(ret) < 0.4: print(f"  segment {s}: month retention min {min(ret):.2f} -> seasonally INSEPARABLE (flag)")
    for m in range(1, 13):
        sel = (month == m) & free; clim = np.nanmean(Y[sel][:, cs], 0)
        Yd[np.ix_(month == m, cs)] -= clim[None, :]
deepc = np.isfinite(dep) & (dep >= 30) & (dep <= 45) & (pdist <= 15)
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111])
XY = km(clat, clon)

# ---------- episode structures ----------
class Ep: pass
EPS = []
for _, r in okep.iterrows():
    e = Ep(); e.seg = r.seg; e.key = r.key; e.ipeak = int(r.ipeak); e.fp = zf[f"fp_{r.key}"]; e.lo = zf[f"lo_{r.key}"]
    e.ion = int(gi[pd.Timestamp(r.fp_onset)]); e.doy = int(r.doy)
    d = np.sqrt(((XY[:, None, :] - XY[e.fp][None, :, :]) ** 2).sum(-1)).min(1)
    e.ring = np.where(deepc & (d >= 100) & (d <= 300))[0]
    e.ring_ok = len(e.ring) >= 10 and (np.ptp(XY[e.ring, 1]) >= 100 if len(e.ring) else False)
    sr = sta_rows[sta_rows.ci.isin(e.fp)].groupby("tag").n.sum(); sr = sr / sr.sum()
    e.stas = [t for t in sr[sr >= 0.05].index if t in SITES.columns]
    EPS.append(e)
print(f"episodes ok: {len(EPS)} | primary N+C: {sum(e.seg in ('N','C') for e in EPS)} | ring valid: {sum(e.ring_ok for e in EPS)}")

def aligned(X, e, shift=0, use_local=True, cellset=None):
    """footprint-mean of X as a function of LAGS; each cell aligned on its own local onset (+shift), else on fp onset."""
    cs = e.fp if cellset is None else cellset
    on = (e.lo if (use_local and cellset is None) else np.full(len(cs), e.ion)) + shift
    idx = on[:, None] + LAGS[None, :]; bad = (idx < 0) | (idx >= G); idx = np.clip(idx, 0, G - 1)
    v = X[idx, cs[:, None]]; v[bad] = np.nan
    return np.nanmean(v, 0) if np.isfinite(v).any() else np.full(len(LAGS), np.nan)
L0 = -LAGS[0]
def contrast(y, w_on=ON, w_pre=PRE, minfrac=0.6):
    a = y[L0 + w_on[0]:L0 + w_on[1] + 1]; b = y[L0 + w_pre[0]:L0 + w_pre[1] + 1]
    if np.isfinite(a).mean() < minfrac or np.isfinite(b).mean() < minfrac: return np.nan
    return np.nanmean(a) - np.nanmean(b)
def run(X, eps, shifts, use_local=True, ring=False):
    out = []
    for e in eps:
        sh = shifts.get(e.seg, 0)
        if ring:
            if not e.ring_ok: out.append(np.nan); continue
            out.append(contrast(aligned(X, e, sh, use_local, e.ring)))
        else: out.append(contrast(aligned(X, e, sh, use_local)))
    return np.array(out)
def stats(s):
    s = s[np.isfinite(s)]; return int((s < 0).sum()), len(s), float(np.mean(s)) if len(s) else np.nan
def nullA(): return {s: int(rng.randint(90, span - 90)) for s in SEGLAT}
def nullB(): return {s: int(rng.choice([-1, 1]) * rng.randint(1, 15) * 365 + rng.randint(-15, 16)) for s in SEGLAT}
def wrapped_eps(eps, shifts):
    """episodes with onsets shifted (wrapped into [E0, E1))"""
    out = []
    for e in eps:
        f = Ep(); f.__dict__.update(e.__dict__); sh = shifts.get(e.seg, 0)
        f.lo = wrap(e.lo + sh); f.ion = wrap(e.ion + sh); out.append(f)
    return out
def run_w(X, eps, shifts, use_local=True, ring=False): return run(X, wrapped_eps(eps, shifts), {}, use_local, ring)

PRIM = [e for e in EPS if e.seg in ("N", "C")]; SEC = [e for e in EPS if e.seg == "S"]
def two_sided(null, obs): return float((np.abs(null) >= abs(obs) - 1e-12).mean())
def full_test(X, eps, label, ndraw=NDRAW, use_local=True, do_ring=True):
    s_obs = run(X, eps, {}, use_local); K, N, A = stats(s_obs)
    r_obs = run(X, eps, {}, use_local, ring=True) if do_ring else np.full(len(eps), np.nan)
    p2_obs = np.nanmean(s_obs - r_obs) if do_ring and np.isfinite(r_obs).any() else np.nan
    dev = {"A": [], "B": []}; amp = {"A": [], "B": []}; p2n = {"A": [], "B": []}
    for nm, gen in (("A", nullA), ("B", nullB)):
        for _ in range(ndraw):
            sh = gen(); s = run_w(X, eps, sh, use_local); k, n, a = stats(s); dev[nm].append(k - n / 2); amp[nm].append(a)
            if do_ring: r = run_w(X, eps, sh, use_local, ring=True); p2n[nm].append(np.nanmean(s - r))
    res = dict(K=K, N=N, mean=A, p2=p2_obs)
    for nm in ("A", "B"):
        res[f"pK_{nm}"] = two_sided(np.array(dev[nm]), K - N / 2); res[f"pmean_{nm}"] = two_sided(np.nan_to_num(np.array(amp[nm])), A)
        res[f"sd_{nm}"] = float(np.nanstd(amp[nm])); res[f"pP2_{nm}"] = two_sided(np.nan_to_num(np.array(p2n[nm])), p2_obs) if do_ring else np.nan
    print(f"{label}: K={K}/{N} (decrease) mean {A:+.3f}% | Null A: pK {res['pK_A']:.3f} pmean {res['pmean_A']:.3f} (sd {res['sd_A']:.3f}) | Null B: pK {res['pK_B']:.3f} pmean {res['pmean_B']:.3f}"
          + (f" | P2 (fp-ring) {p2_obs:+.3f}% pA {res['pP2_A']:.3f} pB {res['pP2_B']:.3f}" if do_ring else ""))
    return res, s_obs, r_obs

# ---------- step-3 stress acceptance (raw-amplitude windows, local orientation) ----------
print("\nSTEP 3 stress acceptance (footprint-median resolved shear change, ON [0,+45] - PRE [-90,-30] of A_raw, local strike/dip):")
neg_frac = {}
for k, off in enumerate(zg["offsets"]):
    Ek = zg["E"][k].astype(float); MU = zg["mu"][k]; LAM = zg["lam"][k]; negs = []
    for e in PRIM:
        on = np.clip(np.arange(e.ion, e.ion + 46), 0, G - 1); pre = np.clip(np.arange(e.ion - 90, e.ion - 29), 0, G - 1)
        dA = A_raw[on].mean(0) - A_raw[pre].mean(0)
        eps_ = np.einsum("obk,k->ob", Ek[e.fp], dA)[None]; sig = st.hooke(eps_, MU[e.fp][None], LAM[e.fp][None])
        sh = np.array([float(st.coulomb_stress(sig[:, j:j + 1], float(zg["strike"][c]), float(zg["dip"][c]), 90.0, friction=0.4)[1][0, 0]) for j, c in enumerate(e.fp)])
        negs.append(np.median(sh) < 0)
    neg_frac[float(off)] = np.mean(negs); print(f"  {off:4.1f} km: footprint-median shear change negative in {100*np.mean(negs):.0f}% of {len(PRIM)} N+C episodes")
stress_ok = all(v >= 0.8 for v in neg_frac.values()); print(f"  -> stress hypotheses {'TESTABLE' if stress_ok else ('MARGINAL' if min(neg_frac.values()) >= 0.6 else 'UNTESTABLE (slip-rate only)')}")

# ---------- primaries ----------
print("\nPRIMARY P1/P2 (N+C, local-onset aligned, deseasoned):")
R1, s1, r1 = full_test(Yd, PRIM, "  P1/P2 N+C")
print("  leave-one-segment-out:")
RN, _, _ = full_test(Yd, [e for e in PRIM if e.seg == "N"], "    N only", ndraw=1000)
RC, _, _ = full_test(Yd, [e for e in PRIM if e.seg == "C"], "    C only", ndraw=1000)
print("  sensitivities:")
RS, _, _ = full_test(Yd, SEC, "    S (secondary)", ndraw=1000)
RE, _, _ = full_test(Yd, PRIM, "    episode-onset aligned", ndraw=1000, use_local=False)
RR, _, _ = full_test(Y, PRIM, "    raw (undeseasoned)", ndraw=1000)
STRICT = [e for e in PRIM if okep.set_index("key").loc[e.key, "ok_strict"]]
RT, _, _ = full_test(Yd, STRICT, f"    strict catalog ({len(STRICT)})", ndraw=1000, do_ring=False) if len(STRICT) >= 3 else (None, None, None)
# controls
print("  controls (same contrast at true episodes vs Null A):")
ctl = {}
for nm, X in (("cc", CC), ("rows/cell-day", NR)):
    s = run(X, PRIM, {}); K, N, A = stats(s); dev = []
    for _ in range(NCTL):
        k, n, a = stats(run_w(X, PRIM, nullA())); dev.append(a)
    ctl[nm] = two_sided(np.nan_to_num(np.array(dev)), A); print(f"    {nm}: mean contrast {A:+.4f} ({K}/{N} neg) p={ctl[nm]:.3f}")
def sites_contrast(eps, shifts):
    out = []
    for e in wrapped_eps(eps, shifts):
        if not e.stas: out.append(np.nan); continue
        y = SITES[e.stas].mean(axis=1).values; idx = e.ion + LAGS; bad = (idx < 0) | (idx >= G); v = y[np.clip(idx, 0, G - 1)]; v[bad] = np.nan
        out.append(contrast(v))
    return np.array(out)
s = sites_contrast(PRIM, {}); K, N, A = stats(s); dev = [stats(sites_contrast(PRIM, nullA()))[2] for _ in range(NCTL)]
ctl["SITES"] = two_sided(np.nan_to_num(np.array(dev)), A); print(f"    SITES (footprint stations, near-station dv/v): mean contrast {A:+.4f} ({K}/{N} neg) p={ctl['SITES']:.3f}")
doyR = {s_: float(np.abs(np.mean(np.exp(2j * np.pi * np.array([e.doy for e in EPS if e.seg == s_]) / 365.25)))) for s_ in SEGLAT if any(e.seg == s_ for e in EPS)}
print(f"    day-of-year R: {doyR}")

# ---------- P3: strain functional ----------
lagset = np.arange(-90, 151, 15); rows = []
for e in PRIM:
    for L in lagset:
        t = e.ion + L
        if t < 30 or t >= G: continue
        dy = np.nanmean(Yd[t, e.fp] - Yd[t - 30, e.fp])
        if not np.isfinite(dy): continue
        de = np.einsum("obk,k->ob", E10[e.fp], dAt30[t]).mean(0) * 1e8
        rows.append([dy, *de, np.nanmean(inc30[t, e.fp]) * 1e3, e.seg, L])
P3 = pd.DataFrame(rows, columns=["dy", "exx", "eyy", "ezz", "exy", "exz", "eyz", "slip", "seg", "lag"])
def r2(df):
    Xm = np.column_stack([np.ones(len(df)), df[["exx", "eyy", "ezz", "exy", "exz", "eyz", "slip"]].values]); y = df.dy.values
    b = np.linalg.lstsq(Xm, y, rcond=None)[0]; return 1 - ((y - Xm @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum(), b
R2_obs, beta = r2(P3)
null_r2 = []
for _ in range(1000):
    sh = nullA(); rr = []
    for e in wrapped_eps(PRIM, sh):
        for L in lagset:
            t = e.ion + L
            if t < 30 or t >= G: continue
            dy = np.nanmean(Yd[t, e.fp] - Yd[t - 30, e.fp])
            if not np.isfinite(dy): continue
            de = np.einsum("obk,k->ob", E10[e.fp], dAt30[e.ion - sh.get(e.seg, 0) + L if False else t]).mean(0) * 1e8
            rr.append([dy, *de, np.nanmean(inc30[t, e.fp]) * 1e3])
    d = pd.DataFrame(rr, columns=["dy", "exx", "eyy", "ezz", "exy", "exz", "eyz", "slip"]); null_r2.append(r2(d)[0])
pP3 = float((np.array(null_r2) >= R2_obs).mean())
print(f"\nP3 strain functional (N+C, {len(P3)} episode-lag samples): R^2 = {R2_obs:.3f}, Null A p = {pP3:.3f} (null median {np.median(null_r2):.3f}); "
      f"coefficients (per 1e-8 strain, per mm slip): {np.round(beta[1:], 4)}")

# ---------- exploratory lag composite with FWER ----------
comp = np.nanmean([aligned(Yd, e) for e in PRIM], 0)
def lagstat(eps, sh):
    ys = np.array([aligned(Yd, e) for e in wrapped_eps(eps, sh)]); c = np.nanmean(ys, 0)
    base = np.nanmean(c[L0 - 90:L0 - 30]); return c - base
obs_l = lagstat(PRIM, {}); mx = []
for _ in range(500): mx.append(np.nanmax(np.abs(lagstat(PRIM, nullA()))))
fwer = float((np.array(mx) >= np.nanmax(np.abs(obs_l))).mean()); ipk = int(np.nanargmax(np.abs(obs_l)))
print(f"EXPLORATORY lag composite (N+C, rel. PRE mean): extreme {obs_l[ipk]:+.3f}% at lag {LAGS[ipk]:+d} d, max-statistic FWER p = {fwer:.3f}")

# ---------- verdict ----------
det = (R1["pK_A"] < 0.01 and R1["pK_B"] < 0.01 and np.sign(RN["mean"]) == np.sign(RC["mean"]) == np.sign(R1["mean"])
       and RN["pK_A"] < 0.05 and RC["pK_A"] < 0.05 and all(v >= 0.05 for v in ctl.values()) and max(doyR.get("N", 0), doyR.get("C", 0)) < 0.3)
bound = float(np.nanpercentile(np.abs(np.array([R1["sd_A"]])), 50)) * 1.96
verdict = "DETECTION" if det else ("BOUND" if (R1["pK_A"] >= 0.05 and R1["pmean_A"] >= 0.05) else "AMBIGUOUS")
print(f"\nVERDICT (pre-registered): {verdict} | P1 mean {R1['mean']:+.3f}% (Null-A sd {R1['sd_A']:.3f}%) -> 95% bound on |footprint response| ~ {abs(R1['mean']) + 1.96 * R1['sd_A']:.3f}% (estimator units, uncalibrated)")
pd.DataFrame(dict(key=[e.key for e in PRIM], seg=[e.seg for e in PRIM], contrast=s1, ring=r1, n_fp=[len(e.fp) for e in PRIM],
                  n_ring=[len(e.ring) for e in PRIM])).to_csv(f"{SD}/primary_episode_table.csv", index=False)
np.savez_compressed(f"{SD}/primary_test.npz", R1=R1, RN=RN, RC=RC, RS=RS, RE=RE, RR=RR, ctl=ctl, doyR=doyR, P3_R2=R2_obs, P3_p=pP3, P3_beta=beta,
                    comp=comp, lags=LAGS, obs_lag=obs_l, fwer=fwer, verdict=verdict, stress_neg_frac=neg_frac, product=PRODUCT)
# figure
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
ys = np.array([aligned(Yd, e) for e in PRIM]); c = np.nanmean(ys, 0); se = np.nanstd(ys, 0) / np.sqrt(np.isfinite(ys).sum(0))
ax[0].fill_between(LAGS, c - se, c + se, color="#bfdbfe"); ax[0].plot(LAGS, c, color="#1d4ed8", lw=1.6); ax[0].axhline(0, color="k", lw=.5)
ax[0].axvspan(*PRE, color="#94a3b8", alpha=.15); ax[0].axvspan(*ON, color="#f59e0b", alpha=.15)
sl = np.nanmean([aligned(inc30 * 1e3, e) for e in PRIM], 0); a2 = ax[0].twinx(); a2.plot(LAGS, sl, color="#f59e0b", lw=1.2); a2.set_ylabel("footprint slip in window (mm)", color="#b45309", fontsize=8)
ax[0].set_title(f"N+C composite, local-onset aligned, deseasoned (N={len(PRIM)})\nP1 K={R1['K']}/{R1['N']} mean {R1['mean']:+.3f}%  pA {R1['pK_A']:.2f} pB {R1['pK_B']:.2f}", fontsize=9)
ax[0].set_xlabel("days from local slip onset"); ax[0].set_ylabel("δβ/β (%)")
ax[1].bar(range(len(PRIM)), s1, color=["#1d4ed8" if e.seg == "N" else "#059669" for e in PRIM]); ax[1].axhline(0, color="k", lw=.5)
ax[1].set_xticks(range(len(PRIM))); ax[1].set_xticklabels([e.key[2:9] for e in PRIM], rotation=90, fontsize=6); ax[1].set_title("per-episode footprint contrast ON−PRE (blue N, green C)", fontsize=9); ax[1].set_ylabel("%")
ax[2].plot(LAGS, obs_l, color="#1d4ed8"); ax[2].axhline(0, color="k", lw=.5); ax[2].axhline(np.nanpercentile(mx, 95), color="r", ls="--", lw=.8); ax[2].axhline(-np.nanpercentile(mx, 95), color="r", ls="--", lw=.8)
ax[2].set_title(f"exploratory lag scan rel. PRE (red = Null-A max-statistic 95%)\nextreme {obs_l[ipk]:+.3f}% at {LAGS[ipk]:+d} d, FWER p={fwer:.2f}", fontsize=9); ax[2].set_xlabel("lag (d)")
fig.suptitle(f"Modeled slow-slip episodes vs interface δβ/β (2-4 s joint inversion) — {verdict}", fontsize=10)
fig.savefig("figures/sse_primary_composite.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_primary_composite.png, primary_test.npz, primary_episode_table.csv")
