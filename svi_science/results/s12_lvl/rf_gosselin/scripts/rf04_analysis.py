#!/usr/bin/env python
"""rf04_analysis: Gosselin et al. 2020 Fig. 2 analysis on our receiver functions.

Steps per station (paper Materials & Methods):
 1. RFs of the chosen period with SNR(Z) > 5 dB; amplitudes normalised by the P-mode RF at zero lag.
 2. time -> depth (Ps moveout, 1-D SVI-average model), 0-60 km, dz 0.25 km.
 3. harmonic decomposition k = 0,1,2 jointly on SV + SH (Eq. 1), least squares (SVD), all RFs of the period;
    predicted SV and SH per event; CC(observed, predicted) over 0-60 km; MAD 2.5 outliers removed.
 4. LVL top (negative) and bottom (positive) Ps depths picked on the k=0 term in the slab window; per-RF amplitude
    (mean over +-0.5 km) -> dlnVs with the linearised Aki-Richards P->SV coefficient at slab-normal incidence
    (Slab2 strike/dip at the station), dVs = dlnVs * 4.6 km/s; cubic-in-time trend removed; MAD 2.5 outliers.
 5. ETS bins pre [mid-W, mid) / post [mid, mid+W), W = 15..150 d; Welch t test; scrambled control = RF dates shuffled
    among events (without replacement), NBOOT draws, 95 % envelope.
Outputs under rf_gosselin/{tables,figures}/<tag>/
Usage: python -P rf04_analysis.py --tag repro --t0 1993-01-01 --t1 2019-06-01 --ets SVI_ALBH --band 05 --stations CN.PGC,...
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
W = f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin"
RFD = f"{ROOT}/data/rf_teleseis/rf"
BINS = [15, 30, 45, 60, 90, 120, 150]
ZMAX, DZ = 60.0, 0.25
# 1-D southern Vancouver Island average (after Ramachandran et al. 2006), Vp/Vs 1.75 except LVL-free average
LAY_Z = np.array([0, 5, 15, 25, 35, 45]); LAY_VP = np.array([5.8, 6.3, 6.7, 6.9, 7.2, 7.8]); VPVS = 1.75
BETA_SCALE = 4.6


def vel_at(z):
    i = np.searchsorted(LAY_Z, z, side="right") - 1
    return LAY_VP[i], LAY_VP[i] / VPVS


def tps_of_z(p, zz):
    a, b = vel_at(zz)
    dt = (np.sqrt(np.maximum(1 / b ** 2 - p ** 2, 0)) - np.sqrt(np.maximum(1 / a ** 2 - p ** 2, 0))) * DZ
    return np.concatenate([[0], np.cumsum(dt)[:-1]])


def slab_geom(lat, lon, r_km=25.0):
    f = f"{ROOT}/data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz"
    S = pd.read_csv(f, header=None, names=["lon", "lat", "dep"]).dropna()
    S["lon"] = np.where(S.lon > 180, S.lon - 360, S.lon)
    x = (S.lon - lon) * 111.19 * np.cos(np.deg2rad(lat)); y = (S.lat - lat) * 111.19
    m = np.hypot(x, y) < r_km
    G = np.c_[np.ones(m.sum()), x[m], y[m]]
    c = np.linalg.lstsq(G, -S.dep[m].values, rcond=None)[0]  # depth positive down; slab2 dep negative
    depth = c[0]; gx, gy = c[1], c[2]  # d(depth)/dx, d(depth)/dy
    dip = np.degrees(np.arctan(np.hypot(gx, gy)))
    dipdir = (np.degrees(np.arctan2(gx, gy)) + 360) % 360  # azimuth of steepest descent
    return float(depth), float(dip), float(dipdir)


def load(key, band, t0, t1):
    z = np.load(f"{RFD}/{key}_rf.npz"); M = pd.read_csv(f"{RFD}/{key}_meta.csv")
    lag = z["lag"]; i0 = int(np.argmin(np.abs(lag)))
    P0 = z["P" + band][:, i0]
    V = z["V" + band] / P0[:, None]; H = z["H" + band] / P0[:, None]
    M["t"] = pd.to_datetime(M["tP"]).astype("datetime64[ns]")
    k = (M.snr > 5).values & (M.t >= pd.Timestamp(t0)).values & (M.t < pd.Timestamp(t1)).values & (P0 > 0)
    return lag, V[k], H[k], M[k].reset_index(drop=True)


def to_depth(lag, X, slow):
    zz = np.arange(0, ZMAX + DZ / 2, DZ)
    out = np.empty((len(X), len(zz)))
    for i, (x, p) in enumerate(zip(X, slow)):
        out[i] = np.interp(tps_of_z(p, zz), lag, x)
    return zz, out


def harmonics(SV, SH, baz):
    """Eq. 1: SV rows [1, cos f, sin f, cos 2f, sin 2f]; SH rows [0, cos(f+pi/2), sin(f+pi/2), cos(2f+pi/2), sin(2f+pi/2)]
    (rfpy dcomp_fix_azim with azim = 0). Returns coefficients (5, nz) and predictions."""
    f = np.deg2rad(baz)
    Gv = np.c_[np.ones_like(f), np.cos(f), np.sin(f), np.cos(2 * f), np.sin(2 * f)]
    Gh = np.c_[np.zeros_like(f), np.cos(f + np.pi / 2), np.sin(f + np.pi / 2),
               np.cos(2 * (f + np.pi / 4)), np.sin(2 * (f + np.pi / 4))]
    G = np.r_[Gv, Gh]
    C = np.linalg.lstsq(G, np.r_[SV, SH], rcond=1e-3)[0]
    return C, Gv @ C, Gh @ C


def rowcc(a, b):
    a = a - a.mean(1, keepdims=True); b = b - b.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))


def mad_keep(x, thr=2.5):
    med = np.nanmedian(x); mad = 1.4826 * np.nanmedian(np.abs(x - med))
    return np.abs(x - med) <= thr * mad


def ps_coef(p_par, a, b):
    """Linearised Aki-Richards P->SV coefficient per unit dlnVs (density term dropped, as in the paper)."""
    ci = np.sqrt(np.maximum(1 - (p_par * a) ** 2, 0)); cj = np.sqrt(np.maximum(1 - (p_par * b) ** 2, 0))
    return (p_par * a / (2 * cj)) * (4 * b * b * ci * cj / (a * b) - 4 * b * b * p_par ** 2)


def slab_parallel_slowness(p, baz, dip, dipdir, a):
    """Magnitude of the interface-parallel slowness of the incident P (velocity a above the interface)."""
    prop = np.deg2rad((baz + 180) % 360)  # propagation azimuth
    pz = np.sqrt(np.maximum(1 / a ** 2 - p ** 2, 0))  # upward
    s = np.c_[p * np.sin(prop), p * np.cos(prop), pz * np.ones_like(p)]  # E, N, Up
    d = np.deg2rad(dip); dd = np.deg2rad(dipdir)
    # upward normal of a plane dipping toward dipdir tilts toward dipdir: n = (sin d sin dd, sin d cos dd, cos d)
    n = np.array([np.sin(d) * np.sin(dd), np.sin(d) * np.cos(dd), np.cos(d)])
    sn = s @ n
    return np.sqrt(np.maximum((s * s).sum(1) - sn ** 2, 0))


def ets_mids(name, t0, t1):
    E = pd.read_csv(f"{W}/tables/ets_mid_{name}.csv")
    m = pd.to_datetime(E["mid"])
    return np.sort(m[(m >= pd.Timestamp(t0)) & (m < pd.Timestamp(t1))].values.astype("datetime64[ns]"))


def bin_masks(t, mids, Wd):
    """pre/post membership (events in any episode's pre/post bin); rel = days from the nearest midpoint."""
    tt = t.values.astype("datetime64[ns]")
    rel = (tt[:, None] - mids[None, :]) / np.timedelta64(1, "D")
    j = np.argmin(np.abs(rel), axis=1)
    r = rel[np.arange(len(tt)), j]
    return (r >= -Wd) & (r < 0), (r >= 0) & (r < Wd), r


def delta_stats(x, pre, post):
    a, b = x[pre & np.isfinite(x)], x[post & np.isfinite(x)]
    if len(a) < 3 or len(b) < 3:
        return dict(n_pre=len(a), n_post=len(b), d=np.nan, p=np.nan, m_pre=np.nan, m_post=np.nan,
                    se_pre=np.nan, se_post=np.nan)
    t = stats.ttest_ind(b, a, equal_var=False)
    return dict(n_pre=len(a), n_post=len(b), d=b.mean() - a.mean(), p=float(t.pvalue), m_pre=a.mean(),
                m_post=b.mean(), se_pre=a.std(ddof=1) / np.sqrt(len(a)), se_post=b.std(ddof=1) / np.sqrt(len(b)))


def scrambled(x, masks, nboot, rng, batch=500):
    """Scrambled control: values shuffled among event dates (without replacement); returns (nboot, nbins)."""
    ok = np.isfinite(x); xv = x[ok]
    pre = np.array([m[0][ok] for m in masks]); post = np.array([m[1][ok] for m in masks])
    out = np.empty((nboot, len(masks)))
    for s in range(0, nboot, batch):
        nb = min(batch, nboot - s)
        perm = np.argsort(rng.random((nb, len(xv))), axis=1)
        X = xv[perm]
        for k in range(len(masks)):
            npre, npost = pre[k].sum(), post[k].sum()
            if npre < 3 or npost < 3:
                out[s:s + nb, k] = np.nan; continue
            out[s:s + nb, k] = X[:, post[k]].mean(1) - X[:, pre[k]].mean(1)
    return out


def analyse(key, a):
    rng = np.random.default_rng(20261013)
    lag, V, H, M = load(key, a.band, a.t0, a.t1)
    if len(M) < 30:
        print(key, "too few RFs", len(M)); return None
    zz, Vz = to_depth(lag, V, M.slow.values)
    _, Hz = to_depth(lag, H, M.slow.values)
    C, Vp_, Hp_ = harmonics(Vz, Hz, M.baz.values)
    M["cc_sv"] = rowcc(Vz, Vp_); M["cc_sh"] = rowcc(Hz, Hp_)
    # LVL picks on k=0 term, in slab window
    S = pd.read_csv(f"{ROOT}/data/rf_teleseis/stations.csv").set_index("sta").loc[key.split(".")[1]]
    sdep, dip, dipdir = slab_geom(S.lat, S.lon)
    A0 = C[0]
    win = (zz >= max(15, sdep - 15)) & (zz <= min(ZMAX - 2, sdep + 12))
    zi = np.where(win)[0]
    itop = zi[np.argmin(A0[zi])]
    below = np.where((zz >= zz[itop] + 2) & (zz <= min(zz[itop] + 15, ZMAX - 1)))[0]
    ibot = below[np.argmax(A0[below])] if len(below) else itop
    ztop, zbot = zz[itop], zz[ibot]
    def amp(Xz, iz):
        m = np.abs(zz - zz[iz]) <= 0.5
        return Xz[:, m].mean(1)
    M["amp_top"] = amp(Vz, itop); M["amp_bot"] = amp(Vz, ibot)
    vp_a, vs_a = vel_at(max(zz[itop] - 2, 0))
    pp = slab_parallel_slowness(M.slow.values, M.baz.values, dip, dipdir, vp_a)
    coef = ps_coef(pp, vp_a, vs_a)
    coef = np.where(np.abs(coef) < 0.02, np.nan, coef)
    tyr = ((M.t - M.t.min()) / pd.Timedelta(days=365.25)).values
    for nm in ("top", "bot"):
        dl = M[f"amp_{nm}"].values / coef
        dv = dl * BETA_SCALE
        ok = np.isfinite(dv)
        c3 = np.polyfit(tyr[ok], dv[ok], 3)
        M[f"dvs_{nm}_raw"] = dv
        M[f"dvs_{nm}"] = dv - np.polyval(c3, tyr)
        keep = mad_keep(M[f"dvs_{nm}"].values)
        M.loc[~keep, f"dvs_{nm}"] = np.nan
    for nm in ("cc_sv", "cc_sh"):
        keep = mad_keep(M[nm].values)
        M[f"{nm}_clean"] = np.where(keep, M[nm].values, np.nan)
    mids = ets_mids(a.ets, a.t0, a.t1)
    masks = [bin_masks(M.t, mids, Wd)[:2] for Wd in BINS]
    rel = bin_masks(M.t, mids, 15)[2]
    M["rel_days"] = rel
    res = dict(station=key, tag=a.tag, band=a.band, t0=a.t0, t1=a.t1, ets=a.ets, n_rf=len(M), n_ets=int(len(mids)),
               slab_depth_km=sdep, dip=dip, dipdir=dipdir, z_top=float(ztop), z_bot=float(zbot),
               A0_top=float(A0[itop]), A0_bot=float(A0[ibot]), median_cc_sv=float(np.nanmedian(M.cc_sv_clean)),
               median_cc_sh=float(np.nanmedian(M.cc_sh_clean)), quantities={})
    boots = {}
    for q in ("cc_sv_clean", "cc_sh_clean", "dvs_top", "dvs_bot"):
        x = M[q].values.astype(float)
        rows = []
        B = scrambled(x, masks, a.nboot, rng)
        boots[q] = B
        for k, Wd in enumerate(BINS):
            d = delta_stats(x, *masks[k])
            b = B[:, k]
            d.update(W=Wd, lo95=float(np.nanpercentile(b, 2.5)), hi95=float(np.nanpercentile(b, 97.5)),
                     pct_below=float(np.nanmean(b <= d["d"])) if np.isfinite(d["d"]) else np.nan,
                     pct_above=float(np.nanmean(b >= d["d"])) if np.isfinite(d["d"]) else np.nan)
            rows.append(d)
        res["quantities"][q] = rows
    od = f"{W}/tables/{a.tag}"; fd = f"{W}/figures/{a.tag}"
    os.makedirs(od, exist_ok=True); os.makedirs(fd, exist_ok=True)
    M.to_csv(f"{od}/{key}_events.csv", index=False)
    np.savez_compressed(f"{od}/{key}_harmonics.npz", z=zz, C=C, boot_cc_sv=boots["cc_sv_clean"][:500],
                        boot_cc_sh=boots["cc_sh_clean"][:500], boot_dvs_top=boots["dvs_top"][:500],
                        boot_dvs_bot=boots["dvs_bot"][:500])
    json.dump(res, open(f"{od}/{key}_result.json", "w"), indent=1, default=float)
    fig2(key, M, res, boots, zz, C, fd)
    return res


def fig2(key, M, res, boots, zz, C, fd):
    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(3, 3)
    q = [("cc_sv_clean", "A  SV: CC(obs, harmonic pred.)"), ("cc_sh_clean", "B  SH: CC"),
         ("dvs_top", "C  dVs, LVL top (km/s)"), ("dvs_bot", "D  dVs, LVL bottom (km/s)")]
    for i, (nm, lab) in enumerate(q):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        r = M.rel_days.values; x = M[nm].values
        m = np.abs(r) < 15
        ax.scatter(r[m], x[m], s=10, c=M.snr.values[m], cmap="viridis", vmin=5, vmax=20, lw=0)
        d = res["quantities"][nm][0]
        for side, xs in (("pre", [-15, 0]), ("post", [0, 15])):
            mu, se = d[f"m_{side}"], d[f"se_{side}"]
            ax.plot(xs, [mu, mu], color="k", lw=2)
            ax.fill_between(xs, mu - se, mu + se, color="tab:red" if side == "post" else "tab:blue", alpha=0.3)
        ax.axvline(0, color="0.5", ls=":")
        ax.set_title(f"{lab}\nΔ = {d['d']:+.3f}, Welch p = {d['p']:.3f}, n = {d['n_pre']}/{d['n_post']}", fontsize=9)
        ax.set_xlabel("days from ETS midpoint (15-day bins)")
    for j, (nm, lab) in enumerate([("cc_sv_clean", "E  ΔCC (SV black, SH grey) vs bin width"),
                                   ("dvs_top", "F  ΔdVs top (black), bottom (red) vs bin width")]):
        ax = fig.add_subplot(gs[2, j])
        B = boots[nm][:500]
        for b in B:
            ax.plot(BINS, b, color="0.85", lw=0.4)
        rows = res["quantities"][nm]
        ax.plot(BINS, [r["d"] for r in rows], "k-o", lw=2)
        ax.plot(BINS, [r["lo95"] for r in rows], "k--", lw=1); ax.plot(BINS, [r["hi95"] for r in rows], "k--", lw=1)
        other = "cc_sh_clean" if nm == "cc_sv_clean" else "dvs_bot"
        ax.plot(BINS, [r["d"] for r in res["quantities"][other]], "-o", color="tab:red" if other == "dvs_bot" else "0.4", lw=1.5)
        ax.axhline(0, color="0.5", lw=0.5); ax.set_xlabel("bin width (days)"); ax.set_title(lab, fontsize=9)
    ax = fig.add_subplot(gs[:2, 2])
    for k, (c, nm) in enumerate(zip(C, ["A", "B1", "B2", "C1", "C2"])):
        ax.plot(c * 4 + k, zz, color="k"); ax.fill_betweenx(zz, k, c * 4 + k, where=c > 0, color="tab:red", alpha=.5)
        ax.fill_betweenx(zz, k, c * 4 + k, where=c < 0, color="tab:blue", alpha=.5)
    ax.axhline(res["z_top"], color="tab:blue", ls=":"); ax.axhline(res["z_bot"], color="tab:red", ls=":")
    ax.set_ylim(ZMAX, 0); ax.set_xticks(range(5)); ax.set_xticklabels(["A", "B1", "B2", "C1", "C2"])
    ax.set_ylabel("depth (km)"); ax.set_title(f"harmonics; slab2 {res['slab_depth_km']:.0f} km, dip {res['dip']:.0f}°", fontsize=9)
    fig.suptitle(f"{key}  {res['tag']}  {res['t0']}–{res['t1']}  band 0.05–0.{res['band'][1]} Hz  ETS {res['ets']} (n={res['n_ets']})  RFs {res['n_rf']}\n"
                 "scrambled control: RF dates shuffled among events (grey = 500 of the draws, dashed = 95 %)", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{fd}/fig2_{key}.png", dpi=110); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="repro"); ap.add_argument("--t0", default="1993-01-01")
    ap.add_argument("--t1", default="2019-06-01"); ap.add_argument("--ets", default="SVI_ALBH")
    ap.add_argument("--band", default="05"); ap.add_argument("--nboot", type=int, default=20000)
    ap.add_argument("--stations", default="CN.PGC,CN.SNB,CN.VGZ,CN.LZB,CN.PFB")
    a = ap.parse_args()
    out = []
    for key in a.stations.split(","):
        if not os.path.exists(f"{RFD}/{key}_rf.npz"):
            print("missing", key); continue
        r = analyse(key, a)
        if r is None:
            continue
        q = r["quantities"]
        out.append(dict(station=key, n_rf=r["n_rf"], n_ets=r["n_ets"], z_top=r["z_top"], z_bot=r["z_bot"],
                        **{f"{nm}_{f}": q[nm][0][f] for nm in q for f in ("d", "p", "pct_below", "pct_above", "n_pre", "n_post")},
                        **{f"{nm}_d150": q[nm][-1]["d"] for nm in q}))
        print(pd.DataFrame([out[-1]]).T.to_string(), flush=True)
    pd.DataFrame(out).to_csv(f"{W}/tables/{a.tag}/summary.csv", index=False)


if __name__ == "__main__":
    main()
