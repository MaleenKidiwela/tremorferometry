#!/usr/bin/env python
"""E4 step 5 — THOMAS-2025 ANALYSIS HARNESS, catalog-agnostic.

Runs on ANY catalog with (time, family_id, lat, lon, depth, magnitude); validated here on Bostock's
own Sept-2012 catalog so the analysis can be separated from the catalog as the bottleneck, and
re-run unchanged on our own catalogs later.

WHAT THOMAS ET AL. (2025, Sci. Adv. adz2832) ACTUALLY SPECIFY (PMC12629176, fetched 2026-09-11):
the paper's front and tidal analysis is QUALITATIVE -- "vertically oriented groups of events closely
spaced in time" read off a space-time plot. Its stated numbers are:
  * main front propagates SOUTHWARD, "~5 km/day" (also quoted "~10 km/day"), ~80 km in 11 days;
  * secondary fronts initiate 1-10 km from the main front, propagate 10-40x faster (example
    ~12 km/hour), last 2-6 hours, median stress drop 8 kPa (~Mw 5.1); four are marked in their
    Fig. 2A and family 141 hosted six of them between 12 and 16 September;
  * early LFEs are "not significantly influenced by the tides" while later secondary fronts are
    "almost always tidally triggered"; ~2 kPa shear stress change is enough.
No detection algorithm, no binning, no completeness cut, no rate-decay functional form and no tidal
method are given. The operational criteria below are therefore OURS, chosen to be as literal a
reading of those numbers as an algorithm allows, and they are printed with every run.

OPERATIONAL CRITERIA (ours; --print-criteria dumps them)
  burst      a family's events with successive gaps < BURST_GAP (600 s) and >= BURST_MIN (4) events
  main front per family, the time of its FIRST event in the episode; robust (Theil-Sen) fit of
             along-strike position vs first-activation time over families with >= MIN_FAM_EV events
             -> speed in km/day. Also reported: the 6-h event-weighted activity centroid track and
             the 6-h leading-edge (10th percentile of active-family position in the propagation
             direction).
  secondary  a sliding WIN (6 h) window; take every family burst STARTING in it; require
  front        >= SF_MIN_FAM (3) distinct families, along-strike span >= SF_SPAN (5 km),
               |Spearman rho(position, burst start time)| >= SF_RHO (0.7),
               |speed| between SF_LO (10x) and SF_HI (40x) the fitted main-front speed,
               duration (last minus first burst start) between 2 and 6 h,
               group centroid 1-10 km BEHIND the main front (on the already-ruptured side).
             Overlapping detections are merged (keep the one with most families).
             A time-shuffled null (each family's burst start times replaced by uniform draws in the
             episode, counts preserved) gives the chance count.
  re-rupture per family: intervals between successive burst starts; fraction of bursts that follow
             another burst of the SAME family within RERUP_MAX (5 d), and the distribution.
  rate decay per family, from its first event: counts in log-spaced bins; Omori fit
             rate(t) = K/(t+c)^p (c = 0.05 d) by weighted least squares on log-log, and an
             exponential fit rate = R0 exp(-t/tau); both reported.
  tides      *** PLACEHOLDER ***  No Cascadia tidal-stress model and no tide package is on disk
             (the only tidal code, catalogs/ducellier_catalog/src/, reads a California-only .mat).
             The phase here is from the EQUILIBRIUM solid-Earth tide: the degree-2 tidal potential
             of Sun and Moon from low-precision ephemerides, converted to an equilibrium tide
             height, band-passed to 0.8-2.2 cpd and phase-extracted by Hilbert transform. This is a
             proxy for the VERTICAL tide, not for the shear stress resolved on the plate interface,
             and its amplitude is meaningless. Schuster test p-values and phase histograms are
             computed from it and are labelled PLACEHOLDER everywhere.

DECIMATION LEG (Merlin vet, prescription 4): every metric is re-run after a magnitude cut
Mc in {none, 1.7, 1.8, 1.9, 2.0, 2.1}, so the output is the COMPLETENESS EACH METRIC NEEDS.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
RES = f"{ROOT}/lfe_catalog_repair/results"
SLAB2 = f"{ROOT}/data/cas_slab2_input_04-18.csv"
R_KM_DEG = 111.195

CRIT = dict(BURST_GAP_S=600.0, BURST_MIN=4, MIN_FAM_EV=20, BIN_H=(1, 6),
            SF_WIN_H=6.0, SF_STEP_H=1.0, SF_MIN_FAM=3, SF_SPAN_KM=5.0, SF_RHO=0.7,
            SF_SPEED_LO=10.0, SF_SPEED_HI=40.0, SF_DUR_H=(2.0, 6.0), SF_BEHIND_KM=(1.0, 10.0),
            SF_ABS_KMDAY=(50.0, 400.0), VMAIN_FALLBACK_KMDAY=5.0, SF_RATE_ANOM=3.0,
            RERUP_MAX_D=5.0, OMORI_C_D=0.05, TIDE="PLACEHOLDER equilibrium solid-Earth tide")


# ---------------------------------------------------------------- geometry
MESH = "/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt"


def slab_interp():
    """Interface depth interpolator. The gridded MESH first; the scattered Slab2 INPUT file only as
    a fallback -- cas_slab2_input_04-18.csv holds the earthquake/receiver-function observations that
    were INPUT to Slab2, not the fitted surface, and its local gradient is noise (it gave an
    along-strike azimuth of 86 deg at southern Vancouver Island, i.e. E-W, against the correct
    -43 deg = N317E from the mesh)."""
    try:
        sys.path.insert(0, f"{ROOT}/src")
        from sse_strain import load_fault_ascii
        from sse_strain.megathrust import slab_depth_interpolator
        return slab_depth_interpolator(load_fault_ascii(MESH)), "cascadia_mesh.txt"
    except Exception as exc:
        from scipy.interpolate import LinearNDInterpolator
        s2 = pd.read_csv(SLAB2, low_memory=False)
        s2 = s2[pd.to_numeric(s2.depth, errors="coerce").notna()]
        return (LinearNDInterpolator(np.c_[s2.lon.astype(float), s2.lat.astype(float)],
                                     s2.depth.astype(float).to_numpy()),
                f"slab2 scattered input (FALLBACK: {exc})")


def strike_frame(lat0, lon0):
    """Along-strike (north-positive) and down-dip unit vectors from the interface depth gradient."""
    f, src = slab_interp()
    h = 0.05
    cs = np.cos(np.radians(lat0))
    dzdx = float(f(lon0 + h, lat0) - f(lon0 - h, lat0)) / (2 * h * R_KM_DEG * cs)
    dzdy = float(f(lon0, lat0 + h) - f(lon0, lat0 - h)) / (2 * h * R_KM_DEG)
    g = np.array([dzdx, dzdy])
    if not np.all(np.isfinite(g)) or np.linalg.norm(g) == 0:
        # fall back on the mean Cascadia strike (N320E) if the local gradient is undefined
        s = np.array([np.sin(np.radians(-40.0)), np.cos(np.radians(-40.0))])
        return s, np.array([s[1], -s[0]]), "fallback N320E"
    d = g / np.linalg.norm(g)
    s = np.array([-d[1], d[0]])
    if s[1] < 0:
        s = -s
    return s, d, (f"{src}: interface depth {float(f(lon0, lat0)):.1f} km, down-dip azimuth "
                  f"{np.degrees(np.arctan2(d[0], d[1])):.1f} deg, along-strike azimuth "
                  f"{np.degrees(np.arctan2(s[0], s[1])):.1f} deg from N "
                  f"(+ = NW along strike, - = SE = Thomas's 'southward')")


# ---------------------------------------------------------------- tides (PLACEHOLDER)
def equilibrium_tide(t_utc_s, lat, lon):
    """Degree-2 equilibrium tide height (m) from low-precision Sun+Moon ephemerides.

    PLACEHOLDER for a real interface shear-stress model -- see the module docstring.
    """
    jd = t_utc_s / 86400.0 + 2440587.5
    T = (jd - 2451545.0) / 36525.0
    d2r = np.pi / 180.0
    # Moon (Meeus low precision)
    Lp = (218.316 + 481267.8813 * T) * d2r
    M = (134.963 + 477198.8676 * T) * d2r
    F = (93.272 + 483202.0175 * T) * d2r
    Mm = (357.529 + 35999.0503 * T) * d2r
    lam_m = Lp + 6.289 * d2r * np.sin(M)
    beta_m = 5.128 * d2r * np.sin(F)
    r_m = 385001.0 - 20905.0 * np.cos(M)
    # Sun
    lam_s = ((280.460 + 36000.771 * T) * d2r + (1.915 * np.sin(Mm) + 0.020 * np.sin(2 * Mm)) * d2r)
    r_s = 149598023.0 * (1 - 0.0167 * np.cos(Mm))
    eps = 23.439 * d2r
    def eq(lam, beta):
        ra = np.arctan2(np.sin(lam) * np.cos(eps) - np.tan(beta) * np.sin(eps), np.cos(lam))
        dec = np.arcsin(np.sin(beta) * np.cos(eps) + np.cos(beta) * np.sin(eps) * np.sin(lam))
        return ra, dec
    ra_m, dec_m = eq(lam_m, beta_m)
    ra_s, dec_s = eq(lam_s, 0.0)
    gmst = (280.46061837 + 360.98564736629 * (jd - 2451545.0)) * d2r
    phi = lat * d2r
    a_e, g_e = 6371.0, 9.81
    out = 0.0
    for ra, dec, r, gm in ((ra_m, dec_m, r_m, 4.9028e3), (ra_s, dec_s, r_s, 1.32712e11)):
        H = gmst + lon * d2r - ra
        cz = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(H)
        out = out + gm * a_e ** 2 / r ** 3 * 0.5 * (3 * cz ** 2 - 1) / g_e * 1e3
    return out                      # metres (scale is meaningless; only the phase is used)


def tidal_phase(t_s, lat, lon, t0, t1):
    """Phase (deg, 0 = tidal maximum) of the band-passed equilibrium tide at each event time."""
    from scipy.signal import butter, sosfiltfilt, hilbert
    fs = 1.0 / 600.0                                     # 10-min grid
    tg = np.arange(t0 - 86400, t1 + 86400, 600.0)
    z = equilibrium_tide(tg, lat, lon)
    sos = butter(4, [0.8 / 86400 / (fs / 2), 2.2 / 86400 / (fs / 2)], btype="band", output="sos")
    zb = sosfiltfilt(sos, z - z.mean())
    ph = np.angle(hilbert(zb))
    return np.degrees(np.interp(t_s, tg, np.unwrap(ph))) % 360.0, tg, zb


def schuster(phases_deg):
    """Schuster test: p = exp(-R^2/n) for n phases; R = length of the vector sum."""
    n = len(phases_deg)
    if n < 5:
        return np.nan, np.nan
    a = np.radians(phases_deg)
    R = np.hypot(np.cos(a).sum(), np.sin(a).sum())
    return float(np.exp(-R * R / n)), float(R / n)


# ---------------------------------------------------------------- analysis
def bursts(t, gap, nmin):
    """Burst start/end/size for one family's sorted event times."""
    if len(t) == 0:
        return np.zeros((0, 3))
    cut = np.where(np.diff(t) > gap)[0]
    st = np.r_[0, cut + 1]
    en = np.r_[cut, len(t) - 1]
    n = en - st + 1
    k = n >= nmin
    return np.c_[t[st[k]], t[en[k]], n[k]]


def theil_sen(x, y, maxpair=200000, rng=None):
    rng = rng or np.random.default_rng(0)
    n = len(x)
    if n < 3:
        return np.nan, np.nan
    if n * (n - 1) // 2 <= maxpair:
        i, j = np.triu_indices(n, 1)
    else:
        i = rng.integers(0, n, maxpair); j = rng.integers(0, n, maxpair)
        k = i != j; i, j = i[k], j[k]
    dx = x[j] - x[i]
    ok = np.abs(dx) > 1e-9
    sl = np.median((y[j][ok] - y[i][ok]) / dx[ok])
    return float(sl), float(np.median(y - sl * x))


def analyse(ev, fam, crit, tag, outdir, vmain_ref=None, rng=None):
    """Every Thomas metric for one (possibly decimated) catalog. Returns a dict of scalars."""
    rng = rng or np.random.default_rng(7)
    out = dict(tag=tag, n_events=len(ev), n_families=ev.family_id.nunique())
    if len(ev) < 50:
        return out
    t0, t1 = ev.t.min(), ev.t.max()
    out["days"] = (t1 - t0) / 86400.0

    # ---- (i) rates
    rates = {}
    for bh in crit["BIN_H"]:
        b = np.floor((ev.t - t0) / (bh * 3600)).astype(int)
        r = ev.groupby([ev.family_id, b]).size().rename("n").reset_index()
        r.columns = ["family_id", "bin", "n"]
        r["bin_h"] = bh
        rates[bh] = r
    R = pd.concat(rates.values(), ignore_index=True)
    R.to_csv(f"{outdir}/e4_thomas_rates_{tag}.csv", index=False)
    out["peak_rate_1h_per_family"] = float(R[R.bin_h == 1].n.max())

    # ---- (ii) main front
    g = ev.groupby("family_id")
    first = g.t.min()
    nev = g.size()
    use = nev[nev >= crit["MIN_FAM_EV"]].index
    F = fam.loc[fam.family_id.isin(use)].set_index("family_id")
    x = F.x_strike.to_numpy()
    tt = (first.loc[F.index].to_numpy() - t0) / 86400.0
    sl, ic = theil_sen(tt, x, rng=rng)                    # km per day of first activation
    out["n_fam_front"] = int(len(x))
    out["front_speed_kmday"] = sl
    out["front_span_km"] = float(np.ptp(x)) if len(x) else np.nan
    out["front_duration_d"] = float(np.ptp(tt)) if len(tt) else np.nan
    rho = pd.Series(tt).corr(pd.Series(x), method="spearman")
    out["front_rho"] = float(rho)
    # least-squares for comparison + 6-h tracks
    if len(x) >= 3:
        out["front_speed_ls_kmday"] = float(np.polyfit(tt, x, 1)[0])
    fam_x = fam.set_index("family_id").x_strike
    ev = ev.assign(x_strike=ev.family_id.map(fam_x))
    b6 = np.floor((ev.t - t0) / (6 * 3600)).astype(int)
    trk = ev.groupby(b6).agg(t_mid=("t", "mean"), x_mean=("x_strike", "mean"),
                             x_lead=("x_strike", lambda v: np.percentile(v, 10 if sl < 0 else 90)),
                             n=("t", "size")).reset_index()
    trk["day"] = (trk.t_mid - t0) / 86400.0
    trk.to_csv(f"{outdir}/e4_thomas_fronttrack_{tag}.csv", index=False)
    k = trk.n >= 10
    if k.sum() >= 3:
        out["front_speed_centroid_kmday"] = theil_sen(trk.day[k].to_numpy(),
                                                      trk.x_mean[k].to_numpy(), rng=rng)[0]
        out["front_speed_lead_kmday"] = theil_sen(trk.day[k].to_numpy(),
                                                  trk.x_lead[k].to_numpy(), rng=rng)[0]

    # ---- bursts (shared by secondary fronts and re-rupture)
    B = []
    for fid, gg in ev.groupby("family_id"):
        bb = bursts(np.sort(gg.t.to_numpy()), crit["BURST_GAP_S"], crit["BURST_MIN"])
        for r in bb:
            B.append(dict(family_id=fid, t_start=r[0], t_end=r[1], n=int(r[2])))
    B = pd.DataFrame(B)
    out["n_bursts"] = len(B)
    if len(B):
        B["x_strike"] = B.family_id.map(fam_x)
        B.to_csv(f"{outdir}/e4_thomas_bursts_{tag}.csv", index=False)

    # ---- (iii) secondary fronts
    # the main-front speed the secondary-front ratio test is referenced to: the fitted value if it
    # is a real front (|v| >= 1 km/day AND |rho| >= 0.3), otherwise Thomas's own ~5 km/day, stated.
    vfit = vmain_ref if vmain_ref is not None else out.get("front_speed_centroid_kmday", sl)
    if not (np.isfinite(vfit) and abs(vfit) >= 1.0):
        vmain = -crit["VMAIN_FALLBACK_KMDAY"] if (sl <= 0) else crit["VMAIN_FALLBACK_KMDAY"]
        out["vmain_used_kmday"] = vmain
        out["vmain_source"] = "Thomas fallback 5 km/day"
    else:
        vmain = vfit
        out["vmain_used_kmday"] = vmain
        out["vmain_source"] = "fitted"
    span_d = (t1 - t0) / 86400.0
    bgr = (B.groupby("family_id").size() / span_d).to_dict() if len(B) else {}
    sf = detect_secondary(B, vmain, sl, ic, t0, crit, bg_rate=bgr)
    sf["tag"] = tag
    sf.to_csv(f"{outdir}/e4_thomas_secondary_{tag}.csv", index=False)
    out["n_secondary"] = len(sf)
    if len(sf):
        out["sf_speed_kmh_median"] = float(np.median(np.abs(sf.speed_kmday)) / 24.0)
        out["sf_speed_ratio_median"] = float(np.median(np.abs(sf.speed_kmday / vmain)))
        out["sf_duration_h_median"] = float(np.median(sf.duration_h))
        out["sf_behind_km_median"] = float(np.median(sf.behind_km))
    # time-shuffled null: each family's burst starts redrawn uniformly in the episode
    if len(B):
        # NULL (fixed 2026-09-11): a CIRCULAR per-family time shift. Each family's whole burst
        # sequence is shifted by its own random amount and wrapped inside the episode, so the
        # within-family burst clustering, the burst count and the spatial pattern are preserved
        # exactly while the inter-family timing that a secondary front consists of is destroyed.
        # (The first version redrew every burst start uniformly, which SPREADS each family's
        # clustered bursts and therefore INFLATES the number of distinct families per 6-h window --
        # it gave a null ABOVE the observed count, which is a broken null, not a negative result.)
        nulls = []
        span = t1 - t0
        for s_ in range(5):
            Bn = B.copy()
            sh = {f: rng.uniform(0, span) for f in Bn.family_id.unique()}
            Bn["t_start"] = t0 + ((Bn.t_start - t0 + Bn.family_id.map(sh)) % span)
            nulls.append(len(detect_secondary(Bn, vmain, sl, ic, t0, crit, bg_rate=bgr)))
        out["n_secondary_null"] = float(np.mean(nulls))
        out["n_secondary_null_sd"] = float(np.std(nulls))

    # ---- (iv) re-rupture
    if len(B):
        B = B.sort_values(["family_id", "t_start"])
        dt = B.groupby("family_id").t_start.diff() / 86400.0
        out["rerupture_frac_le5d"] = float((dt <= crit["RERUP_MAX_D"]).mean())
        out["rerupture_interval_d_median"] = float(dt.median())
        out["rerupture_n_families"] = int((B.groupby("family_id").size() > 1).sum())

    # ---- (v) per-family rate decay from first activation
    dec = []
    for fid, gg in ev.groupby("family_id"):
        tt2 = np.sort(gg.t.to_numpy())
        if len(tt2) < crit["MIN_FAM_EV"]:
            continue
        a = (tt2 - tt2[0]) / 86400.0
        edges = np.logspace(np.log10(1.0 / 48), np.log10(max(a.max(), 0.1)), 9)
        cnt, _ = np.histogram(a, bins=edges)
        w = np.diff(edges)
        mid = np.sqrt(edges[:-1] * edges[1:])
        k = cnt > 0
        if k.sum() < 4:
            continue
        p = np.polyfit(np.log10(mid[k] + crit["OMORI_C_D"]), np.log10(cnt[k] / w[k]), 1)
        rate = cnt[k] / w[k]
        pe = np.polyfit(mid[k], np.log(rate), 1)
        dec.append(dict(family_id=fid, n=len(tt2), omori_p=-p[0],
                        tau_d=(-1.0 / pe[0]) if pe[0] < 0 else np.nan,
                        r_first_half=float((a <= a.max() / 2).mean())))
    DEC = pd.DataFrame(dec)
    DEC.to_csv(f"{outdir}/e4_thomas_decay_{tag}.csv", index=False)
    if len(DEC):
        out["omori_p_median"] = float(DEC.omori_p.median())
        out["omori_p_iqr"] = float(DEC.omori_p.quantile(.75) - DEC.omori_p.quantile(.25))
        out["frac_events_first_half_median"] = float(DEC.r_first_half.median())
        out["n_fam_decay"] = len(DEC)

    # ---- (vi) tides (PLACEHOLDER)
    lat0, lon0 = float(fam.lat.mean()), float(fam.lon.mean())
    ph, _, _ = tidal_phase(ev.t.to_numpy(), lat0, lon0, t0, t1)
    ev = ev.assign(tide_phase=ph)
    p_all, R_all = schuster(ph)
    out["tide_schuster_p_ALL_PLACEHOLDER"] = p_all
    out["tide_R_ALL_PLACEHOLDER"] = R_all
    thirds = np.floor((ev.t - t0) / ((t1 - t0) / 3 + 1e-9)).clip(0, 2).astype(int)
    for i in range(3):
        p_, R_ = schuster(ev.tide_phase[thirds == i].to_numpy())
        out[f"tide_R_third{i+1}_PLACEHOLDER"] = R_
        out[f"tide_p_third{i+1}_PLACEHOLDER"] = p_
    # per family, by time since its own first activation (Thomas: sensitivity GROWS through the SSE)
    # R (the mean resultant length) is biased upward at small n, so early and late are compared at
    # a MATCHED sample size: n = min(n_early, n_late), drawn without replacement, 20 repeats.
    tf = []
    for fid, gg in ev.groupby("family_id"):
        if len(gg) < 40:
            continue
        a = ((gg.t - gg.t.min()) / 86400.0).to_numpy()
        pe = gg.tide_phase.to_numpy()[a <= 1.0]
        pl = gg.tide_phase.to_numpy()[a > 1.0]
        nm = min(len(pe), len(pl))
        if nm < 15:
            continue
        for lab, arr in (("early", pe), ("late", pl)):
            Rs, ps = [], []
            for _ in range(20):
                sub = rng.choice(arr, nm, replace=False)
                p_, R_ = schuster(sub)
                Rs.append(R_); ps.append(p_)
            tf.append(dict(family_id=fid, when=lab, n=nm, R=float(np.mean(Rs)),
                           p=float(np.median(ps))))
    TF = pd.DataFrame(tf)
    TF.to_csv(f"{outdir}/e4_thomas_tides_{tag}.csv", index=False)
    if len(TF):
        for lab in ("early", "late"):
            out[f"tide_R_{lab}_median_PLACEHOLDER"] = float(TF[TF.when == lab].R.median())
            out[f"tide_fracp01_{lab}_PLACEHOLDER"] = float((TF[TF.when == lab].p < 0.01).mean())
    ev[["t", "family_id", "mag", "x_strike", "tide_phase"]].to_parquet(
        f"{outdir}/e4_thomas_events_{tag}.parquet", index=False)
    return out


def detect_secondary(B, vmain, sl, ic, t0, crit, bg_rate=None):
    """Sliding-window detection of secondary fronts. See the module docstring for the criteria."""
    if not len(B) or not np.isfinite(vmain) or vmain == 0:
        return pd.DataFrame(columns=["t_start", "t_end", "n_fam", "span_km", "speed_kmday",
                                     "duration_h", "behind_km", "rho"])
    B = B.sort_values("t_start").reset_index(drop=True)
    ts = B.t_start.to_numpy(); xs = B.x_strike.to_numpy(); fid = B.family_id.to_numpy()
    win = crit["SF_WIN_H"] * 3600.0
    step = crit["SF_STEP_H"] * 3600.0
    cands = []
    for w0 in np.arange(ts.min(), ts.max(), step):
        k = (ts >= w0) & (ts < w0 + win)
        if k.sum() < crit["SF_MIN_FAM"]:
            continue
        tw, xw, fw = ts[k], xs[k], fid[k]
        if len(np.unique(fw)) < crit["SF_MIN_FAM"]:
            continue
        if np.ptp(xw) < crit["SF_SPAN_KM"]:
            continue
        dur = (tw.max() - tw.min()) / 3600.0
        if not (crit["SF_DUR_H"][0] <= dur <= crit["SF_DUR_H"][1]):
            continue
        rho = pd.Series(tw).corr(pd.Series(xw), method="spearman")
        if not np.isfinite(rho) or abs(rho) < crit["SF_RHO"]:
            continue
        # Thomas's secondary fronts are "vertically oriented groups of events closely spaced in
        # time": require the participating families' burst rate in the window to exceed their own
        # episode-average burst rate by SF_RATE_ANOM.
        if bg_rate is not None:
            exp = sum(bg_rate.get(f_, 0.0) for f_ in np.unique(fw)) * (win / 86400.0)
            if exp > 0 and len(tw) < crit["SF_RATE_ANOM"] * exp:
                continue
        sp = theil_sen((tw - tw.min()) / 86400.0, xw)[0]        # km/day
        if not np.isfinite(sp):
            continue
        ratio = abs(sp / vmain)
        # Thomas gives the secondary-front speed BOTH as a ratio (10-40x the main front) and in
        # absolute terms (~12 km/h, i.e. ~288 km/day).  Require either, and record which.
        ok_ratio = crit["SF_SPEED_LO"] <= ratio <= crit["SF_SPEED_HI"]
        ok_abs = crit["SF_ABS_KMDAY"][0] <= abs(sp) <= crit["SF_ABS_KMDAY"][1]
        if not (ok_ratio or ok_abs):
            continue
        # position of the MAIN front at this time, from the first-activation fit
        xm = ic + sl * ((tw.mean() - t0) / 86400.0)
        behind = (np.mean(xw) - xm) * (1 if sl < 0 else -1)     # + = on the already-ruptured side
        if not (crit["SF_BEHIND_KM"][0] <= behind <= crit["SF_BEHIND_KM"][1]):
            continue
        cands.append(dict(t_start=tw.min(), t_end=tw.max(), n_fam=int(len(np.unique(fw))),
                          span_km=float(np.ptp(xw)), speed_kmday=float(sp), duration_h=float(dur),
                          behind_km=float(behind), rho=float(rho), speed_ratio=float(ratio),
                          by_ratio=bool(ok_ratio), by_abs=bool(ok_abs)))
    if not cands:
        return pd.DataFrame(columns=["t_start", "t_end", "n_fam", "span_km", "speed_kmday",
                                     "duration_h", "behind_km", "rho"])
    C = pd.DataFrame(cands).sort_values(["n_fam", "span_km"], ascending=False)
    keep, taken = [], []
    for r in C.itertuples():
        if any(not (r.t_end < a or r.t_start > b) for a, b in taken):
            continue
        keep.append(r._asdict()); taken.append((r.t_start, r.t_end))
    return pd.DataFrame(keep).drop(columns=["Index"], errors="ignore").sort_values("t_start")


# ---------------------------------------------------------------- catalog loading
def load_catalog(path, t0, t1, cols):
    if path.endswith(".parquet"):
        c = pd.read_parquet(path)
    else:
        c = pd.read_csv(path)
    if "time" not in c.columns and {"year", "month", "day", "hour", "min", "sec"} <= set(c.columns):
        t = pd.to_datetime(dict(year=c.year, month=c.month, day=c.day, hour=c.hour,
                                minute=c["min"], second=c.sec.astype(int))) \
            + pd.to_timedelta(c.sec % 1.0, unit="s")
        c["time"] = t
    ren = dict(kv.split("=") for kv in cols.split(",")) if cols else {}
    c = c.rename(columns=ren)
    if "mag" not in c.columns and "magnitude" in c.columns:
        c = c.rename(columns={"magnitude": "mag"})
    if "depth" in c.columns and "dep" not in c.columns:
        c = c.rename(columns={"depth": "dep"})
    c["time"] = pd.to_datetime(c["time"]).astype("datetime64[ns]")
    c = c[(c.time >= pd.Timestamp(t0)) & (c.time < pd.Timestamp(t1))].copy()
    if "family_id" not in c.columns:
        key = c.groupby(["lat", "lon", "dep"]).size().reset_index(name="n").reset_index()
        c = c.merge(key[["index", "lat", "lon", "dep"]].rename(columns={"index": "family_id"}),
                    on=["lat", "lon", "dep"], how="left")
    c["t"] = c.time.astype("int64") / 1e9
    return c.sort_values("t").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default=f"{RES}/external/bostock/"
                                         "bostock2015_sloweq_FULL_2003-2013.csv")
    ap.add_argument("--t0", default="2012-08-15")
    ap.add_argument("--t1", default="2012-10-16")
    ap.add_argument("--cols", default="", help="renames, e.g. 'origin_time=time,depth_km=dep'")
    ap.add_argument("--tag", default="bostock2012")
    ap.add_argument("--outdir", default=RES)
    ap.add_argument("--mc", default="none,1.7,1.8,1.9,2.0,2.1")
    a = ap.parse_args()

    print("[criteria] " + json.dumps(CRIT), flush=True)
    ev = load_catalog(a.catalog, a.t0, a.t1, a.cols)
    print(f"[cat] {len(ev)} events, {ev.family_id.nunique()} families, "
          f"{ev.time.min()} .. {ev.time.max()}", flush=True)
    fam = ev.groupby("family_id").agg(lat=("lat", "first"), lon=("lon", "first"),
                                      dep=("dep", "first"), n=("t", "size")).reset_index()
    s, d, how = strike_frame(float(fam.lat.mean()), float(fam.lon.mean()))
    cs = np.cos(np.radians(float(fam.lat.mean())))
    X = (fam.lon - fam.lon.mean()).to_numpy() * R_KM_DEG * cs
    Y = (fam.lat - fam.lat.mean()).to_numpy() * R_KM_DEG
    fam["x_strike"] = X * s[0] + Y * s[1]
    fam["x_dip"] = X * d[0] + Y * d[1]
    print(f"[geom] along-strike unit vector from {how}; family along-strike span "
          f"{fam.x_strike.max() - fam.x_strike.min():.1f} km, down-dip span "
          f"{fam.x_dip.max() - fam.x_dip.min():.1f} km", flush=True)
    fam.to_csv(f"{a.outdir}/e4_thomas_families_{a.tag}.csv", index=False)

    rows = []
    full = None
    for mcs in a.mc.split(","):
        mc = None if mcs == "none" else float(mcs)
        sub = ev if mc is None else ev[ev.mag >= mc]
        tag = f"{a.tag}_mcall" if mc is None else f"{a.tag}_mc{mcs}"
        nf50 = int((sub.groupby('family_id').size() >= 50).sum())
        print(f"\n=== {tag}: {len(sub)} events ({len(sub)/len(ev)*100:.1f} %), "
              f"{sub.family_id.nunique()} families, {nf50} with >= 50 events", flush=True)
        # every decimation uses the SAME main-front speed reference -- the full catalog's 6-h
        # activity-centroid track speed -- so the 10-40x secondary-front band does not move with Mc
        r = analyse(sub.copy(), fam, CRIT, tag, a.outdir,
                    vmain_ref=(full.get("front_speed_centroid_kmday") if full else None))
        r["mc"] = mcs
        r["kept_frac"] = len(sub) / len(ev)
        r["n_fam_ge50"] = nf50
        if full is None:
            full = r
        rows.append(r)
        print("   " + ", ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}"
                                for k, v in r.items()
                                if k not in ("tag",)), flush=True)
    S = pd.DataFrame(rows)
    S.to_csv(f"{a.outdir}/e4_thomas_summary_{a.tag}.csv", index=False)
    with open(f"{a.outdir}/e4_thomas_criteria_{a.tag}.json", "w") as f:
        json.dump(CRIT, f, indent=2)
    print(f"\n[write] {a.outdir}/e4_thomas_summary_{a.tag}.csv")


if __name__ == "__main__":
    main()
