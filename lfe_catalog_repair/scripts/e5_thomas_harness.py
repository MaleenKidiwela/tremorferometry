#!/usr/bin/env python
"""E5 step 4 -- the Thomas-2025 harness with the estimators FIXED (Merlin E4-audit F4/F5/F6, prescription 4).

Runs on any catalog with (time, lat, lon, dep, mag) -- validated here on Bostock's own Sept-2012 catalog
(zero download).  Reuses the geometry, tide placeholder, burst and Theil-Sen helpers of
scripts/e4_thomas_harness.py; everything below is new and e5_-prefixed (no E4 output is overwritten).

(a) MAIN FRONT   per family the FIRST BURST (>= 4 events with gaps < 10 min) after 2012-08-30; Theil-Sen of
                 along-strike position vs first-burst time over families whose first burst falls inside the
                 PROPAGATION WINDOW; the window is derived from the 6-h leading-edge track (10th percentile of
                 event position, bins with >= 10 events): from the bin where the leading edge has moved >= 5 km
                 south of its initial position to the bin where it first comes within 5 km of its southern
                 extreme.  Also the leading-edge and centroid track speeds INSIDE that window.  Reference
                 7.3 km/day (80 km / 11 d); accept 5-10.  Bootstrap CIs by resampling families / bins.
(b) SECONDARY    front position = the interpolated 6-h leading-edge track (NOT a first-activation fit);
    FRONTS       3-h sliding window, 30-min step, EVENT-level times; >= 8 events from >= 4 families spanning
                 >= 5 km; Theil-Sen along-strike speed over BETWEEN-family event pairs 1-10 km/h with a 200-draw
                 bootstrap CI excluding 0; |Spearman(position, time)| >= 0.6; overlapping windows merged (keep
                 the most events); "distance behind the front" recorded as a DIAGNOSTIC, not a gate.  Nulls:
                 (i) circular per-family time shift of the event times, (ii) within-UTC-day shuffle of family
                 labels across events.  Validation targets from the PMC text of Thomas et al. 2025: a front on
                 14 September at ~12 km/h (Fig. 2C) and six fronts hosted by family 141 between 12 and 16
                 September; the four Fig.-2A arrows carry no dates in the text.
(c) TIDES        one phase per BURST START; R = mean resultant length; significance from a TIME-SHIFT null
                 (tide shifted by a random 1-30 d, >= 200 draws); Schuster p dropped.  The tide is the E4
                 PLACEHOLDER equilibrium solid-Earth tide -- NOT interface shear stress; no claim of tidal
                 sensitivity is made without a Cascadia interface tidal-stress model (external dependency).
(d) RE-RUPTURE   re-activation = a burst that starts >= 12 h after the family's previous event; the interval
                 distribution is reported.
(e) DECIMATION   every metric re-run at Mc = none, 1.7, 1.8, 1.9, 2.0, 2.1 -> the completeness each metric NEEDS.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e4_thomas_harness as H  # noqa: E402  (geometry, tide placeholder, bursts, theil_sen, load_catalog)

RES = H.RES
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
CRIT = dict(T_START="2012-08-30", BURST_GAP_S=600.0, BURST_MIN=4, LEAD_BIN_H=6.0, LEAD_PCT=10, LEAD_MIN_EV=10,
            WIN_DEPART_KM=5.0, WIN_ARRIVE_KM=5.0,
            SF_WIN_H=3.0, SF_STEP_H=0.5, SF_MIN_EV=8, SF_MIN_FAM=4, SF_SPAN_KM=5.0, SF_RHO=0.6,
            SF_SPEED_KMH=(1.0, 10.0), SF_BOOT=200, SF_NULL_DRAWS=20,
            TIDE_SHIFT_D=(1.0, 30.0), TIDE_NULL_DRAWS=200, RERUP_QUIET_H=12.0,
            REF_FRONT_KMDAY=7.3, ACCEPT_FRONT_KMDAY=(5.0, 10.0), DIAG_DAY="2012-09-14",
            TIDE="PLACEHOLDER equilibrium solid-Earth tide (E4); not interface shear stress")


# ------------------------------------------------------------------ helpers
def ts_between(t, x, fam, rng, maxpair=60000):
    """Theil-Sen slope dx/dt over BETWEEN-family event pairs (same-family pairs have dx = 0 and would
    drag the median to zero).  t in hours, x in km -> km/h.  NaN if < 3 usable pairs."""
    n = len(t)
    if n < 3:
        return np.nan
    if n * (n - 1) // 2 <= maxpair:
        i, j = np.triu_indices(n, 1)
    else:
        i = rng.integers(0, n, maxpair); j = rng.integers(0, n, maxpair)
        k = i != j; i, j = i[k], j[k]
    k = (fam[i] != fam[j]) & (np.abs(t[j] - t[i]) > 1e-4)
    if k.sum() < 3:
        return np.nan
    return float(np.median((x[j][k] - x[i][k]) / (t[j][k] - t[i][k])))


def spearman(a, b):
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return float(pd.Series(a).corr(pd.Series(b), method="spearman"))


def lead_track(ev, crit, south_sign=-1):
    """6-h track: leading edge (10th pct of event position toward the propagation direction), centroid."""
    t0 = ev.t.min()
    b = np.floor((ev.t - t0) / (crit["LEAD_BIN_H"] * 3600)).astype(int)
    pct = crit["LEAD_PCT"] if south_sign < 0 else 100 - crit["LEAD_PCT"]
    trk = ev.groupby(b).agg(t_mid=("t", "mean"), x_mean=("x_strike", "mean"),
                            x_lead=("x_strike", lambda v: np.percentile(v, pct)),
                            n=("t", "size"), n_fam=("family_id", "nunique")).reset_index(drop=True)
    return trk[trk.n >= crit["LEAD_MIN_EV"]].reset_index(drop=True)


def propagation_window(trk, crit):
    """[t_a, t_b]: leading edge departs >= WIN_DEPART_KM south of its initial value .. first within
    WIN_ARRIVE_KM of its southern extreme."""
    if len(trk) < 4:
        return np.nan, np.nan, "track too short"
    x = trk.x_lead.to_numpy(); t = trk.t_mid.to_numpy(); n = trk.n.to_numpy()
    # ONSET (fix 2026-09-11, run 1): the pre-front days (08-30..09-05, 8-100 events/day scattered along strike)
    # put the 10th-percentile edge far south already, so "the first bins" is not the northern start.  The
    # onset is the first 6-h bin holding >= 20 % of the peak bin count (main-front arrival); the start
    # position is the leading edge over that bin and the next.
    on = int(np.where(n >= 0.2 * n.max())[0][0])
    x_init = float(np.max(x[on:on + 2]))
    dep = np.where((np.arange(len(x)) >= on) & (x <= x_init - crit["WIN_DEPART_KM"]))[0]
    if not len(dep):
        return np.nan, np.nan, "no departure"
    ia = dep[0]
    x_min = x[ia:].min()
    arr = np.where(x[ia:] <= x_min + crit["WIN_ARRIVE_KM"])[0]
    ib = ia + arr[0]
    return float(t[ia]), float(t[ib]), f"onset {pd.Timestamp(t[on], unit='s').round('h')} lead {x_init:+.1f} -> {x_min:+.1f} km"


def front_position(trk, t):
    """Interpolated leading-edge position (km) at time t (clamped outside the track)."""
    return np.interp(t, trk.t_mid.to_numpy(), trk.x_lead.to_numpy())


def detect_secondary(ev, trk, crit, rng, do_boot=True, diag=None):
    """3-h / 30-min sliding window over EVENT times; see module docstring.  Returns (detections, funnel).
    diag = (t_a, t_b, list): per-window statistics (before the gates) appended for windows in [t_a, t_b]."""
    ts = ev.t.to_numpy(); xs = ev.x_strike.to_numpy(); fs = ev.family_id.to_numpy()
    o = np.argsort(ts); ts, xs, fs = ts[o], xs[o], fs[o]
    win = crit["SF_WIN_H"] * 3600; step = crit["SF_STEP_H"] * 3600
    funnel = dict(windows=0, ge_events=0, ge_fam=0, ge_span=0, rho=0, speed=0, ci=0)
    cands = []
    for w0 in np.arange(ts.min(), ts.max() - win + step, step):
        funnel["windows"] += 1
        i0, i1 = np.searchsorted(ts, w0), np.searchsorted(ts, w0 + win)
        if i1 - i0 < crit["SF_MIN_EV"]:
            continue
        funnel["ge_events"] += 1
        tw, xw, fw = ts[i0:i1], xs[i0:i1], fs[i0:i1]
        if diag is not None and diag[0] <= w0 <= diag[1]:
            r_ = spearman(xw, tw); s_ = ts_between((tw - tw.min()) / 3600.0, xw, fw, rng)
            diag[2].append(dict(w0=str(pd.Timestamp(w0, unit="s").round("min")), n_ev=int(len(tw)),
                                n_fam=int(len(np.unique(fw))), span_km=float(np.ptp(xw)), rho=r_, speed_kmh=s_,
                                x_front_km=float(front_position(trk, w0 + win / 2))))
        if len(np.unique(fw)) < crit["SF_MIN_FAM"]:
            continue
        funnel["ge_fam"] += 1
        if np.ptp(xw) < crit["SF_SPAN_KM"]:
            continue
        funnel["ge_span"] += 1
        rho = spearman(xw, tw)
        if not np.isfinite(rho) or abs(rho) < crit["SF_RHO"]:
            continue
        funnel["rho"] += 1
        th = (tw - tw.min()) / 3600.0
        sp = ts_between(th, xw, fw, rng)
        if not (np.isfinite(sp) and crit["SF_SPEED_KMH"][0] <= abs(sp) <= crit["SF_SPEED_KMH"][1]):
            continue
        funnel["speed"] += 1
        lo = hi = np.nan
        if do_boot:
            bs = []
            for _ in range(crit["SF_BOOT"]):
                ii = rng.integers(0, len(tw), len(tw))
                bs.append(ts_between(th[ii], xw[ii], fw[ii], rng, maxpair=20000))
            bs = np.array([v for v in bs if np.isfinite(v)])
            if len(bs) < 0.8 * crit["SF_BOOT"]:
                continue
            lo, hi = np.percentile(bs, [2.5, 97.5])
            if lo <= 0 <= hi:
                continue
        funnel["ci"] += 1
        tmid = 0.5 * (tw.min() + tw.max())
        xf = float(front_position(trk, tmid))
        # DIAGNOSTICS (not gates): scatter of the events about the Theil-Sen line (a coherent front has a
        # few km; a zone-wide trend tens of km) and the share of events from the single busiest family
        x0 = float(np.median(xw - sp * th)); rr = np.abs(xw - (x0 + sp * th)); resid = float(np.median(rr))
        top = float(pd.Series(fw).value_counts().iloc[0] / len(fw))
        core = rr <= 5.0            # events within 5 km of the fitted line = the front itself; the rest = background
        cands.append(dict(t_start=float(tw.min()), t_end=float(tw.max()), t_mid=float(tmid),
                          resid_km=resid, top_fam_frac=top, core_frac=float(core.mean()),
                          core_span_km=float(np.ptp(xw[core])) if core.sum() > 1 else 0.0,
                          # CRASH GUARD 2026-09-11 (E5 step 5): the day-shuffle null on OUR 2012 catalog produced a window with NO event within
                          # 5 km of the fitted line -> min() of an empty array raised.  Diagnostics only (not gates); NaN when the core is empty.
                          core_n_fam=int(len(np.unique(fw[core]))), core_x_min=float(xw[core].min()) if core.any() else np.nan,
                          core_x_max=float(xw[core].max()) if core.any() else np.nan,
                          n_ev=int(len(tw)), n_fam=int(len(np.unique(fw))), span_km=float(np.ptp(xw)),
                          x_min=float(xw.min()), x_max=float(xw.max()), x_mean=float(xw.mean()),
                          speed_kmh=float(sp), speed_ci_lo=float(lo), speed_ci_hi=float(hi), rho=float(rho),
                          x_front_km=xf, behind_km=float(xw.mean() - xf),   # + = north of (behind) the front
                          families=";".join(str(int(v)) for v in np.unique(fw))))
    if not cands:
        return pd.DataFrame(columns=["t_start", "t_end", "t_mid", "n_ev", "n_fam", "span_km", "speed_kmh",
                                     "rho", "behind_km", "families"]), funnel
    C = pd.DataFrame(cands).sort_values(["n_ev", "rho"], ascending=False)
    keep, taken = [], []
    for r in C.itertuples(index=False):
        if any(not (r.t_end < a or r.t_start > b) for a, b in taken):
            continue
        keep.append(r._asdict()); taken.append((r.t_start, r.t_end))
    return pd.DataFrame(keep).sort_values("t_start").reset_index(drop=True), funnel


def null_circular(ev, rng):
    t0, t1 = ev.t.min(), ev.t.max(); span = t1 - t0
    sh = {f: rng.uniform(0, span) for f in ev.family_id.unique()}
    e2 = ev.copy(); e2["t"] = t0 + ((e2.t - t0 + e2.family_id.map(sh)) % span)
    return e2.sort_values("t").reset_index(drop=True)


def null_dayshuffle(ev, rng, fam_x):
    """Within each UTC day, permute the family labels among that day's events (times kept)."""
    e2 = ev.copy()
    day = np.floor(e2.t / 86400.0).astype(int)
    fam = e2.family_id.to_numpy().copy()
    for d in np.unique(day):
        k = np.where(day == d)[0]
        fam[k] = fam[rng.permutation(k)]
    e2["family_id"] = fam
    e2["x_strike"] = e2.family_id.map(fam_x)
    return e2


def resultant(ph_deg):
    a = np.radians(ph_deg)
    return float(np.hypot(np.cos(a).mean(), np.sin(a).mean())) if len(a) else np.nan


# ------------------------------------------------------------------ one catalog
def analyse(ev, fam, crit, tag, outdir, rng, fam_x, x141=np.nan):
    out = dict(tag=tag, n_events=len(ev), n_families=int(ev.family_id.nunique()))
    if len(ev) < 50:
        return out
    t_start = pd.Timestamp(crit["T_START"]).value / 1e9
    ev = ev[ev.t >= t_start].copy()
    t0, t1 = ev.t.min(), ev.t.max()

    # ---- (a) main front
    trk = lead_track(ev, crit)
    trk["day"] = (trk.t_mid - t_start) / 86400.0
    trk.to_csv(f"{outdir}/e5_thomas_leadtrack_{tag}.csv", index=False)
    ta, tb, how = propagation_window(trk, crit)
    out["win_how"] = how
    if np.isfinite(ta):
        out["win_start"] = str(pd.Timestamp(ta, unit="s").round("min"))
        out["win_end"] = str(pd.Timestamp(tb, unit="s").round("min"))
        out["win_days"] = (tb - ta) / 86400.0
    B = []
    for fid, g in ev.groupby("family_id"):
        bb = H.bursts(np.sort(g.t.to_numpy()), crit["BURST_GAP_S"], crit["BURST_MIN"])
        for r in bb:
            B.append(dict(family_id=fid, t_start=r[0], t_end=r[1], n=int(r[2])))
    B = pd.DataFrame(B)
    out["n_bursts"] = len(B)
    if len(B):
        B["x_strike"] = B.family_id.map(fam_x)
        B = B.sort_values(["family_id", "t_start"]).reset_index(drop=True)
        B.to_csv(f"{outdir}/e5_thomas_bursts_{tag}.csv", index=False)
        fb = B.groupby("family_id").t_start.min()
        FB = pd.DataFrame(dict(t_fb=fb, x=fb.index.map(fam_x)))
        FB["day"] = (FB.t_fb - t_start) / 86400.0
        FB.to_csv(f"{outdir}/e5_thomas_firstburst_{tag}.csv")
        out["n_fam_firstburst"] = len(FB)
        if np.isfinite(ta):
            W = FB[(FB.t_fb >= ta) & (FB.t_fb <= tb)]
            out["n_fam_firstburst_inwin"] = len(W)
            if len(W) >= 5:
                sl, ic = H.theil_sen(W.day.to_numpy(), W.x.to_numpy(), rng=rng)
                out["front_firstburst_kmday"] = sl
                out["front_firstburst_rho"] = spearman(W.day.to_numpy(), W.x.to_numpy())
                bs = [H.theil_sen(*(lambda ii: (W.day.to_numpy()[ii], W.x.to_numpy()[ii]))(
                    rng.integers(0, len(W), len(W))), rng=rng)[0] for _ in range(300)]
                out["front_firstburst_ci_lo"], out["front_firstburst_ci_hi"] = np.percentile(bs, [2.5, 97.5])
            k = (trk.t_mid >= ta) & (trk.t_mid <= tb)
            if k.sum() >= 3:
                d_, xl_, xc_ = trk.day[k].to_numpy(), trk.x_lead[k].to_numpy(), trk.x_mean[k].to_numpy()
                out["front_lead_kmday"] = H.theil_sen(d_, xl_, rng=rng)[0]
                out["front_centroid_kmday"] = H.theil_sen(d_, xc_, rng=rng)[0]
                bs = [H.theil_sen(*(lambda ii: (d_[ii], xl_[ii]))(rng.integers(0, len(d_), len(d_))),
                                  rng=rng)[0] for _ in range(300)]
                out["front_lead_ci_lo"], out["front_lead_ci_hi"] = np.percentile(bs, [2.5, 97.5])
                out["front_lead_travel_km"] = float(xl_[0] - xl_[-1])
                out["front_lead_travel_days"] = float(d_[-1] - d_[0])

    # ---- (b) secondary fronts
    if len(trk) >= 3:
        dg = (pd.Timestamp(crit["DIAG_DAY"]).value / 1e9, pd.Timestamp(crit["DIAG_DAY"]).value / 1e9 + 86400 - 1, [])
        sf, fun = detect_secondary(ev, trk, crit, rng, diag=dg)
        pd.DataFrame(dg[2]).to_csv(f"{outdir}/e5_thomas_windows_{crit['DIAG_DAY']}_{tag}.csv", index=False)
        sf["tag"] = tag
        sf["date"] = [str(pd.Timestamp(v, unit="s").round("min")) for v in sf.t_start] if len(sf) else []
        sf.to_csv(f"{outdir}/e5_thomas_secondary_{tag}.csv", index=False)
        out["n_secondary"] = len(sf)
        out["sf_funnel"] = json.dumps(fun)
        if len(sf):
            out["sf_speed_kmh_median"] = float(np.median(np.abs(sf.speed_kmh)))
            out["sf_span_km_median"] = float(np.median(sf.span_km)); out["sf_resid_km_median"] = float(np.median(sf.resid_km))
            out["sf_n_span_le30km"] = int((sf.span_km <= 30).sum()); out["sf_n_resid_le5km"] = int((sf.resid_km <= 5).sum())
            out["sf_behind_km_median"] = float(np.median(sf.behind_km))
            d12 = pd.Timestamp("2012-09-12").value / 1e9; d17 = pd.Timestamp("2012-09-17").value / 1e9
            d14 = pd.Timestamp("2012-09-14").value / 1e9; d15 = pd.Timestamp("2012-09-15").value / 1e9
            out["sf_14sep"] = int(((sf.t_mid >= d14) & (sf.t_mid < d15)).sum())
            out["sf_14sep_speeds_kmh"] = ";".join(f"{v:.1f}" for v in sf.speed_kmh[(sf.t_mid >= d14) & (sf.t_mid < d15)])
            if np.isfinite(x141):
                k141 = (sf.t_mid >= d12) & (sf.t_mid < d17) & (sf.x_min - 5 <= x141) & (sf.x_max + 5 >= x141)
                out["sf_fam141_12to16sep"] = int(k141.sum())
                k141c = (sf.t_mid >= d12) & (sf.t_mid < d17) & (sf.core_x_min - 5 <= x141) & (sf.core_x_max + 5 >= x141)
                out["sf_fam141_12to16sep_core"] = int(k141c.sum())
            out["sf_core_span_km_median"] = float(np.median(sf.core_span_km)); out["sf_core_frac_median"] = float(np.median(sf.core_frac))
        nulls = {"circ": [], "dayshuf": []}
        for _ in range(crit["SF_NULL_DRAWS"]):
            e1 = null_circular(ev, rng); nulls["circ"].append(len(detect_secondary(e1, lead_track(e1, crit), crit, rng)[0]))
            e2 = null_dayshuffle(ev, rng, fam_x); nulls["dayshuf"].append(len(detect_secondary(e2, lead_track(e2, crit), crit, rng)[0]))
        for k_, v in nulls.items():
            out[f"n_secondary_null_{k_}"] = float(np.mean(v)); out[f"n_secondary_null_{k_}_sd"] = float(np.std(v))
            out[f"n_secondary_null_{k_}_p95"] = float(np.percentile(v, 95))

    # ---- (c) tides (PLACEHOLDER), burst-level, time-shift null
    if len(B) >= 10:
        lat0, lon0 = float(fam.lat.mean()), float(fam.lon.mean())
        ph, tg, zb = H.tidal_phase(B.t_start.to_numpy(), lat0, lon0, t0 - 40 * 86400, t1 + 40 * 86400)
        B["tide_phase_PLACEHOLDER"] = ph
        thirds = np.floor((B.t_start - t0) / ((t1 - t0) / 3 + 1e-9)).clip(0, 2).astype(int).to_numpy()
        ph_unw = np.unwrap(np.angle(__import__("scipy.signal", fromlist=["hilbert"]).hilbert(zb)))
        def R_at(shift_s):
            p_ = np.degrees(np.interp(B.t_start.to_numpy() + shift_s, tg, ph_unw)) % 360.0
            return [resultant(p_)] + [resultant(p_[thirds == i]) for i in range(3)]
        obs = R_at(0.0)
        nul = np.array([R_at(rng.uniform(*crit["TIDE_SHIFT_D"]) * 86400 * rng.choice([-1, 1]))
                        for _ in range(crit["TIDE_NULL_DRAWS"])])
        for i, lab in enumerate(["all", "third1", "third2", "third3"]):
            out[f"tideR_burst_{lab}_PLACEHOLDER"] = obs[i]
            out[f"tideR_burst_{lab}_nullp_PLACEHOLDER"] = float(np.mean(nul[:, i] >= obs[i]))
            out[f"tideR_burst_{lab}_null95_PLACEHOLDER"] = float(np.percentile(nul[:, i], 95))
        out["n_bursts_third"] = ";".join(str(int((thirds == i).sum())) for i in range(3))

    # ---- (d) re-rupture by burst after >= 12 h quiescence
    if len(B):
        rr = []
        for fid, g in ev.groupby("family_id"):
            te = np.sort(g.t.to_numpy())
            bb = B[B.family_id == fid]
            for ts_ in bb.t_start.to_numpy():
                j = np.searchsorted(te, ts_ - 1e-6)
                if j == 0:
                    continue
                gap = (ts_ - te[j - 1]) / 3600.0
                rr.append(dict(family_id=fid, t_start=ts_, quiet_h=gap))
        RR = pd.DataFrame(rr)
        if len(RR):
            RR["reactivation"] = RR.quiet_h >= crit["RERUP_QUIET_H"]
            RR.to_csv(f"{outdir}/e5_thomas_rerupture_{tag}.csv", index=False)
            ra = RR[RR.reactivation]
            out["n_reactivations"] = int(len(ra))
            out["frac_bursts_reactivation"] = float(RR.reactivation.mean())
            out["n_fam_reactivated"] = int(ra.family_id.nunique())
            out["frac_fam_reactivated"] = float(ra.family_id.nunique() / ev.family_id.nunique())
            if len(ra):
                q = ra.quiet_h.quantile([.1, .25, .5, .75, .9]).to_numpy() / 24.0
                out["reactivation_interval_d_p10"], out["reactivation_interval_d_p25"], out["reactivation_interval_d_p50"], \
                    out["reactivation_interval_d_p75"], out["reactivation_interval_d_p90"] = q
                out["reactivation_frac_le5d"] = float((ra.quiet_h <= 120).mean())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default=f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv")
    ap.add_argument("--t0", default="2012-08-15"); ap.add_argument("--t1", default="2012-10-16")
    ap.add_argument("--cols", default=""); ap.add_argument("--tag", default="bostock2012")
    ap.add_argument("--outdir", default=RES); ap.add_argument("--mc", default="none,1.7,1.8,1.9,2.0,2.1")
    ap.add_argument("--null-draws", type=int, default=CRIT["SF_NULL_DRAWS"])
    a = ap.parse_args()
    CRIT["SF_NULL_DRAWS"] = a.null_draws
    rng = np.random.default_rng(11)
    print("[criteria] " + json.dumps(CRIT), flush=True)
    ev = H.load_catalog(a.catalog, a.t0, a.t1, a.cols)
    fam = ev.groupby("family_id").agg(lat=("lat", "first"), lon=("lon", "first"), dep=("dep", "first"),
                                      n=("t", "size")).reset_index()
    s, d, how = H.strike_frame(float(fam.lat.mean()), float(fam.lon.mean()))
    cs = np.cos(np.radians(float(fam.lat.mean())))
    X = (fam.lon - fam.lon.mean()).to_numpy() * H.R_KM_DEG * cs
    Y = (fam.lat - fam.lat.mean()).to_numpy() * H.R_KM_DEG
    fam["x_strike"] = X * s[0] + Y * s[1]; fam["x_dip"] = X * d[0] + Y * d[1]
    # SeisBench / Bostock family numbers by (lat, lon, dep) rounding -> family 141 position
    m = pd.read_csv(f"{SD}/metadata.csv")[["source_family_id", "source_latitude_deg", "source_longitude_deg",
                                            "source_depth_km"]].drop_duplicates()
    m["k"] = list(zip(m.source_latitude_deg.round(2), m.source_longitude_deg.round(2), m.source_depth_km.round(1)))
    fam["k"] = list(zip(fam.lat.round(2), fam.lon.round(2), fam.dep.round(1)))
    fam = fam.merge(m[["k", "source_family_id"]].rename(columns={"source_family_id": "bostock_id"}), on="k", how="left").drop(columns="k")
    x141 = float(fam.loc[fam.bostock_id == 141, "x_strike"].iloc[0]) if (fam.bostock_id == 141).any() else np.nan
    print(f"[geom] {how}; span {fam.x_strike.max()-fam.x_strike.min():.1f} km; family 141 at x = {x141:+.1f} km "
          f"({int(fam.bostock_id.notna().sum())}/{len(fam)} families matched to SeisBench ids)", flush=True)
    fam.to_csv(f"{a.outdir}/e5_thomas_families_{a.tag}.csv", index=False)
    fam_x = fam.set_index("family_id").x_strike
    ev = ev.assign(x_strike=ev.family_id.map(fam_x))
    rows = []
    for mcs in a.mc.split(","):
        mc = None if mcs == "none" else float(mcs)
        sub = ev if mc is None else ev[ev.mag >= mc]
        tag = f"{a.tag}_mcall" if mc is None else f"{a.tag}_mc{mcs}"
        print(f"\n=== {tag}: {len(sub)} events ({len(sub)/len(ev)*100:.1f} %), {sub.family_id.nunique()} families", flush=True)
        r = analyse(sub.copy(), fam, CRIT, tag, a.outdir, rng, fam_x, x141)
        r["mc"] = mcs; r["kept_frac"] = len(sub) / len(ev)
        rows.append(r)
        print("   " + ", ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items() if k != "tag"), flush=True)
    pd.DataFrame(rows).to_csv(f"{a.outdir}/e5_thomas_summary_{a.tag}.csv", index=False)
    json.dump(CRIT, open(f"{a.outdir}/e5_thomas_criteria_{a.tag}.json", "w"), indent=2)
    print(f"\n[write] {a.outdir}/e5_thomas_summary_{a.tag}.csv")


if __name__ == "__main__":
    main()
