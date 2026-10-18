#!/usr/bin/env python
"""Track D §3b-2 — a PURITY-MATCHED main-front estimator (post-registration addition, coordinator
approval 2026-09-17).

WHY.  The E5 harness's main front is a Theil-Sen fit of each family's FIRST BURST against along-strike
position, over a propagation window taken from the 10th-percentile leading-edge track of all events.
At Bostock's purity that works.  At ours (~50 % precision, §3a) the false-alarm background is spread
along the whole strike, so the 10th percentile is already at the southern end before the front arrives:
the window comes out degenerate and the estimator returns nothing in 28 of 31 SSE windows (§3b).

THE ESTIMATOR (v2).  Two changes, both aimed at that failure and nothing else.

  (1) PURITY GATE.  A family's activation is not simply its first burst; the burst must be on a day
      whose event count for that family is inconsistent with that family's OWN measured false-alarm
      floor lambda(f, tier, epoch, season) (§1c).  A family-day is ACTIVATED when
          n_f(d) >= N_MIN           and      P(N >= n_f(d) | Poisson(lambda_f)) < alpha / n_days
      (Bonferroni over the days of the window).  The activation time is the start of the first burst on
      the first activated day, or that day's first event if no burst qualifies.
  (2) WINDOW FROM THE ACTIVATIONS.  The propagation window is the 5th-95th percentile of the activation
      times themselves, not the leading-edge track, so background families cannot set it.

Everything else is the harness's: burst = gap < 600 s with >= 4 events; Theil-Sen; >= 5 families;
300-draw family bootstrap; the circular and day-shuffle nulls.

NULLS.  The harness's two nulls are run unchanged, and a third is added because one of them turns out
to be degenerate for THIS statistic:
  circ    per-family circular time shift -- destroys the relative timing between families while keeping
          each family's own rate and burstiness.  This is the informative null here.
  dayshuf the harness's within-UTC-day permutation of family labels.  It permutes labels among the same
          day's events, so the MULTISET of labels on each day -- and therefore every family's DAILY
          COUNT -- is exactly preserved.  A day-resolution activation statistic is then almost invariant
          under it, and indeed the null mean reproduces the observed speed to within a few per cent in
          nearly every window.  It is reported for continuity and must NOT be read as a significance
          test of this estimator; the evidence for that is in the `null_dayshuf_mean` column itself.
  xshuf   POSITION SHUFFLE (new): keep every activation time, permute the along-strike coordinate across
          families.  This is the proper analogue -- it asks directly whether activation ORDER is related
          to POSITION.

FREE PARAMETERS AND HOW THEY WERE FIXED (no tuning on the extension years).
  N_MIN = 4           inherited from the harness's BURST_MIN — not fitted.
  window 5-95 %       fixed a priori as a symmetric robust range — not fitted.
  MIN_FAM = 5         inherited from the harness — not fitted.
  alpha               the ONLY fitted parameter: chosen from the pre-declared grid {0.05, 0.01, 0.001}
                      as the value that puts ALL FOUR 2012 referee catalogs within 20 % of Bostock's
                      2012 first-burst speed (-4.695 km/day).  The four referees are Bostock's own
                      catalog, Bostock thinned to our per-family completeness plus our measured floor
                      (20 seeds), and our own 2012 cat_core and cat_ext at the frozen operating point.
                      The extension years are not looked at until alpha is frozen.

Also fixed a priori: ONE along-strike frame for every window, computed once from all 130 Bostock family
hypocentres, so speeds are comparable between SSEs (the harness recomputes the frame per catalog).

Outputs: results/trackD/thomas_front_v2_referee.csv (the design), thomas_front_v2.csv (all windows),
thomas_front_v2_activations.csv, figures/d3b2_front_v2_summary.png.
"""
from __future__ import annotations
import argparse, glob, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import poisson
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_thomas_base as H  # noqa: E402  (bursts, theil_sen, strike_frame)
import d_tiers  # noqa: E402

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
OUT = f"{EXT}/results"
TD = f"{OUT}/trackD"
TH = f"{TD}/thomas"
THIN = f"{TD}/thinned"
BOSTOCK = f"{IN}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
HF_CUT = 3.6384
BOSTOCK_2012_SPEED = -4.694908          # e5_thomas_summary_bostock2012.csv, mc = none
CRIT = dict(BURST_GAP_S=600.0, N_MIN=4, MIN_FAM=5, WIN_LO=5.0, WIN_HI=95.0, BOOT=300,
            NULL_DRAWS=20, SPEED_ABS_MAX=20.0)
ALPHA_GRID = [0.05, 0.01, 0.001]


# ------------------------------------------------------------------ geometry (one frame, always)
def strike_map():
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv")
    lat0, lon0 = float(F.lat.mean()), float(F.lon.mean())
    s, d, how = H.strike_frame(lat0, lon0)
    cs = np.cos(np.radians(lat0))
    X = (F.lon - lon0).to_numpy() * H.R_KM_DEG * cs
    Y = (F.lat - lat0).to_numpy() * H.R_KM_DEG
    return dict(zip(F.bkey.astype(int), X * s[0] + Y * s[1])), how


# ------------------------------------------------------------------ the estimator
def activations(ev, lam, crit, alpha):
    """-> DataFrame(family_id, t_act, n_day, lam, p_day) for every family that activates."""
    out = []
    day = np.floor(ev.t.to_numpy() / 86400.0).astype(int)
    n_days = max(int(day.max() - day.min() + 1), 1)
    ev = ev.assign(_day=day)
    for f_, g in ev.groupby("family_id"):
        mu = float(lam.get(int(f_), 0.0))
        mu = max(mu, 1e-4)                       # a zero floor would make every day "significant"
        cnt = g.groupby("_day").size().sort_index()
        for d_, n_ in cnt.items():
            if n_ < crit["N_MIN"]:
                continue
            p = float(poisson.sf(n_ - 1, mu))
            if p >= alpha / n_days:
                continue
            tt = np.sort(g[g._day == d_].t.to_numpy())
            bb = H.bursts(tt, crit["BURST_GAP_S"], crit["N_MIN"])
            t_act = float(bb[0, 0]) if len(bb) else float(tt[0])
            out.append(dict(family_id=int(f_), t_act=t_act, n_day=int(n_), lam=mu, p_day=p))
            break                                 # FIRST activated day only
    return pd.DataFrame(out)


def fit_front(A, fam_x, crit, rng, boot=True):
    r = dict(n_fam_activated=len(A))
    if len(A) < crit["MIN_FAM"]:
        r["note"] = "fewer than MIN_FAM activations"
        return r
    t0 = A.t_act.min()
    lo, hi = np.percentile(A.t_act, [crit["WIN_LO"], crit["WIN_HI"]])
    W = A[(A.t_act >= lo) & (A.t_act <= hi)].copy()
    r.update(win_start=str(pd.Timestamp(lo, unit="s").round("min")),
             win_end=str(pd.Timestamp(hi, unit="s").round("min")),
             win_days=float((hi - lo) / 86400.0), n_fam_in_win=len(W))
    if len(W) < crit["MIN_FAM"]:
        r["note"] = "fewer than MIN_FAM inside the activation window"
        return r
    dd = (W.t_act.to_numpy() - t0) / 86400.0
    xx = np.array([fam_x.get(int(f_), np.nan) for f_ in W.family_id])
    k = np.isfinite(xx)
    dd, xx = dd[k], xx[k]
    if len(dd) < crit["MIN_FAM"] or np.ptp(dd) <= 0:
        r["note"] = "degenerate"
        return r
    sl, ic = H.theil_sen(dd, xx, rng=rng)
    bs = np.array([])
    if boot:   # the nulls do not need a CI, only the point estimate -- skipping it is 300x cheaper
        bs = [H.theil_sen(*(lambda i: (dd[i], xx[i]))(rng.integers(0, len(dd), len(dd))), rng=rng)[0]
              for _ in range(crit["BOOT"])]
        bs = np.array([v for v in bs if np.isfinite(v)])
    r.update(speed_kmday=float(sl),
             ci_lo=float(np.percentile(bs, 2.5)) if len(bs) else np.nan,
             ci_hi=float(np.percentile(bs, 97.5)) if len(bs) else np.nan,
             rho=float(pd.Series(dd).corr(pd.Series(xx), method="spearman")),
             span_km=float(np.ptp(xx)), n_fit=int(len(dd)),
             direction=("SE (Thomas's southward)" if sl < 0 else "NW"))
    return r


def plausible(r, crit):
    """pre-declared: >= MIN_FAM families fitted, a bootstrap CI that excludes zero, and a speed inside
    a loose physical bound (4x Thomas's 7.3 km/day reference)."""
    return bool(r.get("n_fit", 0) >= crit["MIN_FAM"]
                and np.isfinite(r.get("ci_lo", np.nan)) and np.isfinite(r.get("ci_hi", np.nan))
                and not (r["ci_lo"] <= 0 <= r["ci_hi"])
                and abs(r.get("speed_kmday", np.inf)) <= crit["SPEED_ABS_MAX"])


def null_circular(ev, rng):
    t0, t1 = ev.t.min(), ev.t.max(); span = max(t1 - t0, 1.0)
    sh = {f: rng.uniform(0, span) for f in ev.family_id.unique()}
    e2 = ev.copy(); e2["t"] = t0 + ((e2.t - t0 + e2.family_id.map(sh)) % span)
    return e2.sort_values("t").reset_index(drop=True)


def null_dayshuffle(ev, rng):
    e2 = ev.copy()
    day = np.floor(e2.t / 86400.0).astype(int)
    fam = e2.family_id.to_numpy().copy()
    for d_ in np.unique(day):
        k = np.where(day == d_)[0]
        fam[k] = fam[rng.permutation(k)]
    e2["family_id"] = fam
    return e2


def run_one(ev, lam, fam_x, crit, alpha, rng, with_nulls=True):
    A = activations(ev, lam, crit, alpha)
    r = fit_front(A, fam_x, crit, rng)
    r["alpha"] = alpha
    r["plausible"] = plausible(r, crit)
    if with_nulls and r.get("n_fit", 0) >= crit["MIN_FAM"]:
        A0 = A
        for nm in ("circ", "dayshuf", "xshuf"):
            sp = []
            for _ in range(crit["NULL_DRAWS"]):
                if nm == "xshuf":
                    # POSITION SHUFFLE: keep every activation time, permute the along-strike
                    # coordinate across families.  This is the null that actually tests "is the
                    # activation ORDER related to POSITION" for a first-activation estimator.
                    ks = list(fam_x)
                    vs = rng.permutation([fam_x[k] for k in ks])
                    rr = fit_front(A0, dict(zip(ks, vs)), crit, rng, boot=False)
                else:
                    e2 = (null_circular if nm == "circ" else null_dayshuffle)(ev, rng)
                    rr = fit_front(activations(e2, lam, crit, alpha), fam_x, crit, rng, boot=False)
                if np.isfinite(rr.get("speed_kmday", np.nan)):
                    sp.append(rr["speed_kmday"])
            sp = np.array(sp)
            r[f"null_{nm}_n"] = int(len(sp))
            r[f"null_{nm}_mean"] = float(sp.mean()) if len(sp) else np.nan
            r[f"null_{nm}_sd"] = float(sp.std()) if len(sp) else np.nan
            r[f"null_{nm}_p"] = (float((np.abs(sp) >= abs(r["speed_kmday"])).sum() + 1) / (len(sp) + 1)
                                 if len(sp) else np.nan)
    return r, A


# ------------------------------------------------------------------ inputs
def lam_2012(T):
    """per-family 2012 quiet-day detection rate at threshold T with the frozen HF veto — the same floor
    `d3_thinned_reference.py` injects into the thinned catalogs."""
    det = pd.read_parquet(f"{OUT}/e5r_detections_a1_mergeB_L12.parquet")
    det = det[det.L == 12.0]
    Hf = pd.read_parquet(f"{OUT}/a1_hf_ratio_a1_mergeB_L12.parquet")
    Hf = Hf[Hf.L == 12.0][["doy", "fam", "t", "hfr_PGC"]]
    det = det.merge(Hf, on=["doy", "fam", "t"], how="left")
    v = det.hfr_PGC.to_numpy(); det = det[~(np.isfinite(v) & (v > HF_CUT))]
    bud = pd.read_csv(f"{OUT}/e5_budget_a1_mergeB_L12.csv"); bud = bud[bud.L == 12.0]
    bud = bud[bud["T"] == sorted(bud["T"].unique())[0]]
    nd = bud[bud.doy.isin(range(228, 243))].drop_duplicates(["fam", "doy"]).groupby("fam").size()
    q = det[(det.net_mad >= T) & det.doy.isin(range(228, 243))]
    return {int(f): float(v_) for f, v_ in
            (q.groupby("fam").size().reindex(nd.index).fillna(0) / nd).items()}


def ev_from_catalog(path, t0, t1, fam_x):
    c = pd.read_csv(path) if path.endswith(".csv") else pd.read_parquet(path)
    if "time" not in c.columns and {"year", "month", "day", "hour", "min", "sec"} <= set(c.columns):
        c["time"] = (pd.to_datetime(dict(year=c.year, month=c.month, day=c.day, hour=c.hour,
                                         minute=c["min"], second=c.sec.astype(int)))
                     + pd.to_timedelta(c.sec % 1.0, unit="s"))
    c["time"] = pd.to_datetime(c["time"]).astype("datetime64[ns]")
    c = c[(c.time >= pd.Timestamp(t0)) & (c.time < pd.Timestamp(t1))].copy()
    if "family_id" not in c.columns:
        F = pd.read_csv(f"{IN}/e5r_families_2012.csv")
        key = dict(zip(zip(F.lat.round(3), F.lon.round(3), F.dep.round(2)), F.bkey.astype(int)))
        dep = c.dep if "dep" in c.columns else c.depth
        c["family_id"] = [key.get((round(a_, 3), round(b_, 3), round(z_, 2)), -1)
                          for a_, b_, z_ in zip(c.lat, c.lon, dep)]
    c = c[c.family_id.astype(int) >= 0]
    c["family_id"] = c.family_id.astype(int)
    c["t"] = c.time.astype("int64") / 1e9
    c = c[c.family_id.isin(fam_x)]
    return c[["t", "time", "family_id"]].sort_values("t").reset_index(drop=True)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="both", choices=["referee", "production", "both"])
    ap.add_argument("--alpha", type=float, default=np.nan, help="skip the design stage, use this alpha")
    ap.add_argument("--null-draws", type=int, default=40)
    ap.add_argument("--thin-seeds", type=int, default=20)
    ap.add_argument("--catalog", default=f"{TD}/lfe_extension_catalog_v1.parquet")
    ap.add_argument("--out", default=f"{TD}/thomas_front_v2.csv")
    a = ap.parse_args()
    CRIT["NULL_DRAWS"] = a.null_draws
    rng = np.random.default_rng(17)
    fam_x, how = strike_map()
    print(f"[geom] one frame for every window: {how}")
    alpha = a.alpha

    # ---------------- STAGE 1: design on the 2012 referee only
    if a.stage in ("referee", "both") and not np.isfinite(alpha):
        sys.path.insert(0, f"{EXT}/scripts")
        import d3_run_sse as D3
        O12 = D3.ref2012_ours()
        O12["cat_core_hp"] = O12.cat_core & (O12.n_sites >= 5)
        lam_core, lam_ext = lam_2012(8.75), lam_2012(8.25)
        legs = []
        eb = ev_from_catalog(BOSTOCK, "2012-08-15", "2012-10-16", fam_x)
        legs.append(("bostock2012", eb, {f: 0.0 for f in fam_x}))
        for tier, col, lm in (("ours2012_core", "cat_core", lam_core),
                              ("ours2012_ext", "cat_ext", lam_ext)):
            g = O12[O12[col]][["t", "family_id"]].copy()
            g["family_id"] = g.family_id.astype(int)
            g = g[g.family_id.isin(fam_x)]
            g["time"] = pd.to_datetime(g.t, unit="s")
            legs.append((tier, g.sort_values("t").reset_index(drop=True), lm))
        thin = sorted(glob.glob(f"{THIN}/e6_ref_a1_mergeB_L12_T8.75_floor_s*.csv"))[:a.thin_seeds]
        rows = []
        for al in ALPHA_GRID:
            for nm, ev, lm in legs:
                r, _ = run_one(ev, lm, fam_x, CRIT, al, rng, with_nulls=False)
                r.update(leg=nm, n_events=len(ev))
                rows.append(r)
            sp = []
            for f in thin:
                ev = ev_from_catalog(f, "2012-08-15", "2012-10-16", fam_x)
                r, _ = run_one(ev, lam_core, fam_x, CRIT, al, rng, with_nulls=False)
                if np.isfinite(r.get("speed_kmday", np.nan)):
                    sp.append(r["speed_kmday"])
            rows.append(dict(leg="bostock_thinned_to_our_completeness", alpha=al, n_events=np.nan,
                             speed_kmday=float(np.median(sp)) if sp else np.nan,
                             ci_lo=float(np.percentile(sp, 2.5)) if sp else np.nan,
                             ci_hi=float(np.percentile(sp, 97.5)) if sp else np.nan,
                             n_fit=len(sp), n_fam_activated=np.nan, plausible=bool(len(sp) >= 10)))
        R = pd.DataFrame(rows)
        R["rel_diff_vs_bostock2012"] = (R.speed_kmday - BOSTOCK_2012_SPEED).abs() / abs(BOSTOCK_2012_SPEED)
        R["within_20pct"] = R.rel_diff_vs_bostock2012 <= 0.20
        R.to_csv(f"{TD}/thomas_front_v2_referee.csv", index=False)
        pd.set_option("display.width", 220)
        print("\n[DESIGN — 2012 referee only]")
        print(R[["alpha", "leg", "n_events", "n_fam_activated", "n_fit", "win_days", "speed_kmday",
                 "ci_lo", "ci_hi", "rel_diff_vs_bostock2012", "within_20pct"]].round(3).to_string(index=False))
        ok = R.groupby("alpha").within_20pct.all()
        cand = [al for al in ALPHA_GRID if bool(ok.get(al, False))]
        if not cand:
            best = R.groupby("alpha").rel_diff_vs_bostock2012.max().idxmin()
            print(f"\n[DESIGN] NO alpha puts all four referee legs within 20 % — the estimator is NOT "
                  f"frozen. Worst-leg deviation is minimised at alpha = {best} "
                  f"({R.groupby('alpha').rel_diff_vs_bostock2012.max()[best]:.3f}).")
            json.dump(dict(frozen=False, reason="no alpha passes all four referee legs within 20 %",
                           grid=ALPHA_GRID, crit=CRIT),
                      open(f"{TD}/thomas_front_v2_design.json", "w"), indent=1, default=str)
            if a.stage == "both":
                raise SystemExit("[d3b2] estimator not frozen -> production stage skipped by protocol")
            return
        alpha = cand[0]
        print(f"\n[DESIGN] alpha FROZEN at {alpha} (first value of the pre-declared grid "
              f"{ALPHA_GRID} that puts all four referee legs within 20 % of Bostock 2012)")
        json.dump(dict(frozen=True, alpha=alpha, grid=ALPHA_GRID, crit=CRIT,
                       bostock_2012_speed_kmday=BOSTOCK_2012_SPEED, strike_frame=how,
                       rule="first alpha in the grid with all four 2012 referee legs within 20 %"),
                  open(f"{TD}/thomas_front_v2_design.json", "w"), indent=1, default=str)

    if a.stage == "referee":
        return

    # ---------------- STAGE 2: the frozen estimator on every SSE window
    E = pd.read_parquet(a.catalog)
    E["family_id"] = E.family_id.astype(int)
    W = pd.read_csv(f"{TD}/sse_windows.csv")
    L = pd.read_csv(f"{TD}/floor_family_epoch_season.csv")
    EPOCHS, TIERS, TSRC = d_tiers.load()
    rows, acts = [], []
    for r in W.itertuples():
        for tier, col in (("cat_core", "cat_core"), ("cat_ext", "cat_ext")):
            sub = E[(E.date >= r.t0) & (E.date <= r.t1) & E[col]]
            if len(sub) < 50:
                rows.append(dict(sse_id=r.sse_id, tier=tier, t0=r.t0, t1=r.t1, major=bool(r.major),
                                 n_events=len(sub), alpha=alpha, note="fewer than 50 events"))
                continue
            eid = sub.epoch_id.mode().iloc[0]; sea = sub.season.mode().iloc[0]
            tname = "T_core" if tier == "cat_core" else "T_ext"
            lm = L[(L.tier == tname) & (L.epoch_id == eid) & (L.season == sea)]
            lam = dict(zip(lm.family_id.astype(int), lm.lambda_fa.astype(float)))
            ev = sub[["t", "family_id"]].copy()
            ev["time"] = pd.to_datetime(ev.t, unit="s")
            ev = ev[ev.family_id.isin(fam_x)].sort_values("t").reset_index(drop=True)
            res, A = run_one(ev, lam, fam_x, CRIT, alpha, rng)
            res.update(sse_id=r.sse_id, tier=tier, t0=r.t0, t1=r.t1, major=bool(r.major),
                       n_events=len(ev), epoch_id=eid, season=sea,
                       lam_median=float(np.median(list(lam.values()))) if lam else np.nan)
            rows.append(res)
            if len(A):
                A = A.assign(sse_id=r.sse_id, tier=tier)
                A["x_strike"] = A.family_id.map(fam_x)
                acts.append(A)
            print(f"  {r.sse_id} {tier}: {len(ev)} ev, {res.get('n_fam_activated', 0)} activated, "
                  f"speed {res.get('speed_kmday', float('nan')):.2f} km/day, "
                  f"plausible={res.get('plausible')}", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(a.out, index=False)
    if acts:
        pd.concat(acts, ignore_index=True).to_csv(a.out.replace(".csv", "_activations.csv"),
                                                  index=False)
    pd.set_option("display.width", 250)
    for tier in ("cat_core", "cat_ext"):
        g = R[(R.tier == tier) & R.speed_kmday.notna()]
        print(f"\n[{tier}] {len(g)} of {int((R.tier == tier).sum())} windows fitted, "
              f"{int(g.plausible.sum())} plausible")
        if len(g):
            print(g[["sse_id", "n_events", "n_fam_activated", "n_fit", "win_days", "speed_kmday",
                     "ci_lo", "ci_hi", "rho", "null_circ_p", "null_xshuf_p", "null_dayshuf_mean",
                     "plausible"]].round(3).to_string(index=False))
    print(f"\n[write] {a.out} (+ _activations.csv; design in thomas_front_v2_design.json)")


if __name__ == "__main__":
    main()
