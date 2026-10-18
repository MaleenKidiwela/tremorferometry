#!/usr/bin/env python
"""E6 step 4(iv) -- SECONDARY-FRONT statistics that replace the COUNT (Merlin E5r-audit F6: the count is not density-monotonic).
Same window machinery as e5_thomas_harness.detect_secondary (3-h windows, 30-min step, >= 8 events from >= 4 families spanning >= 5 km,
between-family Theil-Sen speed in km/h, Spearman rho, bootstrap CI), but the products are:
 (a) in Thomas's NAMED windows -- 14 September (all 3-h windows of that UTC day) and 12-16 September within +-15 km along strike of
     family 141 -- the speed, rho and bootstrap CI of the window with the largest |rho| plus the number of windows passing all gates;
 (b) a DENSITY-ROBUST statistic: the fraction of qualifying 3-h windows (2012-09-06..09-24) whose between-family Theil-Sen speed is
     1-10 km/h with a bootstrap CI excluding 0; for Bostock's catalog also COMPLETENESS-MATCHED: in every window his events are
     subsampled to our count in that window (20 seeds); the thinned reference catalogs (step 4iii) give the same fraction per seed.
Catalogs: Bostock full, ours at the operating point, the thinned references.  Along-strike frame = Bostock's (e5_thomas_families_bostock2012.csv,
merged by lat/lon/dep -- our families inherit his hypocentres).  Writes results/e6_secondary_stats_<tag>_T<T>.csv."""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e4_thomas_harness as H  # noqa: E402
from e5_thomas_harness import ts_between, spearman  # noqa: E402
RES = H.RES
WIN_H, STEP_H, MIN_EV, MIN_FAM, SPAN_KM, RHO, SPEED = 3.0, 0.5, 8, 4, 5.0, 0.6, (1.0, 10.0)


def load(path):
    ev = H.load_catalog(path, "2012-08-15", "2012-10-16", "")
    fam = pd.read_csv(f"{RES}/e5_thomas_families_bostock2012.csv")
    key = fam.set_index([fam.lat.round(3), fam.lon.round(3), fam.dep.round(2)])
    k = list(zip(ev.lat.round(3), ev.lon.round(3), ev.dep.round(2)))
    ev["x_strike"] = pd.Series(k).map(key.x_strike.to_dict()).to_numpy()
    ev["fid"] = pd.Series(k).map(key.family_id.to_dict()).to_numpy()          # Bostock's family numbering (bkey) for every catalog
    ev = ev.dropna(subset=["x_strike"])
    return ev.sort_values("t").reset_index(drop=True), fam


def window_stats(tw, xw, fw, rng, boot=100, maxpair=5000):
    th = (tw - tw.min()) / 3600.0
    rho = spearman(xw, tw); sp = ts_between(th, xw, fw, rng, maxpair=maxpair)
    lo = hi = np.nan
    if np.isfinite(sp):
        bs = []
        for _ in range(boot):
            ii = rng.integers(0, len(tw), len(tw)); bs.append(ts_between(th[ii], xw[ii], fw[ii], rng, maxpair=maxpair))
        bs = np.array([v for v in bs if np.isfinite(v)])
        if len(bs) >= 0.8 * boot:
            lo, hi = np.percentile(bs, [2.5, 97.5])
    return rho, sp, lo, hi


def scan(ev, t_a, t_b, rng, x_restrict=None, counts=None, seed_counts=None):
    """3-h windows with start in [t_a, t_b); counts: optional dict w0 -> n to subsample to (completeness matching).
    Returns per-window rows (all windows meeting the count/family/span gates)."""
    e = ev if x_restrict is None else ev[np.abs(ev.x_strike - x_restrict[0]) <= x_restrict[1]]
    ts = e.t.to_numpy(); xs = e.x_strike.to_numpy(); fs = e.fid.to_numpy()
    rows = []
    for w0 in np.arange(t_a, t_b, STEP_H * 3600):
        i0, i1 = np.searchsorted(ts, w0), np.searchsorted(ts, w0 + WIN_H * 3600)
        tw, xw, fw = ts[i0:i1], xs[i0:i1], fs[i0:i1]
        n_raw = len(tw)
        if counts is not None:
            n_t = int(counts.get(float(w0), 0))
            if n_t < len(tw):
                ii = np.sort(rng.choice(len(tw), n_t, replace=False)); tw, xw, fw = tw[ii], xw[ii], fw[ii]
        row = dict(w0=float(w0), n_raw=n_raw, n=len(tw), n_fam=int(len(np.unique(fw))) if len(fw) else 0, span_km=float(np.ptp(xw)) if len(xw) else 0.0, qualifies=False)
        if len(tw) >= MIN_EV and row["n_fam"] >= MIN_FAM and row["span_km"] >= SPAN_KM:
            rho, sp, lo, hi = window_stats(tw, xw, fw, rng)
            row.update(qualifies=True, rho=rho, speed_kmh=sp, ci_lo=lo, ci_hi=hi,
                       passes=bool(np.isfinite(rho) and abs(rho) >= RHO and np.isfinite(sp) and SPEED[0] <= abs(sp) <= SPEED[1] and np.isfinite(lo) and not (lo <= 0 <= hi)),
                       speed_ok_ci=bool(np.isfinite(sp) and SPEED[0] <= abs(sp) <= SPEED[1] and np.isfinite(lo) and not (lo <= 0 <= hi)))
        rows.append(row)
    return pd.DataFrame(rows)


def summarize(W, label):
    q = W[W.qualifies]
    out = dict(catalog=label, n_windows=int(len(W)), n_qualifying=int(len(q)), frac_speed_1to10_ci=float(q.speed_ok_ci.mean()) if len(q) else np.nan,
               frac_pass_all=float(q.passes.mean()) if len(q) else np.nan, n_pass_all=int(q.passes.sum()) if len(q) else 0)
    return out


def named(W, label, window):
    q = W[W.qualifies & W.rho.notna()]
    if not len(q):
        return dict(catalog=label, window=window, n_qualifying=0)
    b = q.iloc[int(np.argmax(np.abs(q.rho.to_numpy())))]
    return dict(catalog=label, window=window, n_qualifying=int(len(q)), n_pass_all=int(q.passes.sum()), best_w0=str(pd.Timestamp(b.w0, unit="s").round("min")),
                best_n=int(b.n), best_rho=float(b.rho), best_speed_kmh=float(b.speed_kmh), best_ci_lo=float(b.ci_lo), best_ci_hi=float(b.ci_hi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True); ap.add_argument("--T", type=float, required=True)
    ap.add_argument("--ours-catalog", required=True); ap.add_argument("--seeds", type=int, default=20)
    a = ap.parse_args()
    rng = np.random.default_rng(21)
    his, fam = load(f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv")
    ours, _ = load(a.ours_catalog)
    x141 = float(fam.loc[fam.bostock_id == 141, "x_strike"].iloc[0])
    T = lambda s: pd.Timestamp(s).value / 1e9
    named_rows, dens_rows = [], []
    refs = sorted(glob.glob(f"{RES}/e6_thinned/e6_ref_{a.tag}_T{a.T:g}_floor_s*.csv"))
    cats = [("bostock_full", his), ("ours", ours)] + [(f"ref_s{i}", load(p)[0]) for i, p in enumerate(refs)]
    # (a) named windows
    for lab, ev in cats:
        W14 = scan(ev, T("2012-09-14"), T("2012-09-15"), rng); named_rows.append(named(W14, lab, "14 Sept (all 3-h windows)"))
        W141 = scan(ev, T("2012-09-12"), T("2012-09-17"), rng, x_restrict=(x141, 15.0)); named_rows.append(named(W141, lab, "12-16 Sept, +-15 km of family 141"))
        if lab == "ours":
            W14.to_csv(f"{RES}/e6_secondary_windows_14sep_ours_{a.tag}_T{a.T:g}.csv", index=False)
    # (b) density-robust fraction over 09-06..09-24
    ta, tb = T("2012-09-06"), T("2012-09-24")
    Wo = scan(ours, ta, tb, rng); dens_rows.append(summarize(Wo, "ours"))
    Wh = scan(his, ta, tb, rng); dens_rows.append(summarize(Wh, "bostock_full"))
    counts = dict(zip(Wo.w0, Wo.n))
    fr = []
    for s in range(a.seeds):
        Wm = scan(his, ta, tb, np.random.default_rng(500 + s), counts=counts); fr.append(summarize(Wm, f"bostock_matched_s{s}"))
    fm = pd.DataFrame(fr)
    dens_rows.append(dict(catalog="bostock_matched_to_ours (20 seeds)", n_windows=int(fm.n_windows.iloc[0]), n_qualifying=float(fm.n_qualifying.mean()),
                          frac_speed_1to10_ci=float(fm.frac_speed_1to10_ci.mean()), frac_sd=float(fm.frac_speed_1to10_ci.std()),
                          frac_pass_all=float(fm.frac_pass_all.mean()), n_pass_all=float(fm.n_pass_all.mean())))
    rr = [summarize(scan(ev, ta, tb, np.random.default_rng(700 + i)), lab) for i, (lab, ev) in enumerate(cats[2:])]
    if rr:
        fr_ = pd.DataFrame(rr)
        dens_rows.append(dict(catalog=f"thinned_reference ({len(fr_)} seeds)", n_windows=int(fr_.n_windows.iloc[0]), n_qualifying=float(fr_.n_qualifying.mean()),
                              frac_speed_1to10_ci=float(fr_.frac_speed_1to10_ci.mean()), frac_sd=float(fr_.frac_speed_1to10_ci.std()),
                              frac_pass_all=float(fr_.frac_pass_all.mean()), n_pass_all=float(fr_.n_pass_all.mean())))
    N = pd.DataFrame(named_rows); D = pd.DataFrame(dens_rows)
    N.to_csv(f"{RES}/e6_secondary_named_{a.tag}_T{a.T:g}.csv", index=False); D.to_csv(f"{RES}/e6_secondary_stats_{a.tag}_T{a.T:g}.csv", index=False)
    pd.set_option("display.width", 300)
    print("\n[named windows]"); print(N[N.catalog.isin(["bostock_full", "ours"]) | N.catalog.str.startswith("ref_s0")].round(3).to_string(index=False))
    refN = N[N.catalog.str.startswith("ref_")]
    for w, g in refN.groupby("window"):
        print(f"  thinned refs, {w}: best |rho| mean {g.best_rho.abs().mean():.2f} +- {g.best_rho.abs().std():.2f}, best speed {g.best_speed_kmh.abs().mean():.1f} +- {g.best_speed_kmh.abs().std():.1f} km/h, "
              f"windows passing all gates {g.n_pass_all.mean():.1f} +- {g.n_pass_all.std():.1f}")
    print("\n[density-robust fraction, 3-h windows 09-06..09-24]"); print(D.round(3).to_string(index=False))
    print(f"[write] results/e6_secondary_stats_{a.tag}_T{a.T:g}.csv, e6_secondary_named_{a.tag}_T{a.T:g}.csv")


if __name__ == "__main__":
    main()
