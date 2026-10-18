#!/usr/bin/env python
"""STAGE 1 -- CALIBRATE THE CONFIRMATION RULE ON THE NORTH, where truth exists.  Prefix `r_`.

The Washington confirmation rules used so far (>= 3 days in the confirmation window, quiet-day
FA <= 0.1 per template-day, and a re-located epicentre INTERIOR to the grid and within 10 km of
the discovery location) were never measured against anything.  `results/SCREEN_VALIDATION_NORTH.md`
showed the discovery screen they belong to keeps 2 of the 58 Bostock families the same recipe
rediscovered blind, and that a TRUE family's re-stacked location moves a median 35.5 km -- so the
10 km clause and the interior clause cannot both be right.

This module measures the confirmation, end to end, on the 2005 POLARIS bed:

    H1 = the first 12 ETS days (2005-08-29..2005-09-09)   -- the DISCOVERY half
    H2 = the last  12 ETS days (2005-09-10..2005-09-21)   -- the CONFIRMATION half
    Q  = the 20 literature-quiet days (2005-06-25..2005-07-14), outside both halves

  1. every grown template is RE-STACKED from its H1 detections only (k_grow.build_templates
     verbatim) and located                                        -> the DISCOVERY location;
  2. those H1-only templates are RE-FIRED, unchanged, at 8 x MAD on the 12 H2 days and on the
     20 quiet days                                    -> H2 detection days, quiet-day FA;
  3. the H1-only template is re-stacked on its H2 re-fire detections and re-located
                                                      -> the CONFIRMATION location;
  4. for reference the symmetric H2-only re-stack (from the ORIGINAL detections restricted to
     H2) is located too, so the distance can be read both ways.

The distributions of (2)-(3) over the 177 templates of the 58 rediscovered Bostock families set
the rule at the 85th-percentile-inclusive value; the same rule is then applied to the 1,201
unlocatable/unmatched templates (the noise population) to give its false-confirmation rate.

PRE-STATED (fixed before any number was looked at, and the same bar agent Q used):
  * a clause threshold is the value that keeps >= 85 % of the TRUE families;
  * a clause is KEPT only if at that threshold it also rejects >= 50 % of the noise population;
    otherwise it is reported as powerless and dropped from the rule.

  K_YEAR=2005 K_DAY0=... K_DAY1=... python scripts/r_north.py --stage tpl --half 1
                                     ...  --stage fire --target h2|quiet
                                     ...  --stage reloc
                                     ...  --stage calib
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import k_common as K                                                        # noqa: E402
import k_grow as KG                                                         # noqa: E402
import k_locate as KL                                                       # noqa: E402
import r_common as R                                                        # noqa: E402

S15 = ["GLBC", "TWKB", "MGCB", "TSJB", "TWBB", "PGC", "SILB", "SSIB", "LCBC",
       "B003", "B005", "JRBC", "LZB", "VGZ", "SNB"]
TPLF = "results/repro/grow15/tpl_it0.npz"
DETF = "results/repro/grow15/det_it0.parquet"
CLS = "results/repro/classify15.csv"
OUT = os.environ.get("R_NORTH_OUT", f"{R.OUT}/north")
N_H1 = 12                       # the split: first 12 ETS days / last 12 ETS days


def halves():
    d = K.bed_days()
    return d[:N_H1], d[N_H1:]


def grid():
    """k_locate's own 2005 grid -- the one loc15_it0.csv and q_sv_data used."""
    return (np.arange(47.9, 49.31, 0.03), np.arange(-124.9, -122.79, 0.04),
            np.arange(18.0, 56.1, 2.0))


def is_edge(la, lo, dp, lats, lons, deps, tol=1e-6):
    return bool(abs(la - lats[0]) < tol or abs(la - lats[-1]) < tol
                or abs(lo - lons[0]) < tol or abs(lo - lons[-1]) < tol
                or abs(dp - deps[0]) < tol or abs(dp - deps[-1]) < tol)


def locate(TPL, nd, ndays, nstack, tag):
    """k_locate verbatim on the north grid, with the 50-realisation shuffled-station control."""
    import p_common as P_
    sl = P_.slab2_interp()
    P = {g: {s: (TPL[g][s][1] + KL.TPL_PRE, TPL[g][s][2]) for s in TPL[g]} for g in TPL}
    tstas = sorted({s for g in P for s in P[g]})
    lats, lons, deps = grid()
    R.log(f"[{tag}] grid {lats.size}x{lons.size}x{deps.size}, {len(tstas)} stations")
    NODE, TTg = KL.tt_table(tstas, lats, lons, deps)
    si = {s: i for i, s in enumerate(tstas)}
    rng = np.random.default_rng(3)
    rows = []
    for g in sorted(P):
        ss = [s for s, (o, h) in P[g].items() if h >= KL.MIN_PEAK and s in si]
        base = dict(grp=g, n_sta=len(ss), n_det=int(nd.get(g, 0)), n_days=int(ndays.get(g, 0)),
                    n_stack=int(nstack.get(g, 0)))
        if len(ss) < KL.MIN_STA:
            rows.append({**base, "located": False})
            continue
        obs = np.array([P[g][s][0] for s in ss])
        idx = np.array([si[s] for s in ss])
        i, r, _ = KL.fit(obs, TTg, idx)
        rs = np.array([KL.fit(obs, TTg, idx[rng.permutation(len(idx))])[1]
                       for _ in range(KL.N_SHUF)])
        la, lo, dp = NODE[i]
        try:
            sz = float(np.asarray(sl(np.array([la]), np.array([lo]))).ravel()[0])
        except Exception:
            sz = np.nan
        rows.append({**base, "located": True, "lat": float(la), "lon": float(lo),
                     "dep": float(dp), "rms": r, "rms_shuf_med": float(np.median(rs)),
                     "rms_shuf_p05": float(np.quantile(rs, .05)),
                     "ctrl_pass": bool(r < np.quantile(rs, .05)), "slab2": sz,
                     "dz": (float(dp) - sz) if np.isfinite(sz) else np.nan,
                     "edge": is_edge(la, lo, dp, lats, lons, deps)})
    return pd.DataFrame(rows)


def restack_and_locate(D, tag, dates=None):
    """D = detections; restrict to `dates` if given; k_grow.build_templates verbatim; locate."""
    if dates is not None:
        D = D[D.date.isin(set(dates))]
    groups, skipped = {}, []
    nd = D.groupby("grp").size()
    ndays = D.groupby("grp").date.nunique()
    for g, gg in D.groupby("grp"):
        t = gg.nlargest(KG.MAX_STACK, "cc").t.to_numpy()
        if len(t) >= 3:
            groups[int(g)] = np.sort(t)
        else:
            skipped.append(int(g))
    R.log(f"[{tag}] {len(groups)} groups to re-stack ({len(D)} detections), "
          f"{len(skipped)} with < 3 detections")
    dc = K.DayCache(S15)
    TPL, info = KG.build_templates(dc, S15, groups)
    R.save_tpl(f"{OUT}/{tag}_tpl.npz", TPL)
    nstack = {g: len(v) for g, v in groups.items()}
    L = locate(TPL, nd, ndays, nstack, tag)
    for g in skipped:
        L = pd.concat([L, pd.DataFrame([dict(grp=g, n_sta=0, n_det=int(nd.get(g, 0)),
                                             n_days=int(ndays.get(g, 0)), n_stack=0,
                                             located=False)])], ignore_index=True)
    L = L.sort_values("grp")
    L.to_csv(f"{OUT}/{tag}.csv", index=False)
    lo = L[L.located == True]                                               # noqa: E712
    R.log(f"[{tag}] {len(L)} templates; located {len(lo)}; "
          f"ctrl_pass {int(lo.ctrl_pass.sum()) if len(lo) else 0}")
    return L


# ---------------------------------------------------------------- stages ----------------------
def _sel(a):
    if not a.grps:
        return None
    return set(int(x) for x in a.grps.split(","))


def stage_tpl(a):
    H1, H2 = halves()
    dates = H1 if a.half == 1 else H2
    D = pd.read_parquet(DETF)
    sel = _sel(a)
    if sel:
        D = D[D.grp.isin(sel)]
    R.log(f"[tpl h{a.half}] dates {dates[0]}..{dates[-1]} ({len(dates)} days)")
    restack_and_locate(D, f"h{a.half}date", dates=dates)


def stage_fire(a):
    H1, H2 = halves()
    TPL = R.load_tpl(f"{OUT}/h1date_tpl.npz", _sel(a))
    dates = H2 if a.target == "h2" else list(K.QUIET_DAYS)
    dc = K.DayCache(S15)
    use = []
    for d in dates:
        n = len(dc.day(d))
        if n >= 4:
            use.append(d)
        R.log(f"[fire {a.target}] {d} {n} stations{'' if n >= 4 else '  SKIPPED'}")
    R.log(f"[fire {a.target}] {len(TPL)} H1-only templates on {len(use)} days")
    D = KG.matched_filter(dc, S15, TPL, use)
    D.to_parquet(f"{OUT}/det_h1on_{a.target}.parquet", index=False)
    json.dump(dict(target=a.target, dates=use, n_days=len(use), n_templates=len(TPL),
                   n_det=int(len(D)),
                   per_tpl_day=float(len(D) / max(len(TPL) * len(use), 1))),
              open(f"{OUT}/det_h1on_{a.target}_summary.json", "w"), indent=1)
    R.log(f"[fire {a.target}] {len(D)} detections, "
          f"{D.grp.nunique() if len(D) else 0} templates firing, "
          f"{len(D) / max(len(TPL) * len(use), 1):.4f} per template-day")


def stage_reloc(a):
    D = pd.read_parquet(f"{OUT}/det_h1on_h2.parquet")
    restack_and_locate(D, "confirm_h2")


# ---------------------------------------------------------------- calibration -----------------
def _fam_best(T, col, how):
    """Per family, the BEST value among its templates (a family passes if any template does)."""
    g = T.groupby("fam")[col]
    return (g.max() if how == "max" else g.min())


def stage_calib(a):
    C = pd.read_csv(CLS)[["grp", "cls", "bkey"]]
    H1 = pd.read_csv(f"{OUT}/h1date.csv").set_index("grp")
    H2r = pd.read_csv(f"{OUT}/h2date.csv").set_index("grp")
    CF = pd.read_csv(f"{OUT}/confirm_h2.csv").set_index("grp")
    d2 = pd.read_parquet(f"{OUT}/det_h1on_h2.parquet")
    dq = pd.read_parquet(f"{OUT}/det_h1on_quiet.parquet")
    s2 = json.load(open(f"{OUT}/det_h1on_h2_summary.json"))
    sq = json.load(open(f"{OUT}/det_h1on_quiet_summary.json"))
    nq = max(sq["n_days"], 1)
    nt = max(s2["n_days"], 1)

    nd2 = d2.groupby("grp").size() if len(d2) else pd.Series(dtype=int)
    ndy2 = d2.groupby("grp").date.nunique() if len(d2) else pd.Series(dtype=int)
    ndq = dq.groupby("grp").size() if len(dq) else pd.Series(dtype=int)

    rows = []
    for r in C.itertuples():
        g = int(r.grp)
        h1 = H1.loc[g] if g in H1.index else None
        h2 = H2r.loc[g] if g in H2r.index else None
        cf = CF.loc[g] if g in CF.index else None

        def get(x, k):
            if x is None:
                return np.nan
            v = x.get(k, np.nan)
            try:
                return float(v)
            except Exception:
                return np.nan

        loc1 = h1 is not None and bool(h1.get("located", False))
        locc = cf is not None and bool(cf.get("located", False))
        loc2 = h2 is not None and bool(h2.get("located", False))
        d_conf = (float(R.epi_km(get(h1, "lat"), get(h1, "lon"), get(cf, "lat"), get(cf, "lon")))
                  if (loc1 and locc) else np.nan)
        d_sym = (float(R.epi_km(get(h1, "lat"), get(h1, "lon"), get(h2, "lat"), get(h2, "lon")))
                 if (loc1 and loc2) else np.nan)
        rows.append(dict(
            grp=g, cls=r.cls, fam=(int(r.bkey) if r.cls == "kept_family" else -1),
            # discovery half
            h1_located=loc1, h1_n_sta=get(h1, "n_sta"), h1_n_det=get(h1, "n_det"),
            h1_lat=get(h1, "lat"), h1_lon=get(h1, "lon"), h1_dep=get(h1, "dep"),
            h1_dz=get(h1, "dz"), h1_edge=(bool(h1.get("edge")) if loc1 else None),
            h1_ctrl=(bool(h1.get("ctrl_pass")) if loc1 else False),
            # confirmation half
            n_det_h2=int(nd2.get(g, 0)), n_days_h2=int(ndy2.get(g, 0)),
            n_det_quiet=int(ndq.get(g, 0)),
            fa_quiet_per_tpl_day=float(ndq.get(g, 0)) / nq,
            rate_h2_per_tpl_day=float(nd2.get(g, 0)) / nt,
            cf_located=locc, cf_n_sta=get(cf, "n_sta"),
            cf_lat=get(cf, "lat"), cf_lon=get(cf, "lon"), cf_dep=get(cf, "dep"),
            cf_dz=get(cf, "dz"), cf_edge=(bool(cf.get("edge")) if locc else None),
            cf_ctrl=(bool(cf.get("ctrl_pass")) if locc else False),
            move_km=d_conf,
            # symmetric reference
            h2_located=loc2, h2_lat=get(h2, "lat"), h2_lon=get(h2, "lon"),
            h2_dep=get(h2, "dep"), h2_dz=get(h2, "dz"), move_km_sym=d_sym))
    T = pd.DataFrame(rows)
    T.to_csv(f"{OUT}/calib_per_template.csv", index=False)

    TRUE = T[T.cls == "kept_family"].copy()
    NOISE = T[T.cls == "unlocatable"].copy()
    nfam = TRUE.fam.nunique()

    # ---- family-level best values (a family passes a clause if ANY template does) -------------
    fd = _fam_best(TRUE, "n_days_h2", "max")                 # more days = better
    ff = _fam_best(TRUE, "fa_quiet_per_tpl_day", "min")      # fewer false alarms = better
    fm = _fam_best(TRUE[np.isfinite(TRUE.move_km)], "move_km", "min") if \
        np.isfinite(TRUE.move_km).any() else pd.Series(dtype=float)
    fms = _fam_best(TRUE[np.isfinite(TRUE.move_km_sym)], "move_km_sym", "min") if \
        np.isfinite(TRUE.move_km_sym).any() else pd.Series(dtype=float)
    fz = _fam_best(TRUE[np.isfinite(TRUE.cf_dz)], "cf_dz".replace("cf_dz", "cf_dz"), "min") \
        if np.isfinite(TRUE.cf_dz).any() else pd.Series(dtype=float)

    # the days clause must be satisfiable by families that have ANY H2 detection at all
    N_days = int(R.pct_inclusive(fd.reindex(sorted(TRUE.fam.unique()), fill_value=0),
                                 0.85, "ge"))
    FA_cap = float(R.pct_inclusive(ff.reindex(sorted(TRUE.fam.unique()), fill_value=np.inf)
                                   .replace(np.inf, np.nan).dropna(), 0.85, "le"))
    D_km = float(R.pct_inclusive(fm, 0.85, "le")) if len(fm) else np.nan
    D_km_sym = float(R.pct_inclusive(fms, 0.85, "le")) if len(fms) else np.nan

    def clause_stats(df, name, fn):
        m = fn(df)
        return dict(clause=name, n=int(len(df)), n_pass=int(m.sum()),
                    frac=float(m.mean()) if len(df) else np.nan)

    def fam_pass(df, fn):
        """fraction of the 58 families with >= 1 passing template (families with no row = fail)"""
        m = fn(df)
        k = df[m].fam.nunique()
        return k, nfam, k / nfam

    c_days = lambda d: d.n_days_h2 >= N_days                                # noqa: E731
    c_fa = lambda d: d.fa_quiet_per_tpl_day <= FA_cap                       # noqa: E731
    c_ratio = lambda d: d.rate_h2_per_tpl_day >= 3 * d.fa_quiet_per_tpl_day  # noqa: E731
    c_dz1 = lambda d: d.h1_dz.abs() <= 19.0                                 # noqa: E731
    c_dz2 = lambda d: d.cf_dz.abs() <= 19.0                                 # noqa: E731
    c_move = lambda d: d.move_km <= D_km                                    # noqa: E731

    rep = dict(
        split=dict(h1=halves()[0], h2=halves()[1], n_quiet_days=nq, n_h2_days=nt),
        populations=dict(n_templates=int(len(T)), n_true_templates=int(len(TRUE)),
                         n_true_families=int(nfam), n_noise_templates=int(len(NOISE))),
        network=dict(h2_det_per_tpl_day=s2["per_tpl_day"], quiet_det_per_tpl_day=sq["per_tpl_day"],
                     ratio=s2["per_tpl_day"] / max(sq["per_tpl_day"], 1e-9)),
        distributions=dict(
            true_family_days_h2={str(q): float(np.quantile(fd.reindex(
                sorted(TRUE.fam.unique()), fill_value=0), q)) for q in (.05, .15, .25, .5, .75)},
            true_family_fa={str(q): float(np.nanquantile(ff, q))
                            for q in (.25, .5, .75, .85, .95)},
            true_family_move_km={str(q): float(np.nanquantile(fm, q))
                                 for q in (.25, .5, .75, .85, .95)} if len(fm) else {},
            true_family_move_km_sym={str(q): float(np.nanquantile(fms, q))
                                     for q in (.25, .5, .75, .85, .95)} if len(fms) else {},
            true_template_move_km={str(q): float(np.nanquantile(TRUE.move_km, q))
                                   for q in (.25, .5, .75, .85, .95)},
            n_true_families_with_both_locations=int(len(fm)),
            true_family_cf_dz={str(q): float(np.nanquantile(fz, q))
                               for q in (.25, .5, .75, .85, .95)} if len(fz) else {}),
        thresholds=dict(N_days=N_days, FA_cap=FA_cap, D_km=D_km, D_km_sym=D_km_sym,
                        DZ_km=19.0),
        clauses={})
    for nm, fn in (("days", c_days), ("fa", c_fa), ("fa_ratio", c_ratio),
                   ("dz_discovery", c_dz1), ("dz_confirm", c_dz2), ("move", c_move)):
        k, n, f = fam_pass(TRUE, fn)
        nz = clause_stats(NOISE, nm, fn)
        rep["clauses"][nm] = dict(true_families_pass=int(k), true_families=int(n),
                                  true_frac=float(f),
                                  true_templates_frac=float(fn(TRUE).mean()),
                                  noise_templates_pass=int(nz["n_pass"]),
                                  noise_frac=float(nz["frac"]),
                                  noise_rejected=float(1 - nz["frac"]))

    # ---- the adopted rule.  PRE-STATED: keep a clause only if it rejects >= 50 % of the noise.
    keep_move = bool(np.isfinite(D_km) and rep["clauses"]["move"]["noise_rejected"] >= 0.50)
    rule = dict(N_days=N_days, FA_cap=FA_cap, FA_ratio=3.0, DZ_km=19.0,
                D_km=(D_km if keep_move else None), use_location_clause=keep_move,
                use_interior_clause=False,
                note=("the location clause is retained only if, at the 85 %-inclusive distance, "
                      "it rejects >= 50 % of the north noise population"))

    def full(d):
        m = c_days(d) & c_fa(d) & c_ratio(d) & c_dz1(d) & c_dz2(d)
        if keep_move:
            m = m & c_move(d)
        return m.fillna(False)

    kk, nn, ff2 = fam_pass(TRUE, full)
    rep["adopted_rule"] = rule
    rep["adopted_rule_performance"] = dict(
        true_families_confirmed=int(kk), true_families=int(nn), true_frac=float(ff2),
        true_templates_confirmed=int(full(TRUE).sum()),
        true_templates=int(len(TRUE)),
        noise_templates_confirmed=int(full(NOISE).sum()),
        noise_templates=int(len(NOISE)),
        noise_confirmation_rate=float(full(NOISE).mean()))
    # the OLD rule, for contrast
    old = ((TRUE.n_days_h2 >= 3) & (TRUE.fa_quiet_per_tpl_day <= 0.1)
           & (TRUE.move_km <= 10.0) & (TRUE.cf_edge == False)).fillna(False)   # noqa: E712
    oldn = ((NOISE.n_days_h2 >= 3) & (NOISE.fa_quiet_per_tpl_day <= 0.1)
            & (NOISE.move_km <= 10.0) & (NOISE.cf_edge == False)).fillna(False)  # noqa: E712
    rep["old_rule_performance"] = dict(
        true_families_confirmed=int(TRUE[old].fam.nunique()), true_families=int(nfam),
        true_templates_confirmed=int(old.sum()),
        noise_templates_confirmed=int(oldn.sum()),
        noise_confirmation_rate=float(oldn.mean()))
    json.dump(rep, open(f"{OUT}/calibration.json", "w"), indent=1, default=float)
    R.write_rule(rule)
    print(json.dumps(rep, indent=1, default=float), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["tpl", "fire", "reloc", "calib"])
    ap.add_argument("--half", type=int, default=1, choices=[1, 2])
    ap.add_argument("--target", default="h2", choices=["h2", "quiet"])
    ap.add_argument("--grps", default="", help="comma list of grp ids (SMOKE only)")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    {"tpl": stage_tpl, "fire": stage_fire, "reloc": stage_reloc, "calib": stage_calib}[a.stage](a)
