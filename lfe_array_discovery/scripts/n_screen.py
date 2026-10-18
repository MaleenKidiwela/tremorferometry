#!/usr/bin/env python
"""STAGE 3 -- the SCREEN, applied to EVERY grown template before anything is counted.  Prefix `n_`.

These clauses are not new inventions: they are the five things the N-Olympic 2010 run and its
2011 confirmation re-fire showed are mandatory (PLAN.md sec.15 RESULT / results/
WA_2011_REFIRE_RESULTS.md).  The 2010 chain produced "50 candidates vs 17 null" and every one of
them died in 2011 because (i) 46/50 sat on a face of the search box, (ii) the templates fired
1.89 per template-day on NON-tremor days against 0.66 on tremor days -- persistent cultural
repeaters, anti-correlated with tremor -- and (iii) their locations moved 42-150 km when re-stacked
on other data.  Each clause below is exactly one of those failure modes turned into a test:

  (a) QUIET-DAY FALSE ALARMS, per template.  detections/template-day on the literature-quiet days
      <= 0.1 AND the tremor-day rate at least 3x the quiet-day rate.  The 2010 candidates fail
      this by construction; an LFE family cannot.
  (b) TREMOR CO-TIMING above the same-day base rate.  There is no 2008 tremor catalogue, so the
      reference is the AUDITED tremor-hour set of `n_hours.py`: the fraction of the template's
      detections falling inside an audited tremor hour, against the fraction obtained by redrawing
      the same number of detections uniformly on the same days (200 realisations).
  (c) LOCATION: >= 5 picked stations, the 50-realisation shuffled-station control, a location
      INTERIOR to the wide search box (46.5-48.8 N / -125.2..-121.8 / 15-60 km), and |dz| <= 8 km
      of Slab2.  The wide box is used FROM THE START so that "interior" means something.
  (d) >= 3 distinct tremor days.
  (e) STABILITY: re-stack from the SECOND HALF of the detections and re-locate; the location must
      reproduce within 10 km (3-D).  Evaluated by `n_relocate.py --half 2` and merged here; it is
      the only clause applied to the survivors of (a)-(d) rather than to all templates, because a
      re-stack of every grown template is a second full pass over the waveforms.

Every count is produced for the REAL chain and for the whole-chain SITE-SHIFT NULL, which has run
through byte-identical code.

  python scripts/n_screen.py --stage screen  --loc ... --det ... --det-quiet ... --out ...
  python scripts/n_screen.py --stage stability --screen ... --reloc2 ... [--reloc1 ...] --out ...
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
import n_common as N                                                        # noqa: E402

FA_MAX = 0.1             # (a) detections per template-day on the literature-quiet days
FA_RATIO = 3.0           # (a) ... and the tremor-day rate must exceed it by this factor
DZ_MAX = 8.0             # (c)
MIN_DAYS = 3             # (d)
STAB_KM = 10.0           # (e)
N_BASE = 200


def hyp_km(la1, lo1, dp1, la2, lo2, dp2):
    import s_common as S
    d = S.epi_km(la1, lo1, la2, lo2)
    return float(np.sqrt(d ** 2 + (dp1 - dp2) ** 2))


def stage_screen(a):
    N.activate()
    R = pd.read_csv(a.loc)
    D = pd.read_parquet(a.det)
    Q = pd.read_parquet(a.det_quiet) if a.det_quiet and os.path.exists(a.det_quiet) \
        else pd.DataFrame(columns=["grp", "t", "cc", "nmad", "date"])
    THR = N.tremor_hour_set()
    # the denominator is the number of days the matched filter was RUN over, not the number it
    # happened to fire on
    tdays = N.tremor_dates()
    qdays = N.quiet_days()
    # tremor hours as integer UTC hour-epochs, so the 200-realisation base rate is vectorised
    THR_EP = {int(pd.Timestamp(f"{d} {h:02d}:00:00").timestamp()) // 3600 for (d, h) in THR}
    n_t, n_q = max(len(tdays), 1), max(len(qdays), 1)
    rng = np.random.default_rng(17)

    nq = Q.groupby("grp").size() if len(Q) else pd.Series(dtype=int)
    BYG = {int(g): gg for g, gg in D.groupby("grp")} if len(D) else {}
    N.log(f"[screen {a.tag}] {len(R)} templates; {len(D)} tremor-day detections on {len(tdays)} "
          f"days; {len(Q)} quiet-day detections on {len(qdays)} days "
          f"(network rate {len(Q) / max(len(R) * n_q, 1):.3f} per template-day)")

    # base rate of an audited tremor hour among the hours of the run's own days
    base_hours = sum(1 for d in tdays for h in range(24) if (d, h) in THR) / (24 * n_t)

    rows = []
    for r in R.itertuples():
        d = BYG.get(int(r.grp), D.iloc[0:0])
        q = int(nq.get(r.grp, 0))
        rec = dict(grp=int(r.grp), n_sta=int(r.n_sta), n_det=len(d),
                   n_days=int(d.date.nunique()) if len(d) else 0,
                   n_det_quiet=q,
                   fa_quiet_per_tpl_day=q / n_q,
                   rate_tremor_per_tpl_day=len(d) / n_t,
                   located=bool(r.located))
        rec["c_a"] = bool(rec["fa_quiet_per_tpl_day"] <= FA_MAX
                          and rec["rate_tremor_per_tpl_day"]
                          >= FA_RATIO * rec["fa_quiet_per_tpl_day"]
                          and rec["n_det"] > 0)
        if len(d):
            tv = d.t.to_numpy(np.float64)

            def frac(ts):
                ep = np.floor_divide(ts.astype(np.int64), 3600)
                return float(np.mean([int(x) in THR_EP for x in ep]))

            obs = frac(tv)
            dsec = np.array([pd.Timestamp(x).timestamp() for x in d.date], dtype=np.float64)
            bs = np.array([frac(dsec + rng.random(len(dsec)) * 86400.0) for _ in range(N_BASE)])
            rec.update(co_timing=obs, co_base=float(bs.mean()), co_p=float((bs >= obs).mean()))
            rec["c_b"] = bool(obs > bs.mean())
        else:
            rec.update(co_timing=np.nan, co_base=np.nan, co_p=np.nan, c_b=False)
        if bool(r.located):
            ed = bool(r.edge)
            rec.update(lat=float(r.lat), lon=float(r.lon), dep=float(r.dep),
                       slab2=float(r.slab2), dz=float(r.dz), rms=float(r.rms),
                       ctrl_pass=bool(r.ctrl_pass), edge=ed)
            rec["c_c1"] = True
            rec["c_c2"] = bool(r.ctrl_pass)
            rec["c_c3"] = (not ed)
            rec["c_c4"] = bool(np.isfinite(r.dz) and abs(float(r.dz)) <= DZ_MAX)
        else:
            rec.update(lat=np.nan, lon=np.nan, dep=np.nan, slab2=np.nan, dz=np.nan,
                       rms=np.nan, ctrl_pass=False, edge=None,
                       c_c1=False, c_c2=False, c_c3=False, c_c4=False)
        rec["c_d"] = bool(rec["n_days"] >= MIN_DAYS)
        rec["pass_abcd"] = bool(rec["c_a"] and rec["c_b"] and rec["c_c1"] and rec["c_c2"]
                                and rec["c_c3"] and rec["c_c4"] and rec["c_d"])
        rows.append(rec)
    V = pd.DataFrame(rows)
    V.to_csv(a.out, index=False)

    clauses = ["c_c1", "c_c2", "c_c3", "c_c4", "c_d", "c_a", "c_b"]
    names = {"c_c1": "(c) located on >= 5 stations",
             "c_c2": "(c) shuffled-station control passes",
             "c_c3": "(c) INTERIOR to the wide grid",
             "c_c4": "(c) |dz| <= 8 km of Slab2",
             "c_d": "(d) >= 3 distinct tremor days",
             "c_a": "(a) quiet-day FA <= 0.1/tpl-day and >= 3x lower than tremor",
             "c_b": "(b) tremor co-timing > same-day base rate"}
    cum = np.ones(len(V), bool)
    tab = []
    for c in clauses:
        cum = cum & V[c].to_numpy(bool)
        tab.append(dict(clause=c, name=names[c], alone=int(V[c].sum()), cumulative=int(cum.sum())))
    rep = dict(tag=a.tag, n_templates=int(len(V)), n_tremor_days=len(tdays),
               n_quiet_days=len(qdays), n_det_tremor=int(len(D)), n_det_quiet=int(len(Q)),
               network_fa_per_tpl_day=float(len(Q) / max(len(R) * n_q, 1)),
               network_rate_tremor_per_tpl_day=float(len(D) / max(len(R) * n_t, 1)),
               tremor_hour_base_rate=float(base_hours),
               clause_table=tab, n_pass_abcd=int(V.pass_abcd.sum()))
    json.dump(rep, open(a.out.replace(".csv", "_summary.json"), "w"), indent=1, default=float)
    N.log(json.dumps(rep, indent=1, default=float))


def stage_stability(a):
    V = pd.read_csv(a.screen)
    R2 = pd.read_csv(a.reloc2) if os.path.exists(a.reloc2) else pd.DataFrame(columns=["grp"])
    R1 = (pd.read_csv(a.reloc1) if a.reloc1 and os.path.exists(a.reloc1)
          else pd.DataFrame(columns=["grp"]))
    RA = (pd.read_csv(a.reloc_all) if a.reloc_all and os.path.exists(a.reloc_all)
          else pd.DataFrame(columns=["grp"]))
    idx = V.set_index("grp")

    def dist(Rx, pref):
        out = {}
        for r in Rx.itertuples():
            g = int(r.grp)
            if g not in idx.index or not bool(getattr(r, "located", False)):
                continue
            row = idx.loc[g]
            if not np.isfinite(row.get("lat", np.nan)):
                continue
            out[g] = dict(
                **{f"{pref}_lat": float(r.lat), f"{pref}_lon": float(r.lon),
                   f"{pref}_dep": float(r.dep),
                   f"{pref}_edge": bool(r.edge) if r.edge == r.edge else None,
                   f"{pref}_km": hyp_km(float(row.lat), float(row.lon), float(row.dep),
                                        float(r.lat), float(r.lon), float(r.dep))})
        return out

    for Rx, pref in ((R2, "h2"), (R1, "h1"), (RA, "all")):
        if not len(Rx):
            continue
        dd = dist(Rx, pref)
        for k in next(iter(dd.values())).keys() if dd else []:
            V[k] = [dd.get(int(g), {}).get(k, np.nan) for g in V.grp]
    if "h2_km" not in V.columns:
        V["h2_km"] = np.nan
    V["c_e"] = [bool(np.isfinite(x) and x <= STAB_KM) for x in V.h2_km]
    if "h1_lat" in V.columns and "h2_lat" in V.columns:
        V["h1_h2_km"] = [hyp_km(r.h1_lat, r.h1_lon, r.h1_dep, r.h2_lat, r.h2_lon, r.h2_dep)
                         if np.isfinite(getattr(r, "h1_lat", np.nan))
                         and np.isfinite(getattr(r, "h2_lat", np.nan)) else np.nan
                         for r in V.itertuples()]
    V["candidate"] = V.pass_abcd & V.c_e
    V.to_csv(a.out, index=False)
    rep = dict(tag=a.tag, n_pass_abcd=int(V.pass_abcd.sum()),
               n_stability_tested=int(np.isfinite(V.h2_km).sum()),
               n_pass_e=int((V.pass_abcd & V.c_e).sum()),
               n_candidate=int(V.candidate.sum()),
               median_h2_km=float(np.nanmedian(V.h2_km[V.pass_abcd])) if V.pass_abcd.any()
               else np.nan)
    json.dump(rep, open(a.out.replace(".csv", "_summary.json"), "w"), indent=1, default=float)
    N.log(json.dumps(rep, indent=1, default=float))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["screen", "stability"], required=True)
    ap.add_argument("--loc", default="")
    ap.add_argument("--det", default="")
    ap.add_argument("--det-quiet", default="")
    ap.add_argument("--screen", default="")
    ap.add_argument("--reloc2", default="")
    ap.add_argument("--reloc1", default="")
    ap.add_argument("--reloc-all", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="real")
    a = ap.parse_args()
    if a.stage == "screen":
        stage_screen(a)
    else:
        stage_stability(a)


if __name__ == "__main__":
    main()
