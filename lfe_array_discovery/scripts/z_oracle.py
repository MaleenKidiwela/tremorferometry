#!/usr/bin/env python
"""STEP 1c' STAGE 6 -- the ORACLE CEILING and the 2005 SCORING REFERENCE (Merlin, prescription 6).

Step 1c scored 2005 against the FROZEN TEMPLATE detector's 2012 per-family recall/purity, simply
transferred: no ceiling existed for this bed, so neither 0/10 nor a future 3/10 was interpretable.
This stage builds the ceiling BEFORE any 1c' discovery result is read.

For each of the 30 HIDDEN families (the only place in the whole study where a hidden family is
used outside the scorer): seed the FROZEN growth chain with the family's OWN catalogued events from
its busiest hour (<= MAX_SEED_EVENTS of them, on the same absolute-time axis the pair seeds live on
-- t_origin + the network-median predicted S delay -- because the 2005 catalog times are quantised
to whole seconds and `align_stack` only searches +- 1.5 s), then run handoff -> quiet-day T_seed ->
growth -> template-scramble null -> survival -> 2 iterations -> dedupe -> locate/gates exactly as
frozen, and measure what came back.

Outputs results/step1c_prime/oracle_reference_2005.csv:
    bkey, stratum, n_events_window, survived, located, dist_km, recall_ref, purity_ref,
    purity_cluster3km_ref, kept_explained, counts
and oracle_reference_2005.json with C_iso (of 10) and C_30 (of 30).

PRE-REGISTERED FALLBACK: C_iso < 4 -> the 1c' pass bar becomes >= 0.5 x C_30 (minimum 5) hidden
families of all 30, instead of >= 3 of the 10 isolated.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import z_common as Z                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

PY = sys.executable
MAX_SEED_EVENTS = 40


def sh(cmd):
    W.log("$ " + " ".join(os.path.basename(c) for c in cmd[:3]) + " ...")
    r = subprocess.run(cmd, cwd=os.path.dirname(HERE))
    if r.returncode:
        raise SystemExit(f"oracle stage failed ({r.returncode}): {' '.join(cmd)}")


def build_seeds(out_dir):
    """Hidden-family seeds on the pair-seed time axis."""
    import p_common as P
    mem, of, order, det = W.sites()
    ST = W.station_table()
    FT = Z.fam_tt_map(det, ST)
    hid, iso, _ = S.holdouts()
    days = {d.strftime("%Y-%m-%d") for d in W.bed_days()}
    B = P.bostock()
    B = B[(B.year == W.STEP1_YEAR) & (B.bkey.isin(hid))].copy()
    B["date"] = B.time.dt.strftime("%Y-%m-%d")
    B["hh"] = B.time.dt.hour.astype(int)
    B = B[B.date.isin(days)]
    rows, mems = [], []
    for bk in sorted(int(x) for x in hid):
        g = B[B.bkey == bk]
        if not len(g):
            continue
        key = g.groupby(["date", "hh"]).size().sort_values(ascending=False)
        d, h = key.index[0]
        ev = np.sort(g[(g.date == d) & (g.hh == h)].t_origin.to_numpy(float))
        if len(ev) < 3:
            continue
        if len(ev) > MAX_SEED_EVENTS:
            ev = ev[np.linspace(0, len(ev) - 1, MAX_SEED_EVENTS).round().astype(int)]
        tt = [FT[(bk, s)] for s in det if (bk, s) in FT]
        shift = float(np.median(tt)) - 1.0 if tt else 0.0
        t = ev + shift
        cid = f"H{bk:04d}"
        rows.append(dict(cid=cid, date=str(d), n_members=int(len(t)), n_in_hour=int(len(t)),
                         n_in_3h=int(len(t)), n_bins=1, t0=float(t.min()), t1=float(t.max()),
                         t_seed=float(np.median(t)), max_net_mad=np.nan, n_pairs=0,
                         bkey=int(bk), hour=int(h), n_events_hour=int(key.iloc[0]),
                         axis_shift_s=shift))
        mems.append(pd.DataFrame(dict(cid=cid, t=t)))
    Sd = pd.DataFrame(rows)
    M = pd.concat(mems, ignore_index=True) if mems else pd.DataFrame(columns=["cid", "t"])
    Sd.to_csv(f"{out_dir}/cl_seeds.csv", index=False)
    M.to_parquet(f"{out_dir}/cl_clusters.parquet", index=False)
    W.log(f"[oracle] {len(Sd)} hidden-family seeds ({int(Sd.n_members.sum()) if len(Sd) else 0} "
          f"catalogued events)")
    return Sd


def oracle_survival(d, tag, Sd):
    """DEFECT FIX (2026-09-20, documented in freeze_step1c_prime v2).

    The frozen survival rule vetoes a template whose hits are more than SURV_MAX_NEAR_SEED within
    +- NEAR_SEED_S of its seed window -- a BURST GUARD, written for a pair seed found in tremor.
    Merlin's oracle protocol seeds each family from its **busiest hour**, i.e. deliberately from
    its biggest burst, so the guard is guaranteed to fire on a correct template: in the first 1c'
    oracle run it vetoed 21 of 28 templates (frac_near_seed 0.24-0.88) even though those templates
    had 20-46 hits on 2-8 distinct days at scramble ratios up to 17.  C_30 = 0 was an artefact of
    applying a burst guard to a seed chosen as a burst, not a property of the bed.

    For the ORACLE ONLY the clause is therefore REPORTED, not applied, and the independent-hit and
    distinct-day counts are taken OUTSIDE the seed window -- strictly more conservative than the
    frozen rule, which counts them everywhere.  Nothing in the discovery path changes, and because
    the ceiling enters the bars as 0.5 x C_30 / 0.4 x C_iso, a higher ceiling can only make the
    1c' bar HARDER: this fix cannot manufacture a PASS.
    """
    sp = f"{d}/{tag}_survival.csv"
    S_ = W.read_csv_safe(sp, W.SURV_COLS)
    if not len(S_):
        return S_
    S_.to_csv(f"{d}/{tag}_survival_frozenrule.csv", index=False)
    ip = f"{d}/{tag}_indep.parquet"
    H = pd.read_parquet(ip) if os.path.exists(ip) else pd.DataFrame(columns=["tid", "t", "date"])
    seed = Sd.set_index("cid")
    rows = []
    for r in S_.itertuples():
        tid = str(r.tid)
        g = H[H.tid.astype(str) == tid] if len(H) else H
        n_out, d_out = 0, 0
        if len(g) and tid in seed.index:
            t0s, t1s = float(seed.t0[tid]), float(seed.t1[tid])
            tv = g.t.to_numpy(float)
            out = ~((tv >= t0s - W.NEAR_SEED_S) & (tv <= t1s + W.NEAR_SEED_S))
            n_out = int(out.sum())
            d_out = int(g[out].date.nunique()) if n_out else 0
        ratio = r.scramble_ratio
        okr = (not np.isfinite(ratio)) or (ratio >= W.SCRAMBLE_RATIO)
        rows.append(dict(tpl=r.tpl, tid=r.tid, n_hits=r.n_hits, n_indep=r.n_indep,
                         n_days=r.n_days, frac_near_seed=r.frac_near_seed,
                         n_scramble=r.n_scramble, scramble_ratio=ratio,
                         n_indep_outside_seed=n_out, n_days_outside_seed=d_out,
                         survives=bool(n_out >= W.SURV_MIN_INDEP and d_out >= W.SURV_MIN_DAYS
                                       and okr)))
    R = pd.DataFrame(rows)
    R.to_csv(sp, index=False)
    W.log(f"[oracle/survival] {tag}: frozen rule {int(S_.survives.astype(bool).sum())} survivors, "
          f"oracle rule (outside the seed window, burst guard reported not applied) "
          f"{int(R.survives.sum())}")
    return R


def run_chain(d, Sd=None):
    G = f"{HERE}/w_grow.py"
    sh([PY, G, "handoff", "--clusters", f"{d}/cl_clusters.parquet", "--seeds",
        f"{d}/cl_seeds.csv", "--out", f"{d}/bank0.npz"])
    sh([PY, G, "quiet", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0"])
    sh([PY, G, "grow", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
        "--thr", f"{d}/p0_thresholds.csv"])
    sh([PY, G, "scramble", "--clusters", f"{d}/cl_clusters.parquet", "--seeds",
        f"{d}/cl_seeds.csv", "--out", f"{d}/bank0scr.npz"])
    sh([PY, G, "grow", "--bank", f"{d}/bank0scr.npz", "--out", f"{d}/p0scr",
        "--thr", f"{d}/p0_thresholds.csv"])
    sh([PY, G, "survive", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
        "--seeds", f"{d}/cl_seeds.csv", "--scramble-hits", f"{d}/p0scr_hits.parquet",
        "--scramble-bank", f"{d}/bank0scr.npz"])
    if Sd is not None:
        oracle_survival(d, "p0", Sd)
    prev, tag = f"{d}/bank0.npz", "p0"
    for it in range(1, W.N_ITERATIONS + 1):
        sh([PY, f"{HERE}/w_iterate.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
            "--survivors", f"{d}/p0_survival.csv", "--out-bank", f"{d}/bank{it}.npz",
            "--out-diag", f"{d}/iter{it}_diag.csv"])
        sh([PY, G, "quiet", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}"])
        sh([PY, G, "grow", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
            "--thr", f"{d}/p{it}_thresholds.csv"])
        sh([PY, G, "survive", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
            "--seeds", f"{d}/cl_seeds.csv"])
        if Sd is not None:
            oracle_survival(d, f"p{it}", Sd)
        prev, tag = f"{d}/bank{it}.npz", f"p{it}"
    sh([PY, f"{HERE}/w_dedupe.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
        "--out", f"{d}/dedupe.csv"])
    sh([PY, f"{HERE}/w_locate_gates.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
        "--dedupe", f"{d}/dedupe.csv", "--survival", f"{d}/{tag}_survival.csv",
        "--out", f"{d}/lg"])
    return prev, tag


def score(d, Sd, out_csv):
    import p_common as P
    import g_rediscover as GR
    import w_score_step1 as SC
    F = P.families()
    fam = {int(b): (float(F.lat[b]), float(F.lon[b])) for b in F.index}
    B = P.bostock()
    B = B[B.year == W.STEP1_YEAR]
    ev = {int(b): np.sort(g.t_origin.to_numpy(float)) for b, g in B.groupby("bkey")}
    hid, iso, _ = S.holdouts()
    t_lo = pd.Timestamp(W.STEP1_DAY0).timestamp()
    t_hi = pd.Timestamp(W.STEP1_DAY1).timestamp() + 86400
    ev = {k: v[(v >= t_lo) & (v <= t_hi)] for k, v in ev.items()}
    CT = SC.cluster_truth(fam, ev, sorted(hid))
    kept = [f for f in fam if f not in hid]

    C = (pd.read_parquet(f"{d}/lg_clusters.parquet")
         if os.path.exists(f"{d}/lg_clusters.parquet")
         else pd.DataFrame(columns=["cid", "t", "lat", "lon"]))
    L = W.read_csv_safe(f"{d}/lg_locations.csv")
    SV = W.read_csv_safe(f"{d}/p{W.N_ITERATIONS}_survival.csv", W.SURV_COLS)
    surv = set(SV[SV.survives.astype(bool)].tid.astype(str)) if len(SV) else set()
    loc = {str(r.tid): r for r in L.itertuples()} if len(L) else {}
    rows = []
    for r in Sd.itertuples():
        bk, cid = int(r.bkey), str(r.cid)
        g = C[C.cid == cid] if len(C) else C
        t = np.sort(g.t.to_numpy(float)) if len(g) else np.zeros(0)
        e = ev.get(bk, np.zeros(0))
        rec = pur = pc3 = kb = np.nan
        dist = np.nan
        if len(t):
            mr, _ = GR.one2one(e, t, W.SCORE_TOL_S)
            mp, _ = GR.one2one(t, e, W.SCORE_TOL_S)
            mc, _ = GR.one2one(t, CT.get(bk, np.zeros(0)), W.SCORE_TOL_S)
            rec = float(mr.mean()) if len(e) else 0.0
            pur = float(mp.mean())
            pc3 = float(mc.mean())
            kb = max((GR.one2one(t, ev.get(k, np.zeros(0)), W.SCORE_TOL_S)[0].mean()
                      for k in kept), default=0.0)
            dist = float(S.epi_km(float(g.lat.iloc[0]), float(g.lon.iloc[0]), *fam[bk]))
        lr = loc.get(cid)
        counts = bool(len(t) and np.isfinite(dist) and dist <= W.SCORE_R_KM and rec > 0
                      and (not np.isfinite(kb) or kb <= W.SCORE_MAX_EXPLAINED))
        rows.append(dict(bkey=bk, stratum=("isolated" if bk in iso else "other"),
                         seed_date=r.date, seed_hour=int(r.hour),
                         n_seed_events=int(r.n_members), n_events_window=int(len(e)),
                         survived=bool(cid in surv), n_cluster=int(len(t)),
                         located=bool(lr is not None and bool(getattr(lr, "located", False))),
                         dist_km=dist, recall_ref=rec, purity_ref=pur,
                         purity_cluster3km_ref=pc3, kept_explained=kb,
                         frac_near_seed=float(SV[SV.tid.astype(str) == cid].frac_near_seed.iloc[0])
                         if len(SV) and (SV.tid.astype(str) == cid).any() else np.nan,
                         counts=counts))
    R = pd.DataFrame(rows)
    R.to_csv(out_csv, index=False)
    c_iso = int(R[(R.stratum == "isolated") & R.counts].shape[0]) if len(R) else 0
    c_30 = int(R.counts.sum()) if len(R) else 0
    fb = c_iso < 4
    bar = max(int(np.ceil(0.5 * c_30)), 5) if fb else W.ACCEPT_REDISCOVERED
    out = dict(
        n_seeded=int(len(R)), n_isolated_seeded=int((R.stratum == "isolated").sum())
        if len(R) else 0,
        C_iso=c_iso, C_30=c_30,
        n_survived=int(R.survived.sum()) if len(R) else 0,
        n_located=int(R.located.sum()) if len(R) else 0,
        median_recall_ref=float(np.nanmedian(R.recall_ref)) if len(R) else None,
        median_purity_ref=float(np.nanmedian(R.purity_ref)) if len(R) else None,
        median_dist_km=float(np.nanmedian(R.dist_km)) if len(R) else None,
        fallback_active=bool(fb),
        pass_bar=dict(mode=("all30_half_C30" if fb else "isolated_ge3"), value=int(bar),
                      note=("C_iso < 4: the isolated stratum is uninformative on this bed, so the "
                            "1c' bar becomes >= 0.5 x C_30 (minimum 5) of all 30 hidden families"
                            if fb else "C_iso >= 4: the isolated stratum is informative")),
        protocol=dict(seed="the family's own catalogued events from its busiest hour",
                      max_seed_events=MAX_SEED_EVENTS,
                      axis="t_origin + network-median predicted S delay - 1 s (the pair-seed axis)",
                      align="align_stack +- 1.5 s (2005 catalog times are 1-s quantised)",
                      chain="handoff -> quiet-day T_seed -> growth -> scramble -> survival -> "
                            f"{W.N_ITERATIONS} iterations -> dedupe -> locate/gates, all frozen",
                      counts_clause=f"centroid <= {W.SCORE_R_KM} km AND recall > 0 AND "
                                    f"kept-explained <= {W.SCORE_MAX_EXPLAINED}",
                      survival_defect_fix="ORACLE ONLY: the near-seed burst guard is reported, "
                                          "not applied (the oracle seed IS the family's busiest "
                                          "burst by protocol, so the guard vetoed 21/28 correct "
                                          "templates); independent hits and distinct days are "
                                          "counted OUTSIDE the seed window, which is stricter "
                                          "than the frozen rule.  The discovery path is "
                                          "unchanged, and a higher ceiling only raises the 1c' "
                                          "bar (0.5 x C_30 / 0.4 x C_iso)."),
        verdict="DONE")
    json.dump(out, open(out_csv.replace(".csv", ".json"), "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=f"{W.RES}/step1c_prime/oracle")
    ap.add_argument("--out", default=f"{W.RES}/step1c_prime/oracle_reference_2005.csv")
    ap.add_argument("--seeds-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    Sd = build_seeds(a.dir)
    if a.seeds_only or not len(Sd):
        if not len(Sd):
            raise SystemExit("no hidden family has >= 3 catalogued events in the bed window")
        return
    run_chain(a.dir, Sd)
    score(a.dir, Sd, a.out)


if __name__ == "__main__":
    main()
