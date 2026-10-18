#!/usr/bin/env python
"""Stage 9 -- the Deliverable-A CANDIDATE TABLE, the claim ladder, and the Chestler comparison.

Assembles `results/step5/CANDIDATES.csv`: one row per surviving, deduped template, carrying

  N and the sites it registers on; the location and its FOUR controls (shuffled-pick ratio,
  synthetic floor, leave-one-site-out, free-depth sanity); the Gate-0 verdict and its parts; the
  Gate-1 verdict PER WITHHELD SITE; the whole-chain site-shift null counts; and the Step-5b
  injection recovery at the matched amplitude.

CLAIM LADDER (verbatim from Merlin Q5(b), applied here and nowhere else)
  candidate            N >= 20, Gate 0 PASS (incl. the HF veto at the WA 100-Hz borehole, and
                       distinct from every already-accepted WA family by the adopted registering-
                       site rule), Gate 1 PASS at >= 1 withheld site.
  gated family         Gate 1 PASS at >= 2 withheld SITES (both nulls: 200 random-time p99 with
                       >= +3 dB excess, and 100 moveout-shuffled p95) AND "located" per Q4 AND the
                       whole-chain site-shift null condition (real >= 5 x null mean, null mean <= 1).
  confirmed recurrent  re-fires in the sealed 2011-08 window (or in >= 2 later ETS in the
                       extension) with Gate 3 PASS on the long record.  NOT EVALUABLE HERE -- the
                       2011 window is sealed and unread; this column is always "not_evaluated".

GATE-1 NON-PASS: FAIL or UNTESTABLE (Merlin F4).  Without a Washington control a non-PASS is
UNTESTABLE by construction.  Once the Step-5b injections exist they provide one: a non-PASS at a
withheld station is re-labelled FAIL only if an INJECTED repeater at matched N (nearest planted N)
and comparable distance PASSED at that same station; otherwise it stays UNTESTABLE.  Both the raw
`g_gates` verdict and the re-labelled one are written.

CHESTLER (descriptive prediction, never used to select, merge or locate): the distribution of
nearest-Chestler distances for the gated list against 1,000 draws of the same number of random
active cells.

    python scripts/x_candidates.py --root results/step5
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
import x_common as X                                                        # noqa: E402

N_RANDOM_DRAWS = 1000


def gate1_relabel(raw_verdict, sta, N, dist_km, ctl):
    """UNTESTABLE unless an injected repeater at matched N and comparable distance PASSED there."""
    if raw_verdict == "PASS":
        return "PASS", "injection control not needed"
    if ctl is None or not len(ctl):
        return "UNTESTABLE", "no WA Gate-1 control (Merlin F4): non-PASS is untestable"
    c = ctl[ctl.sta == sta]
    if not len(c):
        return "UNTESTABLE", f"no injected control at {sta}"
    c = c.assign(dn=(c.N - N).abs())
    nsel = c.N.iloc[int(np.argmin(c.dn.to_numpy()))]
    c = c[c.N == nsel]
    if np.isfinite(dist_km) and c.dist_km.notna().any():
        c = c.assign(dd=(c.dist_km - dist_km).abs()).sort_values("dd").head(3)
    if (c.verdict == "PASS").any():
        return "FAIL", (f"injected control at {sta} PASSED at matched N={int(nsel)} "
                        f"-> the test had power here")
    return "UNTESTABLE", (f"injected control at {sta} did not pass at matched N={int(nsel)} "
                          f"-> the test has no power here")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=X.OUT)
    ap.add_argument("--lg", default="")
    ap.add_argument("--inject", default="")
    a = ap.parse_args()
    import s_common as S
    W = X.activate()
    C = X.cfg()
    root = a.root
    lg = a.lg or f"{root}/lg"
    inj = a.inject or f"{root}/inject"

    L = pd.read_csv(f"{lg}_locations.csv") if os.path.exists(f"{lg}_locations.csv") \
        else pd.DataFrame()
    Gv = pd.read_csv(f"{lg}_gates.csv") if os.path.exists(f"{lg}_gates.csv") else pd.DataFrame()
    last = "2"
    if os.path.exists(f"{root}/last_iteration.txt"):
        last = open(f"{root}/last_iteration.txt").read().strip() or "2"
    svp = f"{root}/p{last}_survival.csv"
    SV = pd.read_csv(svp) if os.path.exists(svp) else pd.DataFrame()
    DD = pd.read_csv(f"{root}/dedupe.csv") if os.path.exists(f"{root}/dedupe.csv") \
        else pd.DataFrame()
    NL = pd.read_csv(f"{root}/null_summary.csv") if os.path.exists(f"{root}/null_summary.csv") \
        else pd.DataFrame()
    RC = pd.read_csv(f"{inj}/recovery.csv") if os.path.exists(f"{inj}/recovery.csv") \
        else pd.DataFrame()
    CT = pd.read_csv(f"{inj}/gate1_controls_wa.csv") \
        if os.path.exists(f"{inj}/gate1_controls_wa.csv") else None

    null_mean = float(NL.n_gate_pass.mean()) if len(NL) else np.nan
    null_counts = list(NL.n_gate_pass) if len(NL) else []
    wsites = list(X.WITHHELD_SITES)
    sta2site = {m: q for q, mm in X.WITHHELD_SITES.items() for m in mm}

    rows = []
    for t in (L.itertuples() if len(L) else []):
        tid = str(t.tid)
        g = Gv[Gv.pid.astype(str) == tid] if len(Gv) else pd.DataFrame()
        sv = SV[SV.tid.astype(str) == tid] if len(SV) else pd.DataFrame()
        full = f"{lg}_gate_full_{tid}.json"
        J = json.load(open(full)) if os.path.exists(full) else {}
        r = dict(tid=tid,
                 N=int(sv.n_indep.iloc[0]) if len(sv) else int(getattr(t, "n_members", 0)),
                 n_days=int(sv.n_days.iloc[0]) if len(sv) else np.nan,
                 scramble_ratio=float(sv.scramble_ratio.iloc[0]) if len(sv) else np.nan,
                 sites=getattr(t, "sites", ""),
                 n_picked=getattr(t, "n_picked", np.nan),
                 lat=getattr(t, "lat", np.nan), lon=getattr(t, "lon", np.nan),
                 dep=getattr(t, "dep", np.nan), rms_s=getattr(t, "rms_s", np.nan),
                 located=bool(getattr(t, "located", False)),
                 ctl_shuffled_ratio=getattr(t, "shuffled_ratio", np.nan),
                 ctl_synth_med_km=getattr(t, "synth_med_km", np.nan),
                 ctl_loso_max_km=getattr(t, "loso_max_km", np.nan),
                 ctl_depth_diff_km=getattr(t, "depth_diff_km", np.nan),
                 n_sites_60km=getattr(t, "n_sites_60km", np.nan),
                 gap_deg=getattr(t, "gap_deg", np.nan),
                 gate0=(str(g.g0.iloc[0]) if len(g) and "g0" in g else None),
                 g0_hf_median=(float(g.g0_hf_median.iloc[0])
                               if len(g) and "g0_hf_median" in g else np.nan),
                 g0_min_dist_km=(float(g.g0_min_dist_km.iloc[0])
                                 if len(g) and "g0_min_dist_km" in g else np.nan),
                 gate1_raw=(str(g.g1.iloc[0]) if len(g) and "g1" in g else None),
                 gate2=(str(g.g2.iloc[0]) if len(g) and "g2" in g else None),
                 verdict_gates=(str(g.verdict.iloc[0]) if len(g) else None),
                 null_counts=";".join(str(int(x)) for x in null_counts),
                 null_mean=null_mean)
        npass_sites = set()
        for p in J.get("gate1", {}).get("per_station", []):
            sta = p.get("sta")
            q = sta2site.get(sta, sta)
            raw = p.get("verdict")
            rel, why = gate1_relabel(raw, sta, int(p.get("n_used", 0)),
                                     float(p.get("dist_km", np.nan)), CT)
            r[f"g1_{sta}_raw"] = raw
            r[f"g1_{sta}"] = rel
            r[f"g1_{sta}_excess_db"] = p.get("excess_db", np.nan)
            r[f"g1_{sta}_n"] = p.get("n_used", 0)
            r[f"g1_{sta}_why"] = why
            if rel == "PASS":
                npass_sites.add(q)
        r["g1_sites_pass"] = len(npass_sites)
        r["g1_sites_pass_list"] = ",".join(sorted(npass_sites))
        if len(RC):
            near = RC.assign(d=lambda z: S.epi_km(z.lat, z.lon, r["lat"], r["lon"]))
            near = near.sort_values("d")
            r["inject_nearest_sid"] = str(near.sid.iloc[0]) if len(near) else None
            r["inject_recovered_frac"] = float(RC.recovered.mean()) if len(RC) else np.nan
            r["inject_median_loc_err_km"] = (float(RC[RC.recovered].loc_err_km.median())
                                             if RC.recovered.any() else np.nan)
        # ---- claim ladder
        cand = bool(r["N"] >= 20 and r["gate0"] == "PASS" and r["g1_sites_pass"] >= 1)
        gated = bool(cand and r["g1_sites_pass"] >= 2 and r["located"]
                     and np.isfinite(null_mean) and null_mean <= 1.0)
        r["claim_candidate"] = cand
        r["claim_gated_family"] = gated
        r["claim_confirmed_recurrent"] = "not_evaluated (2011 window sealed and unread)"
        rows.append(r)
    R = pd.DataFrame(rows)
    os.makedirs(root, exist_ok=True)
    R.to_csv(f"{root}/CANDIDATES.csv", index=False)

    # ---- the whole-chain null condition, applied to the run as a whole
    n_gated = int(R.claim_gated_family.sum()) if len(R) else 0
    null_ok = bool(np.isfinite(null_mean) and null_mean <= 1.0
                   and (null_mean == 0 or n_gated >= 5 * null_mean))

    # ---- Chestler descriptive prediction
    ch = dict(available=False)
    if os.path.exists(X.CHESTLER) and n_gated:
        Ch = pd.read_csv(X.CHESTLER)
        cells = pd.read_csv(f"{X.OUT}/active_cells.csv")
        gl = R[R.claim_gated_family]
        d_obs = np.array([np.min(S.epi_km(Ch.latitude.to_numpy(), Ch.longitude.to_numpy(),
                                          r.lat, r.lon)) for r in gl.itertuples()])
        rng = np.random.default_rng(20260920)
        med_null = []
        for _ in range(N_RANDOM_DRAWS):
            k = rng.choice(len(cells), size=len(gl), replace=False)
            dn = np.array([np.min(S.epi_km(Ch.latitude.to_numpy(), Ch.longitude.to_numpy(),
                                           cells.lat.iloc[i], cells.lon.iloc[i])) for i in k])
            med_null.append(float(np.median(dn)))
        med_null = np.asarray(med_null)
        ch = dict(available=True, n_gated=int(len(gl)), n_chestler=int(len(Ch)),
                  median_nearest_km=float(np.median(d_obs)),
                  per_family_km=[float(x) for x in d_obs],
                  null_median_km=float(np.median(med_null)),
                  null_p05_km=float(np.quantile(med_null, .05)),
                  null_p95_km=float(np.quantile(med_null, .95)),
                  p_closer=float(((med_null <= np.median(d_obs)).sum() + 1)
                                 / (len(med_null) + 1)),
                  note=("DESCRIPTIVE PREDICTION ONLY -- the Chestler coordinates were never used "
                        "to select, merge or locate anything, and this comparison is made only "
                        "after the gated list is frozen"))

    summary = dict(
        n_templates=int(len(R)), n_located=int(R.located.sum()) if len(R) else 0,
        n_candidates=int(R.claim_candidate.sum()) if len(R) else 0,
        n_gated_families=n_gated,
        null_counts=null_counts, null_mean=null_mean, null_condition_met=null_ok,
        withheld_sites=wsites,
        n_injected_recovered=int(RC.recovered.sum()) if len(RC) else 0,
        n_injected=int(len(RC)),
        injection_median_loc_err_km=(float(RC[RC.recovered].loc_err_km.median())
                                     if len(RC) and RC.recovered.any() else None),
        chestler=ch,
        accept_bar="Merlin Step 5: >= 5 gated families, null mean <= 1, real >= 5 x null mean",
        sealed_2011="untouched",
    )
    json.dump(summary, open(f"{root}/step5_summary.json", "w"), indent=1, default=str)
    X.log(f"[candidates] {len(R)} templates -> {summary['n_candidates']} candidates, "
          f"{n_gated} gated families; null {null_counts} (mean {null_mean}); "
          f"injection {summary['n_injected_recovered']}/{summary['n_injected']} recovered")
    if ch.get("available"):
        X.log(f"[chestler] median nearest {ch['median_nearest_km']:.2f} km vs random-cell null "
              f"{ch['null_median_km']:.2f} km [{ch['null_p05_km']:.2f}, {ch['null_p95_km']:.2f}], "
              f"p={ch['p_closer']:.3f}")


if __name__ == "__main__":
    main()
