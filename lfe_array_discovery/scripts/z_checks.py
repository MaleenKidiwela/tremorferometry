#!/usr/bin/env python
"""STEP 1c' -- the ACCEPT / ABANDON gates of the queue, one per stage.

Each subcommand reads the stage's own output and exits

    0   accept (the queue continues)
    20  ABANDON -- Merlin's pre-registered abandon condition fired; the reason is printed and
        written to results/step1c_prime/ABANDON_REASON.txt so the FAILED marker carries it.

No stage result is ever re-interpreted here: the verdicts are the ones the stage scripts computed
against the bars frozen in results/freeze_step1c_prime.json.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402

RUN = f"{W.RES}/step1c_prime"
MIN_SEEDS_HANDED_OFF = 500


def abandon(stage, reason):
    os.makedirs(RUN, exist_ok=True)
    txt = f"ABANDON at {stage}: {reason}"
    open(f"{RUN}/ABANDON_REASON.txt", "w").write(txt + "\n")
    print(txt)
    sys.exit(20)


def ok(stage, msg):
    print(f"[check/{stage}] ACCEPT -- {msg}")
    sys.exit(0)


def c_d8():
    j = json.load(open(f"{RUN}/d8.json"))
    v = j["verdict"]
    if v.startswith("ABANDON"):
        abandon("D8", f"{v}: {j['reason']}")
    if v == "BELOW_ACCEPT":
        # below the accept bar but not an abandon condition: recorded, the run continues and the
        # write-up must carry it (Merlin's accept clause is diagnostic, his abandon clause is not)
        print(f"[check/D8] BELOW_ACCEPT (recorded, not an abandon condition): {j['reason']}")
    ok("D8", f"{v}; peaks within 3 s at {j['n_sites_within_3s']} sites; control "
             f"{j['control_excess_median']}; chain ratio {j['chain_ratio_median']}")


def c_health():
    j = json.load(open(f"{RUN}/station_hour_veto.json"))
    if j["verdict"] == "ABANDON_BED_SHRINKS":
        abandon("station-hour health veto", j["reason"])
    ok("health", f"{j['n_vetoed']}/{j['n_station_hours']} station-hours vetoed; SHVB flagged "
                 f"{j['shvb_veto_frac_on_bad_doys']} of its doy 241-251/255-261 hours")


def c_quiet():
    j = json.load(open(f"{RUN}/quiet_hours.json"))
    if j["verdict"] != "ACCEPT":
        abandon("quiet set", j["reason"])
    ok("quiet", f"{j['n_quiet_hours']} hours, {j['clock_hours_covered']}/24 clock hours, "
                f"{j['n_hours_eq_vetoed']} hours earthquake-vetoed")


def c_ac():
    """Stage-4 accept: no hour with zero pairs for a NON-veto reason; the 2005-07-09 23h burst."""
    I = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f"{RUN}/ac.parts/*_info.csv"))],
                  ignore_index=True) if glob.glob(f"{RUN}/ac.parts/*_info.csv") else pd.DataFrame()
    P = pd.read_parquet(f"{RUN}/ac.parquet")
    n_pairs = P.groupby(["date", "hour"]).size() if len(P) else pd.Series(dtype=int)
    ran = {(os.path.basename(p)[:10], int(os.path.basename(p)[11:13]))
           for p in glob.glob(f"{RUN}/ac.parts/*.parquet")}
    zero = [k for k in sorted(ran) if int(n_pairs.get(k, 0)) == 0]
    bad = []
    for (d, h) in zero:
        r = I[(I.date == d) & (I.hour == h)] if len(I) else I
        if not len(r):
            bad.append((d, h, "no census row"))
        elif int(r.n_hour_sites.iloc[0]) >= W.MIN_SITES and int(r.n_veto_health.iloc[0]) == 0 \
                and int(r.n_veto_flat.iloc[0]) == 0:
            bad.append((d, h, f"{int(r.n_hour_sites.iloc[0])} sites, no veto"))
    # the 2005-07-09 23h quiet-hour burst
    Q = pd.read_csv(f"{RUN}/quiet_hours.csv")
    burst_in = bool(((Q.date == "2005-07-09") & (Q.hour == 23)).any())
    burst_pairs = None
    if burst_in and os.path.exists(f"{RUN}/ac_quiet.parquet"):
        QP = pd.read_parquet(f"{RUN}/ac_quiet.parquet")
        burst_pairs = int(((QP.date == "2005-07-09") & (QP.hour == 23)
                           & (QP.net_mad >= 8.0)).sum()) if len(QP) else 0
    out = dict(n_hours_run=len(ran), n_zero_pair_hours=len(zero),
               n_zero_for_non_veto_reasons=len(bad), examples=bad[:10],
               burst_hour_in_quiet_set=burst_in, burst_pairs_ge8=burst_pairs,
               n_windows_vetoed_median=float(I.n_windows_vetoed.median()) if len(I) else None,
               n_station_hours_health_vetoed=int(I.n_veto_health.sum()) if len(I) else 0,
               n_station_hours_flat_dropped=int(I.n_veto_flat.sum()) if len(I) else 0)
    json.dump(out, open(f"{RUN}/ac_check.json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
    if bad:
        abandon("autocorrelation", f"{len(bad)} hours recorded zero pairs for non-veto reasons "
                                   f"(e.g. {bad[:3]}) -- the engine fixes did not take")
    if burst_in and burst_pairs and burst_pairs > 0:
        abandon("autocorrelation", f"2005-07-09 23h survived the earthquake screen and still "
                                   f"produced {burst_pairs} quiet pairs above net_mad 8")
    ok("ac", f"{len(ran)} hours, {len(zero)} zero-pair hours all explained by vetoes; burst hour "
             f"{'in' if burst_in else 'excluded from'} the quiet set")


def c_tpair():
    j = json.load(open(f"{RUN}/T.json"))
    if j["verdict"] != "ACCEPT":
        abandon("T_pair", j["reason"] + " -- inspect the dominating hour before proceeding")
    ok("T_pair", f"T = {j['T']:.2f} ({j['binds']}); LOO range {j['loo_range']:.2f}")


def c_seeds():
    import v_common as V
    W.patch_v()
    tids, stas, Wt, off, meta = V.load_bank(f"{RUN}/bank0.npz")
    Sd = W.read_csv_safe(f"{RUN}/cl_seeds.csv")
    print(f"[check/seeds] {len(Sd)} seed clusters -> {len(tids)} handed off with >= "
          f"{W.MIN_SITES} sites")
    if len(tids) < MIN_SEEDS_HANDED_OFF:
        abandon("seeds", f"only {len(tids)} seeds handed off with >= {W.MIN_SITES} sites "
                         f"(accept >= {MIN_SEEDS_HANDED_OFF}); the pair-seeded construction has "
                         f"no material to grow")
    ok("seeds", f"{len(tids)} templates in bank0")


def c_oracle():
    j = json.load(open(f"{RUN}/oracle_reference_2005.json"))
    print(json.dumps(j, indent=1, default=str))
    if int(j["C_30"]) < 3:
        abandon("oracle ceiling", f"C_30 = {j['C_30']} (< 3): the frozen chain cannot recover even "
                                  f"perfectly seeded families on this bed -- the control bed is "
                                  f"unusable, per Merlin's 'what would change my mind'")
    ok("oracle", f"C_iso {j['C_iso']}/10, C_30 {j['C_30']}/30, fallback "
                 f"{j['fallback_active']}, bar {j['pass_bar']['value']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["d8", "health", "quiet", "ac", "tpair", "seeds", "oracle"])
    a = ap.parse_args()
    {"d8": c_d8, "health": c_health, "quiet": c_quiet, "ac": c_ac, "tpair": c_tpair,
     "seeds": c_seeds, "oracle": c_oracle}[a.cmd]()


if __name__ == "__main__":
    main()
