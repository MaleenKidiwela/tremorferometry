#!/usr/bin/env python
"""Step 1c -- BED SETUP.  Decides the window, the quiet days, the sites and the withheld set ONCE,
writes results/step1c_bed.json, and is never re-run after scoring starts.

Window.  The published southern-Vancouver-Island ETS of September 2005 (the 14-month sequence
Feb 2003 / Jul 2004 / Sep 2005 / Jan 2007), taken as 2005-09-03..09-16 and padded +- 5 d ->
2005-08-29..2005-09-21.  The WTD (Cascadia) tremor catalogue could NOT be retrieved (see
w_energy_hours.py), so hour selection uses the array-energy flag and the window comes from the
literature.  NO Bostock event time is used for the window, the hours or the quiet days.

Quiet days.  2005-06-25..2005-07-14 (20 days), >= 50 d before the ETS onset and >= 11 months after
the July 2004 episode -- "literature-quiet", recorded as such.

Sites.  1-km single-linkage merge of the DELIVERING stations (w_probe2005.py), then the SAME
written rule as Step 1b: add sites in order of increasing median epicentral distance to the 130
Bostock family HYPOCENTRES (coordinates only) while the median S-moveout spread over those
families stays <= AC_WIN_S (the frozen 10-s common window).

Withheld.  VGZ and SNB (as in every earlier step, minus B926, which does not exist in 2005) plus
TWO dense stations chosen by the coverage rule: the pair whose removal leaves the most families
with >= 5 detection sites within 60 km epicentral.
"""
from __future__ import annotations
import argparse, itertools, json, os, sys

import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

ETS = ("2005-09-03", "2005-09-16")
PAD_D = 5
QUIET = ("2005-06-25", "2005-07-14")
ETS_SOURCE = ("literature: southern Vancouver Island ETS of 2005-09 (14-month sequence "
              "2003-02 / 2004-07 / 2005-09 / 2007-01); WTD catalogue not retrievable")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", default=f"{W.RES}/step1c_stations.csv")
    ap.add_argument("--out", default=f"{W.RES}/step1c_bed.json")
    ap.add_argument("--min-days", type=int, default=2)
    a = ap.parse_args()
    import p_vel1d as V1
    F = S.families()
    D = pd.read_csv(a.probe)
    D = D[(D.n_days >= a.min_days) & np.isfinite(D.epi_med)].copy()
    ST = W.station_table()
    for r in D.itertuples():                      # make sure every delivering station has coords
        if r.sta not in ST.index or not np.isfinite(ST.lat.get(r.sta, np.nan)):
            ST.loc[r.sta, "lat"], ST.loc[r.sta, "lon"] = float(r.lat), float(r.lon)
    stas = sorted(D.sta)
    mem, of, order = W.merge_sites(stas, ST=ST)
    # per-site geometry and travel times (family COORDINATES only)
    geo, TT = {}, {}
    for q in order:
        s0 = mem[q][0]
        la, lo = float(ST.lat[s0]), float(ST.lon[s0])
        e = S.epi_km(F.lat.values, F.lon.values, la, lo)
        geo[q] = dict(lat=la, lon=lo, epi_med=float(np.median(e)),
                      n_fam_epi25=int((e <= 25).sum()), members=mem[q],
                      epi=e)
        TT[q] = np.array([V1.tt(float(x), float(z), "S") for x, z in zip(e, F.dep.values)])

    def spread(qs):
        A = np.stack([TT[q] for q in qs])
        return float(np.median(A.max(0) - A.min(0)))

    withheld_fixed = [x for x in ("VGZ", "SNB") if x in of]
    pool = [q for q in sorted(order, key=lambda q: geo[q]["epi_med"])
            if q not in withheld_fixed and all(m not in withheld_fixed for m in mem[q])]
    # --- the written spread rule
    chosen = []
    for q in pool:
        if not chosen:
            chosen = [q]
            continue
        if spread(chosen + [q]) <= W.AC_WIN_S:
            chosen.append(q)
    # --- withheld dense pair by the coverage rule
    def coverage(qs):
        n = 0
        for i in range(len(F)):
            c = sum(1 for q in qs if geo[q]["epi"][i] <= 60.0)
            n += (c >= 5)
        return n

    dense = [q for q in chosen if geo[q]["epi_med"] <= 45.0]
    best, bestcov = None, -1
    for pair in itertools.combinations(dense, 2):
        keep = [q for q in chosen if q not in pair]
        if len(keep) < W.MIN_SITES:
            continue
        c = coverage(keep)
        if c > bestcov:
            best, bestcov = pair, c
    if best is None:
        best = tuple(sorted(chosen, key=lambda q: geo[q]["epi_med"])[-2:])
    final = [q for q in chosen if q not in best]
    wh = withheld_fixed + [m for q in best for m in mem[q]]
    d0 = (pd.Timestamp(ETS[0]) - pd.Timedelta(days=PAD_D)).strftime("%Y-%m-%d")
    d1 = (pd.Timestamp(ETS[1]) + pd.Timedelta(days=PAD_D)).strftime("%Y-%m-%d")
    qd = [d.strftime("%Y-%m-%d") for d in pd.date_range(QUIET[0], QUIET[1], freq="D")]
    cfg = dict(
        year=2005, day0=d0, day1=d1, ets=list(ETS), pad_days=PAD_D, ets_source=ETS_SOURCE,
        quiet_days=qd, quiet_rule="literature-quiet: >= 50 d before the 2005-09 ETS onset; no "
                                  "tremor catalogue exists for 2005",
        hour_source="array-energy flag (w_energy_hours.py); Merlin Q3 -- hour selection only",
        sites=sorted(final), site_keys=sorted(mem[q][0] for q in final),
        site_members={q: mem[q] for q in final},
        site_id_note="site ids are a data-driven tie-break and are NOT stable across cache roots; "
                     "site_keys (one MEMBER station per site) is what the engine matches on",
        withheld=sorted(set(wh)),
        withheld_dense_sites=list(best), withheld_rule="VGZ+SNB (B926 does not exist in 2005) plus "
        "the dense PAIR whose removal leaves the most families with >= 5 detection sites within "
        "60 km epicentral",
        site_rule=f"1-km single linkage, then add sites by increasing median epicentral distance "
                  f"while the median S-moveout spread stays <= AC_WIN_S = {W.AC_WIN_S} s",
        moveout_spread_s_chosen=spread(chosen), moveout_spread_s_final=spread(final),
        coverage_families_ge5sites_60km=int(coverage(final)),
        n_families=int(len(F)),
        site_table=[dict(site=q, members=mem[q], lat=geo[q]["lat"], lon=geo[q]["lon"],
                         epi_med=round(geo[q]["epi_med"], 1),
                         n_fam_epi25=geo[q]["n_fam_epi25"],
                         detection=bool(q in final)) for q in order],
        stations_download=[dict(net=r.net, sta=r.sta, band=r.band, prov=r.prov,
                                n_days_probed=int(r.n_days), epi_med=float(r.epi_med))
                           for r in D.itertuples()],
    )
    json.dump(cfg, open(a.out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in cfg.items()
                      if k not in ("site_table", "stations_download", "quiet_days")},
                     indent=1, default=str))
    print("\nSITES (detection = enters every sum):")
    for r in cfg["site_table"]:
        print(f"  {'*' if r['detection'] else ' '} {r['site']:<7} {'+'.join(r['members']):<28} "
              f"epi_med {r['epi_med']:5.1f} km  n_fam<=25km {r['n_fam_epi25']:3d}")
    print(f"\nwindow {d0}..{d1} ({(pd.Timestamp(d1)-pd.Timestamp(d0)).days+1} d) + "
          f"{len(qd)} quiet days; {len(D)} stations -> "
          f"{len(D)*((pd.Timestamp(d1)-pd.Timestamp(d0)).days+1+len(qd))} station-days")


if __name__ == "__main__":
    main()
