#!/usr/bin/env python
"""Aggregate the per-station noise-repeater census and DEFINE the quiet-subset 2008 bed.  `p_`.

Census columns (per station, from p_census.py's station-block rows):
  per_hour_quiet / per_hour_tremor  detections of the NETWORK in which this station is above its
                                    own 8 x MAD bar, per block it took part in;
  ratio                             per_hour_tremor / per_hour_quiet;
  top3_share_quiet                  fraction of the quiet-hour detections of the blocks it took
                                    part in where it is among the top 3 contributors;
  noise_db                          P4-0 statistic: the station's median 2-8 Hz hourly power over
                                    the literature-quiet days, in dB relative to the QUIETEST PB
                                    borehole (results/n_puget2008/hour_energy28_2008.parquet).

SUBSET RULE (fixed here, before any subset result exists; the brief's rule):
  the 10 PB boreholes that are in the bed, PLUS the CAFE (XU) / YW stations with ratio >= 3 and a
  top3_share_quiet at most TOP3_BAR, ranked by (ratio desc, top3_share_quiet asc), so that the
  subset has at least MIN_N and at most MAX_N stations.  TOP3_BAR = 2 x the chance share
  (3 / n_reporting_stations); if too few CAFE stations clear both bars to reach MIN_N, the
  next-best by the same ranking are added and the shortfall is RECORDED in the meta file.

Coverage is reported exactly as results/P_STEP2_SITES.md §4 does: the % of the sub-region's
active cells with >= 5 subset stations within 60 km (hypocentral, and epicentral for comparison).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
sys.path.insert(0, HERE)
import n_common as L                                                        # noqa: E402

MIN_N, MAX_N = 12, 20
RATIO_BAR = 3.0


def agg(p, tag):
    if not os.path.exists(p):
        return pd.DataFrame(columns=["sta"])
    B = pd.read_parquet(p)
    g = B.groupby("sta").agg(**{
        f"blocks_{tag}": ("date", "size"),
        f"det_{tag}": ("n_det", "sum"),
        f"above_{tag}": ("n_above", "sum"),
        f"above6_{tag}": ("n_above6", "sum"),
        f"above4_{tag}": ("n_above4", "sum"),
        f"zmean_{tag}": ("z_mean", "mean"),
        f"top3_{tag}": ("n_top3", "sum"),
        f"medcc_{tag}": ("med_cc", "median"),
        f"madcc_{tag}": ("mad_cc", "median"),
    })
    g[f"per_hour_{tag}"] = g[f"above_{tag}"] / g[f"blocks_{tag}"].clip(lower=1)
    g[f"per_hour4_{tag}"] = g[f"above4_{tag}"] / g[f"blocks_{tag}"].clip(lower=1)
    g[f"top3_share_{tag}"] = g[f"top3_{tag}"] / g[f"det_{tag}"].clip(lower=1)
    return g.reset_index()


def noise_db(stas_netsta):
    """P4-0-style quiet-day 2-8 Hz noise, dB relative to the quietest PB borehole."""
    p = f"{L.OUTROOT}/hour_energy28_2008.parquet"
    q = [x for x in open(f"{L.OUTROOT}/quiet_dates_2008.txt").read().split(",") if x]
    E = pd.read_parquet(p)
    E = E[E.date.isin(set(q))]
    hc = [c for c in E.columns if c.startswith("h") and len(c) == 3]
    M = E.melt(id_vars=["netsta"], value_vars=hc, value_name="pw").dropna(subset=["pw"])
    M = M[M.pw > 0]
    g = M.groupby("netsta").pw.median()
    pb = [n for n in g.index if n.startswith("PB.")]
    ref = float(g[pb].min()) if pb else float(g.min())
    return {n: float(10 * np.log10(g[n] / ref)) for n in g.index if n in set(stas_netsta)}, ref


def greedy_coverage(base, cand, ST, AC, room):
    """Pick `room` of `cand` one at a time, each time the one that raises the number of active
    cells with >= 5 stations within 60 km (hypocentral) the most; ties go to the earlier (i.e.
    higher-ratio) candidate."""
    import s_common as S
    cla, clo, cdp = AC.lat.to_numpy(), AC.lon.to_numpy(), AC.dep.to_numpy()

    def near(s):
        d = S.epi_km(float(ST.lat[s]), float(ST.lon[s]), cla, clo)
        return (np.sqrt(d ** 2 + cdp ** 2) <= 60.0).astype(int)

    NE = {s: near(s) for s in set(base) | set(cand)}
    cnt = sum((NE[s] for s in base), np.zeros(len(AC), int))
    out, left = [], list(cand)
    while left and len(out) < room:
        best, bestv = None, -1
        for s in left:
            v = int(((cnt + NE[s]) >= 5).sum())
            if v > bestv:
                best, bestv = s, v
        out.append(best)
        cnt = cnt + NE[best]
        left.remove(best)
    return out


def coverage(stas, ST, AC):
    import s_common as S
    cla, clo, cdp = AC.lat.to_numpy(), AC.lon.to_numpy(), AC.dep.to_numpy()
    n_hyp = np.zeros(len(AC), int)
    n_epi = np.zeros(len(AC), int)
    for s in stas:
        d = S.epi_km(float(ST.lat[s]), float(ST.lon[s]), cla, clo)
        n_epi += (d <= 60.0).astype(int)
        n_hyp += (np.sqrt(d ** 2 + cdp ** 2) <= 60.0).astype(int)
    return dict(n_cells=int(len(AC)),
                pct_cells_ge5_60km_hyp=float(100 * np.mean(n_hyp >= 5)),
                pct_cells_ge5_60km_epi=float(100 * np.mean(n_epi >= 5)),
                median_sites_60km_hyp=float(np.median(n_hyp)),
                median_sites_60km_epi=float(np.median(n_epi)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tremor", required=True)
    ap.add_argument("--quiet", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    L.activate()
    KEEP = L.keepers()[["net", "sta", "netsta"]]
    T, Q = agg(a.tremor, "tremor"), agg(a.quiet, "quiet")
    C = KEEP.merge(T, on="sta", how="left").merge(Q, on="sta", how="left")
    # +0.5 counts (Laplace) on BOTH sides: the 8 x MAD contribution count is sparse -- about
    # 0.6 stations per network detection across 50 stations -- so a raw quotient of two small
    # integers is unstable and undefined at zero.  ratio4 (the 4 x MAD count, ~10x denser) is
    # reported beside it as the stability check; the subset rule uses `ratio`, as briefed.
    C["ratio"] = ((C.above_tremor + 0.5) / C.blocks_tremor) / \
                 ((C.above_quiet + 0.5) / C.blocks_quiet)
    C["ratio4"] = ((C.above4_tremor + 0.5) / C.blocks_tremor) / \
                  ((C.above4_quiet + 0.5) / C.blocks_quiet)
    C["ratio_top3"] = (C.top3_tremor / C.det_tremor) / (C.top3_quiet / C.det_quiet)
    nd, ref = noise_db(set(C.netsta))
    C["noise_db"] = [nd.get(n, np.nan) for n in C.netsta]
    nrep = float(np.nanmedian(C.blocks_quiet.notna().sum()))
    chance = TOP3_CHANCE = 3.0 / max(len(C), 1)
    C["top3_excess_quiet"] = C.top3_share_quiet / chance
    C["noise_repeater"] = (C.ratio < 2) & (C.top3_share_quiet > 2 * chance)
    C = C.sort_values("ratio", na_position="first").reset_index(drop=True)
    C.to_csv(f"{a.outdir}/station_noise_census.csv", index=False)

    by_net = (C.groupby("net")
               .agg(n=("sta", "size"),
                    median_ratio=("ratio", "median"),
                    median_top3_share_quiet=("top3_share_quiet", "median"),
                    median_per_hour_quiet=("per_hour_quiet", "median"),
                    median_noise_db=("noise_db", "median"),
                    n_noise_repeater=("noise_repeater", "sum"))
               .reset_index())
    by_net.to_csv(f"{a.outdir}/census_by_network.csv", index=False)

    # ------------------------------------------------------------------ the subset -------------
    top3_bar = 2 * chance
    pb = sorted(C[C.net == "PB"].sta)
    cand = C[(C.net != "PB") & C.ratio.notna()].copy()
    cand["qual"] = (cand.ratio >= RATIO_BAR) & (cand.top3_share_quiet <= top3_bar)
    cand = cand.sort_values(["qual", "ratio", "top3_share_quiet"],
                            ascending=[False, False, True])
    ST, AC = L.station_table(), L.active_cells()
    room = max(MAX_N - len(pb), 0)
    qual = list(cand[cand.qual].sta)
    short = 0
    if len(qual) <= room:
        pick = qual                       # quality decides; there is room for all of them
    else:
        # more stations clear the QUALITY bars than there are slots -> the slots go to the ones
        # that add the most CELL COVERAGE (Merlin's rule: fill by quality x distinct coverage,
        # never by count).  Quality is already decided by the bars; this only breaks the tie.
        pick = greedy_coverage(pb, qual, ST, AC, room)
    if len(pb) + len(pick) < MIN_N:
        extra = [s for s in cand.sta if s not in set(pick)][:MIN_N - len(pb) - len(pick)]
        short = len(extra)
        pick += extra
    subset = sorted(set(pb) | set(pick))
    cov_sub = coverage(subset, ST, AC)
    cov_all = coverage(list(C.sta), ST, AC)
    cov_pb = coverage(pb, ST, AC)
    ns = dict(zip(C.sta, C.netsta))
    meta = dict(rule=dict(min_n=MIN_N, max_n=MAX_N, ratio_bar=RATIO_BAR,
                          top3_bar=float(top3_bar), top3_chance=float(chance)),
                n_subset=len(subset), subset=[ns[s] for s in subset],
                n_pb=len(pb), pb=[ns[s] for s in pb],
                n_cafe_picked=len(pick), cafe=[ns[s] for s in sorted(pick)],
                n_cafe_qualified=int(cand.qual.sum()),
                n_below_bar_added=short,
                shortfall_note=("" if short == 0 else
                                f"{short} station(s) added below the ratio/top3 bars to reach the "
                                f"{MIN_N}-station floor; they are the best-ranked of the rest"),
                coverage_subset=cov_sub, coverage_51=cov_all, coverage_pb_only=cov_pb,
                noise_ref_netsta_db0=float(ref),
                n_noise_repeater=int(C.noise_repeater.sum()),
                worst=[dict(netsta=r.netsta, ratio=None if not np.isfinite(r.ratio) else
                            float(r.ratio), top3=float(r.top3_share_quiet),
                            per_hour_quiet=float(r.per_hour_quiet))
                       for r in C.head(8).itertuples()],
                best=[dict(netsta=r.netsta, ratio=float(r.ratio),
                           top3=float(r.top3_share_quiet),
                           per_hour_quiet=float(r.per_hour_quiet))
                      for r in C.dropna(subset=["ratio"]).tail(8).iloc[::-1].itertuples()])
    json.dump(meta, open(f"{a.outdir}/subset_meta.json", "w"), indent=1, default=float)
    open(f"{a.outdir}/subset_stations.txt", "w").write(",".join(meta["subset"]))
    open(f"{a.outdir}/boreholes_stations.txt", "w").write(",".join(meta["pb"]))
    L.log(json.dumps({k: v for k, v in meta.items() if k not in ("worst", "best")},
                     indent=1, default=float))
    L.log("\n" + by_net.to_string(index=False))


if __name__ == "__main__":
    main()
