"""Lin-seeded family discovery at PGC + NLLB (network autocorrelation).

Per the §8.2 plan, branch A of the NLLB discovery pipeline:

1. Subset Lin's catalog to the southern V.I. bbox.
2. Cluster Lin detections into proto-families by lat/lon grid (0.05 deg).
3. For each proto-family:
   a. Cut envelope-aligned 2-s window at PGC and NLLB around each Lin OT
      (each station finds its own envelope peak in [OT+5, OT+17] s).
   b. Keep only detections with valid cuts at BOTH stations.
   c. All-pairs max-shifted CC at PGC and at NLLB; network CC = mean.
   d. Complete-linkage cluster at network CC >= 0.80.
   e. Keep clusters with >=3 members spanning >=3 years.
4. Save each surviving family: PGC template, NLLB template, member times,
   centroid lat/lon.

This is the Shelly-Beroza-Brown recipe applied to PGC+NLLB; same idea
as the original E.7e but with NLLB instead of LZB.
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from tremorferometry.repeater import (  # noqa: E402
    all_pairs_cc_max_shifted,
    cluster_matches,
    cut_all_detections,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--lin-csv", default="data/raw_lfe/lin2023_lfe.csv")
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--stations", nargs="+", default=["PGC", "NLLB"])
    p.add_argument("--bbox", nargs=4, type=float,
                   default=[48.0, 50.0, -125.5, -122.5],
                   help="lat_min lat_max lon_min lon_max")
    p.add_argument("--grid-deg", type=float, default=0.05)
    p.add_argument("--min-proto-detections", type=int, default=20)
    p.add_argument("--max-per-proto", type=int, default=2000,
                   help="Sub-sample proto-families with more than this many "
                        "detections to keep CC matrices tractable")
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--search-pre", type=float, default=5.0,
                   help="Envelope search window: [OT+pre, OT+post]")
    p.add_argument("--search-post", type=float, default=17.0)
    p.add_argument("--template-pre", type=float, default=-1.0)
    p.add_argument("--template-post", type=float, default=1.0)
    p.add_argument("--max-shift-samples", type=int, default=20)
    p.add_argument("--cc-threshold", type=float, default=0.80)
    p.add_argument("--min-family-members", type=int, default=3)
    p.add_argument("--min-years", type=int, default=3)
    p.add_argument("--out", default="data/nllb_lin_seeded_families.npz")
    p.add_argument("--report", default="data/nllb_lin_seeded_report.csv")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    print(f"[1/4] Loading Lin catalog and filtering to V.I. bbox...")
    lin = pd.read_csv(args.lin_csv)
    lin["OT"] = pd.to_datetime(lin["OT"])
    lat_min, lat_max, lon_min, lon_max = args.bbox
    mask = ((lin["lat"] >= lat_min) & (lin["lat"] <= lat_max) &
            (lin["lon"] >= lon_min) & (lin["lon"] <= lon_max))
    lin = lin[mask].copy()
    print(f"  {len(lin):,} Lin detections in bbox "
          f"({lat_min}-{lat_max} lat, {lon_min}-{lon_max} lon)")
    print(f"  date range: {lin.OT.min()} .. {lin.OT.max()}")

    print(f"[2/4] Clustering into proto-families at {args.grid_deg} deg grid...")
    lin["lat_bin"] = np.round(lin["lat"] / args.grid_deg) * args.grid_deg
    lin["lon_bin"] = np.round(lin["lon"] / args.grid_deg) * args.grid_deg
    lin["proto"] = (lin["lat_bin"].map(lambda v: f"{v:.3f}") + "_"
                    + lin["lon_bin"].map(lambda v: f"{v:.3f}"))
    proto_counts = lin["proto"].value_counts()
    keep_protos = proto_counts[proto_counts >= args.min_proto_detections].index.tolist()
    print(f"  {len(proto_counts)} proto-families, "
          f"{len(keep_protos)} with >= {args.min_proto_detections} detections")

    bandpass = (args.fmin, args.fmax)
    search_window = (args.search_pre, args.search_post)
    out_window = (args.template_pre, args.template_post)
    n_template = int(round((out_window[1] - out_window[0]) * args.fs))

    families = []          # list of dicts: proto, family_id, members (OT), templates per station, lat, lon
    report_rows = []

    print(f"[3/4] Running PGC+NLLB network autocorrelation per proto-family...")
    n_total_kept = 0
    for k, proto in enumerate(keep_protos):
        sub = lin[lin["proto"] == proto].reset_index(drop=True)
        if len(sub) > args.max_per_proto:
            sub = sub.iloc[rng.choice(len(sub), args.max_per_proto, replace=False)]
            sub = sub.sort_values("OT").reset_index(drop=True)
        n_proto = len(sub)
        t0 = time.time()

        # Cut at each station independently
        Xs = []
        oks = []
        for station in args.stations:
            X, peaks, df_ok = cut_all_detections(
                sub.rename(columns={"OT": "OT"}),
                Path(args.wfdir), station, fs=args.fs, bandpass=bandpass,
                search_window=search_window, out_window=out_window,
            )
            # Map back to original indices in sub (df_ok already filtered)
            ok = np.zeros(n_proto, dtype=bool)
            ok[df_ok.index.values if hasattr(df_ok, "index") else np.arange(len(df_ok))] = True
            # Actually we lost the original indices. Recompute by matching OT.
            ok = sub["OT"].isin(df_ok["OT"]).values
            Xs.append((X, ok))
            oks.append(ok)
        # detections valid at BOTH stations
        both = np.ones(n_proto, dtype=bool)
        for ok in oks:
            both &= ok
        n_both = int(both.sum())
        if n_both < max(args.min_family_members, 3):
            report_rows.append(dict(
                proto=proto, n_proto=n_proto, n_both=n_both,
                n_families=0, status="too few both-valid"))
            print(f"  [{k+1}/{len(keep_protos)}] {proto} n={n_proto} both={n_both}: skip")
            continue

        # Re-assemble per-station X aligned to `both` mask
        # X from cut_all_detections is in sub-original index order via df_ok;
        # safest: rebuild aligned X by going through sub and pulling the row
        # for each both-True index.
        X_per_sta = []
        for (X, ok), station in zip(Xs, args.stations):
            X_aligned = np.zeros((n_both, n_template), dtype=np.float32)
            # X has rows in df_ok order; df_ok rows correspond to sub indices where ok=True
            df_ok_idx = np.flatnonzero(ok)
            both_idx = np.flatnonzero(both)
            # Map each both_idx to its position in df_ok_idx
            pos = np.searchsorted(df_ok_idx, both_idx)
            X_aligned = X[pos]
            X_per_sta.append(X_aligned)

        # Per-station all-pairs CC
        cc_sum = np.zeros((n_both, n_both), dtype=np.float32)
        for X in X_per_sta:
            cc = all_pairs_cc_max_shifted(X, max_shift_samples=args.max_shift_samples)
            cc_sum += cc
        network_cc = cc_sum / len(X_per_sta)

        # Cluster
        labels = cluster_matches(network_cc, threshold=args.cc_threshold)
        n_clusters = int(labels.max() + 1) if labels.size and labels.max() >= 0 else 0

        # Filter clusters
        n_kept = 0
        for c in range(n_clusters):
            members_local = np.flatnonzero(labels == c)
            if members_local.size < args.min_family_members:
                continue
            sub_both_indices = np.flatnonzero(both)
            member_orig_idx = sub_both_indices[members_local]
            member_OT = pd.to_datetime(sub.iloc[member_orig_idx]["OT"])
            n_years = member_OT.dt.year.nunique()
            if n_years < args.min_years:
                continue
            # Templates: mean of L2-normalized members, then re-L2-normalize
            templates = {}
            for X, station in zip(X_per_sta, args.stations):
                T = X[members_local].mean(axis=0)
                T = T / (np.linalg.norm(T) + 1e-12)
                templates[station] = T.astype(np.float32)
            fam = dict(
                proto=proto, family_id=f"{proto}__c{c}",
                n_members=int(members_local.size),
                n_years=int(n_years),
                year_span=int(member_OT.dt.year.max() - member_OT.dt.year.min()),
                lat=float(sub.iloc[member_orig_idx]["lat"].mean()),
                lon=float(sub.iloc[member_orig_idx]["lon"].mean()),
                first_year=int(member_OT.dt.year.min()),
                last_year=int(member_OT.dt.year.max()),
                templates=templates,
                member_OT=member_OT.astype(str).tolist(),
            )
            families.append(fam)
            n_kept += 1
            n_total_kept += 1

        dt_s = time.time() - t0
        report_rows.append(dict(
            proto=proto, n_proto=n_proto, n_both=n_both,
            n_families=n_kept, status="ok", seconds=round(dt_s, 1)))
        print(f"  [{k+1}/{len(keep_protos)}] {proto} n={n_proto} "
              f"both={n_both} -> {n_kept} families  ({dt_s:.1f}s)")

    print(f"[4/4] Saving {n_total_kept} discovered families to {args.out}...")
    # Save as npz: per station + per family_id => template
    npz_arrays = {}
    for f in families:
        for station, T in f["templates"].items():
            npz_arrays[f"{f['family_id']}__{station}"] = T
    np.savez(args.out, **npz_arrays)
    # Save a member list as parquet
    rows = []
    for f in families:
        for ot in f["member_OT"]:
            rows.append(dict(family_id=f["family_id"], OT=ot,
                             lat=f["lat"], lon=f["lon"]))
    members_df = pd.DataFrame(rows)
    members_df.to_parquet(Path(args.out).with_suffix(".members.parquet"), index=False)
    # Save a per-family summary
    summary_rows = []
    for f in families:
        summary_rows.append({k: v for k, v in f.items()
                             if k not in ("templates", "member_OT")})
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(Path(args.out).with_suffix(".summary.csv"), index=False)
    pd.DataFrame(report_rows).to_csv(args.report, index=False)
    print(f"  saved templates -> {args.out}")
    print(f"  saved members -> {Path(args.out).with_suffix('.members.parquet')}")
    print(f"  saved summary -> {Path(args.out).with_suffix('.summary.csv')}")
    print(f"  saved per-proto report -> {args.report}")


if __name__ == "__main__":
    main()
