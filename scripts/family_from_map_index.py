#!/usr/bin/env python
"""Map the borehole-map's family slider label "N/M" to the stored family id — EXACTLY as the map
builds its list (build_borehole_dvv_map.py): preferred CSV (_unif05 > _eps05 > original), certified
families only, cc_max > 0.70, ids sorted lexicographically; the display is 1-BASED
(DATA[s].fam[cur.fi-1]; slider 0 = station overview).

Usage: python scripts/family_from_map_index.py STATION N[/M] [N[/M] ...]
The /M part, when given, is checked against the station's family count as a safety interlock.
"""
import sys, os, glob
import pandas as pd

def preferred_csv(tag):
    for suf in ("_unif05", "_eps05", ""):
        p = f"data/daily_dvv_{tag}_Z_2to4{suf}.csv"
        if os.path.exists(p):
            return p
    raise SystemExit(f"no dv/v CSV for {tag}")

def cert_set(tag):
    stem = tag.lower().replace("p90f40", "")
    p = f"data/{stem}_causality_cert.csv"
    if os.path.exists(p):
        c = pd.read_csv(p)
        return set(c[c.reliable].fam)
    p = f"data/{stem}_fwd_vs_rev_coda.csv"
    if os.path.exists(p):
        c = pd.read_csv(p)
        return set(c[c.ratio > 1.5].fam)
    raise SystemExit(f"no certification table for {tag}")

def map_family_list(station):
    """The map's per-station family list, in the map's exact order."""
    tags = [os.path.basename(f).split("daily_dvv_")[1].split("_Z_2to4")[0]
            for f in glob.glob("data/daily_dvv_*_Z_2to4.csv")]
    tag = next((t for t in tags if t.upper().replace("P90F40", "") == station.upper()), None)
    if tag is None:
        raise SystemExit(f"station {station} not found among {len(tags)} stations")
    # replicate build_borehole_dvv_map.py EXACTLY: preferred file; cycle-skip prefilter on
    # unif05/eps05; cc_max > 0.6 (CCMIN); certified minus user blocklist; family kept only with
    # >= 60 distinct days (MINDAYS); ids sorted lexicographically.
    fp = preferred_csv(tag)
    d = pd.read_csv(fp, usecols=["patch", "date", "dvv", "cc_max"])
    if fp.endswith(("_unif05.csv", "_eps05.csv")):
        d = d[~((d.dvv.abs() > 0.02) & (d.cc_max < 0.80))]
    fams = cert_set(tag)
    bp = "data/family_user_blocklist.csv"
    if os.path.exists(bp):
        b = pd.read_csv(bp)
        stem = tag.lower().replace("p90f40", "")
        fams = fams - set(b[b.stem == stem].fam)
    d = d[(d.cc_max > 0.6) & (d.patch.isin(fams))]
    nd = d.groupby("patch").date.nunique()
    return tag, sorted(nd[nd >= 60].index)

if __name__ == "__main__":
    station, picks = sys.argv[1], sys.argv[2:]
    tag, fams = map_family_list(station)
    print(f"{station} (file tag {tag}): {len(fams)} certified families in map order")
    for p in picks:
        n, m = (p.split("/") + [None])[:2]
        n = int(n)
        if m is not None and int(m) != len(fams):
            print(f"  {p}: MISMATCH — map shows /{m} but this station has {len(fams)} families. "
                  f"Wrong station or stale map — NOT resolving.")
            continue
        if not 1 <= n <= len(fams):
            print(f"  {p}: index out of range 1..{len(fams)}")
            continue
        print(f"  {n}/{len(fams)} -> {fams[n-1]}")
