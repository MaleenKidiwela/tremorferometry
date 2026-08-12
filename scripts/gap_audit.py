#!/usr/bin/env python
"""Gap audit: per-year day coverage of every completed station's dv/v, flag INTERIOR years with
suspiciously low coverage (a truncated/gapped year fakes a dv/v discontinuity — the 2017 B011/B926 bug).
First and last year are naturally partial and never flagged. Prints one FLAG line per (station, year);
if traces are on disk it reports the raw trace-day count too (distinguishes a re-downloadable download
truncation from a genuine station outage). Exit 0 always; grep 'FLAG' for the hits.
Usage: python scripts/gap_audit.py [min_days=200]
"""
import sys, glob, os, json, pandas as pd

MIN = int(sys.argv[1]) if len(sys.argv) > 1 else 200
CAP = 300          # selection hard cap (user decision) — flag any station that exceeded it
flags = []

# over-cap (wasted densify) + weak (<20 cert) checks, per finalized station
for jf in sorted(glob.glob('data/*_3comp_summary.json')):
    d = json.load(open(jf)); s = d['station']; sl = s.lower()
    zc = d.get('z_certified')
    self = f'data/{sl}_disc_p70_2010_2026_m3_sel300.summary.csv'
    if os.path.exists(self):
        nsel = sum(1 for _ in open(self)) - 1
        if nsel > CAP + 10:
            flags.append(f"NOTE {s}: densified {nsel} > cap {CAP} (extra COMPUTE spent; dv/v correctly uses "
                         f"all reliable/causality-certified families — do NOT re-run or recap)")
    if zc is not None and zc < 20:
        flags.append(f"FLAG {s}: only {zc} certified families (<20) — weak/starved")

for f in sorted(glob.glob('data/daily_dvv_*p90f40_Z_2to4.csv')):
    sta = os.path.basename(f).split('daily_dvv_')[1].split('p90f40')[0]
    dv = pd.read_csv(f, usecols=['date']); dv['date'] = pd.to_datetime(dv.date)
    yr = dv.groupby(dv.date.dt.year).date.agg(lambda x: x.dt.normalize().nunique())
    years = sorted(yr.index)
    for y in years[1:-1]:                        # interior years only
        if yr[y] < MIN:
            tdir = f'data/waveforms/PB.{sta}/{y}'
            traw = len(glob.glob(f'{tdir}/*.mseed')) if os.path.isdir(tdir) else None
            hint = f'traces on disk: {traw} day-files' if traw is not None else 'traces deleted (re-download to check/fix)'
            flags.append(f"FLAG {sta} {y}: dv/v {int(yr[y])} days (<{MIN}) — {hint}")

for line in flags:
    print(line)
print(f"[gap_audit] {len(flags)} flagged (station,year) across {len(glob.glob('data/daily_dvv_*p90f40_Z_2to4.csv'))} stations")
