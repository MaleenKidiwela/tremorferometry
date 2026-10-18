#!/usr/bin/env python
"""Append the dated NOTES.md entry for the quiet-subset Puget/S-Olympic 2008 run.  Prefix `p_`."""
from __future__ import annotations

import datetime
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)


def j(p):
    try:
        return json.load(open(p))
    except Exception:                                                       # noqa: BLE001
        return {}


def f(x, n=2):
    try:
        return f"{float(x):.{n}f}"
    except Exception:                                                       # noqa: BLE001
        return "--"


def main():
    p, r, nu, c7, bh = sys.argv[1:6]
    sm = j(f"{p}/subset_meta.json")
    ps51 = j(f"{STUDY}/results/n_puget2008/2008/pair_pairstats.json")
    ps, psn, psb = (j(f"{r}/pair_pairstats.json"), j(f"{nu}/pair_pairstats.json"),
                    j(f"{bh}/pair_pairstats.json"))
    sr, sn = j(f"{r}/screen_summary.json"), j(f"{nu}/screen_summary.json")
    fr, fn = j(f"{r}/screen_final_summary.json"), j(f"{nu}/screen_final_summary.json")
    cr, cn = j(f"{c7}/confirm_real_summary.json"), j(f"{c7}/confirm_null_summary.json")
    bed = j(f"{p}/bed_2008q.json")
    try:
        import pandas as pd
        n_cens = len(pd.read_csv(f"{p}/station_noise_census.csv"))
    except Exception:                                                       # noqa: BLE001
        n_cens = "--"
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    txt = f"""
## 2026-09-21 — agent P: Puget / S-Olympic 2008, the QUIET-SUBSET bed (station quality as the lever)
Agent N's 51-station CAFE+PB run gave a pair-stage tremor/quiet ratio of
{f(ps51.get('tremor_quiet_ratio'))} and grown templates firing equally on quiet and tremor days —
the pair stage was dominated by repeating cultural/site noise.  Royer & Bostock (2014) built
COMPOSITE ARRAYS out of chosen stations, so this run tests station QUALITY, not count.
**Per-station noise-repeater census** (`scripts/p_census.py`): each station's own 3-C mean cc
recomputed at exactly agent N's stored (i, j) pairs and compared with its own median + 8 × MAD
over 4,000 random pairs of the same block (validated — the mean over stations reproduces the
stored network cc to 7e-8).  {sm.get('n_noise_repeater', '--')} of
{n_cens} stations are
noise repeaters (tremor/quiet contribution ratio < 2 AND a top-3 share on quiet hours above
twice chance).
**Quiet subset** (rule fixed before any subset result: the PB boreholes + CAFE/YW with ratio ≥ 3
and low top-3 share, floor 12 / cap 20): **{sm.get('n_subset', '--')} stations**
({sm.get('n_pb', '--')} PB + {sm.get('n_cafe_picked', '--')} CAFE/YW);
coverage {f((sm.get('coverage_subset') or {}).get('pct_cells_ge5_60km_hyp'), 1)} % of the
sub-region's active cells with ≥ 5 subset stations within 60 km hypocentral
(51-station bed: {f((sm.get('coverage_51') or {}).get('pct_cells_ge5_60km_hyp'), 1)} %).
**Pair stage, unchanged recipe:** 51 stations {f(ps51.get('tremor_quiet_ratio'))}× →
**quiet subset {f(ps.get('tremor_quiet_ratio'))}×** → boreholes only
{f(psb.get('tremor_quiet_ratio'))}× (north POLARIS 2005: 11×; N-Olympic 2010: 4.2×).
Caveat recorded in the write-up: the quiet blocks run on fewer stations than the tremor blocks
(CAFE caching on calendar-remote days), which raises the absolute 8 × MAD bar on quiet hours.
**Five-clause screen (agent N's, unchanged) real vs whole-chain site-shift null: CANDIDATES
REAL {fr.get('n_candidate', sr.get('n_pass_abcd', 0))} vs NULL
{fn.get('n_candidate', sn.get('n_pass_abcd', 0))}**; network false-alarm rate
{f(sr.get('network_fa_per_tpl_day'), 3)} per template-day on quiet days against
{f(sr.get('network_rate_tremor_per_tpl_day'), 3)} on tremor days.
**2007-01 confirmation:** REAL {cr.get('n_CONFIRMED', 0)}/{cr.get('n_refired', 0)},
NULL {cn.get('n_CONFIRMED', 0)}/{cn.get('n_refired', 0)}.
Bed as run: {bed.get('n_stations', '--')} stations {bed.get('by_net', {})},
{bed.get('n_tremor_hours', '--')} audited tremor hours, {bed.get('n_quiet_days', '--')}
literature-quiet days.  Write-up `results/WA_2008Q_RESULTS.md`; scripts `p_census.py`,
`p_subset.py`, `p_report.py`, `p_notes.py`, `p_puget2008q_queue.sh`; the only change to an
existing file is the optional `N_STATION_SUBSET` whitelist in `n_common.admitted()` (unset →
agent N's numbers reproduce byte-for-byte).
- {now} P_PUGET2008Q_DONE.
"""
    with open(f"{STUDY}/NOTES.md", "a") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
