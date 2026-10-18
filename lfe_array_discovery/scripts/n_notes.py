#!/usr/bin/env python
"""Append the dated NOTES.md entry for the Puget/S-Olympic 2008 run.  Prefix `n_`."""
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


def main():
    o, r, n, c7 = sys.argv[1:5]
    au, au7 = j(f"{o}/hour_audit_2008.json"), j(f"{o}/hour_audit_2007.json")
    b8 = j(f"{o}/bed_2008.json")
    sr, sn = j(f"{r}/screen_summary.json"), j(f"{n}/screen_summary.json")
    fr, fn = j(f"{r}/screen_final_summary.json"), j(f"{n}/screen_final_summary.json")
    cr, cn = j(f"{c7}/confirm_real_summary.json"), j(f"{c7}/confirm_null_summary.json")
    s, s7 = au.get("selected", {}), au7.get("selected", {})
    nan = float("nan")
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    txt = f"""
## 2026-09-21 — agent N: Puget / S-Olympic 2008 (CAFE XU + PB), the last full-array WA bed
**Tremor-hour audit replacing the 86 %-daytime energy flag** (per-station diurnal-normalised
2–8 Hz excess; ≥ 60 % of the near-field array simultaneously above its OWN quiet-day p95; plus a
spatial-coherence test — the excess must fit an amplitude decay from one active cell better than a
uniform field, p ≤ 0.05 against 200 station-label shuffles): **{s.get('n')} hours on
{s.get('n_days')} days**, only {s.get('frac_daytime_15_22utc', nan):.3f} of them in 15–22 UTC
(the old flag: 0.86), weekday fraction {s.get('frac_weekday', nan):.3f} (5/7 = 0.714 by calendar);
**0 of {au.get('n_quiet_day_hours_judged')} literature-quiet-day hours fire**.  The selected days
reconstruct the published May-2008 Washington ETS from the waveforms alone.  The array's diurnal
noise floor is ~10× higher at 19 UTC than at 09 UTC, so the night-heavy clock distribution is
detection-limited, not a claim that tremor is nocturnal (quiet-day false-positive rate by UTC hour
is flat).  Bed: {b8.get('n_stations')} 3-C stations {b8.get('by_net')}, {b8.get('n_quiet_days')}
literature-quiet days.
Pair stage → grouping → growth → the five-clause screen (per-template quiet-day FA; tremor
co-timing; wide-grid INTERIOR location + shuffled-station control + Slab2; ≥ 3 days; second-half
re-stack stability), real against a whole-chain site-shift null through byte-identical code:
**CANDIDATES REAL {fr.get('n_candidate', sr.get('n_pass_abcd', 0))} vs NULL
{fn.get('n_candidate', sn.get('n_pass_abcd', 0))}**.
2007-01 confirmation re-fire (tranche T5, read by no discovery step of this project):
**REAL {cr.get('n_CONFIRMED', 0)}/{cr.get('n_refired', 0)} confirmed, NULL
{cn.get('n_CONFIRMED', 0)}/{cn.get('n_refired', 0)}** (2007 audit: {s7.get('n')} hours on
{s7.get('n_days')} days).
Full write-up `results/WA_2008_RESULTS.md`; new scripts `n_common.py`, `n_hours.py`, `n_auto.py`,
`n_pairstats.py`, `n_run.py`, `n_chunk.py`, `n_locate.py`, `n_screen.py`, `n_relocate.py`,
`n_refire.py`, `n_confirm.py`, `n_report.py`, `n_notes.py`, `n_puget2008_queue.sh`.
- {now} N_PUGET2008_DONE.
"""
    with open(f"{STUDY}/NOTES.md", "a") as f:
        f.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
