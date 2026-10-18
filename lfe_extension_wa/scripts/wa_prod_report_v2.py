#!/usr/bin/env python
"""Append the post-processing v2 section (user decisions 2026-10-06) to results/WA_PROD_RESULTS.md AFTER the production
queue has written its final report (the queue runs the original scripts/wa_prod_report.py at its very end; that script is
left untouched).  Idempotent: does nothing if the section is already there.  Refuses to run before WA_PROD_DONE exists
unless --force.

  python scripts/wa_prod_report_v2.py            # append once the queue is done
  bash   scripts/wa_prod_report_v2_watch.sh      # background: wait for WA_PROD_DONE, then append
"""
from __future__ import annotations
import os, sys, argparse
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P

REP = os.environ.get("WA_PROD_REPORT", f"{P.RES}/WA_PROD_RESULTS.md")
O2 = f"{P.OUT}/post_v2"
TAG = "<!-- post_v2 section -->"


def section():
    L = [TAG, "## Post-processing v2 (user decisions 2026-10-06)\n"]
    L.append("Full results: [`wa_prod/post_v2/POST_V2_RESULTS.md`](wa_prod/post_v2/POST_V2_RESULTS.md). Day classes per family: "
             "'tremor' = PNSN tremor within 30 km (before 2009-08-06: inside the WECC ETS windows); 'quiet' = the PNSN catalogue "
             "covers the day (2009-08-06 .. 2026-05-27) and shows no tremor within 30 km; 'unknown' = every other day. "
             "Headline selection: P_own_WAVE_PNET >= 0.75 for NW and Arm B on every day class; SW P_own_WAVE_PNET >= 0.5 on tremor "
             "days only; Olympic A position-group rates on P_nbr3_WAVE >= 0.75. The WAVE-model columns and the earlier class rules "
             "(tremor 0.75 / SW 0.5, NW quiet 0.85) are kept as secondary columns.\n")
    try:
        fs = pd.read_csv(f"{O2}/summary_headline_false_share_2012_2026.csv")
        L.append("- Headline 2012-2026 (core network; data / time-reversed control / scrambled control; "
                 "`post_v2/summary_headline_false_share_2012_2026.csv`): " + "; ".join(
                     f"{r.region} {r.day_class} P >= {r.cut}: {int(r.data):,} / {int(r.rev_control):,} / {int(r.scrambled_control):,}"
                     for r in fs.itertuples()) + ".")
    except Exception as e:
        L.append(f"- (headline table not readable: {e})")
    L.append("- Files: `post_v2/catalogue_v2.parquet` (every core-network candidate; `selected_v2` = headline, "
             "`selected_v2_wave_classrules` / `selected_v2_pnet_classrules` = secondary), `post_v2/rates_*`.\n")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if not a.force and not os.path.exists(f"{P.OUT}/WA_PROD_DONE"):
        print("queue not done (no WA_PROD_DONE) -> nothing appended"); return 1
    txt = open(REP).read() if os.path.exists(REP) else ""
    if TAG in txt:
        print("v2 section already present"); return 0
    with open(REP, "a") as f:
        f.write("\n" + section())
    print(f"appended v2 section to {REP}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
