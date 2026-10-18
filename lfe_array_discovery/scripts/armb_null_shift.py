#!/usr/bin/env python
"""Run an Arm-B script (o_chestler.py / uu_stacks.py) under the arm-A SITE-DAY-SHIFT scrambled
control, WITHOUT editing the Arm-B code.  Prefix `armb_`.

Why a wrapper: o_common.activate() (O_NULL_SHIFT=1) shifts only k_common.load_day -- the SITE
series (array beams + singles) that arm A and Arm B's `seeds` stage read.  Arm B's `grow` stage and
the long re-stack read MEMBER stations through o_beam.load_member (o_chestler.BandCache), which
activate() never touches.  Here o_beam.load_member is patched so every member takes its SITE's
offset (sites.csv membership), the same within-partition rule (tremor days shift within the
tremor-day list, quiet days within the quiet-day list) and the same offsets as arm A
(results/o_oly2010/2010_null/null_offsets.json; asserted equal to o_common.null_offsets(seed 11)).
Members of one site therefore shift together; only inter-site coherence is destroyed.
The member's intra-day t0 offset is preserved (t0 + ts(date) - ts(source_date)).

  env O_BED=2010 O_OUT=<dir> O_NULL_SHIFT=1 O_NULL_SEED=11 \
      python scripts/armb_null_shift.py scripts/o_chestler.py --stage grow --out <dir> ...
"""
from __future__ import annotations

import json
import os
import runpy
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import o_common as O                                                        # noqa: E402
import o_beam as OB                                                         # noqa: E402

assert os.environ.get("O_NULL_SHIFT", "") == "1", "armb_null_shift.py requires O_NULL_SHIFT=1"
OFF_JSON = f"{O.OUTROOT}/2010_null/null_offsets.json"
off = json.load(open(OFF_JSON))
calc = O.null_offsets(O.stations())
assert off == calc, f"offsets differ from arm A's file: {OFF_JSON}"

site_of = {}
for r in O.sites().itertuples():
    for m in str(r.members).split(","):
        site_of[m] = r.site
td, qd = O.tremor_dates(), O.quiet_days()
part = {d: ("t", i) for i, d in enumerate(td)}
part.update({d: ("q", i) for i, d in enumerate(qd)})
lists = {"t": td, "q": qd}
_base = OB.load_member


def _shifted_member(sta, date):
    date = str(pd.Timestamp(date).strftime("%Y-%m-%d"))
    site = site_of.get(sta)
    if date not in part or site not in off:
        return _base(sta, date)
    kind, i = part[date]
    L = lists[kind]
    if len(L) < 2:
        return _base(sta, date)
    src = L[(i + off[site]) % len(L)]
    v = _base(sta, src)
    if v is None:
        return None
    dt = pd.Timestamp(date).timestamp() - pd.Timestamp(src).timestamp()
    return float(v[0]) + dt, v[1]


OB.load_member = _shifted_member
O.log(f"[armb-null] member site-shift active: {len(site_of)} members -> {len(off)} sites; "
      f"offsets {OFF_JSON}")

script = sys.argv[1]
sys.argv = sys.argv[1:]
runpy.run_path(script, run_name="__main__")
