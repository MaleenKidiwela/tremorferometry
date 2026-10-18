#!/usr/bin/env python
"""ARRAY-AWARE Puget / S-Olympic bed adapter (2008 discovery, 2007 discovery).  Prefix `y_`.

WHY THIS EXISTS.  The single change that made Washington work on the 2010 N-Olympic bed was
treating each array AS AN ARRAY instead of as independent member stations: agent O measured the
tremor/quiet pair-count ratio at 1.39 member-wise against 15.81 site-averaged
(`results/WA_2010_ARRAY_RESULTS.md` sec.0d), and the bed built that way produced the result of
`results/WA_LOCATION_REBUILD.md` sec.9.  The 2008 Puget / S-Olympic bed has NEVER been run that
way -- agent N (`WA_2008_RESULTS.md`) and agent P (`WA_2008Q_RESULTS.md`) both ran all 51
stations MEMBER-WISE, the configuration that failed on the Olympic bed, and both got a pair ratio
of 0.86-0.89 and zero candidates.  The 2007 window (CAFE at its fullest) has never been run for
discovery at all.

WHAT THIS MODULE IS.  NOT a rewrite.  `o_common.py` and `o_beam.py` are the agent-O array
machinery and are used UNCHANGED; this module rebinds o_common's BED FACTS in process, exactly as
o_common rebinds k_common's.  The bed facts (admitted stations, window, audited tremor hours,
literature-quiet days, active cells) come from `n_common.py`, the agent-N 2008/2007 bed, also
unchanged.  The result is: agent O's site clustering, beam builder, three site treatments and
site-shift null, running on agent N's 2008 and 2007 beds, with agent K's frozen recipe underneath
and no line of any of them edited.

  members       = n_common.keepers() for this bed (probe-admitted, in-bed, health-kept)
  sites         = single-link clustering at 1.5 km partitioned by network (o_common._build_sites)
  hours         = n_common.candidate_hours() -- the waveform-derived audit of n_hours.py
  quiet days    = n_common.quiet_days()      -- literature-quiet calendar days
  active cells  = n_common.active_cells()    -- the Puget/S-Olympic box, NOT the N-Olympic one

Environment
  Y_BED=2008|2007     which bed (default 2008)
  O_NULL_SHIFT=1      whole-chain site-shift null (o_common.activate, seed O_NULL_SEED)
  O_OUT               output directory for this chain
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
sys.path.insert(0, HERE)

BED = os.environ.get("Y_BED", "2008")
if BED not in ("2008", "2007"):
    raise SystemExit(f"Y_BED must be 2008 or 2007, not {BED!r}")

OUTROOT = os.environ.get("Y_OUTROOT", f"{RES}/y_wa_arrays")
BEAMROOT = f"{STUDY}/beamcache_y"

# n_common reads N_BED at import; o_common's import-time guard only accepts 2010/2011, so it is
# imported on its own default and every bed fact is rebound below.  Nothing in either file moves.
# ORDER MATTERS: both modules put `lfe_picker/scripts` on sys.path[0], and that directory holds a
# DIFFERENT o_common.py.  o_common is therefore imported first, from this study's scripts/.
os.environ["N_BED"] = BED
os.environ["O_BED"] = "2010"
import o_common as O                                                        # noqa: E402
import n_common as N                                                        # noqa: E402

if os.path.dirname(os.path.abspath(O.__file__)) != HERE:
    raise SystemExit(f"wrong o_common imported: {O.__file__}")
if os.path.dirname(os.path.abspath(N.__file__)) != HERE:
    raise SystemExit(f"wrong n_common imported: {N.__file__}")

_S = {}


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%H:%M:%S')}]", *a, flush=True)


# =============================================================== the member pool ================
def _members():
    """The bed's member stations: n_common's probe-admitted, in-bed, health-kept set.  Returned
    in the frame o_common expects (net, sta, netsta)."""
    if "M" not in _S:
        K = N.keepers()[["net", "sta", "netsta"]].copy().reset_index(drop=True)
        _S["M"] = K
    return _S["M"]


def _member_az12():
    """{member: (az1, az2)} for members cached as (?HZ, ?H1, ?H2) -- the PB boreholes.  Computed
    over the MEMBER pool (n_common.az12_wa keys on the same pool, so it is reused verbatim)."""
    if "AZ" not in _S:
        _S["AZ"] = N.az12_wa()
    return _S["AZ"]


# =============================================================== the patch ======================
def _patch():
    """Rebind o_common onto this bed.  Idempotent; must run before any o_common call."""
    if _S.get("PATCHED"):
        return O
    O.BED = BED
    O.YEAR, O.DAY0, O.DAY1 = N.YEAR, N.DAY0, N.DAY1
    O.VETO_YEARS = N.VETO_YEARS
    O.OUTROOT = OUTROOT
    O.SITES_CSV = f"{OUTROOT}/sites_{BED}.csv"
    O.STATICS = f"{OUTROOT}/statics_{BED}.json"
    O.BEAMLOG = OUTROOT + "/beam_hours_{bed}.csv"
    O.BEAM = BEAMROOT
    O.OUT = os.environ.get("O_OUT", f"{OUTROOT}/{BED}")
    # bed facts
    O.members_2010 = _members
    O.member_az12 = _member_az12
    O.candidate_hours = N.candidate_hours
    O.quiet_days = N.quiet_days
    O.quiet_hours = N.quiet_hours
    O.active_cells = N.active_cells
    os.makedirs(OUTROOT, exist_ok=True)
    _S["PATCHED"] = True
    return O


_patch()

# ---- the o_common API, re-exported so `import y_common as Y` reads naturally ------------------
sites_all = O.sites_all
sites = O.sites
site_members = O.site_members
arrays = O.arrays
stations = O.stations
station_table = O.station_table
member_table = O.member_table
member_cache_index = O.member_cache_index
cache_index = O.cache_index
beam_path = O.beam_path
window_days = O.window_days
bed_days = O.bed_days
beam_days = O.beam_days
candidate_hours = O.candidate_hours
tremor_dates = O.tremor_dates
quiet_days = O.quiet_days
quiet_hours = O.quiet_hours
health_veto = O.health_veto
station_corrections = O.station_corrections
near_stations = O.near_stations
active_cells = O.active_cells
null_offsets = O.null_offsets
activate = O.activate


def bed_summary():
    """o_common.bed_summary plus what is specific to this bed (the Puget sub-region, the audited
    hours, the literature-quiet days, and the member stations that end up in NO array)."""
    b = O.bed_summary()
    S = sites()
    b["subregion"] = N.SUBREGION
    b["tranche"] = N.TRANCHE
    b["n_member_stations"] = int(len(_members()))
    b["n_quiet_hours"] = int(len(quiet_hours()))
    b["members_by_net"] = {k: int(v) for k, v in _members().groupby("net").size().items()}
    b["arrays_by_net"] = {k: int(v) for k, v in
                          S[S.kind == "ARRAY"].groupby("net").size().items()}
    b["n_active_cells"] = int(len(active_cells()))
    return b
