#!/usr/bin/env python
"""E4 step 0 — THE component-order resolver, shared by every E4/E3 script that reads 3-C data.

WHY THIS FILE EXISTS (bug found by the Merlin E3 re-audit, 2026-09-11).
The SeisBench LFE stacks (`seisbench_lfestacks_cascadia/waveforms.hdf5`) declare
`data_format/component_order = "Z12"` with **1 = N and 2 = E** (verified independently here by the
P-wave polarization test in `e4_polarization_check.py`: the horizontal P principal axis agrees with
`path_backazimuth_deg` to a median ~9 deg under 1=N and ~40 deg under 1=E).  Three E3 scripts built
their data matrix as `[Z] + sorted(non-Z)`, which for CN stations gives BHZ, BHE, BHN = **Z, E, N**
-- the two horizontals swapped relative to the template.  Every 3-C correlation and every injected
waveform in E3 steps 3, 3b and 4 was therefore computed with the horizontals 90 deg rotated.

Fix: resolve the order EXPLICITLY by the last letter of the channel code -- Z first, then N or 1,
then E or 2 -- never by `sorted()`.
"""
from __future__ import annotations

# codes accepted for each slot, in priority order
_SLOTS = (("Z",), ("N", "1"), ("E", "2"))


class ChannelOrderError(RuntimeError):
    pass


def resolve_order(chans, station="?", strict=True, legacy=False):
    """Return the channel codes in the order (Z, N|1, E|2).

    chans   iterable of channel code strings (e.g. {'BHZ','BHN','BHE'} or {'EHZ','EH1','EH2'})
    strict  raise if a slot is missing or ambiguous (both N and 1 present -> STOP, as prescribed)
    legacy  reproduce the OLD (buggy) E3 behaviour, [Z] + sorted(non-Z); kept so the pre-fix
            numbers stay reproducible behind a flag.
    """
    chans = sorted(set(chans))
    if legacy:
        return [c for c in chans if c.endswith("Z")] + [c for c in chans if not c.endswith("Z")]
    out = []
    for slot in _SLOTS:
        hit = [c for c in chans if c[-1:].upper() in slot]
        if len(hit) > 1:
            raise ChannelOrderError(
                f"{station}: slot {slot} is ambiguous -- {hit}. Two different horizontal naming "
                f"conventions on the same station-day; STOP and resolve by hand (E4 step 0 rule).")
        if not hit:
            if strict:
                raise ChannelOrderError(f"{station}: no channel for slot {slot} in {chans}")
            out.append(None)
            continue
        out.append(hit[0])
    return out


def describe(order, station="?"):
    lab = "/".join("-" if c is None else c[-1] for c in order)
    return f"[chanorder] {station}: {order} -> Z/{lab.split('/', 1)[1] if '/' in lab else '?'}"
