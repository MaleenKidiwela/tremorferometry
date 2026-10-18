"""s10_multi_common -- shared paths/constants for the S10 multi-family scale-up
(rules: svi_science/results/s10_dvv/multi/PLAN.md). Pilot machinery is imported from s10_common/s10_windows/s10_measure."""
from __future__ import annotations
import os
import s10_common as S

OUT = f"{S.ROOT}/svi_science/results/s10_dvv/multi"
TB, FIG, WIN, LOGD = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/win", f"{OUT}/logs"
PILOT = S.OUT
PILOT_FAMS = (10, 87, 17)
CORE5 = ["PB.B011", "PB.B009", "PB.B010", "PB.B926", "CN.PGC"]
EXTRA = ["CN.VGZ", "CN.SNB"]
STATIONS = CORE5 + EXTRA
BOREHOLES = ["PB.B011", "PB.B009", "PB.B010", "PB.B926"]
MIN_DAYS, MIN_YEARS = 10, 10
SEED = 20261011


def ensure():
    for d in (TB, FIG, WIN, LOGD):
        os.makedirs(d, exist_ok=True)
