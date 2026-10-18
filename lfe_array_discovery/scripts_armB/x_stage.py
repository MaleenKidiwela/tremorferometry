#!/usr/bin/env python
"""Run one frozen `w_` engine stage ON THE WASHINGTON BED.

    python scripts/x_stage.py w_autocorr --hours tremor --out results/step5/ac --t-record 5.0

`x_common.activate()` rebinds `w_common` onto the WA bed IN PROCESS and then the named `w_` module's
own `main()` is called with the remaining arguments, unchanged.  No `w_` file is edited, and every
stage therefore runs the byte-identical engine that Step 1c ran -- only the bed differs.
"""
from __future__ import annotations

import importlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x_common as X                                                        # noqa: E402

ALLOWED = {"w_autocorr", "w_cluster", "w_grow", "w_iterate", "w_dedupe", "w_locate_gates"}


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"usage: x_stage.py <{'|'.join(sorted(ALLOWED))}> [args...]")
    mod = sys.argv[1]
    X.activate()
    m = importlib.import_module(mod)
    sys.argv = [f"{mod}.py"] + sys.argv[2:]
    m.main()


if __name__ == "__main__":
    main()
