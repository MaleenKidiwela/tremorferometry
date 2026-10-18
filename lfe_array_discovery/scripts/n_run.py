#!/usr/bin/env python
"""Run one of agent K's stage scripts UNCHANGED on the Puget/S-Olympic bed.  Prefix `n_`.

  python scripts/n_run.py k_group --detfile ... --out ...

`n_common.activate()` rebinds k_common's bed facts (station table, borehole azimuths,
net-qualified cache index, optional day-shift null) in this process, then the named k_ module's
own `main()` runs with the remaining arguments.  No k_ file is edited.  (Identical in design to
`l_run.py`, which drove the 2010 run.)
"""
from __future__ import annotations

import importlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import n_common as N                                                        # noqa: E402

if __name__ == "__main__":
    if len(sys.argv) < 2 or not sys.argv[1].startswith("k_"):
        raise SystemExit("usage: n_run.py <k_module> [args...]")
    N.activate()
    name = sys.argv[1]
    sys.argv = [f"{name}.py"] + sys.argv[2:]
    importlib.import_module(name).main()
