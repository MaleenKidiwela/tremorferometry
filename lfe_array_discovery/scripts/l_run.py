#!/usr/bin/env python
"""Run one of agent K's stage scripts UNCHANGED on the Washington bed.  Prefix `l_`.

  python scripts/l_run.py k_grow --nodes ... --summary ... --out ...

`l_common.activate()` rebinds k_common's bed facts (station table, borehole azimuths,
net-qualified cache index, optional day-shift null) in this process, then the named k_ module's
own `main()` runs with the remaining arguments.  No k_ file is edited.
"""
from __future__ import annotations

import importlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import l_common as L                                                        # noqa: E402

if __name__ == "__main__":
    if len(sys.argv) < 2 or not sys.argv[1].startswith("k_"):
        raise SystemExit("usage: l_run.py <k_module> [args...]")
    L.activate()
    name = sys.argv[1]
    sys.argv = [f"{name}.py"] + sys.argv[2:]
    importlib.import_module(name).main()
