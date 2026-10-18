#!/usr/bin/env python
"""Run one of agent K's or agent N's stage scripts UNCHANGED on the ARRAY bed.  Prefix `o_`.

  python scripts/o_run.py k_group  --detfile ... --out ...
  python scripts/o_run.py n_locate --tpl ... --det ... --out ...

`o_common` is registered in sys.modules as `n_common` BEFORE the target module is imported, so
every `n_` script (the five-clause screen, the wide-grid locator with its edge flag, the
re-stack/re-locate, the re-fire and the confirmation) reads the ARRAY bed through it and is not
edited.  `o_common.activate()` then rebinds k_common's bed facts, so the `k_` stages run unchanged
too.  No k_ or n_ file is touched.
"""
from __future__ import annotations

import importlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import o_common as O                                                        # noqa: E402

sys.modules["n_common"] = O

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].split(".")[0][:2] not in ("k_", "n_", "o_"):
        raise SystemExit("usage: o_run.py <k_/n_/o_ module> [args...]")
    name = sys.argv[1]
    if name.startswith("k_"):
        O.activate()
    sys.argv = [f"{name}.py"] + sys.argv[2:]
    importlib.import_module(name).main()
