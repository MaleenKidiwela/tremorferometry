#!/usr/bin/env python
"""Run one of agent K's / N's / O's / U's stage scripts UNCHANGED on the ARRAY-AWARE Puget bed.
Prefix `y_`.  Same device as `o_run.py`, one bed further.

  Y_BED=2008 python scripts/y_run.py o_auto   --mode siteavg --hours tremor --out ...
  Y_BED=2008 python scripts/y_run.py k_group  --detfile ... --out ...
  Y_BED=2008 python scripts/y_run.py n_locate --tpl ... --det ... --out ...
  Y_BED=2008 python scripts/y_run.py uu_stacks --bed wa --det ... --tpl ... --out ...

`y_common` is imported FIRST, which (a) imports this study's `o_common` before anything can put
`lfe_picker/scripts` ahead of it on sys.path and (b) rebinds o_common's bed facts onto the 2008 /
2007 Puget bed.  `o_common` is then registered as `n_common`, so every `n_` stage (the five-clause
screen, the wide-grid locator, the re-stack/re-locate, the re-fire, the confirmation) and
`uu_common.bed("wa")` read the ARRAY bed.  `y_common.activate()` rebinds k_common, so the `k_`
stages run unchanged too.  NO k_, n_, o_ or uu_ file is edited.
"""
from __future__ import annotations

import importlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import y_common as Y                                                        # noqa: E402

sys.modules["n_common"] = Y.O

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].split(".")[0].split("_")[0] not in (
            "k", "n", "o", "uu", "y"):
        raise SystemExit("usage: y_run.py <k_/n_/o_/uu_/y_ module> [args...]")
    name = sys.argv[1]
    if name.startswith("k_"):
        Y.activate()
    sys.argv = [f"{name}.py"] + sys.argv[2:]
    importlib.import_module(name).main()
