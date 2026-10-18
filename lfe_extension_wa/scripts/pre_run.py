#!/usr/bin/env python
"""WA pre-production runner: points the FROZEN modules (wa_mf / pre_mf_cap / wa_hf / wa_cache) at a results/wa_pre/
sub-directory without changing them.  W.RES = --res (must hold templates/ -> symlink or a pre-production template
set), W.SETUP = --setup, download roots cut out of TRACE_ROOTS (pre_common.patch_wa), optional extra cache dirs
searched first.

  python scripts/pre_run.py --res results/wa_pre/cap --setup results/wa_gonogo/setup mf -- --tag cap30 ...
  modules: mf (wa_mf) | mfcap (pre_mf_cap) | hf (wa_hf)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
args = sys.argv[1:]
res = os.path.abspath(args[args.index("--res") + 1])
setup = os.path.abspath(args[args.index("--setup") + 1])
cdirs = args[args.index("--cache-dirs") + 1].split(",") if "--cache-dirs" in args else []
cmd = [x for x in args[:args.index("--")] if x in ("mf", "mfcap", "hf")][0]
rest = args[args.index("--") + 1:]
os.environ["WA_SETUP"] = setup
sys.path.insert(0, HERE)
import pre_common as P
W = P.W
W.RES = res
W.SETUP = setup
P.patch_wa([os.path.abspath(c) for c in cdirs])
if "--inject-window" in args:                     # NAME=csv(date,kind): a custom day list as a window label
    import pandas as pd
    nm, csv = args[args.index("--inject-window") + 1].split("=", 1)
    X = pd.read_csv(csv)
    W.WINDOWS[nm] = [(d, d, k) for d, k in zip(X.date, X.kind)]
assert not any("traces_dl" in r or "traces_sealed" in r for r in W.TRACE_ROOTS), W.TRACE_ROOTS
mod = __import__({"mf": "wa_mf", "mfcap": "pre_mf_cap", "hf": "wa_hf"}[cmd])
sys.argv = [f"{cmd}.py"] + rest
mod.main()
