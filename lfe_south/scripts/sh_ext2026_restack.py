#!/usr/bin/env python
"""SH_EXT2026 wrapper: sh_v3_restack.py (unchanged) restricted to a station subset.
  python scripts/sh_ext2026_restack.py --stations NC.KMR -- --mode det --det X --picks P --days-file D --thr 9 --out O.npz
Everything after '--' is passed to sh_v3_restack.main() verbatim; only its module-level S1_STATIONS is replaced."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
i = sys.argv.index("--")
own, rest = sys.argv[1:i], sys.argv[i + 1:]
stas = own[own.index("--stations") + 1].split(",") if "--stations" in own else None
import sh_v3_restack as R
if stas: R.S1_STATIONS = stas
sys.argv = [R.__file__] + rest
R.main()
