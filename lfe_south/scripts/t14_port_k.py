#!/usr/bin/env python
"""STAGE R2 prep -- port the REPRODUCED Bostock discovery recipe to the southern permanent bed.

The k_*.py scripts of lfe_array_discovery are the recipe as reproduced in the north (58/96 kept,
22/30 sealed).  They are COPIED here unchanged except for two fidelity edits, both forced by the
southern network and both recorded in NOTES.md:
  (1) the cache root and the station list point at lfe_south/cache_duc (the Ducellier window);
  (2) `load_day` accepts 1-C and 2-C station-days.  ~65 % of the southern permanent bed is
      vertical-only; rejecting those stations would run a THINNED network, which is not the
      published method.  Each available component contributes to the per-component cc sum exactly
      as a 3-C station's components do.
Nothing else is touched: window, band, 8 x MAD, 4 s decluster, Bungum-Husebye realignment,
transitive grouping, no merge, no gates.
"""
import os,sys,re,shutil
HERE=os.path.dirname(os.path.abspath(__file__))
SRC=os.path.abspath(f"{HERE}/../../lfe_array_discovery/scripts")
for f in ["k_common.py","k_auto.py","k_group.py","k_grow.py","k_locate.py"]:
    s=open(f"{SRC}/{f}").read()
    s=s.replace("import k_common as K","import ks_common as K")
    s=re.sub(r"\bimport k_common\b","import ks_common",s)
    if f=="k_common.py":
        s=s.replace('''    z = np.load(idx[key])
    x = z["x"]
    order = [str(o) for o in z["order"]]
    if x.shape[0] != 3 or [o[-1] for o in order] != ["Z", "N", "E"]:
        return None
    return float(z["t0"]), np.ascontiguousarray(x, dtype=np.float32)''',
'''    z = np.load(idx[key])
    x = z["x"]
    # SOUTHERN FIDELITY EDIT (agent T, 2026-09-23): accept 1-C and 2-C station-days.  ~65 % of
    # Ducellier's permanent bed is vertical-only (NC EHZ/SHZ); his detector correlates PER CHANNEL
    # and averages, so those stations are full members of the network.  Requiring 3-C here would
    # run a thinned network, which is not the published method.
    if x.shape[0] < 1:
        return None
    return float(z["t0"]), np.ascontiguousarray(x, dtype=np.float32)''')
        s=s.replace('CACHES = [f"{STUDY}/cache20hz", f"{REP}/cache20hz"]',
                    'CACHES = [f"{STUDY}/cache_duc", f"{STUDY}/cache20hz"]')
        s=s.replace('STUDY = os.path.dirname(HERE)',
                    'STUDY = os.path.dirname(HERE)   # = lfe_south when run from lfe_south/scripts')
    open(f"{HERE}/ks_{f[2:]}","w").write(s)
    print("ported ks_"+f[2:])
