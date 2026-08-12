#!/usr/bin/env python
"""Resample a 100 Hz discovery template npz to 40 Hz (resample_poly 2/5, m 200->80) and write both
the 40 Hz npz and a time-reversed copy (last-axis flip) for the noise-floor reversed densify.
Replicates how data/b011_disc_p70_2010_2026_m3_40hz{,_rev}.npz were produced.
Usage: python make_40hz_npz.py <disc_prefix>   (reads <prefix>.npz -> <prefix>_40hz.npz + _40hz_rev.npz)
"""
import sys, numpy as np
from scipy.signal import resample_poly

prefix = sys.argv[1]
src = np.load(prefix + ".npz", allow_pickle=True)
keys = list(src.files)
out40, outrev = {}, {}
for k in keys:
    a = np.asarray(src[k], float)
    b = resample_poly(a, 2, 5, axis=-1)   # 100 -> 40 Hz
    out40[k] = b.astype(np.float32)
    outrev[k] = b[..., ::-1].astype(np.float32)
np.savez(prefix + "_40hz.npz", **out40)
np.savez(prefix + "_40hz_rev.npz", **outrev)
n = len(keys)
sh = np.asarray(out40[keys[0]]).shape if n else ()
print(f"[make_40hz] {n} templates {np.asarray(src[keys[0]]).shape}->{sh} "
      f"-> {prefix}_40hz.npz + {prefix}_40hz_rev.npz")
