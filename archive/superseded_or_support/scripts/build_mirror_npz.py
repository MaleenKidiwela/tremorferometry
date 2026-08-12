#!/usr/bin/env python
"""Build a 'mirror' stacks npz for the mirror-window noise-reference test (Merlin step 2).
Copy the forward Z stacks but REPLACE the coda window [2,4]s samples with the TIME-FLIPPED pre-arrival
window [-3,-1]s (sample nearest the arrival maps to the anchor side). Running the IDENTICAL
dvv_roll30cal --window 2 4 --origin-anchor on this then measures the stretch of the PRE-ARRIVAL NOISE
with the exact same anchor / SVD-Wiener / rolling logic as the real coda dv/v -> a fair free noise reference.
Usage: build_mirror_npz.py <TAG>   e.g. B926p90f40"""
import sys, numpy as np

TAG = sys.argv[1]
d = np.load(f"data/long_window_daily_{TAG}_Z.npz", allow_pickle=True)
t = d["t"]; S = d["stacks"].copy()
coda = np.where((t >= 2) & (t <= 4))[0]
mir = np.where((t >= -3) & (t <= -1))[0]
n = min(len(coda), len(mir))
coda, mir = coda[:n], mir[:n]
# place the pre-arrival noise, time-flipped, into the coda slot (t=-1 -> t=2 near anchor; t=-3 -> t=4)
S[:, coda] = d["stacks"][:, mir][:, ::-1]
out = {"t": t, "fs": d["fs"], "stacks": S, "patches": d["patches"], "dates": d["dates"], "n_det": d["n_det"]}
np.savez(f"data/long_window_daily_{TAG}_MIRROR.npz", **out)
print(f"{TAG}: mirror npz written ({n} coda samples <- flipped pre-arrival) -> data/long_window_daily_{TAG}_MIRROR.npz")
