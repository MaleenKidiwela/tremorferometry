---
name: check-geometry-before-density-claims
description: Before claiming a network is "dense/close", compute per-source hypocentral distances with co-located sensors merged; a nearest-family epicentral list misled the user on 2026-09-19
metadata:
  type: feedback
---

On 2026-09-19 I told the user the north had "ten close stations / ~12 sites at 38–45 km" from a list of each station's distance to
its nearest family (no operating years, no co-location check). Merlin's per-family hypocentral check found ONE new near site (SYMB);
two stations were co-located with existing boreholes, four were 65–78 km away. The user asked whether it was a lie.

**Why:** proximity is the binding term for LFE detection/discovery (P3 ablation); an overstated geometry mis-sets expectations and
the phase plan. **How to apply:** never quote station density for LFE work without (1) hypocentral distance from every source
cell/family, (2) sensors within 1 km merged into sites, (3) operating epochs, (4) counts within 45/55/60 km per source; and say
"station count" vs "close sites" explicitly. Also: no author contact for reference catalogs (user decision 2026-09-19).
