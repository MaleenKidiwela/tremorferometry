---
name: family-veto-mechanisms
description: "User QC vetoes (229 families, 4 date windows) + era-split whitelist (6 stations) — standing mechanisms, 2026-08-31"
metadata: 
  node_type: memory
  type: project
  originSessionId: eeb752df-b1e8-4172-9ab5-f9abf3317b7f
  modified: 2026-08-31T14:33:29.402Z
---

Three standing user-QC mechanisms (2026-08-31), all honored by assemble_res_catalog_daily.py AND
build_borehole_dvv_map.py — full ledger in notes/2026-08-31_family_qc_vetoes.md:
1. `data/family_user_blocklist.csv` — 229 vetoed families (whole stations KMR/SEAS2/LEID included).
2. `data/station_date_blocklist.csv` — date-range vetoes (PETY middle era, B031 post-2014,
   SSO2 post-2025-12, BTB post-2018).
3. `data/era_split_whitelist.txt` — ONLY WISH/SCOT/KHBB/TRIN/LCCR/YACT/PETY get per-era offsets
   (user ruling: visible cross-era steps, plus PETY whose two surviving sensors straddle a vetoed
   gap); other era_table stations are single-era now.

**How to apply:** resolve HTML-map family labels with `scripts/family_from_map_index.py`
(label "N/M" is 1-BASED, /M is a checksum; c-numbers are NOT unique per station — always store full
ids). After ANY veto change: re-run map builder + the full inversion chain (assembler → solver →
transmission_T → transmission_deconv) + regenerate movies/explorer, or products are inconsistent.
PETY: user RESOLVED 2026-08-31 — single c across the vetoed gap (sensors A+C share one offset);
risk knowingly accepted; a step-at-the-gap at PETY would implicate this choice first.
