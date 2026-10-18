---
name: metadata-overlay-qc
description: "Always run the metadata-overlay QC per station; it catches bad-SENSOR eras per-era referencing can't fix"
metadata: 
  node_type: memory
  type: project
  originSessionId: 5cecf374-37e4-4bcc-aea9-79b2b6cdc569
---

For every station's dv/v, run `scripts/plot_dvv_metadata.py --network <NET> --station <STA>`
(fetches FDSN channel epochs, draws sensor/rate changes red, response epochs grey on the dv/v).

**Why:** instrument changes masquerade as velocity. Two failure modes:
1. **Steps at a change** (HDW EHZ→HHZ 2022-09 → the apparent 2023 "rise" was the swap). Fixed by
   **per-era referencing** at the change dates.
2. **A bad SENSOR that is internally unstable for a whole epoch** — per-era referencing centers
   the epoch's mean but CANNOT fix within-epoch garbage. **CPW** is the case study: clean in its
   L-4C eras (1997–2011, 2019–2022) but ±0.1% wild swings during the **ROCK1 sensor era
   (2011–2019, L4-1095 then L4-2310)** — the swings line up exactly with the sensor-swap lines.
   Verdict: trust the L-4C eras, discard 2011–2019.

**How to apply:** get per-era boundaries from the FDSN epochs (NOT guesses); run per-era dv/v +
the metadata overlay; then read the overlay to decide which eras are usable. Same-TYPE sensor
swaps (different serials) can still ruin a record. Also check station end_date up front — many
UW short-period sites are decommissioned; PB borehole (B-codes) are the active replacements.
See [[gpu-family-discovery]], [[pnsn-catalog-cut]], [[proactive-note-taking]].
