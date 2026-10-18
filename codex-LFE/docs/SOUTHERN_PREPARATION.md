# Southern catalog preparation — September 14, 2026

The independent southern pilot now has a training bank for one named published family at four target-network sites: BK.RVIT, PB.B935, PB.B045 and BK.PETY. **No southern detections or evaluation claim yet.**

The published reference release was independently normalized into 61,441 records and 27 named families. Its documented clock is template start, one second before estimated S arrival at the original first station; it is not an origin-time catalog. Original family IDs and metadata coordinates are retained. These published records are separate from our detections.

The initial January 1–5 training selection chose the most frequent named family, `shelly2025_23083010676s`, using reference counts alone. Its 24 events did not produce a passing template at any of the six closest qualifying roster stations. The expanded experiment retained that family, adding January training dates through January 31 while excluding the reserved January 6–10 dates. It contains 95 references on 17 dates. Channel epochs were checked for the whole expanded interval.

| Station | Valid training windows | Split-stack correlation | Training gate |
| --- | ---: | ---: | --- |
| BK.RVIT | 64 | 0.395 | pass |
| PB.B935 | 95 | 0.445 | pass |
| PB.B045 | 95 | 0.423 | pass |
| BK.HAYF | 73 | 0.067 | fail |
| BK.PETY | 73 | 0.312 | pass |
| NC.KHBB | 73 | -0.016 | fail |

Training uses the existing >=20-window and 0.25 split-stack-correlation requirements. Window counts differ because only fully valid, padded waveform packets are accepted. No station support requirement or coherence threshold was lowered. Alternating-event training consistency is not independent event validation, especially within bursts.

Short waveform requests, rather than full continuous days, were acquired. Total saved miniSEED across the initial attempts and expansion is 6,818,816 bytes, including repeated training content. NCEDC first rejected nanosecond/timezone timestamp syntax; a microsecond UTC retry succeeded. Its initial overlapping requests returned truncated coverage, so intervals were merged explicitly before the final retry. Failures and retries are retained. No calibration or test waveform was downloaded.

The original January 6–7 calibration window contains zero published references for this family. It cannot measure known-reference recovery. Next: predeclare a reference-supported calibration block, retain a separate evaluation block, acquire continuous data for the four usable sites, and calibrate the detector with shifted-site controls. This initial family is a regional workflow increment, not southern family completeness.

[Expanded training report](../results/southern_training_expanded_v1/report.json) · [Frozen acquisition design](../results/southern_training_expanded_v1/design_before_acquisition.json) · [Reference provenance](../results/southern_reference_audit_v1/manifest.json).

Software verification: **94 tests passed**, including rejection of duplicate southern source rows, inconsistent clocks and unmapped family IDs, plus subsecond UTC preservation. The saved 61,441-row reference artifact exactly reproduces after the missing-family guard fix.
