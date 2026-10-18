# Four-second packet calibration — September 14, 2026

**Retain the eight-second baseline.** Both arms select 9 MAD under the existing calibration gate. Four-second templates recover 24/336 references with 26 candidates; the baseline recovers 31/336 with 34 candidates. Maximum exposure-adjusted control ratios are 3.85% and 2.94%, respectively. At 8 MAD neither arm passes the 5% criterion. No evaluation data selected the winner.

The experiment crops all 55 frozen training templates to training peak −0.5 through +3.5 seconds. The crop starts 1.5 seconds later than the baseline window, and the station offset is increased by exactly 1.5 seconds to preserve the reference anchor. Source groups, stations, equal site weighting, four-site support, calibration dates, control realizations and matching rules are retained. The four-date scan took 83.8 seconds on four CPU workers.

All 336 reference IDs remain eligible. Shorter windows add 731.5 observed family-seconds around gaps (0.0244% of baseline exposure); exposure ratios use each arm's actual masks. These are therefore identical-reference comparisons, not exactly identical-mask comparisons.

An artifact audit found that pandas' default approximate comparison of epoch timestamps incorrectly accepted these small interval differences. The original decision artifact is preserved because the next block already recorded its hash; `verification.json` explicitly corrects its `identical_waveform_intervals` field. The experiment script now uses exact equality. Counts, exposure-adjusted ratios, reference eligibility, and rejection of the shorter method are unaffected.

[Comparison](../results/method_interleaved_short_packet_v1/comparison.csv), [original decision](../results/method_interleaved_short_packet_v1/decision.json), [artifact verification and correction](../results/method_interleaved_short_packet_v1/verification.json).

The baseline is now frozen at 9 MAD for the September 6–10, 2010 northern catalog block. No further window or threshold search is performed on those dates. Existing software checks passed all 90 tests; template arrays were independently checked for exact cropping and timing offsets.
