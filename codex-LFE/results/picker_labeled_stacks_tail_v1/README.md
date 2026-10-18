# Published-stack picker timing diagnostic

This missing diagnostic uses local jms weights (documented training excludes Cascadia), existing published nominal P/S labels and no waveform downloads. Earlier continuous and injection experiments remain available and were not rerun.

The initial fixed-grid scan excluded most S labels because short records did not supply the next full 60-second window. The tail check evaluates only previously uncovered labels in a final unpadded window placed by record length; already covered results are reused unchanged. Production inference code has not been changed.

| Phase | Tolerance | Covered labels | Matches | Recovery | Family-bootstrap 95% interval |
|---|---:|---:|---:|---:|---:|
| P | ±0.5 s | 1067 | 569 | 53.3% | 49.9–56.6% |
| P | ±1 s | 1067 | 697 | 65.3% | 62.5–68.2% |
| S | ±0.5 s | 1453 | 663 | 45.6% | 42.9–48.6% |
| S | ±1 s | 1453 | 857 | 59.0% | 56.1–62.1% |

All table rows use score 0.3 as a diagnostic, not a selected production threshold. Other predeclared thresholds and tolerances are in the CSV.

Limits: stacks average many events; all P labels are aligned at sample 800; source families overlap the original row-level splits; label uncertainties are absent; unmatched peaks are not proven false LFEs. The original PO station identifiers are preserved rather than silently relabeled as current roster stations. CPU/GPU verification on one stack found identical validity masks and maximum probability difference 0.000858; this is not exhaustive numerical validation.

The timing figure and intervals describe stack behavior only. Natural recordings with independent phase/event review, continuous usable exposure, instrument transfer and independent-site association remain necessary before full-history production.
