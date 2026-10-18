# Calibration packet review — September 14, 2026

Completed the waveform/template diagnostic following the unsuccessful site-weighting experiment. All 336 eligible calibration references were examined at the frozen reference anchor plus each saved station-template offset. No lag search, template fitting, threshold change, evaluation scan, or waveform download was performed.

[Review packet](../results/method_interleaved_packet_review_v1/index.html): 19 multistation panels, selected before waveform loading as the median network-score example in each populated source/status group. Human event/phase labels remain pending. Arrays and per-station correlations are saved alongside the plots.

## Measured waveform similarity

Median fixed-lag, three-component correlation with the frozen template:

| Physical site | Below threshold (298 references) | Recovered (31 references) | Below-threshold shifted medians, −120 / +120 s |
| --- | ---: | ---: | ---: |
| PGC/B011 | 0.184 | 0.369 | −0.002 / 0.002 |
| SNB | 0.159 | 0.306 | 0.002 / 0.009 |
| VGZ | 0.083 | 0.224 | 0.002 / 0.002 |
| B926 | 0.160 | 0.309 | 0.005 / 0.007 |

PGC and B011 instrument correlations are averaged before reporting their shared physical site. All 1,680 unshifted station windows are valid; one VGZ +120-second window is invalid and excluded explicitly. These are normalized waveform correlations, not the block-normalized detector MAD scores.

The missed-reference population retains template similarity above the shifted-window medians at every site. VGZ has the weakest median similarity. For the seven references that reach network threshold but fail simultaneous site support, VGZ has median fixed-lag correlation 0.057, versus 0.288–0.396 at the other sites. This is a descriptive pattern, not proof that VGZ alone caused every failure.

Recovered versus missed groups were selected using detector scores and need not have the same family composition; their separation is not independent validation. Nearby references can overlap, and shifted windows may contain other events. No significance, precision, or false-discovery estimate is inferred from these comparisons. Weak signals, template mismatch, source variability, and residual clock errors remain unresolved explanations.

## Review and next experiment

Plots retain native component coordinates and use a separate common three-component RMS for each template and reference packet. This preserves relative component amplitudes within each packet but removes absolute template/event amplitude differences. One panel was visually checked for layout and timing; event and phase labels were not assigned.

The next bounded sensitivity experiment should compare a shorter training-defined packet window against the existing eight-second template, with equal site weights, identical calibration eligibility, the existing site-support rule, and all shifted-site controls. This is a hypothesis motivated by the concentration of template energy in the packet, not an adopted method. Preserve the baseline and freeze the window rule before scanning. A calibration gain would still require fresh-date validation.

## Reproduction and verification

From `codex-LFE`, run `python -m codex_lfe.packet_review`. The module targets the documented saved inputs and refuses to overwrite its output directory. It records its selection design before waveform loading and saves input hashes and waveform provenance.

Artifact checks confirm exact reference-ID coverage, unique station/shift rows, all 5,040 expected windows, and exact reproduction of correlations from the 95 saved example packets. The existing software suite passes: **90 tests**. See [verification](../results/method_interleaved_packet_review_v1/verification.json).
