# NW dominance cap on two more episodes (V2008 test half, A2011): result

Rules were written before any run: `PREREGISTERED.md` (2026-10-03 01:43Z). Runs 01:45–02:50Z.

## In short
- **Verdict: NW stays uncapped.** The 0.3 cap fails the preregistered rule on both episodes, and only on the recall
  condition:
  - V2008: recall falls 20.5 % (allowed: 20 %).
  - A2011: recall falls 27.8 %.
- The other two conditions pass on both episodes, and by more than they did on V2010:
  - Own precision rises +0.030 (V2008) and +0.083 (A2011).
  - Family-windows above both scrambled controls rise 86 → 114 (V2008) and 37 → 86 (A2011).
- So the pattern from V2010 repeats on both episodes: the cap gives purer, better-separated NW detections, but it
  loses about a fifth to a quarter of the matched repeats. Recall is what matters for following families through time
  and for dv/v stacks.
- One caution about the family-window gain: under the cap, the false-pass level (control 2 vs control 1) also rises,
  8 → 10 (V2008) and 9 → 18 (A2011). So part of the A2011 gain is the bar getting easier to clear. Net of that:
  - V2008: 78 → 104.
  - A2011: 28 → 68.
  - V2010: 77 → 107.

## What was run
- Network: the go/no-go 21-core network. 20 stations / 13 sites can contribute; PB.B926 has no template.
- Templates: the build-half T1 templates (all ≥ 30-event templates).
- Matched filter: `scripts/pre_mf_cap.py` in three weight modes, each with 2 scrambled controls (stations shifted
  in time by a different random amount per site, same seeds as before):
  - uncapped,
  - capped at 0.3,
  - dominant site left out.
- Scoring: HF column + `gonogo_score.py` unchanged; ≥ 8.75 × MAD, one-to-one ±2 s.
- Time-reversed-template control: the same three weight modes with every template channel reversed (reference only).
- Episodes:
  - **V2008**: test half scored, HF cuts from the V2008 build half.
  - **A2011**: all of it scored. It was scored together with V2010, exactly as the go/no-go, so its HF cuts come
    from the V2010 build half.
- Not touched: 2007-01, 2009-05, 2011-04.
- Scripts: `scripts/pre_cap2_queue.sh`, `scripts/pre_cap2.py` (guard, checks, report). Log: `logs/wa_pre_cap2.log`.

## Checks (all passed)
- **Data completeness before the runs** (`guard_summary.json`), on the 20 usable stations:
  - V2008: 1 of 513 station-days missing or unreadable (0.19 %, PB.B018).
  - A2011: 0 of 423.
  - Counting all 21 stations, V2008 is 27 of 539 = 5.0 %. 26 of those are PB.B926, whose raw data is on disk but
    which has no template. This is the same as the go/no-go V2008 supplement, and the rule was fixed in advance to
    count the 20.
- **After each run**: every cached station-day was loaded, or dropped by the frozen rule (> 5 % of the day guarded:
  12 on V2008, 9 on A2011). 0 were lost for any other reason, in all 12 matched-filter runs.
- **HF stage**: every station-day job returned ratios (256/256, 255/255, 291/291, 527/527, 518/518, 767/767).
- **Reproducibility**:
  - The uncapped reruns are identical to the go/no-go matched-filter output: V2008 653,475 candidates, A2011
    894,286, max |Δ net_mad| = 0.
  - The V2010 rows reproduce `cap/cap_comparison.csv` exactly: own precision 0.3177 / 0.3653 / 0.6765 and
    family-windows 81 / 111 / 128.

## Per-episode table (NW, 132 families observable; `cap2_comparison.csv`)

Column meanings:
- **fam-windows ≥ 2× both controls**: the family's raw ≥ 8.75 count on label days is at least 2× the higher of the two
  scrambled controls.
- **ctl vs ctl**: the same test with control 2 against control 1, i.e. how often chance alone passes.
- **reversed control**: ≥ 8.75 × MAD candidates of the time-reversed templates per observable family-day, and how many
  family-windows have data ≥ 2× the reversed control.

| episode | run | own precision | own recall | fam-windows ≥ 2× both controls (ctl vs ctl) | own precision of scrambled control | data / control 1 / control 2 per fam-day | reversed control per fam-day | fam-windows data ≥ 2× reversed | median dominant-site share used |
|---|---|---|---|---|---|---|---|---|---|
| V2008 | uncapped | 0.287 | 0.124 | 86 (8) | 0.004 | 8.96 / 6.38 / 6.43 | 0.67 | 132 | 0.34 |
| V2008 | capped 0.3 | **0.316** (+0.030) | **0.099** (−20.5 %) | **114** (10) | 0.003 | 6.61 / 3.77 / 3.82 | 0.20 | 132 | 0.30 |
| V2008 | dominant site out | 0.682 | 0.012 | 125 (21) | 0.000 | 0.80 / 0.06 / 0.06 | 0.03 | 127 | 0.19 |
| A2011 | uncapped | 0.197 | 0.128 | 37 (9) | 0.003 | 9.40 / 8.39 / 8.40 | 1.40 | 109 | 0.44 |
| A2011 | capped 0.3 | **0.280** (+0.083) | **0.092** (−27.8 %) | **86** (18) | 0.005 | 3.98 / 2.88 / 2.89 | 0.15 | 120 | 0.30 |
| A2011 | dominant site out | 0.487 | 0.015 | 44 (24) | 0.006 | 0.35 / 0.19 / 0.17 | 0.09 | 66 | 0.30 |
| V2010 (earlier, reproduced) | uncapped | 0.318 | 0.124 | 81 (4) | 0.009 | 10.13 / 6.56 / – | – | – | – |
| V2010 (earlier, reproduced) | capped 0.3 | 0.365 | 0.107 | 111 (4) | 0.007 | 7.83 / 3.89 / – | – | – | 0.30 |
| V2010 (earlier, reproduced) | dominant site out | 0.676 | 0.018 | 128 (19) | 0.007 | 1.53 / 0.13 / – | – | – | 0.18 |

Time-reversed control, own-family "precision" (its chance level): 0.000 on V2008 in all three modes; on A2011 0.001
(uncapped), 0.002 (capped) and 0.011 (dominant site out).

## Rule outcome (`cap2_verdict.json`)

| episode | precision +≥ 0.02 | family-windows do not fall | recall falls ≤ 20 % | pass |
|---|---|---|---|---|
| V2008 | +0.030 yes | 86 → 114 yes | −20.5 % **no** | no |
| A2011 | +0.083 yes | 37 → 86 yes | −27.8 % **no** | no |

Decision as preregistered: **NW stays uncapped.**

## Notes (not part of the decision)
- **V2008 recall** misses the bar narrowly: 0.0986 vs the 0.0993 needed. A2011 misses clearly.
- **Two episodes:** the recall loss shows up in both. V2010 was −13.6 % and would have passed; the two new episodes
  would not.
- **Reversed control:** the cap cuts it 3× (V2008) to 9× (A2011), more than it cuts the scrambled controls (about
  1.7–2.9×). One possible reading, not tested here: much of what the reversed templates pick up in the uncapped
  network is dominant-site bursts.
- **Dominant site left out:** again very pure, but recall ≈ 0.012–0.015, and its ctl-vs-ctl false-pass level is the
  highest (21–24). It is not a production option.
