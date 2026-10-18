# Olympic Arm A dominance-cap test: definitions fixed BEFORE the runs (2026-10-02 21:10Z)

Runs: 21-core network (results/wa_gonogo/setup, 20 usable stations), V2010 (39 days; the 2010 truth days), build-half
T1 templates, all >= 30-event templates, --drop-dead 0.05, 2 scrambled controls (same seeds as the go/no-go run, so
the controls differ only through the weights). Variants: frozen (= results/wa_gonogo/mf_core21, re-scored on V2010
only), cap30 (any site's share of the family's site weight <= 0.3, water-filling), lodo (dominant-weight site left out).
Scorer: scripts/gonogo_score.py unchanged, --wins V2010 (test half scored; HF cuts calibrated on the V2010 build half).

Acceptance (fixed by the main session): adopt the cap for production if
  (a) Olympic A nbr3 data/control >= 3, AND
  (b) NW own precision (prec_own_cat_core, V2010, Royer_NW_WA) within 0.02 of the frozen run's.
Else: proceed uncapped; Olympic A gets stacks / group rates only.

Definition of (a), chosen here because no confidence model P exists for the new runs: the ratio of COUNTS of
>= 8.75 x MAD candidates, data vs the mean of the two scrambled controls, counted at the <= 3-km group level:
Olympic A families are joined into groups (connected components of the <= 3-km epicentral graph, lat_use/lon_use),
per group and day the candidates of all its families are merged at +-2 s (one event per cluster), on the test half of
the observable V2010 episode days. This is the analogue of the plan's "P_nbr3 >= 0.75: data 1,043 / control 603"
(ratio 1.7) and of the go/no-go "data / control per family-day 241.9 / 240.8" (ratio 1.00).
Also reported, not used for the decision: the per-family (un-merged) count ratio, and the truth-matched count ratio
(data matched to <= 3-km truth / control matched to <= 3-km truth), which is ~8 already in the frozen run.
