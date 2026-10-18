# Merlin review of the decision-tree page (2026-09-30)
Full text returned to the main session; summary of the spec (implement in this order):
D1 spine layout: chosen card directly under each question on a vertical green spine; failed/other attempts in a side column (collapsed on phones);
   renderer orders cards (chosen first); green stops at the last COMPLETED card; running/planned = dashed blue; 'Logged:' text moves to a `log` field.
D2 status board at the top (one row per catalogue: paper target, ours real/distinct/located, model, extension state, next action, link) replacing the
   extplan Inventory + Readiness blocks; sections grouped A. Reproductions / B. Extension; table of contents; delete duplicate Bostock next cards.
D3 question `why` (because: previous result), `tests` chips (parallel tests that decided it), option `feeds` links between questions (stable ids;
   re-key localStorage), option `also:1` (adopted alongside).
D4 status vocabulary: add `plan`; fix `open` bug (renders 'undefined'); section fields paper/ours/ext drive pills; `updated` date.
D5 density/phone: remove min-width 860 px; word limits (h ≤ 8, r ≤ 25, n ≤ 40 chars) checked at build; fuller legend; 'path only' mode.
Content fixes 2a (1–20), structure fixes S1–S8, missing items 2c, and the restructured extension tree (Q1–Q8 + Z5 side box) as listed in the report.

---------------------------------------------------------------------------------------------------
FULL REPORT (verbatim, Merlin 2026-09-30)

Verdict: the content is ~80 % right but the page's design defeats its own convention: the green path is one card in a horizontal fan of up to 9 equal cards (min-width 860 px, nowrap), planned work is drawn in green, results are parked in "next" rows, and the extension section says the same thing four times. The single decision that matters: switch to a spine layout (goal → Q → chosen card → Q … straight down the middle; failed attempts as side cards or a collapsed "N other attempts") and add a status board at the top that replaces the extplan inventory + readiness blocks. Everything else is card edits.

PART 1 — DESIGN SPEC (ranked)

D1. Spine layout for the decision path. Evidence: .opts{display:flex;flex-wrap:nowrap}, .tree{min-width:860px}, .o{flex:0 1 240px}; drawTree() draws a grey bezier from the question to every card and a green one to the chosen card. Two questions have 8–9 cards (~2,300 px row → horizontal scroll). Chosen cards are usually the last on the right because cards are authored chronologically.
- Each question level becomes a two-column grid: left = the spine (question pill, then the chosen card directly under it, then the next question); right = "other attempts" cards, stacked vertically, small, 60–70 % opacity, each with an "Attempt N · date" tag. The spine is a straight vertical green line (CSS border-left on the spine column when a chosen card exists; teal/grey dashed when no option was chosen).
- Card order enforced by the renderer: chosen first, then part, bad, dead; within a class by date/att.
- The green line must stop at the last completed chosen card. Options with s:"run" or new s:"plan" are never drawn green even if chosen (today Route Steps 5–7 and the Z5 "next request" are run + chosen:1 → green through unfinished work). Draw planned spine segments as dashed blue.
- Keep the SVG only for cross-links (feeds); the spine itself is CSS.
- Card slots: tag row (status word · Attempt N · date), title (≤ 8 words), reason (≤ 25 words), number (one figure, mono), footnote link ("log ↗" from a log field). Move every "Logged: …" sentence out of r into log.

D2. Status board at the top; merge the four extplan blocks. Evidence: extplan has Inventory, Readiness, Merlin review and Route to 2026 blocks that repeat; Inventory Q1 and Readiness Q1 list the same catalogues with different numbers; Route Step 4b duplicates Merlin Q1 Attempt 1 (+52 stations); Route Step 4c duplicates Merlin Q5 (stacks); Route Step 3's chosen card (clear-only templates, 386/61) is the approach Merlin Q1 rejected.
- New top component board: one row per catalogue (SVI, Royer NW, SW, NVI, Olympic A, Olympic B, Oregon, Ducellier, Shelly). Columns: paper target (families / located) · ours real / distinct / located · location model · extension state (done / running / launching / ready / partial / not started) · next action · link to the section. Sourced from section-level fields (D4). Colour the extension cell with the status palette. Replaces Inventory and Readiness. Footer: v1.2 totals (610 / 586 / 551 / 496 with a carrier).
- Group sections under A. Reproductions (Bostock, Royer, Chestler, Ducellier, Shelly, Oregon) and B. Extension to 2026 (one decision tree + one "extra data" box). Table of contents under the lede.
- Extplan becomes: (i) board; (ii) ONE tree "Carrying families onto the permanent network" (outline in Part 3); (iii) a side box "Extra data: Z5 Rainier nodes".
- Bostock "Extension" tree next cards "South (Ducellier) extension" and "Washington extension" duplicate section B → delete; keep "Front estimator v2 (parked)".

D3. Show tests that steer one decision, and results that feed a later question.
- Question field why: one line under the question: "Because: <previous result>". Default = previous chosen card's n.
- Question field tests: array of {h, n, s, log} rendered as small chips between the question and its options, labelled "decided by".
- Option field feeds: array of question ids; renderer draws a dotted teal SVG link when both ends visible and shows chips "→ used in: <tree> Q<n>" / "← from: …"; clicking scrolls. Give trees and questions stable id fields (today ids are positional so localStorage open/shut state attaches to the wrong question after a reorder). Migrate/bump the lfe-shut storage key.
- Option field also:1: "adopted alongside the chosen option" rendered on the spine as a thin second card with an "also adopted" tag. Cases: Bostock Q2 "Realign before every stack", Royer NW Q4 "Speed growth up on the GPU".

D4. Status vocabulary and section metadata.
- Bug: s:"open" on the Shelly next card is not in K → label "undefined" and borderless card. Add statuses plan (planned) and keep run (running); map open → plan. Legend line for each.
- Section fields: paper:{families, located, note}, ours:{real, distinct, located, model}, ext:{s, note, next}. Derive the pill as "ours real / paper".
- Global updated:"2026-09-30" next to the h1.
- end box gets a paper slot rendered as "ours vs paper".

D5. Density and phone.
- Remove min-width:860px; below 720 px the side column collapses into <details> "N other attempts"; SVG cross-links hidden below 720 px (chips remain). Summary line drops below the title on phones.
- Card limits (h ≤ 8 words, r ≤ 25 words, n ≤ 40 chars) warn at build. Worst offenders: Route Step 2b, Route Step 4 Attempt 1, Merlin Q1 Attempt 2, Shelly Q7.
- Legend: add question pill, goal/end boxes, green spine, Attempt tag, dashed-blue planned, dotted feeds link, also adopted; two lines max.
- Numbers: "ours · paper" in the mono slot. Section sum ≤ 2 sentences.
- "path only" mode (existing "Questions only" button) hides side cards.

PART 2 — CONTENT FIXES
2a Factual/stale:
1 royer pill+sum "135 families found" → "241 rows · 230 distinct (NW 141/132, SW 81/79, NVI 19/19); Royer 223".
2 royer NW plateau end "NW_mute (66) stays best" → wrong (NW final = lag15 then union 141).
3 royer NW end "66 real, 65 located" → continue to "Union of 5 runs: 141 real · 140 located · 66 with erh < 5 km; 45 of Royer's 122 within 5 km (p = 0.006)".
4 royer SW end 54/54 → continue to SW union "81 real · 81 located · 58 with erh < 5 km; 29 of Royer's 54 within 5 km (p = 0.009), 46 within 10 km". SW next "Use noise families as detectors" (run) contradicts NW next (dead) → delete from SW.
5 royer NVI end 15 → union 19 (NVI_mute 15 + NVI_lag15 4).
6 chestler pill "53 real families" → "Arm A 84 real (61 at true 2–8 Hz, control 1) · 54 located; Arm B 47 new · 43 located; Chestler 43". Chestler next cards "Count the real families … 53 real · 46 located" and "Implement Chestler's own front end" (run) are stale: front end built (Arm B, 09-29) → results.
7 extplan inventory Olympic "62 real · 37 located (1–4 Hz: 53 · 46)" → "84 real (61 at 2–8 Hz) · 54 located".
8 inventory Arm B "Not located yet" vs "43 located" → v1.2 numbers.
9 inventory Oregon → "30 rows (26 + 4 candidates), 6 located; paper 18".
10 readiness Mendocino → forward scan continuous 2023-12-01..2026-09-28 (12-month overlap = continuity gate within 0.10 of v3b on 100 2024 test days); KMR re-templated per sensor epoch (2025-09-17, 2026-03-12); backward BK1 (KMPB/KMR/B049/B932–934 on 385 held-out days) runs in parallel; BK2 2008–2017 only if BK1 precision ≥ 0.5, after the forward run.
11 readiness Ducellier "running" → "launching" until markers exist.
12 Merlin Q1 Attempt 2 "families at 3× … 28 → 359 (control 10)" → "family-episode windows ≥3× control 28 → 359 (control 1 → 10)"; add costs: NW own precision 0.29 → 0.25, SW recall only 0.010.
13 Merlin Q2 Attempt 3 add: fit on the test half because the build-half fit leaks through the templates; 2011 the only fully independent check; in NW the control still carries 50–65 % of summed confidence (busy-hour context).
14 Merlin Q3 chosen add: "always the most active member" already gives 77–89 %; the gain is in the probabilities (mean P(true) 0.74–0.90 vs 0.56–0.61).
15 Route Step 1 "422 → 416" → v1.2 610 → 586 (11 groups, 35 rows). Step 2 "353 of 422" → 496 of 586 distinct. Step 3 chosen (clear-only) → red side card "Attempt 1 (09-29): clear templates only"; chosen = "Attempt 2 (09-30): keep all 11,358 ≥30-event templates, weighted".
16 extplan pill "catalogue v1 built" → "v1.2: 610 · 586 distinct · 551 located".
17 shelly status run → ok (reproduction done); fix the open card.
18 ducellier Q3 "Blind discovery on the southern permanent bed — not finished" → dead "dropped, no result" (unless known otherwise); its next "Re-run the extension with the fixed detector" → the run/plan option under Q3.
19 royer NVI Q1 two cards both "Attempt 2 (09-27)" → 2a / 2b.
20 no "null" in content.
2b Structural:
S1 results parked in next rows (NW plateau next holds union 141 and full-lag 65 while end says "did not scale"; chestler next holds USGS relocation, real count, fixed depth, and a done item marked run). Results must be questions with a chosen card; next only for run/plan/dead ideas.
S2 NW plateau → Q4 "Full-lag clustering (the paper's ±15 s) — does it change the count?" chosen part (65 real, 38 erh<5, 30 of 122 within 5 km; user 09-29: lag15 = NW final); Q5 "Each run finds sources the others miss — combine?" chosen ok union (141 · 140 located; 15 flagged <2 km; 45/122 within 5 km, p = 0.006); end = union. Fold the 176-seed test as a tests chip on Q1.
S3 SW tree → Q4 "Does full-lag clustering add SW families?" chosen union (81). NVI same (19).
S4 Chestler → Q4 "How many of the 1,479 are real?" side "Attempt 1: 1–4 Hz score (filter slip) 53", chosen "2–8 Hz score, either band: 84 (61 true band; control 1)"; Q5 "Locate them" chosen USGS NW model 54 located (side: GSC/L9; McCrory fixed depth = also/part); Q6 "Chestler's own front end (Arm B)" chosen "135 real; 47 not Arm A copies; 43 located; control 0 of 50". End "131 real · 97 located · Chestler 43".
S5 Shelly Q3 and Q4 are one attempt split in two → one question "Why does it fire on quiet days?" chosen "Attempt 3: dead channels + glitches on 664 family-days → gap guard + muting; held-out recall 0.32–0.34, precision 0.28–0.37". Q5 (iterations, no chosen) correct.
S6 Merlin Q1 Attempt 1 (+52 stations) = Route Step 4b; Merlin Q5 = Route Step 4c → one tree (Part 3).
S7 Z5 question inside Route is not a step → side box.
S8 also cases (Bostock Q2 realign, NW Q4 GPU) look like rejected alternatives.
2c Missing:
- Shelly: Mendocino blind discovery as a dead side card under a question "Are there families Shelly missed?" → "Royer chain on 25 stations, 80-day autocorrelation block: 6 'real', 0 of his 27 rediscovered, 5 look like earthquake repeaters (17–19 impulsive sites, 15–24 km, near the 2021–22 Petrolia/Ferndale sequences); stopped by user 09-30. Use his 27."
- Shelly next: forward scan to 2026 (run), backward check BK1 (run), BK2 (plan, gated), location comparison (plan). Also: v3b 17 stations, recall 0.68 / precision 0.78 (Q after iterations: "Do his 4 missing stations matter?" chosen v3b), and the S1-recipe/catalogue-rules steps (v3: 0.51/0.75) must be on the path.
- Ducellier: the fixed-detector re-run (launching) with first check = his 2008 reproduction; the old run's inflated counts as the red side card.
- Extension tree: the WA trial result that started it (recall 0.16, precision 0.25, control ≈ data; only 3+-site families survive) is the why of Q3/Q4.
- Board row for Oregon extension: "not started; only COR carries 2 families".

PART 3 — RESTRUCTURED OUTLINE, section B (extension)
Board (rows as D2), footer: v1.2 610 / 586 / 551; 496 distinct with ≥1 permanent carrier.
Tree "Carrying families onto the permanent network" (goal: detect every real family through 2026 on permanent stations; validation = array-era answer 2008/2010 + independent 2011, scrambled control alongside)
- Q1 One catalogue, no double counting? why: catalogues overlap in NW/Olympic. chosen: waveform + station-timing test → v1.2 586 distinct. feeds Q2.
- Q2 Which permanent stations see each family? chosen: clear-arrival check on array-era stacks → 496 of 586 (lower bound; Oregon only COR).
- Q3 Do families survive on permanent stations alone? chosen: none (Attempt 1 09-29 part: recall 0.16, precision 0.25, control ≈ data; only 3+-site families). tests chip: +52 stations → 15–18 followable (part). feeds Q4.
- Q4 Why so few? chosen: keep all 11,358 templates weighted (T1) → observable 61 → 274, SW 2 → 73 (precision 0.84, recall 0.010); side: clear-only 386 (bad), +52 stations (part).
- Q5 How do we trust one repeat? chosen: confidence fit on the array-era answer with PNET (2011 AUC 0.87; slope 0.93; control P>0.5 0.7 % vs 2.2 %; caveats item 13); side: data-vs-control lfdr (part, slope 0.5–0.65), PNET in lfdr (bad); tests chip: PNET fine-tune 0.77.
- Q6 Neighbours 1–3 km apart — detect as groups? chosen: family templates + nearest-3-site member probability (item 14); side: group templates (bad).
- Q7 Can post-array stations learn a family? chosen: stack at confident times (NW 0.54, SW 0.33, Arm B 0.42, Olympic 0.85; control 0.14; SW/Arm B coverage 45–65 %); side: all MF peaks (part, NW 0.27), PNET times (bad). feeds from Q5.
- Q8 Families too weak for single repeats? chosen: stacks 1–30 d (2011 1-day AUC 0.84, 20 best families; ETS peaks 2012–14, control flat).
- end (ok): method fixed. next (plan): template post-2014 stations + scan NW/Olympic; SW/Arm B partial → field deployment fallback; completeness + fronts.
Side box "Extra data: Z5 Rainier nodes": day-1 test result (ok), next request (plan).
Other sections needing the same treatment: Royer NW plateau (S2), SW/NVI (S3), Chestler (S4), Shelly (S5 + 2c). Bostock and Oregon structurally sound (only also tags and next-card deletions).
Verify before publishing: render with the open → plan fix (no "undefined"); re-key localStorage after reordering.
If the user prefers the fan layout: keep it but enforce chosen-first ordering, stop green at the last done card, collapse rows > 4 cards, and add the board.
