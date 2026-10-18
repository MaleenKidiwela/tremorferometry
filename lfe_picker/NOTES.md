# NOTES — lfe_picker

## 2026-09-17 plan drafted
PLAN.md written by the orchestrator from the user discussion (array detect-and-locate, cross-station SSL pretraining, family embeddings, injection curriculum; P0-P6 gated; three-gate rule). Sent to Merlin for review before anything is built.

## 2026-09-17 Merlin review received
14 flaws; two premises overturned (√N gain already realised by the template sum; Lin's catalog recovers ~10-15% of Bostock, not 90%). Re-scoped to a discovery engine gated by a CPU oracle (P0-O). Amendments appended to PLAN.md. Awaiting user go for P0-O.

## 2026-09-17 user go (conditional)
User: launch picker creation (Opus agents, Fable orchestrating) once the per-family cumulative figure from Track E (lfe_extension) is delivered. Order: P0-O oracle first (per Merlin), then P0/P1.

## 2026-09-17 P0-O delivered (oracle + referee freeze) — agent P0-O

All four P0-O items and the P0 label table are done; full write-up `results/P0O_RESULTS.md`.
Work lives entirely in `lfe_picker/` (`scripts/o_*`, `scripts/o1..o3_*`); `../lfe_extension`,
`../lfe_catalog_repair`, `../codex-LFE` were read only, nothing deleted.

**O1 referee freeze.** `results/referee_protocol.json` + `results/holdouts.json`, both hashed.
2013 = FINAL referee (one look, spent once here in O2 and recorded); 2011 = one evaluation per stage
gate (one spent); 2010 + 2012 = development; 2008/2009 = training breadth. Time base for EVERY year =
the frozen 2012 per-family offsets, never refit. sha256 of the Bostock CSV `04cb6d76703e9fdc`, offsets
`abfdd198b5ee7a47`, templates `64b45a5bf6a34c1a`, frozen point `5dc7ea70640decfe`, stations
`50a3999c6bbe8e4e`. Holdouts, seed 20260917: 30 families in three NN strata (isolated ≥5 km: 8,10,12,
15,29,62,87,113,115,119; 2.5-5 km: 1,2,4,7,18,28,59,69,88,112; <2.5 km: 17,26,27,43,58,60,78,82,95,125)
and 3 stations B926/NLLB/LRIV. 2008/2009 ETS windows downloaded (851 station-days, 17 GB) —
**PB boreholes only; CN/C8 return HTTP 204 at IRIS for 2008/2009 and the NRCAN FDSN endpoint is
unreachable**, so pre-2010 breadth has no CN surface station.

**O2 P1a baseline: ACCEPT.** 1,989 station-days cached and the frozen configuration run unchanged on
2011 and 2013. 2011 own recall 0.1619 (T_core) / 0.2002 (T_ext) = **1.11× / 1.09×** the 2012 S15
reference (0.1463 / 0.1841) — inside the 0.7-1.3× accept band, nowhere near the <0.5× time-base
failure mode. 2013 (the single allowed look, recorded and closed): 0.1646 / 0.1985. Quiet-day network
false events/day: 2011 6.0 / 21.0, 2013 12.3 / 32.8. Merlin's pre-registered 12/day for T_core is
confirmed to 3 %; his 116/day for T_ext is NOT reproduced (~30/day) — O3 still evaluates at both.

**O3 headroom oracle: DISCOVERY-ONLY, decisively.** Energy-network recall / frozen-template recall =
**0.00** at 12/day and at 116/day (bar for discovery-only is <1.25×). The energy statistic is real but
weak (network AUC 0.785 vs the template's 0.995) and its quiet-day tail is glitch-dominated (peaks p99
14, max 350 summed MAD units, against a positive p90 of 15.5), so at any tolerable false-alarm rate it
recovers nothing. To match the template's T_core recall it needs 7,519 false events/day (626×).
ONE genuinely useful finding: at the seven FAR boreholes (B001/B003/B004/B005/B006/B007/LRIV) the
source-agnostic energy statistic BEATS the family template (AUC ≈0.64 vs ≈0.53, i.e. the template is
near chance there) — that, not the near sites, is where a learned per-station CF has headroom, and it
is exactly where the labels are predicted rather than picked. Single-station recall at 100 peaks/day at
PGC: template cc 0.119, energy 0.001 — so Merlin's P1b bar of 0.30 is a 2.5× jump over the best
statistic on disk and should be treated as a test, not an expectation.

**O4 DEFERRED**, nothing queued: `results/o4_deferred.md` (needs B045/B935/RVIT/PETY, none on disk;
Shelly's 27 families / 61,441 detections ARE on disk).

**P0 labels:** `results/labels_v0.parquet`, 1,479,426 rows — 1,011,174 pos_S (697,325 picked /
313,849 predicted), 433,229 pos_P, 704 certain-negative station-days, 32,644 hard negatives, 1,675
unlabelled ETS station-days; 2008-2013. B004 has zero stack S picks — all 38,925 of its positives are
predicted. Hidden families (421,496 rows) and held-out stations (250,431) are FLAGGED, not dropped.
Trap caught: the PNSN tremor catalog starts 2010-01-22, so a naive "no tremor that day" rule declared
every 2008/2009 day quiet, ETS days included — the builder now refuses certain negatives outside the
catalog's coverage and also requires zero Bostock events that day.

## 2026-09-17 10:34 UTC P0-O delivered (results/P0O_RESULTS.md)
- Referee frozen: 2013 final (one look spent: 0.165/0.199), 2011 one per gate (0.162/0.200 = 1.11x/1.09x the 2012 S15 reference -> P1a ACCEPT), 2010/2012 dev, 2008/09 training-only; frozen 2012 offsets for all years; 30 hidden families (10 per NN stratum, seed 20260917) + B926/NLLB/LRIV held out.
- Quiet-day network-wide false events/day: T_core ~6-12 (Merlin's 12 confirmed), T_ext ~21-33 (Merlin's 116 NOT reproduced).
- ORACLE: DISCOVERY-ONLY. Energy-network recall at matched FA = 0.00x template (template AUC 0.995 vs energy 0.785; matching T_core recall would cost 626x more false events). Headroom exists only at the 7 far boreholes (energy AUC 0.64 vs template 0.53). Single-station PGC recall at 100 peaks/day: cc 0.119 -> Merlin's P1b bar 0.30 is 2.5x the best on-disk statistic.
- Labels v0: 1.48M rows (1.01M pos_S: 697k picked / 314k predicted; 433k pos_P; 704 certain-neg station-days; 32.6k hard negs; ETS non-event = unlabelled); PNSN starts 2010-01-22 -> no certain negatives 2008/09.
- 2008/09 download: PB boreholes only (851 station-days, 17 GB); CN 2008/09 unavailable at IRIS (HTTP 204) and NRCan unreachable. O4 deferred.
- NEXT: P1b baselines + P2 core (discovery scope) and Gate 0-3 machinery, two Opus agents.

## 2026-09-17 P1/P2 build started (agent P1P2)

Workspace rules honoured: everything new is `scripts/p_*`, `scripts/p1b_*`, `scripts/p2_*`,
`results/p1`, `results/p2`, `results/mm`, `figures/p_*`; the `o_*` scripts/results are read-only, as
are `../lfe_extension`, `../lfe_catalog_repair`, `../codex-LFE`.  Nothing deleted.

**Infrastructure.** `p_build_mm.py` -> one memory-mapped float16 store of every cached S15
station-day 2008-2013 (3,735 station-days, 38.7 GB, `results/mm/days_zne_f16.dat` + `day_index`),
ROTATED TO ZNE with the StationXML azimuths.  The boreholes are NOT at 0/90 (B001 335.5/65.5,
B003 239/329, B004 167.5/257.5, B005 90.5/180.5, B006 201/291, B007 36.5/126.5, B010 228.5/318.5,
B011 11/101, B926 299/29; only B009 is 0/90) -- so the frozen chain's `EH1 -> N, EH2 -> E` naming is
a mislabel for the learned model, which is why the rotation is done here.  B006/B007 azimuths were
not in `a4_stationxml_audit.csv` and were fetched from IRIS (`results/p_stationxml_B006_B007.xml`).
The template path in `o_*` is untouched.

**Time-base check (PASS).** Stacked S-window energy at the labels_v0 S times, 2012, 4,000 events per
station: PGC +0.40 s peak / 2.51x pre-window, B010 +0.25/1.99, VGZ +0.35/1.95, SNB +0.50/1.59 --
the time base is sound.  BUT the far boreholes peak late: B005 +0.75, B001 +1.15, B003 +1.95,
**B004 +3.15 s**.  B003/B004 labels are `predicted` (Crosson + a5 correction); the >= 0.5 s label
uncertainty P0-O warned about is really 2-3 s at those two stations.  Recorded, reported by
provenance, and the reason P2 fits a per-SITE CF lag on the dev year before back-projecting.

**Grid (P2).** 1,444 cells of 3 km over 48.0-49.0 N / -124.5..-123.0 E, depth fixed on Slab2; all
130 families inside; nearest cell 1.22 km median.  Slab2 sits **7.5 km SHALLOWER** than Bostock's
family depths (p10 -17.3, p90 +6.6), which leaves a 0.50 s median RMS differential travel-time
mismatch between a family and its own nearest grid cell -- reported, not fudged.

## 2026-09-17 Gates engineer — g_gates built and validated (agent GATES)

Deliverables (all inside `lfe_picker/`, prefix `g_`; `../lfe_extension`, `../lfe_catalog_repair`,
`../codex-LFE` read only, every borrowed routine COPIED): `scripts/g_common.py` (day store, Crosson
travel times, stacking, `align_stack`, native-rate HF ratio), `scripts/g_locate.py` (the a5/d3
estimator + a 2-D Slab2 INTERFACE mode + synthetic floor), `scripts/g_gates.py` (Gates 0-3,
`Proposal` / `GateConfig` / `GateContext` / `evaluate`, CLI), `scripts/g_rediscover.py` (P3 scorer),
`scripts/g_crossfam.py`, `scripts/g_controls.py`, `scripts/g_make_proposals.py`,
`scripts/g_make_g3_proposals.py`, `scripts/g_validate.py`, `scripts/g_report.py`,
`scripts/g_figures.py`, `scripts/a5_env.py`.  Results in `results/gates/`, write-up
`results/GATES_RESULTS.md`, one-page usage doc `results/gates/USAGE.md`.

**Validation headline** (full tables `results/GATES_RESULTS.md` / `results/gates/g_tables.md`;
734 waveform proposals + 180 Gate-3 proposals + 234 control tests):
- **Gate 1** (withheld VGZ/SNB/B926, 200 random-time + 100 moveout-shuffled nulls per station):
  power **0.77 / 0.91 / 0.92** at N = 20 / 50 / 100 (same on Bostock times and on frozen template
  detections: 0.78 / 0.92 / 1.00); median excess +3.6 / +5.5 / +7.2 dB — Track A 3b (76/81 at N≈62,
  +7.6 dB) reproduced. **False pass 1/150** over three negative mechanisms (wrong cell 1/50, time-shift
  0/50, **dv/v families 0/50 even with a ±30 s offset search**).
- **Gate 0**: Merlin's literal rule false-passes **37/50 aliases** (at low-SNR sites the cc is trivially
  below chance). Adopted a corrected rule — cc below chance at ≥2 REGISTERING sites and at >half of
  them — which gives **6/50 alias false-pass**, 100 % power at ≥5 km, 0.64-0.71 at 2.5-5 km, 0.19-0.32
  at <2.5 km (56 % of Bostock families are in that unresolvable regime). Class HF veto fires 6/6 on
  HF-rich proposals, 0/734 spurious; measured per member, not on the coherent stack (the coherent stack
  is depleted 2.6x at 10-19 Hz — `g_hf_stack_vs_member.csv`).
- **Gate 2**: geometry would admit 113/130 families, but the PICK yield binds — testable 0/116 (N=20),
  1/109 (N=50), 9/102 (N=100), 21/62 (N=200), of which 10 pass. Median proposal-to-solution distance
  5.6 km at N=200 while all four controls pass: the controls certify data-driven, not accurate.
- **Gate 3**: validated on a LONG span (catalog v1 2014-2026; on a single 52-d ETS window the ≥30 d
  circular shift does not exist and the gate correctly returns UNTESTABLE). Power 54/60; 0/60 on
  random-time controls; **0/60 on activity-covarying controls** although they pass the shift null
  100 % of the time — the local-vs-distant contrast is what closes Merlin F13.
- **Runtime** 3.7 s median / 9.5 s p90 per proposal for Gates 0+1+2 (warm cache, CPU only; the H100 was
  not needed) — ~1,000 proposals/hour/process, so P5 screening is cheap.
- **P3 scorer calibration**: the frozen template detector's own clusters score only **9/30** hidden
  families (isolated 5/10, mid 1/10, close 3/10), spurious-split rate 0.01. The ≥30 % one-to-one recall
  clause is above the template detector's own per-family recall, so Merlin's "≥2/3 of the isolated
  hidden families" bar is above what any detector matching the templates could reach — read P3 against
  the 0.50 ceiling as well as against the bar.

## 2026-09-17 11:21 UTC Gates 0-3 delivered (results/GATES_RESULTS.md, results/gates/USAGE.md)
- Gate 1 power 0.77/0.91/0.92 at N=20/50/100; false pass 1/150 (dv/v families 0/50 even with +-30 s offset search). Gate 0 (adopted rule: cc below chance at >=2 REGISTERING sites and >half) alias false-pass 12% (Merlin's literal rule 74%). Gate 2: pick yield binds, needs N~150-200; controls certify stability not 5-km accuracy (solution 5.6 km from proposal). Gate 3: 54/60 power, 0/60 false pass incl. activity-covarying (local-vs-distant contrast closes F13); untestable on <60-day spans.
- Runtime 3.7 s median per proposal (CPU); ~1,000 proposals/h/process.
- P3 ceiling: the frozen template's own clusters rediscover only 9/30 hidden families under the >=30% one-to-one clause (template per-family recall ~0.18-0.20) -> read Merlin's 2/3 bar against a 0.50 ceiling.

**P1b done (2026-09-17 ~12:15).** Two baselines trained from scratch on the P0 labels: Lin's
`make_large_unet(fac 2, sr 40, ncomps 3)` re-implemented in torch (15-s windows, 6 log/sign channels,
our 20 Hz data polyphase-upsampled to his 40 Hz) and a PhaseNet-style 20 Hz / 60-s U-Net, plus a
w_pred = 1.0 sensitivity run.  60,000 steps x batch 256 each (~3.5 GPU-hours total).  PNET wins at
every peak rate on both dev years -> selected.  **One 2011 evaluation spent** (ledger updated).
Full write-up `results/P1_RESULTS.md`.  Headline: Merlin's bar (0.30 at <= 100 peaks/station-day at
PGC) NOT MET (0.0218) and **arithmetically unreachable** -- one CF series per station-day caps recall
at 100 x station-days / events = 0.126 at PGC in 2011.  At the template's own total peak budget
(11,600 peaks/station-day = 100 per template-day x 116 families) PNET reaches 0.486 vs the template
cc's 0.119.

**SSL pretraining GATE: FAILED, pretraining dropped.** Cross-station envelope prediction trained on
2014-2016 only, targets >= 15 km from every input station; linear probe on 2012 Bostock windows vs
in-tremor non-event windows of the same station-days, grouped 4-fold: **AUC 0.757** (folds 0.766,
0.749, 0.765, 0.747) against Merlin's 0.85 bar.  Exactly his F4 prediction.  `results/p2/ssl_probe.json`.
(Operational note: the 2014-2016 memmap is on NFS and random 7 kB reads ran the trainer at 1 % CPU;
fixed by pulling 40 whole days into RAM first.)

**P2 delivered (2026-09-17 ~14:30).** Full write-up `results/P2_RESULTS.md`.  CF = the P1b-selected
PhaseNet U-Net retrained 100k x 256 with a station embedding under 30 % dropout and **real-waveform
transplants at random Slab2 cells (6,135,239 windows = 24.0 % of samples, cap 30 %)**; back-projection
= 1,444 Slab2 cells x 0.25 s, Crosson + a5, site-merged log-odds, per-cell per-day MAD, threshold from
quiet-day FA only.  NO learned location head.

*The lever nobody had tried:* Merlin F6 permits the unlabelled ETS-day tremor to be EXCLUDED **or**
used with nnPU.  Excluding it (`p2core`) produces a CF that saturates at ~0.20 across a whole ETS day
and rises only to ~0.27 at a catalogued LFE (1.3x contrast, 2 s FWHM) -> its top network z values come
from QUIET days (glitches/EQ) and at the pre-registered FA it fires less on ETS days than on quiet
days.  **nnPU with pi = 0.10 (set a priori from the ~3 % event duty cycle / ~1/3 catalogue
completeness) raises dev recall at 6 FA/day from 0.0022 to 0.0428 (20x) and the ETS/quiet contrast
from 0.9x to 30x.**  Chosen on dev, before the 2011 look.

**2011 (one pass, budget SPENT):** recall **0.0363 at 6 false events/day** and **0.0674 at 21/day**
vs the frozen template's 0.162 / 0.200 -> **0.22x / 0.34x** (chance-corrected 0.15x / 0.21x).  Hidden
families recovered slightly BETTER than kept (0.0418 vs 0.0341) -> the model learned a class, not the
100 training families; recall rises 8.7x with magnitude -> the detections are real events.
**61 % of detections land > 5 km from every known family, but the median location error is 25.6 km**,
so that number means "the argmax is far from the known patches", NOT "61 % are new sources".
Calibration vs withheld-station realness: NOT within +-0.1 (max 0.29; <= 0.11 over the four bins
holding 96 % of detections); the withheld-station excess rate does rise monotonically 0.12 -> 0.30
against a 0.05 chance floor.

**PRE-REGISTERED MOVEOUT-SHUFFLE NULL: FAIL (0.90-0.93 vs the <= 0.10 bar) -- and the diagnosis
matters.**  The shuffle is real (differential moveout changes by 10.2 s RMS vs a 6.2 s true spread),
but the detection statistic is a MAX OVER 1,444 CELLS, and with that many cells a wrong travel-time
table still finds some cell that matches a tremor burst.  The supplementary **fixed-cell** null
(`scripts/p2_fixedcell_null.py`, maximisation removed: read the site-summed CF at each Bostock event's
OWN family cell) **passes decisively** -- true 4.75 site-MAD units vs shuffled 2.65 vs random-time
0.15; 6.9x enrichment above the shuffled p99, 13.6x above the random-time p99.  **A corrected null
must be pre-registered before P3 and no discovery claim should be made until it is.**

*Housekeeping:* nothing pre-existing was deleted.  The only files removed were four smoke-test
checkpoints this agent had created minutes earlier (`results/p1/models/smoke*.pt`,
`results/p2/models/{tpsmoke,nnpusmoke}.pt`).  `results/mm/` (91 GB: the 2008-2013 ZNE day store plus
the 2014-2016 SSL subsample) and `results/p2/cf_*/` are rebuildable from `scripts/p_build_mm.py` and
`scripts/p2_cf.py` if disk is ever needed.

## 2026-09-19 02:20 UTC orchestrator status after P1/P2 (agent hit the weekly Opus limit while finalising; results complete on disk)
- P1b: PhaseNet-style U-Net beats the Lin re-implementation on both dev years; Merlin's 0.30@100 peaks/day bar NOT MET (0.022) and arithmetically unreachable (ceiling 0.126 = 100 peaks x 22 days / 17,463 events at PGC); at the template's total peak budget PGC single-station recall 0.486 vs template cc 0.119 (4.1x); matches pretrained LFEDetect at 720 peaks/day. Far-borehole labels (B003/B004) are 2-3 s late (predicted labels).
- P2: network recall at matched quiet-day FA 0.036/0.067 vs template 0.162/0.200 (0.22x/0.34x) — discovery-only as licensed; moveout-shuffle null FAILED (0.90 vs <=0.10) because max over 1,444 cells lets a wrong table match some cell; fixed-cell null passes decisively (6.9x above shuffled p99); 61% of detections >5 km from known families BUT median location error 25.6 km -> not new sources; calibration not within 0.1; SSL pretraining failed the probe (AUC 0.757) and was dropped. 2011 spent for P1b and P2; 2013 closed.
- DECISION NEEDED (Merlin): corrected null pre-registration; whether P3 is worth running with 25-km per-event locations; whether discovery should instead be cross-station-correlation clustering triggered by the CF, with Gates 0-3 as the judge.

## 2026-09-19 02:33 UTC Merlin P1/P2 review received; P3 launched
P2 closed as discovery engine; P3 = template-free repeater search with Bostock-times positive control first; 2013 kept sealed (user to confirm). Opus P3 agent launched on Steps 0-2 (stop conditions per Merlin).

## 2026-09-19 P3 engineer — Steps 0 + 1a delivered; **Step 1a STOPPED the phase** (agent P3)

Write-up `results/P3_RESULTS.md`; addenda appended to `results/P1_RESULTS.md` and
`results/P2_RESULTS.md`. All work under the `q_` prefix; `o_*`, `p_*`, `g_*` untouched (the one
override needed — withholding stations from the frozen matched filter — is an in-process
monkey-patch in `scripts/q_run_month.py`). **2011 was NOT re-run** (Step 0 re-scores stored products
of evaluations already spent), 2013 not opened, `referee_ledger.csv` unchanged.

**Step 0 (closure).**
- (a) `q_p1b_chance.py`: every P1b curve chance-corrected with the ±307 s shift. Chance at 11,600
  peaks/station-day is **0.24**, so PGC 2011 goes 0.4856 → **0.3209**. The template cc denominator
  (0.1194) has chance 0.0000 (read at the family's own S). **"4.1×" → 2.69× (2011) / 2.74× (2012).**
  Far boreholes 0.356 → 0.123; B003 @100/day is at chance (−0.0003).
- (b) `q_step0_bc.py`: Merlin's corrected event-level null reproduced exactly — the moveout-shuffled
  detector keeps **0.600** of the true chance-corrected recall at 6/day (0.62/0.71 at 21/100; dev
  2012 0.62–0.75). Uniform-over-grid expectation for ">5 km from every family" = **64.5%** vs the
  observed 60.9% (2011) / 59.1% (2012) — the 61% is withdrawn.
- (c) loud-day collapse: 2011 doy 239–242 → P2 emits **10/10/5/10** detections against 2262/1225/
  1447/1169 Bostock events (template 402/201/195/208); 2012 doy 253/255/257/260 → 1/17/13/17.
- (d) `q_rediscover.py` with RELATIVE clauses (recall ≥0.5×, purity ≥0.8× the truth template's own
  per-family value), plus NO-REFERENCE for families the truth template cannot detect. Reference from
  a frozen-template run with **SNB+B926 withheld** (`q_run_month.py --drop-stations`, 2012-08..10,
  `results/p3/trackQ/`): isolated-hidden median recall_ref **0.196**, purity_ref **0.567** — i.e. the
  old absolute 0.30/0.50 clauses were above the truth templates' own level on 8 of those 10.

**Step 1a — pair-linkability kill test: STOP.** All 27,681 usable Bostock 2012 ETS-day events, 13
stations / 7 sites, 24-s 3-C windows at the grid-centre constant, ±3 s lags, per-(station, day-block)
MAD, site mean, site sum, max over lag; null = per-station circular shift of index×1800 s.
3.83e8 pairs in 36 s on the H100 (fp16 matmul formulation).
- **0 of the 10 isolated hidden families have ANY same-family pair above the null's 1e-6 quantile**
  (bar: ≥30 for ≥5 of 10; STOP if <3). Best family: 3 pairs at the 10× looser 1e-5 quantile.
- Signal is REAL but ~1000× too small: same-family enrichment over far cross-family is
  8.7×/12.6×/22×/38× at the 1e-3/1e-4/1e-5/1e-6 quantiles; only 0.56% of same-family pairs reach
  even 1e-3, where the false-edge ratio is already 11:1.
- **Trap found and fixed:** the first null put 7,472 pairs above 100 with a clean gap below — all of
  them candidate pairs <1 s apart (Bostock lists simultaneous events under several families), whose
  24-s windows overlap ~100% so cc≈1 in BOTH runs. Excluding pairs closer than window+lag (27 s)
  dropped the null 1e-6 quantile from 172.06 to 23.32. Any future pair/graph work must carry this.
- Audits: (i) an in-distribution null (far cross-family pairs) gives thr 19.83 and still 0/10;
  (ii) window ablation — oracle placement + oracle 8-s duration buys only ~3× enrichment, so the
  template-free 24-s grid-centre window is not what fails; (iii) event-vs-own-family-STACK on the
  same windows gives +4.6 z-excess where event-vs-event gives +0.32 — the machinery, time base and
  ZNE rotation are sound, the pair statistic is 14× weaker (the √(signal fraction) penalty);
  (iv) to reach 30 linkable pairs a family would need ~1,670 events in one ETS; the biggest 2012
  family has 561.
- **Verdict: template-free single-ETS discovery is not viable on S15 → Merlin Step 6** (re-scope P5
  to dense post-2019 regions and/or multi-ETS accumulation before stacking). Steps 1b and the
  Step-2 trigger ablation NOT run.

**Step 2, CF completeness kill test (independent of the chain, so it was run).**
`q_cf_completeness.py`: site-merged learned CF (log-odds, lowhalf per station-day) read at the
predicted S of every stored 2012 ETS-day template detection (60,671 with ≥5 sites), AUC against the
own-family one-to-one Bostock match. **7–8 MAD: 0.755 | 8–9: 0.689 | ≥9: 0.645.** Merlin's rule on
the 8–9 bin → "in between, report". Upper bound (CF trained on 2012); `match_own` is itself a lower
bound on realness. Recommendation: do NOT retrain yet — repeat on 2010 (CF saw no labels there)
against a Track-E realness target first.

**Decisions needed from the orchestrator:** (1) accept the Step-1a STOP and re-scope to Merlin
Step 6, or ask for a specific variant first (multi-ETS accumulation on 2010+2012 is the cheapest
test, ~2 GPU-hours, and is the one thing that could move the pair rate); (2) whether the 2010 CF
completeness repeat is worth doing now; (3) 2013 remains sealed (unchanged).

## 2026-09-19 03:12 UTC P3 Step 1a STOP (results/P3_RESULTS.md)
- Pair-linkability on Bostock's own 2012 events (27,681 events, 13 stations/7 sites, 3.8e8 pairs): same-family pairs above the null 1e-6 quantile = 0 for all 10 isolated hidden families (bar >=30 for >=5). Signal real (same-family enrichment 38x at 1e-6) but ~1e3 too small; four audits (in-distribution null, oracle window, event-vs-stack +4.6 z vs event-vs-event +0.32, pair-count arithmetic) confirm. Trap fixed: overlapping-window pairs (<27 s apart) inflated the null.
- Step 0 corrections: P1b gain 4.1x -> 2.69x/2.74x chance-corrected; 61% -> uniform-grid 64.5%; loud-day collapse tabulated; relative rediscovery clauses implemented (truth-template isolated-family recall_ref median 0.196).
- CF completeness AUC (8-9 MAD bin) 0.689 = in-between; do not retrain; repeat on 2010 vs a Track E target if ever pursued.
- VERDICT: template-free single-ETS repeater discovery NOT viable on S15. Next cheap lever: pooled 2010+2012 candidates (quadratic pair gain, ~2 GPU-h) — ordered; then Merlin Step 6 re-scope to post-2019 dense regions / literature seeds.

## 2026-09-19 (later) P3 engineer — coordinator-approved follow-up: POOLED 2010+2012 Step 1a, site ablation, Step-6 sizing

New sections `results/P3_RESULTS.md` § 5 (pooled repeat + ablation) and § 6 (Step-6 re-scope sizing).
Dev years only; 2011/2013 untouched; ledger unchanged. New scripts: `q_cache2010.py`,
`q_pairlink2.py`, `q_pool_report.py`, `q_figs2.py`; figure `figures/q_pairlink_pooled.png`.

**Two things had to be fixed before the pooled run meant anything.**
1. The mm day store has NO B001-B007 for 2010. The raw traces DO exist in `../lfe_extension/traces`
   (doy 222-256), so `q_cache2010.py` built the 209 missing station-days of 20 Hz cache by calling
   `o_cache_month.cache_one` as a library (frozen chain, unmodified). All clean.
2. `q_pairlink2.py` reads the 20 Hz caches instead of the mm store — and those hold RAW counts
   (CN.PGC 2012 peaks at 1.6e7) which **overflow float16 → NaN thresholds**. Fixed by rms-scaling
   each window before the f16 cast (free: every downstream normalisation is per row) + an assert.
   The mm store was safe because it is per-day scaled. **Carry this trap forward.**
   Also: `--topk 50000` clipped the 1e-3/1e-4 null quantiles (rank > 50k at 7.8e8 pairs) -> 1e6.

**Control:** the pooled code on 2012 alone reproduces the original mm-store § 2 result exactly —
same-family pairs 3,994,687 vs 3,986,404, above 1e-6 **11 vs 11**, above 1e-5 **78 vs 78**.

**Pooled result (39,563 events = 11,853 from 2010 + 27,710 from 2012; 782.5M pairs):**
- isolated hidden families above the pooled 1e-6 quantile (22.59): 0,1,1,0,0,0,**3**,0,0,0 —
  **0 of 10 reach the bar of 30.** STOP stands.
- the lever works and is far too small: same-family pairs 3.99M -> 8.12M (2.03x), pairs above 1e-6
  11 -> 49 (4.5x), per-pair rate 2.75e-6 -> 6.04e-6 (2.19x). Enrichment 8.4/16.8/25.7/32.1x at the
  1e-3/1e-4/1e-5/1e-6 quantiles. Alias (families <5 km) still 55-67% of the same-family rate at the
  loose thresholds.
- **arithmetic:** at the pooled isolated-stratum rate (8.6e-6) a family needs ~3.5M same-family pairs
  ~ **2,650 events** for 30 linkable pairs; the median isolated hidden family has 328 pooled events
  over 2 ETS (~164/episode), the largest 759 -> **~16 ETS episodes** (median) / **~10** (prolific).

**Site ablation — the far sites DILUTE (new, and it matters for P5):** keeping every per-site
lag-resolved z let 14 site subsets be scored from one set of matmuls. Isolated-hidden same-family
pairs above threshold: 7 discovery sites 4542/311/30/5 vs the 4 nearest of those 7
7708/**1306**/**221**/4 vs PGCq+B926+VGZ+nearest far borehole 9435/1264/153/8 (1e-3/1e-4/1e-5/1e-6).
Catalogue-wide at 1e-5: 270 (7 sites) vs **1069** (nearest 4) vs 1038 (oracle 4). Per-config
enrichment reaches **113x** (n4_B001+LRIVq+PGCq+VGZ) and 68x (o4_B926+LRIVq+PGCq+VGZ) vs 32x for the
protocol 7. NLLB (82 km), B001 (68), B004 (61) add noise faster than signal. **What P5 needs is
stations CLOSE to the source, not merely more of them.** (n4/o4 use the family position and o4 adds
back the withheld B926 — diagnostics, not detectors.)

**sqrt(N_sites) extrapolation:** same-family shift delta = 0.38 sigma at 7 sites; 7->20 sites cuts
the events-per-family requirement only 3,120 -> 1,730 (1.8x). Site count alone does not close the
gap; with episode pooling it lands at ~7-11 ETS episodes for a prolific family.

**Step-6 sizing (§ 6, nothing run):** roster `data/candidate_stations_post2020.csv`; within 80 km of
the Slab2 30-45 km band — Puget/Olympic 47.0-48.5N: **111 stations, median 30 per cell**;
Oregon 43.5-46.0N: 67 stations, median 12/cell (**not obviously better than S15 — do not pull
first**). SSE days 2020-2026 (>=50 PNSN detections/day in segment): 241 Puget, 193 Oregon, 421 union.
Pull: **Puget-only 26,751 station-days = 548 GB raw + 554 GB cache = 1.1 TB** (union of both: 3.0 TB).
FDSN wall time ~ a few hours at P0-O rates. Holdout without Bostock truth = Gates 0-3 as the only
judge (Gate 3 finally testable on a 2020-2026 span) + station holdout replacing family holdout +
pre-registered N-shift/N-quiet/N-moveout as PRIMARY evidence + cross-episode reproducibility
(8-10 episodes) + literature seeds (Ducellier/Shelly) ONLY after the still-DEFERRED P0-O O4
existence test passes on the new roster.

## 2026-09-19 03:38 UTC P3 pooled test + ablation + Step-6 sizing (P3_RESULTS.md s5-6)
- Pooled 2010+2012 (39,563 events, 782M pairs): isolated hidden families above 1e-6 = 0/1/1/0/0/0/3/0/0/0 -> 0/10 at bar; linkable pairs 11->49 (4.5x), per-pair rate 2.2x; ~30x short. Needs ~2,650 events of one family (~16 ETS episodes median, ~10 most prolific).
- Site ablation: far sites (NLLB 82 km, B001 68 km, B004 61 km) DILUTE; nearest-4 sites give 4x more linkable pairs at 1e-5; per-config enrichment 113x vs 32x. sqrt(N_sites) 7->20 buys only 1.8x; PROXIMITY is the bigger term.
- Step-6 sizing: Puget/Olympic 47.0-48.5N post-2020 = 111 stations, median 30/cell (p10 9); Oregon 67 stations, median 12/cell = not better than S15. Puget SSE-day pull 2020-2026: 26,751 station-days = 1.1 TB incl. cache; holdout = gates + station holdout + pre-registered N-shift/N-quiet/N-moveout + cross-episode reproducibility (8-10 episodes).
- Traps fixed: 20 Hz caches are raw counts -> float16 overflow (rms-scale first); topk cap clipped null quantiles.
- STATE: Step-1a STOP stands. Picker-phase discovery in the north CLOSED. Options: (A) Puget/Olympic 2020-2026 multi-episode repeater search (1.1 TB, gates-judged); (B) literature seeds (Chestler/Royer/Thomas by request); (C) far-south seeds after the O4 existence test.

## 2026-09-19 03:55 UTC P3D launched (user): dense-north 2019-2026 discovery test
Download of the ten close stations started (Opus); Merlin asked to prescribe the transfer protocol, post-2013 truth metric, and the discovery nulls. U-Net = optional trigger arm only.

## 2026-09-19 04:00 UTC P3D step 1 — download plan, SSE windows, P0 launched
- `scripts/r_download_plan.py` (from `o_download_plan.py`): IRIS `level=channel` metadata 2014-01-01..2026-09-19
  for the twelve P3D stations -> `results/p3d/station_epochs.csv` WITH per-channel azimuth/dip
  (`chan_az_dip`, for the ZNE rotation the learned/transfer steps need). All ten 3-C stations are **HH 100 Hz,
  loc '', standard ZNE (Z 0/-90, N 0/0, E 90/0)**; STW/UWFH are EHZ 100 Hz. One epoch per station, no band
  changes, no azimuth changes in 2019-2026.
- Seismometer start dates (NOT the station start dates in the PLAN): CLRS 2014-09-11, SYMB 2016-09-19,
  PTRF 2017-02-17, GOBB HH 2017-05-03, OSQM 2021-03-06, RNWY 2021-04-25, LOPEZ 2021-10-20, TURTL 2023-02-12,
  BUCKS 2023-09-17, **ODUC 2024-05-15** (the 2023-03-17 station start is the strong-motion/GNSS gear; ODUC's
  HH seismometer is a 2024 instrument, so it covers ~2.3 yr, not 3.5).
- SSE windows: PNSN `catalogs/pnsn_tremor_cascadia_full.csv`, box 48.0-49.2N / -124.6..-122.8E, days with >=50
  epicentres, clustered with <=20-day gaps, span >=10 d, padded +-5 d -> **22 episodes / 866 padded days**
  (`results/p3d/sse_windows.csv`). The six largest follow the northern ~14-month ETS recurrence exactly:
  2019-11, 2021-01, 2022-03, 2023-04, 2024-07, 2025-08. CAVEAT: the catalog ends **2026-05-27**, so S22
  (2026-04-01..2026-06-01) is truncated and any 2026-summer episode is invisible.
- Queue `results/p3d/download_queue.csv`: P0 7,953 pending station-days (866 d x in-epoch stations), P1 2,259
  (one quiet month/yr: 2019-09, 2020-12, 2021-12, 2022-11, 2023-03, 2024-12, 2025-11, 2026-06), P2 15,721.
  7,847 out-of-epoch days are written to the manifest as `no_epoch` at plan time (the worker only walks pending).
- `scripts/r_download_worker.py` + `scripts/r_download_run.sh`: <=8 concurrent IRIS requests, resumable,
  atomic manifest append, progress json every 5 min, nohup, logs in `logs/p3d_download_*.log`.

## 2026-09-19 04:13 UTC Merlin P3D review: premise corrected (one new near site, not twelve); P3D agent launched on Steps 0-5

## 2026-09-19 — P3D engineer (prefix `s_`): Step 0 frozen, Steps 1–2 machinery

Spec: `results/merlin_review_p3d_2026-09-19.md`. Everything new is `scripts/s_*.py` →
`results/p3d/` + `figures/s_*`; `q_/p_/g_/o_/r_` are read-only libraries, `e5r_stations.csv`
untouched, 2011/2013 never opened (a hard guard lives in `s_common.cache_index`).

**Step 0 (done, hashed).** `s_step0_roster.py` →
`results/p3d/{stations_p3d.csv, sites_p3d.csv, premise_family_sites.csv, pred_q2iii.csv,
pred_q2iii_recall.csv, pred_site_terms.csv, freeze.json}`.
- Roster = the frozen 15 in `STAS_ORDER` + the ten new 3-C stations + STW/UWFH as Z-only helpers
  (flagged `summable=False`, never summed). Site merge corrected per Merlin flaw 3:
  **B926q = B926+CLRS (25 m)**, **B001q = B001+OSQM (17 m)**, plus the frozen PGCq/LRIVq clusters.
  17 sites total. LRIV carries `band_p3d = HH` and is a transfer TARGET, not a template station.
- **Premise correction reproduced**: median sites within 45 km **1 → 2**, within 55 km **2 → 5**;
  82/130 families gain ≥1 *genuinely new* site ≤45 km. Counting CLRS/OSQM as new sites (the bug
  Merlin's flaw 3 warns about) would have said 93/130 and 0 → 2. Of the ten isolated hidden
  families only **12, 29, 62** gain a genuinely new site ≤45 km (10 only at 55 km) — exactly
  Merlin's list. SYMB is the only new close site: median hypocentral **44.7 km**, 66 families
  within 45 km, better than anything in the frozen roster (PGCq 52.7, VGZ 52.7).
- **Geometry rule** (median hypocentral ≤60 km, withheld out) = **PGCq, VGZ, SYMB, BUCKS**.
  B926q lands at 60.20 km and is withheld anyway. Merlin's illustrative drop list omits B003 (70.1),
  LRIVq (65.0), SNB (68.0), PTRF (66.8), ODUC (73.5) — the literal rule drops those too.
- **Q2(iii) prediction, hashed before any transfer stack was opened.** S-SNR = a + b·log₁₀(hypocentral
  km) per station class on `e5r_template_qc_2012.csv` at matched stack size; per-station residuals
  (site terms) written out as the honest uncertainty. Predicted SYMB median S-SNR(N60) 8.52 dB,
  w = 1.13 — i.e. **one more PGC-class site, not a new class of site**. Recall pushed through the E5
  model (`e5_prediction_model.py` read as a library of formulae): predicted gain over families with
  ≥1 new site ≤45 km = **D1/D2 +17.2 %, D3 +27.9 %** (mean +20.8 %).
- Freeze sha256 recorded in `freeze.json` (hidden families, withheld SITE set
  {B926q, GOBB, LOPEZ, TURTL}, episode split, the ≤30-expected-false-edge discovery threshold rule,
  the prediction table and its hash).

**Engine.** `s_transfer.py` implements Merlin Q1 verbatim (frozen `cut_windows`/`align_stack`
copied, ±20 s placement at Crosson S with zero correction, S pick on the aligned stack, 12-s
template from S−2 s, acceptance (a)–(d), LOO registration ρ, N-matched weights). Speed: one WIDE
(±50 s) cut per event, every offset a numpy slice; the 200 random-time nulls are cut in one
day-sorted pass. ~3 s per leg.

**Trap found in acceptance (b).** The wrong-cell moveout null at a SINGLE station degenerates: a
cell whose predicted S is (S−P) ≈ 6–7 s EARLY puts the family's own **P arrival** inside the
S−1..S+5 window, so the "wrong" cell scores as large an excess as the truth. Measured at PGC,
family 10: excess +16.2 dB at Δ=0 and **+17.2 dB at Δ=−5 s**, falling to ≤2 dB for |Δ| ≥ 10 s and
to −15 dB at Δ=+8 s. So (b) as written fails at the gold-standard reference station too. Both are
reported: the literal `pass_b`, and `pass_b_noP` — the same null with offsets within 4 s of −(S−P)
removed, which both PGC and SYMB pass comfortably. The literal (b) is never used alone to reject a
new station that the reference also fails (Merlin's clause (c): control fails → UNTESTABLE).

## 2026-09-19 04:44 UTC corrections and decisions
- ORCHESTRATOR ERROR acknowledged: 'ten close stations / ~12 sites' was computed as distance to the nearest family only, without co-location or per-family hypocentral geometry; Merlin's check: one genuinely new near site (SYMB) for the median family, BUCKS from 2023; CLRS=B926 site, OSQM=B001 site; GOBB/LOPEZ/TURTL/RNWY far. Lesson: run the per-family hypocentral + co-location check BEFORE quoting network density.
- USER: no author contact (emails closed) -> seeded routes for Washington/Oregon closed; south depends on discovery (Puget, if P3D licenses the method) and on Ducellier/Shelly seeds already on disk (after the existence test).
- Puget/Oregon geometry check launched (per-cell hypocentral sites <=45/55/60 km, co-located merged, vs the north reference) before any density claim is repeated.

## 2026-09-19 04:55 UTC P3D step 1 DONE — P0 (SSE windows) downloaded, QC'd, cached; P1/P2 fill running
- **P0 complete**: 7,953 queued station-days -> **7,384 data / 569 nodata / 0 errors**, **166.86 GB** raw
  miniSEED in **49 min** (~9,700 station-days/h, 57.5 MB/s sustained, 8 concurrent IRIS requests). Not one CN
  request failed at IRIS — NRCan was never needed. Layout `traces/<NET>.<STA>/<YYYY>/<DOY>.mseed`, `.mseed.nd`
  markers for the 569 verified-empty days. Median coverage 1.0000 and median 1 fragment/channel at EVERY station.
- **Delivery verdict (metadata "open" is NOT delivery) — `results/p3d/download_qc_P0.csv`, all twelve YES:**
  ODUC 270/270, OSQM 598/598, RNWY 598/598, TURTL 419/419, SYMB 864/866, CLRS 863/866, STW 838/866,
  BUCKS 308/325, LOPEZ 508/547, UWFH 794/866, PTRF 782/866, **GOBB 542/866 (0.63)**.
- **The odd one is CN.GOBB**: open metadata, whole EPISODES absent from IRIS — S05, S07, S14, S15, S16, S17
  return 0-1 days each (that is 2020-02, 2020-09..11, 2023-04, 2023-06, 2023-08..10, 2024-01), and S04 7/28,
  S18 36/59. Two of those (S07, S14 2023-04) are major northern ETS. GOBB is therefore a PART-TIME site for P3D;
  do not assume it in any per-episode site count. Other whole-episode holes: UWFH 0/29 S09 and 0/37 S10
  (2021-06, 2021-10); LOPEZ 0/32 S11 (2022-03, a major ETS) despite being live since 2021-10; PTRF 42/76 S07,
  43/82 S21, 31/37 S19; BUCKS 6/23 S16; STW 12/32 S17.
- **Cache**: all 7,384 P0 days minus 29 low-coverage days (<0.5) = **7,355 station-days at 20 Hz / 1-8 Hz**
  in `cache20hz/<NET.STA>/<YYYY>/<DOY>.npz` (**158 GB**), built through the FROZEN chain
  (`o_cache_month.cache_one`, called as a library by `scripts/r_cache_p0.py`; the two Z-only helpers go through
  `r_cache_p0.cache_one_z`, which calls the same `o_cache_month.bp` / decimation with strict=False and stores a
  (1, N) array). 3-C days are (3, 1728001) float32, order HHZ/HHN/HHE; helpers (1, 1728001), EHZ. Zero errors
  (`results/p3d/cache_p0.parquet`). NOTE: the sibling P3D science agent's `s_cache.py` calls
  `r_cache_p0.dispatch`, so its live caching during the download used the identical chain — my run re-read and
  verified every one of those files and built the 19 it had not reached.
- **P1/P2 fill running under nohup**: 17,980 station-days, `logs/p3d_download_P1-P2_20260919T045736Z.log`,
  pid in `logs/p3d_download_P1-P2.pid`. Expect ~2 h and ~400 GB. Disk after P0: 71 T free.
- Figures `figures/r_download_coverage.png` (P0 only) and `figures/r_download_coverage_all.png` (the fill).
- Commands: check `cat results/p3d/download_progress.json`; `tail -f logs/p3d_download_P1-P2_*.log`;
  restart (resumable, never re-fetches) `bash scripts/r_download_run.sh P1,P2 8`;
  re-QC `python scripts/r_download_qc.py --priority P1` (or `--all`); cache more
  `python scripts/r_cache_p0.py --priority P1,P2 --workers 14`.

## 2026-09-19 04:59 UTC P3D download P0 complete (results/p3d/download_qc_P0.csv)
All 12 stations deliver: 7,384 station-days / 167 GB in 49 min, 0 errors; cached 158 GB. 22 SSE windows (6 major ETS: 2019-11, 2021-01, 2022-03, 2023-04, 2024-07, 2025-08). Caveats: GOBB part-time (0.63; whole ETS episodes missing); LOPEZ missing S11 (2022-03); ODUC seismometer only from 2024-05-15; PNSN catalog ends 2026-05-27. P1/P2 fill running (~2 h, ~400 GB).

## 2026-09-19 05:06 UTC ALL-NETWORK INVENTORY (results/p3d_puget_geometry/, figures/t_stations_2019plus_map.png)
- 547 3-C stations / 208 sites in the box 2019+. Per Slab2 cell (all networks): sites within 45 km median North 2 (p90 5), Puget 3 (p90 9) from 2023, Oregon 2; within 60 km median North 7, Puget 13 (p90 20), Oregon 7. Roster-only Puget is 7 at 60 km -> non-roster stations double Puget's density.
- KEY FIND: four PQ-network broadbands on S. Vancouver Island, unused, each within 45 km of 36-66 Bostock families: LHLYB (2023-11), JRDRB (2021-04), ALBHB (2022-01), DAOB (2024-02). PTRF 13 families, ODUC 7. These change the northern discovery geometry far more than SYMB. Delivery being probed; added to P3D as pre-registered network variant V2.
- Puget tremor-weighted top adds: ODEER (2025-06), B018, EAGLE, SKOKO, WINDI, PAN4H, GNW, B013/B943, HDW, TKEY.

### Step 0 acceptance, Step 1, Step 2 (same day)

**Step 0 acceptance PASS.** SYMB SSE delivery 864/866 = **99.8 %** (bar 80 %). 74 station-days
masked at coverage_frac < 0.9 and excluded from every later step by a mask enforced inside
`s_common.cache_index`. GOBB delivers only 63 % — a genuine data gap; GOBB is a withheld site.

**Step 1 (co-located machinery).** CLRS vs the ROTATED B926 2012 templates: **44/52 control-OK
families match at cc ≥ 0.348, 0.846 (95 % CI 0.72–0.93) → PASS**. The bar validates itself: CLRS's
own cross-family chance level, measured from its transfer stacks against the other families'
rotated B926 templates (3,192 pairs), is p95 = **0.3435** vs B926's tabulated 0.348. OSQM vs B001 is
3/5 = 0.60 with 95 % CI 0.15–0.95 and P(X ≤ 3 | p = 0.7) = 0.47 → **UNTESTABLE, low power**, not a
machinery failure (B001 has 47 usable 2012 templates and OSQM starts 2021-03, so only 5 families
have a usable control leg). Proceeded on the powered leg; recorded.

**Step 2 (transfer).** 54 (family, SYMB) templates accepted → the ≥30 gate is met. Bank
`bank_transfer_v1/{bank_full,bank_production}.npz` = 251 / 184 (family, station) templates, hidden
families included / excluded, keys listed.

**THE number: ρ_SYMB / ρ_PGC = 0.783** (n = 64 families, same events, same N, leave-one-out
single-event registration). ρ_BUCKS/ρ_PGC = 0.296. No new station reaches PGC. Merlin's scale:
≥3 × → the pairwise discovery arm has a chance; ≤1.2 × → "SYMB is not the prize". On the same events
SYMB's stack S-SNR is **2.73 dB BELOW** PGC's where the hashed class+distance model predicted
+0.36 dB — a 3.1 dB miss, exactly the size of the per-station SITE TERMS the same fit measures
(−3.7 dB LRIV … +2.7 dB PGC, 6.5 dB spread). Distance and class predict the ORDERING of the new
stations; they cannot predict a new site's noise.

**Coordinator scope amendment (PQ), hashed as `amendment_v2` in freeze.json BEFORE any PQ stack was
opened** (`scripts/s_step0_v2.py`; PQ delivery state at hash time recorded as 0 manifest rows and 0
cached days for all six V2 stations). Two pre-registered network variants:
V1 = the frozen roster; V2 = V1 + PQ.LHLYB/JRDRB/ALBHB/DAOB (+ UW.OHOH/ODEER if delivered).
Premise with and without PQ: median sites ≤45 km **1 → 2 (V1) → 4 (V2)**; ≤55 km **2 → 5 → 8**;
families gaining ≥1 genuinely new site ≤45 km **82/130 (V1) → 100/130 (V2)**. The ≤60 km geometry
rule grows from {PGCq, VGZ, SYMB, BUCKS} to {+LHLYB, JRDRB, ALBHB, DAOB}. Predicted recall gain
**V1 +27.9 %, V2 +78.8 %**. Isolated hidden families now gain PQ sites ≤45 km: 12 (+3), 29 (+3),
62 (+4), 15 (+2), 10 (+1), 87 (+1); 8, 113, 115, 119 still gain nothing.
PQ time split as instructed: JRDRB/ALBHB calibrate S08–S13, test S14–S22; LHLYB/DAOB build on S18,
test S21/S22. Step-3 DECISION will be taken on V2; V1 is reported alongside.

## 2026-09-19 — T-study (prefix `t_`): hypocentral network density over the LFE source zones

Read-only geometry, no waveforms. Everything in `results/p3d_puget_geometry/` (scripts `t_*.py`
live there too) + `figures/t_puget_geometry.png`, `figures/t_stations_2019plus_map.png`. Full write-up:
`results/p3d_puget_geometry/SUMMARY.md`.

**Pass A (original brief, project roster, Puget vs Oregon vs North).** Slab2 cells every 5 km at
28-45 km depth; sites = stations merged within 1 km; 3-C epochs from
`nonoverlapping_triplet_epochs.csv` + IRIS `level=channel`. 2024, tremor-active cells, median sites
within 45 km **hypocentral**: North 2, Puget 3 (mean 4.0, p90 9), Oregon 3 (mean 2.9). The only band
that is materially better than the north's 1-2 is **Puget 47.0-47.5 N — median 5, gap 68 deg** (built
by ~10 UW 3-C stations installed 2020-22: WINDI, SKOKO, WYNO, ANGUS, PAN4H, HERD, OHOH, MONTE, DOORS,
EAGLE). Oregon's only match is 45.5-46.0 N; south of 44.5 N Oregon is worse than the north.

**Pass B (revised brief, EVERY network 2019+ in 46.5-49.5 N / -125.5..-121.5 E).** IRIS
`network=*`: 546 stations with a 3-C set >= 20 Hz, of which **237 are real 3-C seismometers, all
open** (303 accelerometer-only, 9 Z-only, flagged separately) -> **207 sites**, 158 alive in 2025.
- **Bostock families, sites <= 45 km: S15 median 1 (2% with >= 3) -> ALL open 3-C in 2025 median 4
  (73% with >= 3).** The discovery network already exists; we are not using it.
- The whole northern gain is **four open NRCan `PQ` broadbands nobody in this project has ever
  touched** — JRDRB (2021), ALBHB (2022), LHLYB (2023), DAOB (2024), 100 Hz, on southern Vancouver
  Island. S15+SYMB+BUCKS+PQ4 = median 3 sites <= 45 km, **72% of families with >= 3**, against
  **Bostock's own 2003-05 21-station array at 77%** (median 5). Nothing else in the north comes close
  (PTRF +2 families, ODUC/RNWY marginal).
- 16 families (deep, 40-47 km, under the Strait of Georgia) have **no** site within 45 km even with
  every station: bkeys 16, 17, 83, 88, 110, 112-114, 116, 118, 119, 123, 124, 127-129.
- Olympics/Puget tremor-weighted ranking (all UW, built 2021-25): ODUC, OHOH, ODEER, B013+B943,
  EQUIL, DOSE, SKOKO, EAGLE, PAN4H, WINDI (+ B018, WYNO, DOORS, GNW, HDW, TKEY).
- Deliverables: `stations_all_2019plus.csv`, `sites.csv`, `family_geometry_by_year.csv`,
  `cell_geometry_by_year.csv`, `station_ranking.csv`, `bostock_array_comparison*.csv`.

**Traps hit.** (i) A cached channel CSV reads an empty location code back as NaN, and
`groupby(['loc','band'])` then silently DROPS those channels — 546 3-C stations collapsed to 129
until `loc` was re-filled with "". (ii) `data/cas_slab2_input_04-18.csv` is Slab2's scattered
control-point INPUT, not the gridded model: fine for point depth queries (the project recipe), but
contours from it need smoothing + a hull mask or they are noise. (iii) Geometry is necessary, not
sufficient — site noise is untested here.

## 2026-09-19 05:11 UTC geometry study complete (results/p3d_puget_geometry/SUMMARY.md)
- 237 real 3-C seismometers in box 2019+ (207 sites, 158 alive 2025; 80 in roster). Sites <=45 km per Bostock family (median / % with >=3): S15 1 / 2%; +SYMB/BUCKS 1 / 19%; +PQ4 (JRDRB 2021, ALBHB 2022, LHLYB 2023, DAOB 2024; NRCan Public Safety Geoscience) 3 / 72%; all open 3-C 4 / 73%; Bostock's 2003-05 array 5 / 77%. PQ4 essentially restores his geometry. 16 deep Strait-of-Georgia families have no site <=45 km under any network.
- Puget: tremor-active cells median 3 sites <=45 km in 2025 (gap 100 deg); best band 47.0-47.5N median 5 (p90 10, gap 68). Oregon no better than the north except 45.5-46N.

### Step 3 — the pair-rate prediction (V1 network)

`s_pairlink3.py` (P3's `q_pairlink2` statistic on 2019–2026 hp candidates) + `s_step3_delta.py`.
2,834 hp events of KEPT families on S01–S13, 4.01 M pairs after the 27-s overlap exclusion,
94,138 same-family; index×1,800 s shift null. BUCKS drops out of the calibration dense set (it
starts 2023-09), so the dense set is PGCq+SYMB+VGZ.

δ fitted from the four tail quantiles (P3 §5.5): **old7 1.25 σ, dense 1.88 σ → ratio 1.50** →
**predicted δ_full = 0.38 × 1.50 = 0.57 σ** against the 1.5 σ licence bar (ratio ≥ 4).
**Step 6 NOT LICENSED on V1; the prediction is the result.**

Two ablation results worth keeping:
* **PGCq ALONE is the best configuration (δ = 2.44 σ).** Adding VGZ → 1.72, adding SYMB+VGZ → 1.88,
  the full seven-site sum → 1.25. P3 §5.4's dilution finding reproduced at single-site granularity
  on a different decade and a different candidate set.
* **SYMB contributes about what its registration says**: PGCq+VGZ 1.72 → 1.88 σ (+9 %), old7
  1.25 → 1.58 σ (+26 %). Real, and an order of magnitude short of the factor 4 needed.

Also recorded, from the Step-5 production runs: **at a FIXED MAD threshold the dense network detects
FEWER events than old13** (2024-08, T_core 9.25: 6,294 → 4,641 → 3,687 for old13 / dense /
dense-production). Adding stations raises the network sum's own MAD faster than the signal, which is
exactly why the comparison has to be at matched network-wide quiet-day FA and not at matched MAD.

### Step 4a — the CF timing kill test: the arm SURVIVES

`s_step4_cf.py` on 5,275 hp events of the four TEST episodes, every station through the
UNKNOWN-station embedding (`sid = nsta = 15`; verified `p2_cf.cf_day_sid` takes `sid` explicitly and
`sid_of` is consulted only under `--sid known`, so nothing is indexed by S15 position). Days read
from the 20 Hz caches via `s_common` (ZNE), not `p_common.DayStore` (S15-filtered). Only the window
grid covering each event's ±5 s is evaluated, by `cf_day_sid`'s own tiling rule, so the values equal
the full-day CF in the interior — a ~20× saving that made the test affordable.

Median |CF peak − Crosson S|: **SYMB 0.95 s, PGCq 1.05 s, VGZ 1.20 s, BUCKS 2.30 s**;
far sites 1.5–2.85 s. A CF with no timing information peaks uniformly in ±5 s → chance = **2.50 s**,
so B004 (2.85), ODUC (2.80), B001q/B003 (2.40) are AT CHANCE. **1 of the 4 ≤60 km sites exceeds the
1.5 s bar → arm NOT closed.** But 0.95–1.20 s is above Merlin's own 0.75 s "5 km plausible" line:
clause (c) (≤10 km) is reachable, clause (d) (≥3/10 isolated hidden within 5 km) is not, for the
truth-scarcity reason of flaw 2. The most interesting number: **SYMB, a station the CF has never
seen, has the best residual of all 16 sites** — the strongest evidence here that the CF generalises
across stations. 4b (back-projection) is licensed and NOT RUN this session.

### Step 7 — Arm B pre-registered (write-only)

`results/p3d/arm_b_preregistration.md`, sha256 in `freeze.json → arm_b`. Seeds are the top 1 % of
tremor-hour windows by the MEDIAN per-site z (not the sum — Step 3 showed the sum is PGCq alone);
matched filter over the pooled test episodes; two iterations; acceptance needs Gate 1 at ≥2 WITHHELD
sites AND the N-shift claim rule AND cross-episode re-firing; and a KILL TEST on the KEPT families
comes first. Its written expectation: Arm B will rediscover LARGE families and will not settle the
isolated-hidden score.

### Step 3 on S18 with the PQ stations (V2) — and why its "licence" is an ARTIFACT

PQ delivery (manifest, in-epoch SSE days): **ALBHB 431 data / 85 nodata, DAOB 217/4, LHLYB 157/72
but only 19 of the 59 days of its own pre-registered BUILD episode S18, JRDRB 87 data / 553 nodata —
JRDRB essentially does not deliver waveforms.** So the V2 decision run on S18 is V1 + ALBHB + DAOB
(LHLYB and JRDRB dropped by the availability rule at 0.7 of the episode days).

Run on S18, 1,384 candidates, 956,444 pairs, 24,019 same-family, IDENTICAL events in both runs:
* V1 config: δ_old7 **0.544**, δ_dense (PGCq+SYMB+VGZ+BUCKS) **1.557** → ratio 2.86 → 1.09 σ → not licensed
* V2 config: δ_old7 **0.304**, δ_dense (+ALBHB+DAOB) **1.611** → ratio 5.30 → 2.01 σ → "licensed"

**Do not read that as a licence.** The two runs use the SAME events and the SAME seven old sites; the
only difference in the `old7` column is which circular shift each station gets in the null (the shift
is index × 1,800 s and the index changes when ALBHB/DAOB are added to the station list). So
δ_old7 = 0.544 and 0.304 are **two draws of the same quantity from two null realisations**, a ±40 %
swing that moves the pre-registered ratio from 2.86 to 5.30 straight across the bar of 4 with no
change in the data. At this sample size old7's tail is populated by 102/18/3/0 pairs at
1e-3/1e-4/1e-5/1e-6 — two or three usable quantiles — so the denominator is not estimable.
Three further null realisations (shift 1,200 / 2,400 / 3,600 s) were run to measure that spread.

The number that IS stable, because it is measured on identical events with identical stations:
**δ_dense V1 1.557 → V2 1.611, i.e. ALBHB + DAOB add +3.5 %.** The well-powered measurement remains
the 13-episode calibration run (2,834 events, 4.0 M pairs, all four quantiles populated in both
configs): ratio **1.50 → predicted δ_full 0.57 σ → STEP 6 NOT LICENSED**.

### Step 5 — certified recall gain at matched FA: bar NOT MET

Three configurations run through a monkey-patched FROZEN `o_run_month` (`s_run_month.py`;
`o_run_month.py`/`o_common.py` byte-identical). Key discovery in the patch: `o_run_month` derives
its station list FROM THE TEMPLATE BANK, so handing it the transferred bank brings the new stations
into the network automatically — the patch only has to teach the supporting tables (station table,
epochs, the CO-LOCATED site merge, the N-matched weights, a P3D quiet-month amp reference, and an
optional un-drop of LRIV).

**At a FIXED MAD the dense network detects FEWER events** (2024-08, T = 9.25: 6,294 / 4,641 / 3,687
for old13 / dense / denseprod) — the production-statistic version of the Step-3 dilution result, and
the reason the FA-only rule is not optional. Matched FA on 61 quiet days (2024-12 + 2025-11):
old13 T 9.25 at **18.28 quiet-day false events/day** → T_dense **8.959**, T_denseprod **8.801**.

S18 at matched FA, 112 families: old13 8,872 vs dense 8,469 detections; 6,503 matched ±2 s,
**1,966 marginal**, **2,369 reverse-marginal**. a3b certification at the withheld site (B926+CLRS),
138/216 legs certifiable: **marginal purity median 0.361, reverse-marginal 0.306** — the dense-only
detections are no more real than the ones it loses.

**Certified gain −5.1 %** (bootstrap 2.5 % LB −8.5 %), **Spearman(gain, Δsites ≤45 km) = −0.013
(p 0.95)**, against the hashed prediction **+27.9 %**. All three clauses of the bar fail.
By geometry: Δ45 = 0 → −3.8 %, 1 → −6.1 %, 2 → −4.5 % — a uniform small LOSS, flat in geometry.
This is NOT Merlin's pre-registered failure signature (uniform GAIN / >2× predicted → site-noise →
suspect the tier derivation); it is the direct consequence of ρ_new < ρ_PGC.

**The transferable lesson: a distance-and-class model is not a substitute for measuring a site.**
The hashed prediction put SYMB 0.36 dB above PGC; it measured 2.73 dB below, and the per-station
site terms of the same fit span 6.5 dB.

### The null-realisation spread — a protocol correction

Five draws of the IDENTICAL pair statistic on IDENTICAL S18 events (1,384 candidates, 956,444
pairs), differing only in which circular shift each station gets in the null (index × SHIFT s):

| run | δ_old7 | δ_dense | ratio | verdict |
|---|---|---|---|---|
| V1, shift 1800 | 0.544 | 1.557 | 2.86 | not licensed |
| V2 (+ALBHB+DAOB), shift 1800 | 0.304 | 1.611 | **5.30** | "licensed" |
| V1, shift 1200 | 0.775 | 1.334 | 1.72 | not licensed |
| V1, shift 2400 | 0.452 | 1.351 | 2.99 | not licensed |
| V1, shift 3600 | 0.336 | 1.803 | **5.36** | "licensed" |

δ_old7 mean 0.482, range 0.304–0.775 (**±49 %**); δ_dense mean 1.531 (±15 %); the ratio spans
**1.72–5.36** across a bar of 4 and "licenses" in 2 of 5 draws with no change in the data. On S18 the
`old7` tail holds 102/18/3/0 same-family pairs at 1e-3/1e-4/1e-5/1e-6, so two or three of the four
fit points exist — the denominator is not estimable. **The δ-ratio licence rule as written needs a
null-realisation AVERAGE and a minimum populated-quantile count; that is a correction to the
protocol, not to the data.** The stable number is δ_dense 1.557 → 1.611 when ALBHB+DAOB are added:
**+3.5 %**, the same order as SYMB's +9 % in the calibration run.

Final: **Step 6 NOT licensed** (well-powered 13-episode calibration run: ratio 1.50 → 0.57 σ).

## 2026-09-19 05:35 UTC P3D RESULT (results/P3D_RESULTS.md)
- Registration rho vs PGC on the same events at matched N: SYMB 0.78, CLRS 0.77, GOBB 0.57, TURTL 0.56, RNWY 0.49, LRIV-HH 0.37, BUCKS 0.30, PTRF 0.28, OSQM 0.17, ODUC 0.16. SYMB's stack S-SNR is 2.7 dB BELOW PGC where the distance model predicted +0.4 dB -> site term (6.5 dB spread across stations). Proximity is necessary, not sufficient; PGC is an exceptionally quiet site.
- Step 3 prediction: delta_old7 1.25 sigma, delta_dense 1.88 -> predicted delta_full 0.57 sigma vs 1.5 bar -> Step 6 discovery NOT LICENSED. PGC-cluster alone is the best configuration; every other site dilutes.
- Certified recall gain from the added stations: -5.1% (LB -8.5%) vs predicted +27.9% -> bar NOT MET; uniform small loss.
- CF 4a survives (SYMB residual 0.95 s, PGC 1.05, VGZ 1.20); 4b licensed, not run.
- PQ: JRDRB delivers 87 days vs 553 nodata; LHLYB 19/59 of its build episode -> V2 = V1+ALBHB+DAOB, adds +3.5% to delta. V2 'licensed' flag is unstable across null realisations (ratios 1.7-5.4) -> NOT reported as a licence. Protocol fix needed: null-realisation average + min populated quantiles.
- Remaining: 4b, pooled Step 5, Arm B (pre-registered, not run).

## 2026-09-19 05:35 UTC P3D roster EXTENDED to 18 — four PQ + two UW Olympic broadbands
Coordinator priority, from `results/p3d_puget_geometry/station_ranking.csv` (PQ sites sit within
45 km of 36-66 Bostock families each). **Probed BEFORE queueing** (the GOBB lesson), then downloaded.

**NRCan is NOT unreachable — the recorded URL was wrong.** `https://earthquakescanada.nrcan.gc.ca/fdsnws`
answers **301**; the real host is `https://www.earthquakescanada.nrcan.gc.ca/fdsnws`, where
`station/1/version` and `dataselect/1/query` both answer 200. obspy still cannot be pointed at it
(`_discover_services` -> FDSNNoServiceException, `event/1` is 404), so NRCan must be queried over plain
HTTP — `r_probe_delivery.nrcan_query()` does that. Supersedes the "NRCan unreachable from this pod" note.
**But NRCan does not serve PQ at all** (station service returns an empty table for net=PQ; dataselect 204/404
on every probed day), and it returns **204 on the CN.GOBB gap days too** — so the GOBB holes are real
archive gaps, not IRIS-side routing. Everything here comes from IRIS.

- Probe bug caught and fixed: the first pass "tried NRCan" on 20 empty days and got `nodata/resumed` for all
  of them — `fetch_day` was short-circuiting on the `.nd` marker the IRIS attempt had just written, so NRCan
  was never actually asked. The fallback is now an independent HTTP request; the honest result is 14x http 404
  + 6x http 204, i.e. NRCan really has nothing.
- **The FDSN availability service is retired** (service.iris.edu and service.earthscope.org both return
  "This service has been retired"), so `scripts/r_probe_extent.py` maps the real extent by asking dataselect
  for 10 minutes on the 5th/15th/25th of every month of the epoch. Results now live in
  `station_epochs.csv` as `first_data_month,last_data_month,months_with_data,months_probed`.
- **StationXML epochs lie.** Declared vs MEASURED first data: ALBHB 2022-01 vs **2022-09**; LHLYB 2023-11 vs
  **2024-08**; JRDRB 2021-04 vs **2024-10** (and only 16/66 months hold anything); OHOH 2022-04 vs 2022-05;
  DAOB and ODEER are honest. The queue is still built from the DECLARED epoch so delivery is recorded
  empirically in the manifest.
- **TRAP - PQ VERTICALS ARE DECLARED POSITIVE-DOWN.** All four PQ stations report `HHZ az=0.0 dip=+90.0`
  (Guralp CMG-3ESP/Taurus), while every CN and UW station here reports `dip=-90.0`. Confirmed against raw IRIS
  StationXML text, not just our parser. PQ Z is therefore POLARITY-INVERTED relative to the rest of the
  network as declared. Anything polarity-sensitive — template transfer, cross-station stacking, back-projection,
  P-polarity — must flip PQ Z (or honour the dip) or the PQ sites will cancel the others. `chan_az_dip` in
  `station_epochs.csv` carries it; the frozen 20 Hz cache stores RAW traces and does NOT correct it.
- **P0 for the new six**: 2,375 queued -> **1,752 data / 623 nodata / 0 errors, 57.1 GB in 17 min**.
  All-station P0 now **9,136 station-days / 224.0 GB**, all 9,100 cacheable days at 20 Hz (0 errors).
- Concurrency: to stay at the <=8 concurrent IRIS cap while keeping the fill alive, P1/P2 was stopped
  (graceful, resumable) and both workers restarted at **4 + 4**.

**New-station P0 delivery (`results/p3d/download_qc_P0.csv`)**

| station | verdict | P0 days | data | nodata | frac | GB | measured first data | usable SSE episodes |
|---|---|---|---|---|---|---|---|---|
| PQ.ALBHB | YES | 532 | 447 | 85 | 0.8402 | 14.596 | 2022-09 | 10 |
| PQ.DAOB | YES | 294 | 290 | 4 | 0.9864 | 11.995 | 2024-02 | 5 |
| PQ.JRDRB | PARTIAL | 598 | 144 | 454 | 0.2408 | 6.653 | 2024-10 | 2 |
| PQ.LHLYB | YES | 302 | 230 | 72 | 0.7616 | 7.403 | 2024-08 | 4 |
| UW.OHOH | YES | 505 | 497 | 8 | 0.9842 | 13.779 | 2022-05 | 11 |
| UW.ODEER | YES | 144 | 144 | 0 | 1.0 | 2.707 | 2025-07 | 2 |

Whole-episode holes: **JRDRB is effectively a 2025+ station** — S09-S20 are ALL 0/N (twelve consecutive
episodes empty), only S21 (2025-08) and S22 (2026-04) are usable; LHLYB S17 0/32 and S18 19/59;
ALBHB S11 0/32 (a major ETS) and S12 7/57; OHOH S11 0/5 (only a 5-day sliver of S11 is inside its epoch);
DAOB S17 20/24; ODEER clean on both episodes it covers.

## 2026-09-19 05:39 UTC downloader final report (new six stations)
- All six deliver from IRIS; NRCan serves no PQ (real NRCan host is www.earthquakescanada.nrcan.gc.ca/fdsnws; GOBB gaps are real archive gaps). Measured extents (StationXML lies): ALBHB from 2022-09, LHLYB 2024-08, JRDRB effectively 2025+ (S21/S22 only), DAOB 2024-02, OHOH 2022-05, ODEER 2025-07. All-station P0: 9,136 station-days / 224 GB cached; P1/P2 fill ETA ~1.6 h.
- TRAP: PQ HHZ dip = +90 (positive DOWN) vs -90 everywhere else -> PQ Z polarity-inverted in the raw cache. P3D agent asked to redo the PQ legs with the dip honoured.
- FDSN availability service retired; r_probe_extent.py month-samples dataselect instead.

## 2026-09-19b — PQ component-dip correction (addendum; results/P3D_RESULTS.md § ADDENDUM)

All four PQ stations declare **HHZ dip +90 (positive DOWN)**; every CN/UW/PB station is −90. The
frozen cache stores raw traces, so PQ Z was inverted relative to the network.

**Audit — exactly ONE delivered leg used PQ data.** Step 2 (transfer/ρ/bank), Step 4a and Step 5 all
ran on the V1 roster (`--variant v1` default), so no PQ station is in any of them; only
`pairlink3_v2_s18_*` read ALBHB + DAOB. The PQ transfer stacks, PQ ρ and PQ CF residuals the
correction asks about did not exist and were built here for the first time, with the dip honoured.

**General fix**: `s_common` parses `chan_az_dip` for every station and applies
`ZSIGN = −1 if dip_Z > 0 else +1` in `load_day_zne` → only {ALBHB, DAOB, JRDRB, LHLYB} flip, the
other 29 stations are bit-identical. `S_P3D_NO_DIP=1` reproduces the old behaviour through the same
code path, so before/after is measured, not argued. Control: SYMB 4a re-run gives 0.95 s, identical.

**Before/after**
* **Pair statistic / δ_dense V2: NO-OP.** Flipping one component of *both* windows cancels in the
  3-C inner product (measured: cc −0.085813289610 both ways, 12 digits). Full V2 S18 re-run:
  δ_old7 0.304, δ_dense 1.611, ratio 5.30 — identical to the raw-Z run. The dip would have mattered
  for anything summing Z coherently ACROSS stations (Step 4b back-projection, not run).
* **PQ registration ρ (new)**: **LHLYB ρ/ρ_PGC = 1.086** (n=8 families) — the first new station to
  match PGC — and **DAOB 0.826** (n=33), both above SYMB's 0.783 and far above BUCKS's 0.296.
  ALBHB untestable (measured extent starts 2022-09 → 129/130 families insufficient on S12–S13).
  ρ is invariant to the flip (0.302307 vs 0.302294; 0.275114 vs 0.275092).
* **PQ CF 4a**: LHLYB 1.10 s (raw 1.05), DAOB 1.20 (1.20), ALBHB 1.40 (1.25), JRDRB 1.65 (1.65) —
  0–3 samples of the 20 Hz grid; the CF is near polarity-blind for timing. Under the V2 ≤60 km rule
  2 of 8 sites exceed 1.5 s → **CF arm still SURVIVES, with more margin than on V1**.

**Measured extents now drive the split**: ALBHB 2022-09 → calibrate S12–S13 (not S08–S13);
JRDRB 2024-10 → build S21, test S22 (no calibration episode at all); LHLYB 2024-08 → S18 build is
only 19 of 59 days; DAOB 2024-02 unchanged. Hashed as `freeze.json → amendment_v2b` (keeps
`amendment_v2` verbatim; adds the z_sign map, the extents and the corrected StationXML coordinates).

**Licence statistic with the protocol fix** (`s_step3_licence.py`: ≥3 populated quantiles with ≥3
pairs; mean over null realisations AND its bootstrap LB; **≥3 admissible realisations** before the
rule may be adjudicated — without that clause the bootstrap over one realisation returns itself):

| set | admissible | ratio | predicted δ_full | verdict |
|---|---|---|---|---|
| **V2 (+ALBHB+DAOB), dip-fixed, S18** | **1 of 4** | — | — | **UNADJUDICABLE** |
| V1, S18 | 3 of 4 | 2.40 (LB 1.72) | 0.91 σ (LB 0.65) | NOT LICENSED |
| **V1, calibration S01–S13** | **3 of 3** | **1.52** (LB 1.50) | **0.579 σ** (LB 0.569) | **NOT LICENSED** |

Two extra calibration null realisations (shift 1,200 / 2,400 s) closed clause (iii) and gave the
cleanest contrast in the whole study: δ_old7 across realisations is **1.253 / 1.242 / 1.189 on the
calibration set (±2.6 %)** against **0.775 / 0.544 / 0.452 / 0.336 / 0.304 on S18 (±49 %)**, while
δ_dense is stable in both (calibration 1.863–1.876, ±0.3 %). It is the DENOMINATOR that decides
whether the rule can be adjudicated at all, and it is estimable only at the calibration set's size
(2,834 events / 4.01 M pairs), not at S18's (1,384 / 0.96 M).

The single "LICENSED" draw of § 8.3 is exactly the artifact the fixed rule rejects. **Step 6 remains
NOT licensed** on the only well-powered measurement. A real V2 decision needs `old7` to keep ≥3
populated quantiles in ≥3 realisations — 2,834 candidates does it, 1,384 does not.

**Two bookkeeping faults found and fixed (recorded, not quietly repaired):** `s_step2_transfer`
wrote its bank to a fixed directory regardless of `--tag`, so the PQ run clobbered the V1 bank —
**Step 5 unaffected**, all ten bank loads in `logs/s5b_*.log` record `+251 templates`; path now
honours `--tag`, both banks regenerated. And the first PQ CF attempt crashed on an empty frame and
its wrapper's `mv` renamed the V1 4a artifacts; regenerated, and `s_step4_cf.py` gained
`--variant` / `--out-tag` and an empty-frame guard.

## 2026-09-19 07:27 UTC P3D addendum (PQ dip-corrected; results/P3D_RESULTS.md A1-A6)
- Only the V2 pair-link leg had used PQ data; PQ transfer/rho/4a built for the first time with the dip honoured (ZSIGN from chan_az_dip).
- PQ registration rho vs PGC: LHLYB 1.09 (n=8 families) — first new station to match PGC; DAOB 0.83; ALBHB untestable (data from 2022-09); JRDRB no calibration episode (2024-10+). Pair statistic invariant to the flip (same-station inner product).
- Licence rule fixed (>=3 populated quantiles, mean over >=3 admissible null realisations): V2 on S18 UNADJUDICABLE (1/4 admissible); V1 S18 ratio 2.40 -> 0.91 sigma NOT licensed; V1 calibration 0.57 sigma NOT licensed. A real V2 decision needs pooled S18+S21+S22 (spends the discovery look) or 2026 data.
- CF 4a survives on V2 too (PQ residuals 1.1-1.65 s).

## 2026-09-19 07:43 UTC USER DECISION: spend the discovery look
Run Step 6 (template-free pair search) AND pre-registered Arm B (seed-and-grow) on the pooled PQ-era episodes S18+S21+S22 with the V2 network, dip honoured; truth = hp events of the 30 hidden families; nulls N-shift/N-quiet/N-moveout; gates 0-3; relative scorer; CF trigger as an ablation arm. One look; outcome closes or opens northern discovery.

## 2026-09-19c — discovery look: PRE-REGISTERED and HELD (nothing run)

Coordinator authorised spending the discovery look on pooled S18+S21+S22 with the V2 network.
`results/p3d/discovery_freeze.json` is written and hashed
(`ce8b50f026a22fbc812d51a0949dacfdb56f7febb0f18b5bf09bebffe7790da3`) **before any pass**, with
status `AWAITING_COORDINATOR_CONFIRMATION`. Sizing done; two blocking findings stop the run.

**Finding 1 — the full V2 network empties S18.** Requiring all 8 V2 ≤60 km sites cached leaves
**112 of 203 pooled days and ZERO of S18's 59**: LHLYB delivers 18/59 S18 days (2024-07 is genuinely
`nodata`) and JRDRB 0/59. The instruction asks for a pooled S18+S21+S22 run *and* a per-episode
readout for S18; with the full network the S18 readout cannot exist. Four branches, all costed:

| branch | sites | days | episodes | windows | pairs |
|---|---|---|---|---|---|
| A all 8 | 8 | 112 | S21 50 + S22 62, **S18 absent** | 228,600 | 2.6e10 |
| **B drop JRDRB (recommended)** | **7** | **129** | **S18 17 + S21 50 + S22 62** | **381,000** | **7.3e10** |
| C drop JRDRB+VGZ | 6 | 161 | S18 17 + S21 82 + S22 62 | 761,400 | 2.9e11 |
| D strict 3 sites | 3 | 203 | S18 59 + S21 82 + S22 62 | 1,003,800 | 5.0e11 |

B keeps an S18 readout, keeps LHLYB (44.77 km median hypocentral, 66 families ≤45 km — ties SYMB as
the best-placed site in the roster), drops only the station already documented as non-delivering,
and at 7.3e10 pairs is ~2× Merlin's own 3.4e10 sizing and tractable. **D is not absurd**: Step 3
found PGCq ALONE is the best pair configuration, so a 3-site sum over all 203 days may beat an
8-site sum over 112.

**Finding 2 — the scoreable isolated families are largely the ones PQ does NOT help.** The bar is
2/3 of the scoreable isolated stratum = bkeys 8, 10, 15, 87, 115, 119, whose new-PQ-sites-≤45 km are
**0, 1, 2, 1, 0, 0**. The three isolated families that gain **5** new PQ sites each (12, 29, 62)
carry **8, 8 and 4** T1 events and are all unscoreable. Medians: scoreable 44.5 T1 / 0.5 new PQ
sites; unscoreable 6.0 T1 / 5.0 new PQ sites. (Spearman(n_T1, d45_V2) = −0.17, p = 0.63, n = 10 —
**not significant**; the structural fact is what is reported, not a correlation.) This is Merlin's
flaw 2, sharpened by the PQ amendment: the look would test V2 mostly on families V2 does not improve.

Also recorded: the PNSN catalog ends 2026-05-27 23:52, so **S22 contributes tremor hours from 57 of
its 62 padded days**; the last 5 days are cached but carry no tremor-hour label and yield no
candidates.

Not run, because choosing a branch after seeing these numbers is the forking path Q4 warns against,
and the look is unrecoverable once spent.

## 2026-09-19d — discovery look CONFIRMED and under way

Coordinator settled the fork: **PRIMARY = NET-B** (7 sites, drop JRDRB; S18 17 + S21 50 + S22 62 d,
381,000 windows, 7.3e10 pairs). **SECONDARY = NET-D**, defined by the stated ρ rule as
PGC cluster + the two best-registering others by Step-2 ρ = **LHLYB (1.086) and DAOB (0.826)**.
Deviation recorded in the freeze: the instruction said "over all 203 days", but LHLYB is cached on
162 of 203, so **NET-D realises 162 days** (S18 18 + S21 82 + S22 62), 762,600 windows, 2.9e11 pairs.
The two 203-day alternatives (PGCq+BUCKS+DAOB, PGCq+SYMB+DAOB) would violate the ρ rule by
substituting a worse-registering site, so the ρ rule is honoured and the shortfall stated.
Bars, secondary readouts (identity-only on 12/29/62; Δsites-weighted fraction) and the structural
fact are all in `discovery_freeze.json`, re-hashed `31da03e3…` **before any pass**.

**Engine.** `s_discover_pairs.py` = the P3 pair statistic on uniform tremor-hour windows, with the
per-site accumulator **FUSED** into one tensor (the old code held one [n_i, n_j·121] fp32 array per
site and summed them — 8× the memory traffic for no arithmetic reason). Validated bit-exactly
against the per-site path: max |difference| **0.000e+00**, argmax lag identical on 100 % of pairs.
`s_discover_chain.py` carries edges → kNN → components → per-site align_stack → S excess → 12-s
templates.

**Throughput probes** (NET-B, same 129 days, candidates subsampled):

| probe | windows | pairs | runtime |
|---|---|---|---|
| 30 k | 29,896 | 4.47e8 | 162 s |
| 120 k | 119,587 | 7.15e9 | (measuring) |

The 30 k probe averages only 232 windows/day, so its GEMMs are ~160× smaller in area than the full
run's and its rate is a floor, not a projection. Threshold rule behaves as designed: at 4.47e8 pairs
it returned the 6.7e-8 quantile for ≤30 expected false edges.

### Engine fix and the real throughput

First full NET-B attempt died with **CUDA OOM**: some tremor days carry ~11 k windows, and the j-side
unfold costs n_j·121·1440 fp16 *per station* (~3.9 GB each, ~39 GB for ten) plus an 11 GB fp32 output
tile. Fixed by **j-blocking** (a day is split into sub-blocks of ≤3,000 windows) and a smaller i-chunk.
Regression-checked on the 30 k probe — pre- vs post-blocking gives **identical** n_pairs
(446,870,460), **identical** threshold (29.682254791259766) and identical quantile, and runs
**2.4× faster** (162 s → 67 s): the blocking is exact, not an approximation.

Measured full-scale rate on NET-B: **2.8e7 pairs/s** (blocks 21→31: 2.16e9 pairs in 77 s), so the
7.26e10-pair pass costs **~45 min**, not the 4 h the 30 k probe's tiny GEMMs implied. NET-D's
2.9e11 pairs would be ~3 h per pass — above the coordinator's ≤2 GPU-h/pass condition, to be
decided once ARM A is scored.

Arm inputs built: **ARM B** 3,809 seeds = top 1 % of the 380,943 scored windows by the MEDIAN
per-site S-energy z (median, not sum — Step 3 showed the sum is PGCq alone); **ARM C** full-day
frozen CF at the seven ≤60 km sites, top 720 peaks/station-day, declustered at 6 s and merged
across sites at 3 s.

### Discovery run state (2026-09-19d, in flight)

Queued and running unattended, all products landing in `results/p3d/disc_*`:

| stage | state |
|---|---|
| `discovery_freeze.json` | **written + hashed before any pass** (`31da03e3…`) |
| fused accumulator | **validated bit-identical** (max diff 0.000e+00, 100 % argmax agreement) |
| j-blocking | **regression-clean** (same pairs, same threshold, 2.4× faster) |
| ARM-B seeds | done — 3,809 of 380,943 windows (top 1 % by MEDIAN per-site S-energy z) |
| ARM-C CF candidates | done — 282,019 peaks → **145,661** merged candidates |
| ARM A NET-B null | running, ~2.8e7 pairs/s, ~45 min/pass |
| ARM A NET-B true + chain | queued behind it (`/tmp/.../armA.sh`) |
| ARM B, ARM C null+true+chain | queued (`/tmp/.../armBC.sh`) |
| scorer | `s_discover_score.py` written (Gate 1 at the withheld sites via the a3b statistic, identity vs the 2012 PGC templates, PRIMARY + both SECONDARY readouts, all clauses read from the freeze) |
| figures | `s_discovery_figs.py` written; the map renders now, the space-time panel fills as arms land |
| N-shift claim-rule nulls (5 seeds) | **deliberately deferred** — they are needed only for gate-PASS clusters matching no Bostock family, so they are run conditionally on ARM A producing one, rather than spending 5 × 45 min unconditionally |
| NET-D | decision pending ARM A: 2.9e11 pairs ≈ 3 h/pass, above the coordinator's ≤2 GPU-h/pass condition |

## 2026-09-19e — THE DISCOVERY LOOK: all three arms, zero clusters, and why the bar was unreachable

Full write-up `results/P3D_DISCOVERY_RESULTS.md`; pre-registration `discovery_freeze.json`
(`31da03e3…`, hashed before any pass); ledger row `P3D-DISC`.

| | ARM A pairwise | ARM B seed-and-grow | ARM C CF-triggered |
|---|---|---|---|
| pairs | 7.21e10 | 9.15e8 | 1.05e10 |
| threshold | 64.795 (q 4.16e-10) | 17.084 (q 1e-4) | 36.048 (q 2.86e-9) |
| edges | 455 | 327,763 | 873 |
| components ≥20 | 2 | 3 | 2 |
| **pass S excess ≥3 dB at ≥3 sites** | **0/2** | **0/3** | **0/2** |
| PRIMARY scoreable isolated | **0/6** | **0/6** | **0/6** |

Closest approach in the whole study: ARM C cluster 10640 (62 members) reaches **+3.27 dB at ALBHB
and +3.61 dB at LHLYB** — two sites, both PQ, against a clause needing three. ARM C not retained
(0 vs 0, needed ≥2 more than A). N-shift claim-rule nulls **not triggered** — no gate-PASS cluster
existed, so the 5 × 45 min was correctly not spent.

**Kill test: the construction is STARVED, not BROKEN.** ARM B links kept families at
**11.7 × chance (5 same-family of 15 labelled-labelled edges, p = 4.4e-5)** — same sign and order as
P3's 8–38 × on catalogue candidates. The physics reproduces template-free, on a different decade and
a different network; there is simply ~10²–10³ too little of it.

**THE finding — the PRIMARY bar was unreachable before any detector ran.** Only **23 %** of hp events
lie within ±2 s of a candidate window, and the shortfall is not uniform: **bkeys 8, 10 and 15 have
ZERO representable truth** because their events live in S18 and NET-B keeps only **17 of S18's 59
days**. The bar needs 4 of 6; at most **3** were reachable. Across the branches:

| configuration | days | S18 | families reachable |
|---|---|---|---|
| **NET-B (chosen)** | 129 | 17 | **3/6** |
| NET-D (ρ rule) | 162 | 18 | 5/6 |
| max-coverage 3 (PGCq+BUCKS+DAOB) | 201 | 57 | **6/6** |

The branch chosen on geometry and ρ was the worst of the three on the axis that decided the test,
and the branch rejected for violating the ρ rule was the only one with full coverage. Not knowable
at decision time — it needed the realised candidate grid. **So the correct reading is not "the dense
PQ network cannot discover families" but "this configuration could not test that."**

Two empty-case bugs found and fixed while scoring (the scorer and the figure loader both assumed a
non-empty cluster table); a chain that yields nothing is a result, not an error, and now keeps its
schema. Also fixed a misleading console line that printed "30 expected false edges" even under the
fixed-quantile rule — the JSON was always right.

## 2026-09-19 09:50 UTC DISCOVERY LOOK SPENT (results/P3D_DISCOVERY_RESULTS.md)
- Arms A (7.2e10 pairs), B (seed-and-grow), C (CF trigger): 0 clusters pass S excess >=3 dB at >=3 sites; PRIMARY 0/6; secondaries 0. Closest: ARM C cluster 10640 (+3.3 dB ALBHB, +3.6 dB LHLYB — two PQ sites). Kill test: Arm B links KEPT families at 11.7x chance (p 4e-5) template-free -> physics present, ~1e2-1e3 short.
- DESIGN FLAW (orchestrator's branch choice): NET-B required LHLYB+ALBHB cached -> only 17/59 S18 days -> bkeys 8, 10, 15 had ZERO representable truth; at most 3/6 scoreable families reachable vs a 4/6 bar. Max-coverage 3-site config (PGCq+BUCKS+DAOB, 201 days) would reach 6/6; NET-D 5/6. Honest reading: this configuration could not test the dense-network question. Ledger row P3D-DISC: look spent.
- Option: post-hoc EXPLORATORY run of the max-coverage 3-site configuration (~3 h/pass), gates + nulls as judge, labelled exploratory; or stop.

## 2026-09-19f — EXPLORATORY re-run on NET-M (max-coverage) and NET-D

`discovery_freeze.json → amendment_exploratory_2026-09-19` (sha256 `c9034088…`), written before any
pass. **Ledger P3D-DISC stays SPENT/CLOSED — this is post-hoc and licenses nothing.** Full write-up
in `results/P3D_DISCOVERY_RESULTS.md` § EXPLORATORY.

**Reachable truth (the reason for the re-run), stated before results:**

| config | sites | days | PRIMARY families reachable | isolated truth representable |
|---|---|---|---|---|
| NET-B (pre-registered) | 7 | 129 | **3/6** | 28 % |
| NET-D (ρ rule) | 3 | 162 | 5/6 | 42 % |
| **NET-M (max-coverage)** | 3 | 201 | **6/6** | **95 %** |

**Results: rediscovery 0 in every stratum, configuration and arm.** NET-M ARM A 4.86e11 pairs →
1,469 edges → 5 components → 0 pass S-excess; ARM C 310 edges → 0 pass; **ARM B percolated** at the
pre-registered 1e-4 cut (one 404,165-window component = 41 % of candidates) — a result, not a bug:
a 3-site sum has a heavier tail and the quantile is too permissive. A z≥16 sensitivity
de-percolates into 9 components, still 0 pass. NET-D 2.90e11 pairs → 2,167 edges → 7 components →
**1 passes the S-excess clause**.

**That one cluster (NET-D 51599) is an artefact.** 124 windows but only **33 independent time
episodes** (91/123 gaps ≤27 s at the 6-s hop — adjacency joins transitively through a shared distant
partner, which the 27-s *pair* exclusion does not prevent). On those 33 episodes, alignment-free:
**Gate 1 fails at all 5 withheld sites** (0.48–1.97 dB vs null p99 1.98–2.92), Gate 0 identity top
|cc| **0.231** vs the 0.524 chance level, and the PGC stack shows **S-SNR 1.70 dB / excess 0.59 dB**
— no S arrival at a site inside the detection network. Its chain excess came from `align_stack` on
adjacency-inflated members. **No gate-PASS cluster anywhere → the 5-seed N-shift claim-rule nulls
were not triggered.**

**The real gain: the kill test now PASSES.** On NET-M ARM B the construction links KEPT families
**196 times of 476 both-labelled edges (41.2 %) = 14.5 × chance** — against 5 edges on NET-B. Truth
representability 28 % → 95 % converted a starved kill test into a working one. **So linkage physics
and truth coverage are both confirmed; what fails is the step from linked pairs to a clustered,
gated family.** That relocates the problem from "can we see them" to "can we assemble them".

Structural note: the S-excess clause needs ≥+3 dB at **≥3 sites**, and NET-M/NET-D have exactly 3 —
so it demands unanimity there. Several clusters cleared +3 dB at one or two sites (LHLYB 9.77, DAOB
8.73, DAOB 6.87/LHLYB 3.83/PGCq 3.60) and died on the third.

Two engine faults fixed, both of which turned a result into a crash: `components()` died on
percolation (giant component above the size cap) and the chain wrote a column-less CSV when nothing
survived; both now emit typed empty products and report the component-size distribution.

## 2026-09-20 23:59 UTC EXPLORATORY DISCOVERY RE-RUN (results/P3D_DISCOVERY_RESULTS.md, EXPLORATORY section)
- Reachable truth: NET-B 3/6 (28%), NET-D 5/6 (42%), NET-M (PGCq+BUCKS+DAOB, 201 d) 6/6 (95%).
- Results: NET-M A 5 comps >=20, 0 pass S-excess; NET-M B PERCOLATED at the 1e-4 cut (one 404k-window component; 3-site sum tail heavier); NET-M C 0; NET-D A 1 passing cluster (51599) = ARTEFACT (124 windows / 33 independent episodes via transitive adjacency; Gate 1 0/5 withheld sites; Gate 0 top |cc| 0.23; PGC stack no S). Rediscovery 0 in all strata, all configs. No gate-PASS cluster -> N-shift nulls not triggered.
- Kill test (the gain): NET-M Arm B links KEPT families 196 times at 14.5x chance (476 labelled edges, 41% same-family) -> linkage physics + truth coverage confirmed; failure is at linked pairs -> clustered, gated family (percolation; unanimity of 3 sites at +3 dB).
- P3D-DISC stays spent/closed; re-run licenses nothing. Agent stopped.

## 2026-09-20 03:28 UTC user: next = Puget/Olympic. P4 section added; Merlin review + P4-0 site-noise inventory launched (agents instructed to launch-and-exit).

## 2026-09-20 — P4 launched; Merlin review in; P4-0b (north kill test) is the gate
- **Merlin P4 review** saved verbatim → `results/merlin_review_p4_2026-09-20.md`; amendments appended to PLAN.md §P4. Headline: P3D never ran the true seed-and-grow (pairwise graph only; seeds were earthquakes/bursts) → the Puget discovery premise is unproven; run the real construction on the north first.
- **P4-0 probe queue running** (`scripts/u_p4_0_queue.sh`, pid 1027582, log `logs/p4_probe_queue.log`, markers `results/p4/P4_0_DONE|FAILED`): 139 candidate 3-C stations (≤60 km hypocentral of 463 active cells 47.0–48.3 N, alive 2023+), 4 tremor + 4 matched quiet days (2023-04-14, 2024-08-21, 2025-05-16, 2025-10-04), 969 station-days ≈ 24 GB; 22 stations at 50/125 Hz go through a rational-resample twin (`chain=rational`, never mixed with frozen-chain products); 4 stations carry `zsign=-1` (dip +90). Outputs: `results/p4/{candidates,station_rank_all,shortlist_v0,noise_inventory}.csv`. Smoke test: B011 −7.5 dB vs PGC with 14.2 dB tremor contrast; MORSE +4 dB but 2.9 dB contrast (noisy, blind).
- **Follow-on for the P4-0 analysis pass (per Merlin F5):** calibrate the proxy on the 14 north stations with measured ρ (LHLYB 1.09 … ODUC 0.16) before the Puget shortlist is used; report Spearman + LOO error; fall back to geometry+contrast if Spearman <0.4.
- **P4-0b agent launched** (scripts `v_*`, queue `scripts/v_p4_0b_queue.sh`, markers `results/p4/P4_0B_DONE|FAILED`): true Arm B on the PQ-era north, given-seed + template-free arms, per Merlin §P4-0b and P4-2 steps 2–6. Decides whether P4-2 (Puget discovery) runs at all.
- Governance unchanged: 2013 sealed; 30 hidden families untouched; Puget TEST 2025+ sealed.

## 2026-09-20 — P4-0a ANALYSIS: the proxy CALIBRATES (LOO ρ_s +0.69) but all three sub-regions STOP
Full write-up `results/P4_0_RESULTS.md`; scripts `u_north.py`, `u_calibrate.py`, `u_fartremor.py`,
`u_hf_north.py`, `u_shortlist_v1.py`.
- **Calibration (Merlin F5) — ACCEPT.** 15 north stations × 4 tremor + 4 matched quiet days inside
  S18/S21 (≥100 epicentres ≤40 km of the Bostock centroid; quiet = none within 150 km ±1 d), run
  through the identical chain. Prescribed model ρ/ρ_PGC ~ dist + quiet-night noise + contrast:
  **LOO Spearman +0.688 (p 0.007, n 14), signs all as expected, LOO median |err| 0.166** (bars 0.6 /
  0.25). Held-out VGZ (derived ρ 0.59) predicted 0.70.
- **…but it is mostly geometry.** `dist only` gives LOO +0.732; the contrast coefficient is ≈0.
  Contrast IS a good univariate proxy — **+0.754 in the native 2–8 Hz band vs +0.556 in the 1–8 Hz
  cache band** — but is collinear with distance. Say "calibrated", not "the site measurement works".
- **Far-tremor control had to be rebuilt.** Reusing north tremor days with no Puget tremor is
  geometrically IMPOSSIBLE (the family field is 25 km north of the band top; zero such days
  2023→2026). Substituted: 4 days with ≥100 epicentres 120–350 km away and none within 120 km, both
  azimuths (VI 170 km ×2, Oregon 212/285 km). **The ≤+1 dB bar is a factor 3 inside its own noise**
  — the quiet-day-to-quiet-day scatter of the same statistic is median 3.4 dB — so both the
  specified bar and an intent bar (local − far ≥ 3 dB) are carried. The intent bar is what catches
  the real failures (B003 +21 dB on a glitch day, LOPEZ, OSQM).
- **10–19 Hz control:** native rate, 106/133 have HF contrast < 2–8 Hz contrast.
- **STOP all three sub-regions**, under both far-bar readings. Binding clause = "≥4 admitted sites
  ≤45 km hypocentral of ≥50 % of cells": achieved fraction **0.00** (0.00–0.07 even at ≥3 sites).
  Admitted sites R1 3 / R3 3 / R2 3 (specified) or 5 / 6 / 4 (intent); the 1-km merge collapses
  B005+B006+B007+LRIV to ONE site.
- **The geometry is fine; the alignment isn't.** All 118 candidate sites give 0.55 / 0.34 / 0.60 —
  the clause needs ~30 sites per sub-region. The geometrically decisive sites are noisy UW lowland
  broadbands that fail the contrast bar (WINDI 62 cells 2/4 days, BUCKS 71 cells −4.8 dB, ANGUS,
  EAGLE, DOORS); the sites that pass are quiet boreholes 45–70 km away.
- **Cheapest thing that could change it:** a 6+6-day re-probe of the ~25 geometrically decisive
  sites (~3 GiB, ~10 min) — "≥3 of 4" is a 75 % clause on a noisy per-day statistic where Merlin
  specified ≥4 of 6.
- **Both flagged judgement calls are moot:** every admitted station passes the PSD power bar
  (≤ +4.8 dB), and no `chain=rational` station is admitted.
- **Bugs fixed:** `na_values=[""]` parsed location code `"00"` as `0.0` → response lookup missed →
  11 stations' far statistics were in raw counts (165 dB "contrasts"); delivery was divided by 8
  instead of by the days in the station's own epoch (89 → 117 of 134 pass).

## 2026-09-20 — P4-0b BUILT AND LAUNCHED: the TRUE Arm B kill test on the north

Merlin's P4 review (F1/F2) showed P3D's "seed-and-grow / ARM B" was only a seed-restricted
PAIRWISE graph (`s_discover_pairs.py:186-190,227`; `s_discover_chain.py:27-46`) whose seeds were
earthquakes and bursts, and that the construction pre-registered in
`results/p3d/arm_b_preregistration.md` steps 2-5 has **never been run**. P4-0b runs it, on the
PQ-era north, where the 130 Bostock families minus the 30 hidden ones are truth. Its result gates
whether Puget discovery (P4-2) happens at all.

**Built (all new, prefix `v_`, nothing existing edited, no frozen json touched):**

| script | what |
|---|---|
| `v_common.py` | sites/roster, per-day availability, eligible families, quiet days, region centre, robust scales, one-to-one matcher, bank IO, freeze helpers |
| `v_engine.py` | the GPU matched filter: frozen-runner cc (rFFT, float64 prefix-sum denominator, 3-C mean) -> per-station-day **lowhalf-MAD z** -> per-site mean **clipped at 8** -> **equal-weight sum over observable sites** -> per-(template,day) med/MAD `net_mad`; hit = `net_mad >= T` AND >= 2 sites at z>=3 AND >= 3 sites observable; decluster 30 s, independent 60 s |
| `v_seeds.py` | GIVEN seeds (3 hp events/family by frozen `net_mad`), FREE seeds (Merlin's corrected rule), and the **seed-z band TRUE hp events occupy** |
| `v_bank.py` | pass-0 banks (20 s from S-6 s), the per-template **N-moveout** wrong-cell banks, and the arm-(i) truth-stack reference bank |
| `v_run.py` | quiet-day FA -> `T_seed` by the frozen FA-only rule -> growth -> survival clauses |
| `v_iterate.py` | `align_stack` x2 (+-1.5 s) -> re-pick S at the highest-S-SNR site -> re-cut 12 s from S-2 s; and the cheap alignment-free S-excess **screen** (reported, never acceptance) |
| `v_dedupe.py` | >=50 % hit overlap on a common (predicted-S-at-PGC) axis, or stack cc above the recomputed per-station cross-template p95 at >= half the registering sites |
| `v_score.py` | arm (i) relative clauses of `q_rediscover.py`; arm (ii) one-to-one +-2 s with ONE nuisance offset charged to a 20-realisation null; burst convergence; verdict |
| `v_freeze.py`, `v_p4_0b_queue.sh` | the freeze writer and the single unattended queue |

**Frozen before launch** (`results/p4/p4_0b_freeze.json`, 10 `v_` script sha256s):
sites PGCq(PGC+B009/10/11), VGZ, SYMB, BUCKS, LHLYB, ALBHB, DAOB; withheld B926/NLLB/LRIV (+ the
four P3D withheld sites) never in any sum; **203/203 pooled S18+S21+S22 days usable** with per-day
availability (5 sites 21 d, 6 sites 53 d, 7 sites 129 d) — days are **not** intersected across
stations, which is what starved NET-B; **20 quiet days** (2024-09-06..2026-06-24, all 7 sites);
**18 eligible kept families** (>=50 hp events in the pooled episodes) -> bar (i) is **>=12/18**.

**Three documented deviations, all in the freeze:**
1. *Region-centre depth.* Slab2 at (48.50, -123.75) returns **21.9 km**, ~13 km above the family
   field; the adopted depth is the **median Bostock family depth 35.3 km** (lat/lon unchanged).
   Using Slab2 would mis-cut every template by ~3 s.
2. *Network comparability.* Merlin's equal-weight site SUM is converted to `net_mad` by the frozen
   runner's per-(template, day) median/MAD step — otherwise one fixed threshold is not comparable
   across days with 5, 6 or 7 sites. Note the site clip at 8 caps the sum at 8 x n_sites, so
   `net_mad` saturates near 20 on a 6-site day.
3. *Arm-(i) reference.* The P3D transfer bank has only 26-184 (family, station) pairs and none at
   PGC/VGZ/ALBHB, so "the family's own transferred template" cannot be built at 7 sites. The
   reference is the family's own **truth-stack** 12-s template (all hp events, `align_stack`, S-2 s)
   grown through the identical machinery — at least as strong, so the relative bar is conservative.
   Absolute recall/purity are reported alongside.
   Arm (ii) matching uses **one** constant offset per (template, family): a region-moveout seed's
   origin axis is displaced from the family's by the single constant `b = t_seed - tau_0` (same for
   every event); the offset search is charged to a 20-realisation within-day null and the strict
   `b = 0` numbers are reported too.

**Smoke test (3 given seeds, then 54, end-to-end on real days):** the engine self-detects each seed
(`net_mad` 19.8-20.5) and, over 5 episode days at the T=8.0 floor, single-event pass-0 templates
reach recall 0.06-0.15 at purity 0.62-1.00; one iteration of G44_0 lifted recall **0.368 -> 0.474**.
Quiet-day FA at `net_mad >= 8` is **0** for every template tested, so `T_seed` = 8.00. Truth-stacks
are clean (family 13: S-SNR 10.5/10.9/9.2 dB at PGC/SYMB/DAOB). Marginal cost measured at
**~810 template-station-days/s** (day setup ~6 s dominates), so the full queue is ~6-9 GPU-h, well
inside Merlin's 18-24 h allowance.

**Launched:** `nohup bash scripts/v_p4_0b_queue.sh > logs/p4_0b_queue.log 2>&1 &` — **pid 1040534**,
markers `results/p4/P4_0B_DONE` / `P4_0B_FAILED`, products under `results/p4/p4_0b/`.
Queue order is Merlin's: GIVEN arm (bank -> quiet T -> grow -> N-moveout -> survival -> iterate x2
-> dedupe -> screen) and its reference **first**, a partial score written as soon as the decisive arm
lands (`score_partial.txt`), then the TEMPLATE-FREE arm, then the single scoring pass.
Check with `tail -30 logs/p4_0b_queue.log`; every step is guarded by its output file, so re-running
the queue resumes. Governance unchanged: 2011/2013 sealed, 30 hidden families never seeded/scored/
matched, no Puget byte read or written, no FDSN request.

## 2026-09-20 — P4-0a re-probe at Merlin's ORIGINAL clause (6+6, ≥4 of 6): STOP still stands
Addendum in `results/P4_0_RESULTS.md`; script `u_reprobe.py`; pre-registration (hashes + day list +
clause, written BEFORE any fetch) `results/p4/reprobe_prereg.json`.
- Added tremor 2023-04-11 (415 epicentres ≤40 km) and 2024-08-17 (671) + matched quiet 2023-03-19,
  2024-08-03, chosen on epicentre count alone; 101 stations / 84 sites (≥10 cells ≤45 km, or 1–2 of
  4); 323 of 338 station-days delivered, 10 GiB, 2.5 min.
- **The 4-day "≥3 of 4" reading WAS too strict**: contrast-bar pass count 31 → **59** of 134.
- **The absolute far bar is confirmed to be noise.** B005 1.58, B006 1.42, LRIV 1.16, SKOKO 1.75,
  HDW 1.12, HERD 2.58 dB — all now FAIL the ≤+1 dB clause purely because the quiet baseline grew
  from 4 to 6 days, while their local−far margins are 6.5–10.4 dB. Under the intent bar nothing was
  lost and 12 stations were gained.
- **Verdict unchanged: STOP ×3.** Best configuration (6-day, intent bar) gives R1 8 sites / R3 10 /
  R2 7 and `frac cells with ≥4 admitted sites ≤45 km` = **0.024 / 0.083 / 0.000** against a bar of
  0.50 (0.218 / 0.122 / 0.014 even at ≥3 sites). The clause needs ~30 sites per sub-region.
- `shortlist_v1.csv` is now the 6-day product; the first pass is `shortlist_v1_4day*.csv`.
- **2026-09-20 05:05Z P4-0b partial (given-seed arm): 0/18 assembled under the frozen clauses → frozen verdict "GIVEN FAIL".** Mechanics (`results/p4/p4_0b/run_given0_survival.csv`, `score_given.csv`): pass-0 single-event templates yield median 6 independent hits (75th pct 11.5, max 34) on the 203-day pooled record, so the ≥20-independent-hit floor removed 49/54 seeds before iteration; the 5 survivors iterated twice, no burst convergence (median 2.7 % of hits within ±10 min of seed), and two of them reached recall at 0.8× the truth-stack reference (G44 0.36 vs ref 0.45; G75 0.36 vs 0.60) but failed the purity clause (G44 0.30 vs 0.34 needed). Truth-stack references themselves score recall 0.11–0.73 / purity 0.26–0.82 here. True hp events have median-over-sites z_med median 2.3, 65 % below 4, 32 % inside Merlin's 4–25 seed band → the template-free seed floor discards two-thirds of real LFEs (Merlin anticipated this). Sensitivity ablation on the pass-0 floor ordered (labelled, separate from the frozen result).

### 2026-09-20b — P4-0b SENSITIVITY ABLATION launched (frozen verdict unchanged)
Frozen given-seed arm: **0/18 assembled — stands**. Ablation (`scripts/v_abl.py`,
`scripts/v_p4_0b_abl_queue.sh`, pre-reg `results/p4/p4_0b_abl_prereg.json`, dir
`results/p4/p4_0b_abl/`, log `logs/p4_0b_abl_queue.log`, markers `P4_0B_ABL_DONE|FAILED`, **pid
1046043**) tests only the ≥20-independent-hit survival floor, set for a multi-year record and
applied to 203 days. Everything else (seeds, pass-0 hits, T rule, two iterations, dedupe, relative
clauses, truth-stack reference) is the frozen one, imported unchanged; `v_iterate.py`'s hard-coded
10-hit stacking floor is exposed via a verbatim copy in `v_abl.iterate_abl` rather than edited, and
the prereg records that all 10 frozen `v_` script hashes still match.
**A1** floor ≥8 / ≥3 days → **18 of 54 seeds survive (18 families)**; **A2** floor ≥5 / ≥2 days →
**21 survive**, vs 5 frozen. **A3 already landed and is the key diagnostic:** own-family purity of
pass-0 templates falls monotonically with hit count — median purity 1.00 (0–4 hits, 21 tpl), 0.71
(5–7), 0.50 (8–12), 0.53 (13–19), **0.36 (≥20, 5 tpl)**; Spearman(purity, n_indep) = **−0.68**. So
the floor and purity trade directly: the templates that clear ≥20 hits are the impure ones, which
is why G44/G75 passed recall and failed purity. **A4:** truth-stack reference medians — T=8.0
recall 0.454 / purity 0.489 / 70 indep; T_ext 8.25 → 0.448 / 0.515 / 60; T_core 8.75 → **0.398 /
0.565 / 48.5**. Write-up (frozen verdict first, ablation second, seed-z band third) on the next
resume, once both DONE markers exist.
- **2026-09-20 05:50Z P4-0b ablation DONE** (`results/p4/p4_0b_abl/`, labelled sensitivity, not the frozen result): pass-0 floor ≥8/≥3 days → 1/18 assembled (5 near-misses); floor ≥5/≥2 days → 1/18 (6 near-misses); frozen 0/18. The survival floor is not what fails the given-seed arm; iterated templates fall short mainly on purity. Frozen free arm still growing (101/203 days at 05:46Z).

## 2026-09-20c — P4-0b RESULT: the true Arm B construction does NOT assemble known families (full text `results/P4_0B_RESULTS.md`)

Both queues done (frozen 3.04 h pid 1040534; ablation 0.64 h pid 1046043).

**FROZEN VERDICT (pre-registered, `results/p4/p4_0b_freeze.json` sha256 `28fe6d5d02a8cc1c…`):**
- **Arm (i) given seeds: 0/18 eligible kept families assembled** (bar ≥12/18) → **FAIL**. All 18
  testable — a measured zero, not an untestable one.
- **Arm (ii) template-free: 1 kept family rediscovered** (family 44 by `F187`, 34 one-to-one ±2 s
  matches = 36 % of its 95 hp events, strict b=0 gives 33, null mean 0.6 / max 2, 0 aliases)
  against a bar of ≥3 → **FAIL**.
- **Pre-registered consequence: STOP P4-2 — do not run Puget discovery.**

**Mechanics.** Pass-0 survival 54 → 5 on a single binding clause (≥20 independent hits: 5 pass;
median indep hits over 203 days is **6**, p75 11.5, max 34). The five survivors iterated cleanly
(per-site S-SNR 9.85 → 11.95 dB, S-excess 9.64 → 12.15 dB, screen 5/5 at 10/10 sites) and hit counts
tripled, but **purity fell monotonically with every iteration**: G44_0 0.471→0.311→0.298,
G75_1 0.281→0.145→0.122, G110_1 0.391→0.215→0.185, G61_0 0.364→0.130. Near-passes: **family 44**
recall 0.358 (needs 0.226, ref 0.453) but purity 0.298 vs **0.341** needed; **family 75** recall
0.360 (needs 0.300) but purity 0.122 vs 0.207. Burst convergence **did not occur** (0/5 families,
near-seed fractions 0.011–0.065) — the decluster/independence guards work.

**Ablation (labelled not-the-frozen-result, `results/p4/p4_0b_abl_prereg.json`, all 10 frozen
hashes verified unchanged).** Floor ≥8/≥3 d → 18 survivors, 11 families iterate, **1/18 assembled**
(family 110 via G110_0), 5 near-misses; floor ≥5/≥2 d → 21 survivors, 12 families, **1/18**,
6 near-misses. **Every near-miss fails on purity, none on recall.** A3: own-family purity vs hit
count 1.00 (0–4 hits) / 0.71 / 0.50 / 0.53 / **0.36 (≥20)**, Spearman **−0.68**. A4: reference
medians recall/purity 0.454/0.489 at T=8.0, 0.448/0.515 at 8.25, 0.398/0.565 at 8.75 — the
0.11–0.73 reference-recall spread is a property of the families, not of the 8.0 floor.
**So the floor was NOT the cause: lowering it doubles the families reaching iteration (5→11→12) and
moves the score only 0/18 → 1/18, because a lower floor admits more LOW-PURITY templates — the very
axis the bar fails on.**

**Arm (ii) seeds.** 883 tremor hours → 529,695 windows scored → 48,730 in band → 48,220 after the
±60 s EQ veto → **1,590 after the caps** → 1,500 by budget. **The cap worked**: max 3/hour,
top-10-hours **2.0 %** (vs the >20 % failure signature), 609 hours / 59 days. The HF veto never
fired (in-band ratios ≈0.1–0.3 vs the 3.5 cut) — inert on this network. Survival 1500 → 2
(≥20 indep 19; ≥3 days 165; ≥2 episodes 224). The second survivor **F1184** (24 indep hits, 18 days,
3 episodes, +7.7…+19.1 dB S excess at 9/10 sites, N-moveout ratio 4.0) matches **no** eligible kept
family (best 1 match, null max 1). Hidden families were NOT matched against it, by governance — what
it is stays open.

**Seed-z band of TRUE hp events (2,469 events, the number that would fix a Puget seed rule):**
percentiles 5/25/**50**/75/95 = −0.51 / 0.53 / **2.32** / 6.09 / 20.28. **64.6 % below z_med 4**,
3.6 % above 25, only **31.8 % inside Merlin's 4–25 band with ≥3 sites at z≥3**, and only 50.8 % have
≥3 sites at z≥3 at all. Caveat: the 20 Hz cache is 1–8 Hz, not native 2–8 Hz, which biases z_med low
by an unmeasured amount. So the pre-registered band discards ~2/3 of real LFEs on this network.

**Scope (strict).** This shows that on 7 PQ-era sites and 203 pooled episode days, single-event seed
→ matched-filter growth → two iterated stacks does not assemble known families to 0.5×/0.8× of their
own truth stacks, and that the failure is purity. It does NOT show the construction fails with
multi-year records, more than two iterations, a different growth statistic, a different template
length/cell, or anything about Puget — all untested, none advocated.
Measured growth rate ≈940 template-station-days/s; quiet-day FA at net_mad ≥8 is 0, so T_seed is the
8.0 recording floor for 96–99 % of templates.
- **2026-09-20 ~08:00Z — P4 CLOSED.** Merlin close-out saved verbatim (`results/merlin_review_p4_closeout_2026-09-20.md`); PLAN.md §"P4 — CLOSED" records the decision. Corrections pass ordered (cluster-truth re-score, F1184 Gate 1, dated correction sections in both results files; CPU only). Discovery line (P3D + P4) closed; parked items are the front estimator v2 and the far-south existence test (download needs user OK).

## 2026-09-20d — Merlin close-out: two mechanism statements RETRACTED, STOP stands (corrections filed)

Merlin's close-out (`results/merlin_review_p4_closeout_2026-09-20.md`): **STOP for Puget is right and
over-determined, but not for the reasons the write-ups gave**, and two of the failures trace to his
own spec. Three CPU-only steps run; corrections filed as dated sections in
`results/P4_0B_RESULTS.md` and `results/P4_0_RESULTS.md` (originals left unedited).

**Step 1 — post-hoc cluster-truth re-score** (`scripts/v_rescore.py`, `results/p4/p4_0b_rescore/`).
Truth for family X = frozen detections of all **non-hidden** families within 3.0 km hypocentral,
deduped at 2 s; templates/hits/reference/clauses unchanged.
**frozen 1/18 (p_real≥0.2), 2/18 (any p); A2 ablation 5/18, 6/18.** Max 6/18 → **Merlin's ≤6/18
closure criterion met, the ≥12/18 failure signature did NOT fire.** His predictions reproduced
exactly: G44_0 any-p purity 0.667 vs REF44 0.782 → relative **0.852**, assembles; G75_1 0.806,
assembles. Also recorded: quiet-day FA at T=8.0 is **13 hits / 1,080 template-days = 0.012 per
template-day vs the 0.1 allowance** — his F4, the 8.0 growth floor was never exercised, so the
pass-0 starvation is partly a floor convention.

**Step 2 — F1184 Gate 1 at the seven never-summed sites** (`scripts/v_f1184_gate1.py`):
**FAIL, 0 of 7 above null p99.** Every observed alignment-free S excess is NEGATIVE (CLRS −1.82,
GOBB −2.42, LOPEZ −2.60, TURTL −1.02, LRIV −1.58, NLLB −0.67, B926 −1.17 dB) against null p99
+1.94…+3.22; tremor co-timing 0.500 vs same-day base 0.522 (median 0 PNSN epicentres within 60 km
±30 min). So F1184's +7.7…+19.1 dB at the DETECTION sites was selection bias (Merlin F3) — it is not
an LFE family. (Engine note: the first Gate-1 implementation thrashed a per-cut DayCache — 25k day
loads — and was rewritten day-major, 2 h → 2 min.)

**Step 3 — corrections.** (a) **RETRACTED the purity-collapse mechanism**: 97.5 % of hp events are
co-detected by another family's template and the eligible families sit 0.12–1.6 km apart
(52→0.12 km, 75→0.17, 21→0.45), so a single-label p≥0.8 one-to-one purity is not a noise measure and
the A3 −0.68 anti-correlation is what an incomplete truth produces as a template gains sensitivity.
(b) **RETRACTED §5's "lower the seed floor"**: only **2.2 %** of the 1,500 seeds are within ±2 s of a
kept hp event and **12.6 %** of any frozen detection (chance 0.21 %/2.82 % → 10.7×/4.5× enriched, so
not noise — but **flat in z_med**: 9.4/13.5/12.9/13.2 %), so no band, floor or budget can raise it.
There is no energy-based seed rule for LFEs on this network — the decisive negative for Puget, which
has only template-free seeds. (c) The P4-0a **50 % coverage clause was unattainable even with all
118 sites admitted** (R1 0.545 / R3 0.340 / R2 0.599 — R3 could never pass); the standing geometry
finding is the **admitted-set** coverage (0.024/0.083/0.000 at ≥4 sites; 0.218/0.122/0.014 at ≥3).
(d) The **intent far-tremor bar** (local − far ≥ 3 dB) is the operative one. (e) The truth-stack
reference is **circular** (trained on the scored events) — harder, not unsound. (f) The N-moveout
null is a **window-exclusion** null, not a moveout scramble. (g) Credited **F187 → family 44** from a
seed **184 s** from any catalogued hp event, 34 matches vs null max 2, 0 aliases, median 21 PNSN
epicentres — the one positive datum; note F187's 1.000 tremor co-timing is uninformative (base rate
0.993), the epicentre count is the discriminator. (h) **Merlin's F10 VERIFIED as a mechanism,
immaterial here**: a synthetic 30 % gap deflates the lowhalf MAD ≈2× (z +94…+125 %), but the
gappiest real station-day in the run (DAOB 2024-08-12, 4.5 % zeros) costs **+5.9 %** in z and only
0.16 % of 631 scanned station-days exceed 1 % zeros, none 10 %; gap samples also cannot become hits
(the valid mask zeroes them). Worth fixing before reuse on gappier data.
Merlin's Q4 plain-language standing conclusion is pasted verbatim as the closing section of
`results/P4_0B_RESULTS.md`.
- 2026-09-20 ~09:30Z user idea: the temporary arrays behind every published family set are public at EarthScope (FDSN verified: CAFE XU 63 stns 2006–2012 WA; Array of Arrays XG 124 stns 2009–2011; ZK/Z3 2008/2009; FACES YW; COLZA XA; FAME XQ). Route: discover WA/OR families on the arrays for the 2007–2011 ETS windows (Bostock autocorrelation), register onto PB/UW in overlap years, extend to 2026 with the frozen chain. Assessment given; not planned/approved yet.
- 2026-09-20 ~14:00Z (from lfe_array_discovery agent E): bug in scripts/v_engine.py DayData.station_z — one self.NF reused across stations; a short station-day raises a shape error. Not fixed on disk (frozen P4-0b product); fixed in lfe_array_discovery w_ subclass. Add to the P4 close-out defects list (C-8 gap-MAD, this NF reuse) before any reuse of v_engine.
