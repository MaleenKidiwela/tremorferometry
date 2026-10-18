# `g_gates` — one page for the P2 / P5 agents

**What it is.** A detector-agnostic family-acceptance pipeline. You hand it a *proposal* — a set of
detection times on the origin axis plus a proposed interface cell — and it returns per-gate statistics,
every null, and PASS / FAIL / UNTESTABLE per gate and overall. It does not care where the proposal came
from (template matched filter, back-projection, the learned model).

## Call it from Python

```python
import sys; sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import g_gates as GG

ctx = GG.GateContext(GG.GateConfig(withheld=("VGZ", "SNB", "B926")))   # build ONCE, reuse for all proposals
p   = GG.Proposal(pid="cand_0001", times=t_epoch_seconds, lat=48.41, lon=-123.55,
                  dep=np.nan,                 # NaN -> Slab2 at (lat, lon)
                  offset_search_s=0.0,        # >0 searches a CONSTANT time offset; the same search is
                                              # applied to every null, so the freedom is paid for
                  meta=dict(drop_bkey=-1))    # bkey to remove from the accepted list (leave-one-out)
r = GG.evaluate(ctx, p, gates=(0, 1, 2, 3), evaluated_days=days)        # days = np.datetime64[D] array
print(r["verdict"], GG.flatten(r))
```

`GateContext` is expensive (travel-time table, location grid, cross-family table, day cache); build it
once per process. `evaluate` on a warm cache costs **~4 s** median, ~20 s p90 (see GATES_RESULTS.md).

## Or from the command line

```bash
python scripts/g_gates.py proposals.parquet --out results/gates/run1 --gates 0,1,2,3
```

`proposals.parquet` columns: `pid`, `t` (epoch seconds, one row per member), `lat`, `lon`,
optional `dep`, `offset_search_s`, `drop_bkey`, `kind`, `truth_bkey`.

## What each verdict means

| verdict | meaning |
|---|---|
| `PASS` | all of Gates 0-2 pass and Gate 3 passes or is untestable |
| `PASS_inter_ETS_candidate` | Gates 0-2 pass, Gate 3 fails — report it, do **not** reject it |
| `FAIL` | at least one of Gates 0-2 failed |
| `UNTESTABLE` | no gate failed but at least one had no power here (too few members, no registering site, geometry, or the positive control itself fails) |

`UNTESTABLE` is **not** a pass. A family may only enter the catalog on `PASS`.

## The seven things you must not get wrong

1. **The withheld stations must really be withheld.** Gate 1 is only evidence if the stations in
   `cfg.withheld` were excluded from the detection run that produced the proposal, for the whole run.
   Feeding it detections that used those stations turns Gate 1 into a self-confirmation.
2. **`t` is the ORIGIN axis**, not a pick time and not a window start. The predicted S at a station is
   `t + ttS(Crosson, cell->station) + a5 station correction`. If your detector's time is a template
   window start, set `offset_search_s` (30 s is the tested value) so the constant offset is searched —
   and read the fitted `offset_fit_s` back out.
3. **N matters more than anything else.** Gate 1 power is ~0.4 at N=20, ~0.75 at N=50, ~1.0 at N=100
   (own-event proposals). Gate 2 is essentially untestable below N≈200 on the Vancouver Island network.
   Do not propose a family with fewer than 20 members; do not expect a location below ~150.
4. **Gate 0 needs an accepted-family list and a leave-one-out rule.** Pass `meta["drop_bkey"]` when the
   proposal is a re-proposal of a family that is itself in the accepted list; otherwise Gate 0 can only
   say "you are yourself". For genuinely new proposals leave it at -1.
5. **Gate 3 needs a long evaluated span.** The circular-shift null requires ≥ 2 × 30 d of calendar; on a
   single ETS window it returns UNTESTABLE by design. Pass the real `evaluated_days` of your run (the
   days the detector actually scanned), not the member days.
6. **Nothing here is a completeness statement.** The gates say a source exists, is new, is where you say
   it is and fires with local tremor. Per-family completeness still comes from the frozen template
   machinery (Track C/D).
7. **Runtime scales with members and with cold cache.** The first proposal on a new day costs ~20 GB of
   NFS reads per station-year; sort proposals by time so the `DayStore` LRU works for you.

## Rediscovery scoring (P3)

```bash
python scripts/g_rediscover.py --clusters clusters.parquet --events bostock.parquet \
       --families e5r_families_2012.csv --hidden 8,10,12,... --windows 2012-08-30,2012-10-05
```
`clusters.parquet`: `cid`, `t`, `lat`, `lon`. `events.parquet`: `bkey`, `t`. Returns per-hidden-family
rediscovery (centroid ≤ 5 km, ≥ 30 % one-to-one recall at ±2 s, ≥ 50 % member purity, ≤ 50 % explained
by a kept family) and the spurious-split rate of the kept families. **Calibration you need to know:**
feeding it the frozen template detector's own per-family detections scores only 9/30 hidden families as
rediscovered — the ≥ 30 % one-to-one recall clause is above the template detector's own recall. Read a
low P3 number against that ceiling, not against 1.0.

## Files

| file | what |
|---|---|
| `scripts/g_gates.py` | the pipeline (Gates 0-3, `Proposal`, `GateConfig`, `GateContext`, `evaluate`) |
| `scripts/g_common.py` | day store, travel times, stacking, band ratios, `align_stack` |
| `scripts/g_locate.py` | the a5/d3 location estimator + the 2-D Slab2 interface mode + synthetic floor |
| `scripts/g_rediscover.py` | the P3 rediscovery scorer |
| `scripts/g_crossfam.py` | per-station cross-family cc p95 (Gate 0 chance level) |
| `scripts/g_controls.py` | the Gate-1 positive-control table |
| `scripts/g_validate.py`, `g_report.py`, `g_figures.py` | the validation suite |
| `results/gates/g_crossfam_cc.csv` | Gate-0 chance levels |
| `results/gates/g_gate1_controls.csv` | Gate-1 control table (station × family × N) |
| `results/gates/g_val_*.csv/json` | validation results, one file per set |
