# S12-LVL, LFE side — can our LFE data see the post-ETS LVL speed-up of Gosselin et al. (2020)?

Started Tue 2026-10-13 ~9:45 am PT. Rules: `PLAN.md` (§1–5 written ~10:12 am PT; §6 for the coordinator's tasks A–D
~10:25 am PT, both before the numbers they govern). Decision log: `NOTES.md`. Results in plain words: `RESULTS.md`.
Companion: receiver-function reproduction in `../rf_gosselin/` (another agent).

## Layout
- `scripts/s12_forward.py` — task 1: LVL forward calibration (MC coda kernels with an LVL tally; ray-theory reflections).
- `scripts/s12_reflect.py`, `s12_reflect_inject.py` — task 2: LVL-base reflection search on the family grand stacks,
  scrambled control, injection threshold.
- `scripts/s12_epoch.py`, `s12_epoch_pool.py` — tasks 3–4: superposed epoch on the front (lapse-resolved dv/v, decay,
  P/S, S−P) against a between-ETS baseline; SVI and WA-NW.
- `scripts/s12_novel_measure.py`, `s12_novel.py` — coordinator tasks A–D (propagating onset, recharge, source properties,
  episode size).
- `scripts/s12_figures.py` — figures. Queues `scripts/s12_queue.sh`, `s12_queue2.sh` (DONE markers in this folder);
  restart = block 15 of `scripts/resume_after_reboot.sh`.
- `tables/` (forward + reflection), `tables/svi/`, `tables/wa_nw/` (epoch + A–D), `figures/`, `logs/`.

## Check progress
`tail logs/queue.log`; markers `QUEUE_DONE` (tasks 1–4) and `QUEUE2_DONE` (A–D).
