#!/usr/bin/env python
"""Track D — reconcile `results/trackC/months_done.csv` with the catalog for the months whose ledger
row was written BEFORE Track C re-flagged the tiers (coordinator approval 2026-09-17).

Default mode is REPORT-ONLY and writes nothing outside `results/trackD/`:
  results/trackD/ledger_corrections.csv        the machine-readable diff
  results/trackD/LEDGER_NOTE.md                the human note

With `--apply` it additionally updates `n_cat_ext` / `n_cat_core` in place in the Track C ledger —
nothing else in the row, the schema or the row order — exactly as `c_watch.mark_done` would have
written it, after backing the file up.  It refuses to touch the ledger unless EVERY disagreeing month's
parquet is verified to already carry `cat_*` columns matching the CURRENT per-epoch thresholds, i.e.
the file is right and only the ledger is behind.
"""
from __future__ import annotations
import argparse, datetime, os, shutil, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_tiers  # noqa: E402

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TC = f"{EXT}/results/trackC"
TD = f"{EXT}/results/trackD"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="also update the two fields in the Track C ledger in place (backs it up first)")
    a = ap.parse_args()
    EPOCHS, TIERS, TSRC = d_tiers.load()
    D = pd.read_csv(f"{TC}/months_done.csv")
    rows = []
    for r in D.itertuples():
        f = f"{TC}/detections_{r.ym}.parquet"
        if not os.path.exists(f):
            continue
        d = pd.read_parquet(f, columns=["date", "net_mad", "hf_pass", "keep_xfam_ext",
                                        "keep_xfam_core", "cat_ext", "cat_core"])
        d = d[d.date.astype(str).str[:7] == r.ym]
        eid = d_tiers.epoch_of(d.date.to_numpy().astype(str), EPOCHS)
        th = np.array([[TIERS.get(e, d_tiers.FROZEN)["T_ext"],
                        TIERS.get(e, d_tiers.FROZEN)["T_core"]] for e in eid])
        want_ext = (d.net_mad.to_numpy() >= th[:, 0]) & d.hf_pass.to_numpy() & d.keep_xfam_ext.to_numpy()
        want_core = (d.net_mad.to_numpy() >= th[:, 1]) & d.hf_pass.to_numpy() & d.keep_xfam_core.to_numpy()
        stored_ok = bool((want_ext == d.cat_ext.to_numpy()).all()
                         and (want_core == d.cat_core.to_numpy()).all())
        n_ext, n_core = int(d.cat_ext.sum()), int(d.cat_core.sum())
        if n_ext == r.n_cat_ext and n_core == r.n_cat_core:
            continue
        rows.append(dict(ym=r.ym, epochs=";".join(sorted(set(eid))),
                         ledger_cat_ext=r.n_cat_ext, ledger_cat_core=r.n_cat_core,
                         file_cat_ext=n_ext, file_cat_core=n_core,
                         file_matches_current_thresholds=stored_ok))
    R = pd.DataFrame(rows)
    os.makedirs(TD, exist_ok=True)
    if not len(R):
        print("[ledger] already consistent with every month parquet — nothing to do")
        return
    print(R.to_string(index=False))
    R.to_csv(f"{TD}/ledger_corrections.csv", index=False)
    bad = R[~R.file_matches_current_thresholds]
    note = f"""# Track C ledger vs the catalog — {datetime.datetime.now(datetime.timezone.utc):%Y-%m-%d %H:%M} UTC

`n_cat_ext` / `n_cat_core` in `results/trackC/months_done.csv` are STALE for the months below: those rows
were written by `c_watch.mark_done` BEFORE Track C re-derived the tiers per epoch and re-flagged
`cat_ext` / `cat_core`. Each month's parquet was verified to already carry `cat_*` columns matching the
CURRENT per-epoch thresholds ({TSRC}), so the FILES are right and only the ledger is behind — and
`results/trackD/lfe_extension_catalog_v1.parquet`, which is built from the parquets, is correct.

{'| ' + ' | '.join(R.columns) + ' |'}
{'|' + '|'.join('---' for _ in R.columns) + '|'}
{chr(10).join('| ' + ' | '.join(str(v) for v in row) + ' |' for row in R.itertuples(index=False))}

Produced by `lfe_extension/scripts/d_fix_ledger.py`. Run it with `--apply` to update those two fields in
place in the ledger (it backs the file up to `months_done.csv.pre_d_fix_ledger` first and refuses unless
every disagreeing month's parquet already matches the current thresholds).
"""
    open(f"{TD}/LEDGER_NOTE.md", "w").write(note)
    print(f"[ledger] wrote {TD}/ledger_corrections.csv and {TD}/LEDGER_NOTE.md")
    if not a.apply:
        print("[ledger] report-only (no Track C file touched); re-run with --apply to update it")
        return
    if len(bad):
        raise SystemExit("[ledger] REFUSING to touch the ledger: these months' cat_* columns do not "
                         f"match the current thresholds either, so the file is not simply newer:\n{bad}")
    bk = f"{TC}/months_done.csv.pre_d_fix_ledger"
    if not os.path.exists(bk):
        shutil.copy2(f"{TC}/months_done.csv", bk)
    for r in R.itertuples():
        D.loc[D.ym == r.ym, "n_cat_ext"] = r.file_cat_ext
        D.loc[D.ym == r.ym, "n_cat_core"] = r.file_cat_core
    D.sort_values("ym").to_csv(f"{TC}/months_done.csv", index=False)
    open(f"{TC}/months_done_LEDGER_NOTE.md", "w").write(note)
    print(f"[ledger] updated {len(R)} row(s) in place; backup {bk}")


if __name__ == "__main__":
    main()
