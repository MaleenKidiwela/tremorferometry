#!/usr/bin/env python
"""Track A-2 follow-up (2) — build the MERGED per-epoch template sets Track C consumes, and write
`results/trackC/templates_by_epoch.csv` (`sta, epoch_start, npz`).

Why merged and not the identity-filtered sets directly: `a4_templates_<STA>_<epoch>.npz` holds only the
families that PASSED the epoch-identity test, which for PGC is 46 of the ~113 families that have a 2012
template at that station.  Handing Track C that file alone would silently delete two thirds of PGC's
family coverage.  So each merged set is

    every 2012 (family, station) template of that station          <- the frozen baseline, unchanged
    REPLACED, family by family, by the new-epoch re-stack          <- where the identity test ACCEPTED it
    (polarity already applied by a4_epoch_report.py)

which strictly dominates the status quo: full coverage, and an epoch-matched waveform wherever we have
evidence that the epoch-matched waveform is the same LFE.  The 8-s template of a replaced entry is the
first 8 s of its own 12-s template, re-normalised per component (the normalisation is per component, so
the prefix of the 12-s template IS the 8-s template).

B926 is NOT written: it has no epoch boundary, it was only the negative control, and its 3-family
"accepted" set would cripple a station that needs no per-epoch treatment at all.

LRIV IS built but is deliberately NOT listed in `templates_by_epoch.csv`: its epoch identity FAILED
(cc median -0.22 against its own 2012 stacks, worse than the no-change control), so its post-2019 rows
must keep Track C's `2012_provisional` label rather than be relabelled `per_epoch:<file>` as if the
transfer had been verified.  The merged file is on disk for whoever decides to use it; `--with-lriv`
adds the row.
"""
from __future__ import annotations
import argparse, os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
FS = 20.0

# station -> (epoch start Track C keys on, the a4 file that covers it)
EPOCHS = {
    "PGC":  ("2017-08-10", "a4_templates_PGC_2018_HH.npz"),
    "NLLB": ("2017-02-23", "a4_templates_NLLB_2018_HH.npz"),
    "VGZ":  ("2017-03-03", "a4_templates_VGZ_2018_HH.npz"),
    "SNB":  ("2017-06-08", "a4_templates_SNB_2018_HH.npz"),
    "LRIV": ("2019-05-30", "a4_templates_LRIV_2021_HH.npz"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=f"{EXT}/results")
    ap.add_argument("--with-lriv", action="store_true",
                    help="also list LRIV in templates_by_epoch.csv (its identity test FAILED)")
    a = ap.parse_args()
    Z12 = np.load(f"{REP}/results/e5r_templates_2012.npz", allow_pickle=True)
    T12, S12 = Z12["tpl"].item(), Z12["sp"].item()
    ID = pd.read_csv(f"{EXT}/results/a4_epoch_identity.csv")
    rows, idx = [], []
    for sta, (estart, fn) in EPOCHS.items():
        src = f"{EXT}/results/{fn}"
        if not os.path.exists(src):
            print(f"[{sta}] {fn} missing — skipped")
            continue
        Za = np.load(src, allow_pickle=True)
        Ta, Sa = Za["tpl"].item(), Za["sp"].item()
        tpl = {k: {L: v.copy() for L, v in d.items()} for k, d in T12.items() if k[1] == sta}
        sp = {k: S12[k] for k in tpl}
        n12 = len(tpl)
        n_repl = 0
        for k, d in Ta.items():
            t12 = d[12.0]
            seg8 = t12[:, :int(8 * FS)].astype(np.float32)
            nn = np.linalg.norm(seg8, axis=1)
            if np.any(nn <= 0):
                continue
            tpl[k] = {8.0: (seg8 / nn[:, None]).astype(np.float32), 12.0: t12.astype(np.float32)}
            sp[k] = float(Sa[k])
            n_repl += 1
        out = f"{a.out_dir}/a4_templates_merged_{sta}_{estart}.npz"
        np.savez_compressed(out, tpl=np.array(tpl, dtype=object), sp=np.array(sp, dtype=object))
        q = ID[(ID.station == sta) & (ID.key == sta)]
        rows.append(dict(sta=sta, epoch_start=estart, npz=os.path.relpath(out, EXT),
                         n_templates=len(tpl), n_from_2012=len(tpl) - n_repl, n_new_epoch=n_repl,
                         n_tested=int(len(q)), n_ctl_ok=int((q.verdict != "floor_fail").sum()) if len(q) else 0,
                         source_a4=fn))
        if sta != "LRIV" or a.with_lriv:
            idx.append(dict(sta=sta, epoch_start=estart, npz=os.path.relpath(out, EXT)))
        print(f"[{sta} @ {estart}] {len(tpl)} templates = {len(tpl)-n_repl} frozen-2012 + {n_repl} "
              f"new-epoch re-stacks -> {os.path.basename(out)}", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(f"{EXT}/results/a4_templates_by_epoch_provenance.csv", index=False)
    os.makedirs(f"{EXT}/results/trackC", exist_ok=True)
    pd.DataFrame(idx).to_csv(f"{EXT}/results/trackC/templates_by_epoch.csv", index=False)
    pd.set_option("display.width", 200)
    print("\n" + R.to_string(index=False))
    print("\n[write] results/trackC/templates_by_epoch.csv, results/a4_templates_by_epoch_provenance.csv")


if __name__ == "__main__":
    main()
