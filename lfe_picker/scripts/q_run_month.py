#!/usr/bin/env python
"""q_ wrapper around the FROZEN matched filter (scripts/o_run_month.py) with two overrides that the
frozen scripts do not have and that P3 needs.  o_run_month.py and o_common.py are NOT edited; both
overrides are monkey-patches applied in this process only, and both are recorded in the run json.

  --drop-stations SNB,B926      withhold stations from the network for the WHOLE run (Merlin's P3
                                setting: discovery and Gate 1 both run without SNB and B926).  The
                                frozen template npz and the frozen operating-point json are read
                                unchanged; the (family, station) entries of the withheld stations are
                                dropped from the loaded tables, exactly as if those stations had no
                                template.  The frozen cross-family station-sharing matrix is rebuilt
                                on the surviving station set (o_common.share_ok_matrix(stas=...)),
                                because station sharing is a property of the network actually used.
  --bank <npz>                  a DISCOVERED template bank in the frozen npz format, used instead of
                                the frozen 2012 templates (o_run_month's own --templates already does
                                this; the option is mirrored here so a discovered bank and a station
                                override can be combined in one command).

Everything else — placement, weights, moveout, per-family-day MAD, declustering, the FA-only tier
rule, the HF veto — is the frozen code path, untouched.

    python scripts/q_run_month.py 2012-09 --drop-stations SNB,B926 --tag noSNBB926 \
        --out-dir results/p3/trackQ --cache-dirs /home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz
"""
from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o_common as C
import o_run_month as R

_DROP = set()
_orig_load = C.load_frozen_tables
_orig_share = C.share_ok_matrix


def _load(fz, tpl_path, L):
    tpl, sp, F, wgt, mv, ampref = _orig_load(fz, tpl_path, L)
    if not _DROP:
        return tpl, sp, F, wgt, mv, ampref
    tpl = {k: v for k, v in tpl.items() if k[1] not in _DROP}
    sp = {k: v for k, v in sp.items() if k[1] not in _DROP}
    wgt = {k: v for k, v in wgt.items() if k[1] not in _DROP}
    mv = {k: v for k, v in mv.items() if k[1] not in _DROP}
    print(f"[q_run_month] withheld stations {sorted(_DROP)}: "
          f"{len(tpl)} (family, station) templates remain", flush=True)
    return tpl, sp, F, wgt, mv, ampref


def _share(stas=None):
    if _DROP and stas is None:
        stas = [s for s in C.STAS_ORDER if s not in _DROP]
    return _orig_share(stas)


def main():
    argv = sys.argv[1:]
    out = []
    i = 0
    bank = ""
    while i < len(argv):
        if argv[i] == "--drop-stations":
            _DROP.update(x for x in argv[i + 1].split(",") if x)
            i += 2
        elif argv[i] == "--bank":
            bank = argv[i + 1]
            i += 2
        else:
            out.append(argv[i])
            i += 1
    if bank:
        out += ["--templates", bank]
    C.load_frozen_tables = _load
    C.share_ok_matrix = _share
    R.C.load_frozen_tables = _load
    R.C.share_ok_matrix = _share
    sys.argv = [sys.argv[0]] + out
    R.main()


if __name__ == "__main__":
    main()
