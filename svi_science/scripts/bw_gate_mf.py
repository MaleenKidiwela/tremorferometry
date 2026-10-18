#!/usr/bin/env python
"""bw_gate_mf -- run the FORWARD SVI feature matched filter (lfe_extension/scripts/svi_ref_mf.py, unchanged) with
the backward station sets (CN4 / B1 / B2 / B3) and, with --backward, the pre-2010 station epochs, download ledger
and trace root (see bw_gate_common.patch_forward).  Output root = $SVI_REF_DIR (set here to
svi_science/results/backward_gates/mf unless already set), so nothing is written into lfe_extension/results.

  python bw_gate_mf.py [--backward] -- <svi_ref_mf.py arguments>
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("SVI_REF_DIR", "/home/jovyan/tremorferometry/svi_science/results/backward_gates/mf")
import bw_gate_common as G

argv = sys.argv[1:]
backward = "--backward" in argv
argv = [a for a in argv if a not in ("--backward", "--")]
R, C = G.patch_forward(backward=backward)
import svi_ref_mf
sys.argv = ["svi_ref_mf.py"] + argv
svi_ref_mf.main()
