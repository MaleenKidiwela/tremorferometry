#!/usr/bin/env python
"""D3 stage 4 control — how well CAN this station geometry locate, given perfect model knowledge?

Forward-models S and P times for a synthetic source in the same 1-D model, then inverts with
`d3_stage4_locate.locate_from_tables`. Because forward and inverse share the model, this measures
ONLY the geometry + pick-noise limit (the "formal" precision); it says nothing about the 1-D-model
systematic, which is why the located depths are also reported without the Slab2 prior in the main
run. Output: lfe_catalog_repair/results/d3_stage4_synthetic_location_test.txt
"""
import sys
import numpy as np
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import d3_stage4_locate as L
import d3_vel1d as V

STAS = ['PGC', 'SYMB', 'B011', 'VGZ', 'SNB', 'CLRS', 'LRIV', 'PTRF', 'OCP']
LAT0, LON0 = 48.5, -123.8


def kmen(la, lo, la0, lo0):
    return ((lo - lo0) * 111.2 * np.cos(np.radians(la0)), (la - la0) * 111.2)


def main(pick_sigma=(0.1, 0.2, 0.4)):
    out = []
    t = time.time()
    xs, zs, TP, TS = L.build_tt_table()
    out.append(f"tt table {time.time()-t:.1f} s {TP.shape}")
    E, Nn, gz = L.grid_axes(25, 0.5, 10, 60, 0.5)
    TTS, TTP = L.precompute_tt(E, Nn, gz, STAS, LAT0, LON0, (xs, zs, TP), (xs, zs, TS))
    out.append(f"grid {E.size} nodes x {gz.size} depths; stations {STAS}")
    rng = np.random.default_rng(0)
    for (te, tn, tz) in [(6.0, -4.0, 38.0), (0, 0, 32.0), (-12, 9, 45.0), (18, -15, 42.0)]:
        tS, tP = [], []
        for s_ in STAS:
            e, n = kmen(*L.STA_LL[s_], LAT0, LON0)
            d = np.hypot(te - e, tn - n)
            tS.append(V.tt(d, tz, 'S'))
            tP.append(V.tt(d, tz, 'P'))
        tS, tP = np.array(tS), np.array(tP)
        r = L.locate_from_tables(tS, tP, list(range(len(STAS))), TTS, TTP, E, Nn, gz,
                                 tz, LAT0, LON0)
        ge, gn = kmen(r[0], r[1], LAT0, LON0)
        out.append(f"true ({te},{tn},{tz}) noise-free -> ({ge:.2f},{gn:.2f},{r[2]:.1f}) "
                   f"err_h {np.hypot(ge-te, gn-tn):.2f} km err_z {r[2]-tz:+.1f} km")
        for sg in pick_sigma:
            for tag, pp in (("S+P", tP), ("S only", np.full(len(STAS), np.nan))):
                errs = []
                for _ in range(25):
                    r2 = L.locate_from_tables(tS + rng.normal(0, sg, len(STAS)),
                                              pp + rng.normal(0, sg, len(STAS)),
                                              list(range(len(STAS))), TTS, TTP, E, Nn, gz,
                                              tz, LAT0, LON0)
                    g2 = kmen(r2[0], r2[1], LAT0, LON0)
                    errs.append((np.hypot(g2[0] - te, g2[1] - tn), r2[2] - tz))
                errs = np.array(errs)
                out.append(f"   {tag:6s} pick sigma {sg:.2f} s: err_h {errs[:,0].mean():.2f}"
                           f"+-{errs[:,0].std():.2f} km, err_z {errs[:,1].mean():+.2f}"
                           f"+-{errs[:,1].std():.2f} km")
    txt = "\n".join(out)
    print(txt)
    Path("lfe_catalog_repair/results/d3_stage4_synthetic_location_test.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
