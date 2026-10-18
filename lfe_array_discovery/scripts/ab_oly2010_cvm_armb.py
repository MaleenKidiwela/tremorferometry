#!/usr/bin/env python
"""Locate the 47 new distinct Arm-B Olympic 2010 families (armb_new_distinct_28.csv) exactly as
ab_oly2010_cvm.py located Arm A: Hyp2000, USGS CVM v1.7 NW profile, Vp/Vs 1.823, R 1.5, free depth,
picks = ab_hypo.make_picks at FS 20, station codes / site coordinates = ab_oly2010_cvm.sites().
2026-09-29.

Arm A located SITE stacks (array beams + single stations, results/u_locate/wa_real_stacks.npz).
Arm B stacks (results/u_locate/wa_armb_stacks.npz) are per MEMBER station, so they are averaged
into the same sites first, as scripts/armb_site_score.py and catalogue v1.1 do: each member's 3-C
stack normalised by its max|.|, then the plain mean over the members of a site (arrays -> x<prefix>,
B005/6/7 -> xPBB005, B013/B943 -> xPBB013).  Differences from Arm A that remain:
  * the Arm-A beams are delay-and-sum (plane wave + static terms, o_beam.py) over all array
    members; the Arm-B site stack is a plain mean over the 1-3 members Arm B used;
  * the Arm-A beams are 1-8 Hz, the Arm-B member stacks are 2-8 Hz (o_chestler.BandCache);
  * ab_hypo.py now writes "ZTR 35 F" (fixed 2026-09-29); Arm A ran with "ZTR 35" (trial depth
    ignored; the SW_mute ZTR test showed no material difference).
--armA re-runs Arm A's real site stacks through this code (control: reproduces loc_USGS_NW.csv
up to the ZTR fix).

  python scripts/ab_oly2010_cvm_armb.py [--armA]
"""
import argparse, os, re, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_hypo_cvm                                                          # noqa: E402,F401
import ab_hypo as H                                                         # noqa: E402
import ab_common as A                                                       # noqa: E402
import ab_oly2010_cvm as O                                                  # noqa: E402

RES = O.RES
D = f"{RES}/ab_royer/OLY2010_cvm"
MODEL, VPVS, RMIN = "USGS_NW", 1.823, 1.5


def armb_site(s):                     # = scripts/armb_site_score.py site()
    if s in ("B005", "B006", "B007"):
        return "xPBB005"
    if s in ("B013", "B943"):
        return "xPBB013"
    p = s[:2]
    if p in ("BH", "BS", "DR", "GC", "LC", "PA", "PL", "TB") and re.fullmatch(r"[A-Z]{2}\d+", s):
        return "x" + p
    return s


def armb_site_stacks(code, grps):
    z = np.load(f"{O.UL}/wa_armb_stacks.npz", allow_pickle=True)
    want, raw = set(grps), {}
    for k in z.files:
        if k == "meta":
            continue
        g, s = k.split("_", 1)
        g = int(g)
        if g not in want:
            continue
        x = np.asarray(z[k], np.float64)
        m = np.abs(x).max()
        if np.isfinite(m) and m > 0:
            raw.setdefault(g, {}).setdefault(armb_site(s), []).append(x / m)
    return {g: {code[s]: np.mean(L, axis=0) for s, L in d.items()} for g, d in raw.items()}


def locate(ST, wd, pop):
    stas = sorted({s for d in ST.values() for s in d})
    P = H.make_picks(ST)
    ids = H.write_hyp_inputs(wd, P, stas, VPVS, RMIN, model=MODEL)
    ztr = [l.strip() for l in open(f"{wd}/run.hyp") if l.startswith("ZTR")]
    Hy, lg = H.run_hyp(wd)
    M = pd.DataFrame(dict(grp=sorted(ST))).merge(Hy, on="grp", how="left")
    M["located"] = M.lat.notna()
    M["n_input"] = M.grp.isin(ids)
    M["n_sites_stack"] = M.grp.map({g: len(d) for g, d in ST.items()})
    M["model"], M["vpvs"], M["rmin"], M["population"], M["ztr_line"] = MODEL, VPVS, RMIN, pop, ztr[0]
    print(f"[{pop}] {len(ST)} families, {len(P)} site stacks, located {int(M.located.sum())}, {ztr}", flush=True)
    return M, P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--armA", action="store_true")
    a = ap.parse_args()
    A.FS = 20.0
    code = O.sites()
    if a.armA:
        ST = O.load_stacks(f"{O.UL}/wa_real_stacks.npz", code)
        M, _ = locate(ST, f"{D}/run_armA_USGS_NW_ztrF", "real")
        M.to_csv(f"{D}/loc_USGS_NW_armA_ztrF.csv", index=False)
        return
    grps = pd.read_csv(f"{D}/armb_new_distinct_28.csv").grp.astype(int).tolist()
    ST = armb_site_stacks(code, grps)
    M, P = locate(ST, f"{D}/run_armb_USGS_NW", "armb_new")
    M.to_csv(f"{D}/loc_USGS_NW_armb.csv", index=False)
    P.to_csv(f"{D}/picks_armb.csv", index=False)


if __name__ == "__main__":
    main()
