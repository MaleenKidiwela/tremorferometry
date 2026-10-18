#!/usr/bin/env python
"""MENDO (2026-09-29): write results/ab_royer/MENDO/MENDO_SUMMARY.md from the real run and the scrambled
control (MENDO_control/, each site's day series shifted by whole days, AB_NULL=1).  Safe to run while the
control is still missing (its column then reads 'pending')."""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_common as A  # noqa: E402

R = f"{A.RES}/MENDO"
C = f"{A.RES}/MENDO_control"


def stats(run):
    o = {}
    p = f"{run}/auto_window/_params.json"
    if os.path.exists(p):
        j = json.load(open(p)); o["auto_hours"] = j["n_hours"]; o["auto_det"] = j["n_det_total"]
    p = f"{run}/clusters/summary.json"
    if os.path.exists(p):
        j = json.load(open(p)); o["n_det_culled"] = j.get("n_detections"); o["tpl0"] = j.get("n_initial_templates")
    p = f"{run}/grow/iters.csv"
    if os.path.exists(p):
        it = pd.read_csv(p).iloc[-1]; o["grown"] = int(it.n_families); o["grow_det"] = int(it.n_detections)
        o["iters"] = int(it["iter"])
    p = f"{run}/hypo_28/review/auto_score.csv"
    if os.path.exists(p):
        s = pd.read_csv(p); o["scored"] = len(s); o["real"] = int((s.n_sites_3km >= 3).sum())
        o["two_site"] = int((s.n_sites_3km == 2).sum())
        p2 = f"{run}/hypo/loc_headline.csv"
        if os.path.exists(p2):
            h = pd.read_csv(p2)
            rl = h[h.grp.isin(s[s.n_sites_3km >= 3].grp)]
            o["located_all"] = len(h); o["located_real"] = len(rl)
            o["located_real_erh5"] = int((rl.erh < 5).sum())
            o["dep_real_med"] = float(rl.dep.median()) if len(rl) else np.nan
    return o


def main():
    r, c = stats(R), stats(C)
    keys = [("auto_hours", "autocorrelation hours"), ("auto_det", "autocorrelation detections"),
            ("n_det_culled", "culled detection catalogue"), ("tpl0", "initial templates (UPGMA)"),
            ("iters", "growth iterations"), ("grown", "families grown (final iteration)"),
            ("grow_det", "detections (final iteration)"), ("scored", "families scored (2-8 Hz stacks)"),
            ("two_site", "families with exactly 2 sites"), ("real", "REAL families (>= 3 sites, < 3 km = 1 site)"),
            ("located_all", "families located (Hyp2000 headline, all)"),
            ("located_real", "real families located"), ("located_real_erh5", "real located, ERH < 5 km"),
            ("dep_real_med", "real located, median depth (km)")]
    L = ["# MENDO blind Royer-chain discovery: summary", "",
         f"Written {pd.Timestamp.utcnow():%Y-%m-%d %H:%M}Z by scripts/mendo_summary.py. Queue "
         "scripts/mendo_discovery_queue.sh, log logs/mendo_discovery_queue.log.", "",
         "Scrambled control = the same chain with each site's day series shifted by whole days "
         "(AB_NULL=1; shift within each year's list of scan days; stations < 3 km apart shift together).", "",
         "| quantity | real run | scrambled control |", "|---|---|---|"]
    for k, lab in keys:
        f = lambda d: ("pending" if k not in d else f"{d[k]:.1f}" if isinstance(d[k], float) else str(d[k]))
        L.append(f"| {lab} | {f(r)} | {f(c)} |")
    p = f"{R}/compare_shelly/families_vs_shelly.csv"
    if os.path.exists(p):
        V = pd.read_csv(p)
        loc = V.lat.notna() if "lat" in V else pd.Series(False, index=V.index)
        L += ["", "## Real families vs Shelly et al. (2025)", "",
              f"- real families: {len(V)} (located {int(loc.sum())})",
              f"- (a) epicentre within 10 km of a Shelly family: {int(V.a_loc.sum())}",
              f"- (b) waveform duplicate of a Shelly v3b it0 template: {int(V.b_wave.sum())}",
              f"- (c) shared detection times with his catalogue: {int(V.c_times.sum())}",
              f"- known (any of a/b/c): {int(V.known.sum())}",
              f"- **NEW (fail all three): {int(V.new.sum())}** (located {int((V.new & loc).sum())})", ""]
        if V.new.any():
            cols = [x for x in ["grp", "n_sites_3km", "lat", "lon", "dep", "erh", "near_fam", "near_km",
                                "c_peak_frac", "n_det"] if x in V]
            L += ["New families (approximate area = headline Hyp2000 location):", "",
                  V[V.new][cols].round(3).to_markdown(index=False)]
    else:
        L += ["", "Comparison with Shelly: pending (scripts/mendo_compare.py)."]
    open(f"{R}/MENDO_SUMMARY.md", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
