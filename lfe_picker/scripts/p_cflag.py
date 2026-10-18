#!/usr/bin/env python
"""Per-station CF lag on the DEV year: where does the learned CF actually peak relative to the
label, and how much of that is a label error rather than a model error?

The back-projection does NOT read the labels: it reads tt(cell, station) from the Slab2 grid.  What
it needs to know is where the CF actually peaks relative to THAT prediction.  Measured here: for
every Bostock event, the CF peak lag around  t_origin + tt(nearest grid cell of its family, station),
per station, on a development year.  Two systematics are folded into that one number -- the label
error at the far boreholes (P0-O: `predicted` S times, 2-3 s off at B003/B004) and the Slab2-vs-
Bostock depth offset (7.5 km median).  The per-SITE median, referenced to the network median (only
the DIFFERENTIAL delay is identifiable -- a common shift just renames the origin time), is written to
results/p2/site_delay.json and added to the back-projection travel times.  The lag vs the LABEL is
reported alongside as a diagnostic.

  python scripts/p_cflag.py --cf results/p2/cf_p2core_2012 --year 2012 --out results/p2/site_delay.json
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

FS = P.FS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cf", required=True)
    ap.add_argument("--year", type=int, default=2012)
    ap.add_argument("--half", type=float, default=3.0)
    ap.add_argument("--max-ev", type=int, default=4000)
    ap.add_argument("--out", default=f"{P.RES}/p2/site_delay.json")
    a = ap.parse_args()
    assert a.year in P.DEV_YEARS, "the CF lag may only be fitted on a development year"
    L = pd.read_parquet(f"{P.RES}/labels_v0.parquet")
    ps = L[(L.kind == "pos_S") & (L.year == a.year) & (~L.hidden_family)]
    prov_of = ps.groupby("sta").provenance.apply(lambda v: (v == "picked").mean()).to_dict()
    st = P.DayStore()
    z = np.load(f"{P.RES}/p2/grid_tt.npz")
    gs = [str(s) for s in z["stas"]]
    FG = pd.read_csv(f"{P.RES}/p2/family_on_grid.csv").set_index("bkey")
    ttlab = P.fam_tt_labels()
    B = P.bostock()
    B = B[(B.year == a.year) & (~B.hidden)]
    n = int(a.half * FS)
    rows = []
    for sta in P.S15:
        if sta not in gs:
            continue
        k = gs.index(sta)
        lags, lagl, amps = [], [], []
        for doy, g in B.groupby("doy"):
            f = f"{a.cf}/{sta}_{a.year}_{doy}.npy"
            if not os.path.exists(f) or not st.has(sta, a.year, doy):
                continue
            cf = np.load(f).astype(np.float32)
            i = st.idx(sta, a.year, doy)
            g = g.sample(min(a.max_ev // 6, len(g)), random_state=0)
            for r in g.itertuples():
                c = int(FG.cell[r.bkey])
                j = int(round((r.t_origin + z["tt"][c, k] - st.t0[i]) * FS))
                if j - n >= 0 and j + n < len(cf):
                    w = cf[j - n:j + n + 1]
                    lags.append((np.argmax(w) - n) / FS)
                    amps.append(float(w.max()))
                tl = ttlab.get((int(r.bkey), sta))
                if tl is not None:
                    j2 = int(round((r.t_origin + tl[0] - st.t0[i]) * FS))
                    if j2 - n >= 0 and j2 + n < len(cf):
                        lagl.append((np.argmax(cf[j2 - n:j2 + n + 1]) - n) / FS)
        if not lags:
            continue
        lags = np.array(lags)
        rows.append(dict(sta=sta, n=len(lags), lag_med=float(np.median(lags)),
                         lag_p25=float(np.quantile(lags, .25)), lag_p75=float(np.quantile(lags, .75)),
                         lag_vs_label_med=float(np.median(lagl)) if lagl else np.nan,
                         cf_peak_med=float(np.median(amps)),
                         prov=("picked" if prov_of.get(sta, 0) > .5 else "predicted")))
    T = pd.DataFrame(rows)
    T["site"] = T.sta.map(P.SITE_OF)
    site = T.groupby("site").apply(
        lambda d: float(np.average(d.lag_med, weights=d.n)), include_groups=False).to_dict()
    # only the DIFFERENTIAL delay matters (a common shift just moves the origin time)
    ref = float(np.median(list(site.values())))
    site = {k: round(v - ref, 3) for k, v in site.items()}
    T.to_csv(f"{P.RES}/p2/cf_lag_{a.year}.csv", index=False)
    json.dump(site, open(a.out, "w"), indent=1)
    print(T.sort_values("lag_med").to_string(index=False))
    print("\n[site delay, s, differential]", site)


if __name__ == "__main__":
    main()
