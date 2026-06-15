"""Family-level filter (the right place): cluster ALL candidates (baseline), then score each
FAMILY by its members' mean P(LFE) and drop low-scoring families. Keeps real families intact
(no member fragmentation). Sweep thresholds; report LFE retention + purity vs candidate-level filter.

Usage: python lfe_features/fam_level_filter.py --sta B011 --slat 48.65 --slon -123.448
"""
import argparse, os, sys, math
import numpy as np, pandas as pd
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
from compare_discovery import arr_lin, arr_eq, frac_near


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", required=True); ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    a = ap.parse_args(); s = a.sta.lower()
    lin = arr_lin(a.slat, a.slon); eq = arr_eq(f"data/eq_catalog_{s}.csv", a.slat, a.slon); bl = arr_eq(f"data/blast_catalog_{s}.csv", a.slat, a.slon)

    m = pd.read_parquet(f"data/{s}_disc_baseline.members.parquet")
    m["e"] = pd.to_datetime(m.time, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9
    cand = pd.read_parquet(f"data/{s}_cand_baseline.parquet")
    cand["e"] = pd.to_datetime(cand.time, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9
    pmap = dict(zip(np.round(cand.e.values, 3), cand.p_lfe.values))
    m["p"] = np.round(m.e.values, 3); m["p"] = m["p"].map(pmap)

    rows = []
    for fam, g in m.groupby("family_id"):
        t = g.e.values; pl = g.p.dropna().values
        fr = {"LFE": frac_near(t, lin), "EQ": frac_near(t, eq), "BLAST": frac_near(t, bl)}
        cls = max(fr, key=fr.get); cls = cls if fr[cls] >= 0.2 else "UNK"
        rows.append(dict(fam=fam, n=len(t), cls=cls, meanp=np.nanmean(pl) if len(pl) else np.nan,
                         medp=np.nanmedian(pl) if len(pl) else np.nan))
    d = pd.DataFrame(rows).dropna(subset=["meanp"])
    nLFE0 = (d.cls == "LFE").sum()
    print(f"baseline (clustered, family-scored): {len(d)} families, {nLFE0} Lin-confirmed LFE, "
          f"purity {100*nLFE0/len(d):.1f}%\n")
    print("FAMILY-LEVEL filter (keep families with mean member P(LFE) >= thr):")
    print(f"{'thr':>5}{'families':>10}{'LFE':>7}{'retain%':>9}{'purity%':>9}{'contam_removed%':>16}")
    cont0 = len(d) - nLFE0
    for thr in [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
        k = d[d.meanp >= thr]
        nL = (k.cls == "LFE").sum(); cont = len(k) - nL
        print(f"{thr:>5.1f}{len(k):>10}{nL:>7}{100*nL/max(nLFE0,1):>9.0f}{100*nL/max(len(k),1):>9.1f}"
              f"{100*(1-cont/max(cont0,1)):>16.0f}")
    print("\n(compare: candidate-level filter gave LFE retain 19%, purity 38.5%)")
    d.sort_values("meanp", ascending=False).to_csv(f"lfe_features/data/fam_level_{s}.csv", index=False)
    print("FAM LEVEL DONE")


if __name__ == "__main__":
    main()
