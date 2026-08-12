"""Compare BASELINE vs FILTERED discovery families by ground truth.

Classify each discovered family by whether its member times coincide (+/-6 s) with Lin LFEs /
ANSS earthquakes / ANSS blasts. Show whether the tremor-picker filter removes contaminant
(EQ/blast/uncatalogued) families while retaining the Lin-confirmed LFE families.

Usage: python lfe_features/compare_discovery.py --sta B011 --slat 48.65 --slon -123.448
"""
import argparse, os, sys, math
import numpy as np, pandas as pd
from collections import Counter
TOL = 6.0


def arr_lin(slat, slon):
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dk = np.sqrt(((lin.lat - slat) * 111) ** 2 + ((lin.lon - slon) * 111 * math.cos(math.radians(slat))) ** 2)
    lin = lin[dk < 35]
    return np.sort(pd.to_datetime(lin.OT, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9 + 11.0)


def arr_eq(csv, slat, slon):
    d = pd.read_csv(csv); ep = np.sqrt(((d.lat - slat) * 111) ** 2 + ((d.lon - slon) * 111 * math.cos(math.radians(slat))) ** 2)
    return np.sort(pd.to_datetime(d.ot, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9 +
                   np.sqrt(ep ** 2 + np.clip(d.depth, 0, 60) ** 2) / 3.6)


def frac_near(times, arr):
    if len(arr) == 0 or len(times) == 0:
        return 0.0
    k = np.clip(np.searchsorted(arr, times), 1, len(arr) - 1)
    d = np.minimum(np.abs(times - arr[k - 1]), np.abs(times - arr[k]))
    return float((d <= TOL).mean())


def classify(memb_parquet, lin, eq, bl):
    m = pd.read_parquet(memb_parquet)
    m["e"] = pd.to_datetime(m.time, utc=True).values.astype("datetime64[ns]").astype("int64") / 1e9
    rows = []
    for fam, g in m.groupby("family_id"):
        t = g.e.values
        fl, fe, fb = frac_near(t, lin), frac_near(t, eq), frac_near(t, bl)
        fr = {"LFE": fl, "EQ": fe, "BLAST": fb}
        cls = max(fr, key=fr.get)
        cls = cls if fr[cls] >= 0.2 else "UNK"
        rows.append(dict(fam=fam, n=len(t), lfe=fl, eq=fe, blast=fb, cls=cls))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", required=True); ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    a = ap.parse_args(); s = a.sta.lower()
    lin = arr_lin(a.slat, a.slon); eq = arr_eq(f"data/eq_catalog_{s}.csv", a.slat, a.slon); bl = arr_eq(f"data/blast_catalog_{s}.csv", a.slat, a.slon)
    print(f"ground truth: Lin {len(lin)}, EQ {len(eq)}, BLAST {len(bl)} arrivals (within 35km)\n")

    out = {}
    for tag, mp in [("BASELINE", f"data/{s}_disc_baseline.members.parquet"),
                    ("FILTERED", f"data/{s}_disc_filtered.members.parquet")]:
        d = classify(mp, lin, eq, bl); out[tag] = d
        c = Counter(d.cls)
        nfam = len(d); lfe = c["LFE"]
        print(f"=== {tag}: {nfam} families ===")
        print(f"   by type: LFE={c['LFE']}  EQ={c['EQ']}  BLAST={c['BLAST']}  UNK={c['UNK']}")
        print(f"   LFE-confirmed PURITY = {100*lfe/nfam:.1f}%  (contaminant+unk families = {nfam-lfe})")
        print()
    b, f = out["BASELINE"], out["FILTERED"]
    nb, nf = (b.cls == "LFE").sum(), (f.cls == "LFE").sum()
    cb, cf = (b.cls != "LFE").sum(), (f.cls != "LFE").sum()
    print("=== VERDICT ===")
    print(f"   Lin-confirmed LFE families:   baseline {nb} -> filtered {nf}   (retained {100*nf/max(nb,1):.0f}%)")
    print(f"   contaminant/UNK families:     baseline {cb} -> filtered {cf}   (removed {100*(1-cf/max(cb,1)):.0f}%)")
    print(f"   purity (LFE fraction):        baseline {100*nb/len(b):.1f}%  ->  filtered {100*nf/len(f):.1f}%")
    f.sort_values("lfe", ascending=False).to_csv(f"lfe_features/data/disc_filtered_classified_{s}.csv", index=False)
    print("COMPARE DONE")


if __name__ == "__main__":
    main()
