#!/usr/bin/env python
"""Avenue A6: score every earlier located Washington family table (real AND its null) against the
digitized Royer & Bostock (2014) Hyp2000 templates with ab_truth_score.score (unchanged).
No row is filtered except to 'located' (a row without lat/lon cannot be scored); the dedupe-kept
subset is reported as a second line, the full located set is the headline.
-> results/a4_seed/A6_prior_runs_vs_truth.{csv,md}"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ab_truth_score as TS
R = os.path.join(os.path.dirname(HERE), "results")
OUT = f"{R}/a4_seed"
TABLES = [  # (bed, chain, path, locator)
    ("2010 N-Olympic arrays (agent O / sf_final)", "real", "sf_final/table_oly2010_real.csv", "sf grid"),
    ("2010 N-Olympic arrays (agent O / sf_final)", "null", "sf_final/table_oly2010_null.csv", "sf grid"),
    ("2010 N-Olympic, u_locate L9 relocation", "real", "u_locate/wa/wa_real_ranked.csv", "L9_xcdiff"),
    ("2010 N-Olympic, u_locate L9 relocation", "null", "u_locate/wa/wa_null_ranked.csv", "L9_xcdiff"),
    ("2010 Arm B Chestler chain (sf_final)", "real", "sf_final/armb_real_loc.csv", "sf grid"),
    ("2010 Arm B Chestler chain (sf_final)", "null", "sf_final/armb_null_loc.csv", "sf grid"),
    ("2010 Arm B, u_locate L9", "real", "u_locate/wa/wa_armB_ranked.csv", "L9_xcdiff"),
    ("2008 Puget/S-Olympic (agent N)", "real", "sf_final/puget2008_real_accept.csv", "sf grid"),
    ("2008 Puget/S-Olympic (agent N)", "null", "sf_final/puget2008_null_accept.csv", "sf grid"),
    ("2010 WA (agent L, l_wa2010)", "real", "l_wa2010/loc_it0.csv", "L grid"),
    ("2010 WA (agent L, l_wa2010)", "null", "l_wa2010_null/loc_it0.csv", "L grid"),
    ("2010 refire (m_2011)", "real", "m_2011/loc2010_real.csv", "L grid"),
    ("2010 refire (m_2011)", "null", "m_2011/loc2010_null.csv", "L grid"),
    ("2011 refire (m_2011)", "real", "m_2011/loc2011_real.csv", "L grid"),
    ("2011 refire (m_2011)", "null", "m_2011/loc2011_null.csv", "L grid"),
]
truth = pd.read_csv(TS.TRUTH); truth = truth[truth.locator == "Hyp2000"]
rows = []
for bed, chain, p, loc in TABLES:
    fp = f"{R}/{p}"
    if not os.path.exists(fp):
        print("missing", fp); continue
    T = pd.read_csv(fp)
    pass

    T = T[T.get("located", True).astype(bool) & T.lat.notna()].rename(columns={"dep": "depth"})
    subsets = [("all located", T)]
    for c in ("dedupe_kept", "kept_dedupe"):
        if c in T:
            subsets.append(("dedupe-kept", T[T[c].astype(bool)])); break
    for sname, S in subsets:
        res = TS.score(truth, S[["lon", "lat", "depth"]].reset_index(drop=True), ["NW", "SW"], 1000, 0)
        for r in res.to_dict("records"):
            rows.append(dict(bed=bed, chain=chain, locator=loc, subset=sname, n_located=len(S), **r))
D = pd.DataFrame(rows)
D.to_csv(f"{OUT}/A6_prior_runs_vs_truth.csv", index=False)
cols = ["bed", "chain", "subset", "subarray", "n_ours_in_region", "median_match_km",
        "frac_3km", "null_mean_3km", "p_3km", "frac_5km", "null_mean_5km", "null_p95_5km", "p_5km",
        "frac_10km", "null_mean_10km", "p_10km", "n_depth_pairs", "dz_median_km"]
def md(df):
    h = "| " + " | ".join(df.columns) + " |\n|" + "---|" * len(df.columns) + "\n"
    return h + "\n".join("| " + " | ".join(f"{v:.3g}" if isinstance(v, float) else str(v) for v in r) + " |"
                         for r in df.itertuples(index=False))
txt = f"""# A6 -- earlier located Washington family tables vs the Royer & Bostock (2014) Hyp2000 templates

Script `scripts/a6_prior_vs_truth.py`; scorer `scripts/ab_truth_score.py` (unchanged: one-to-one
Hungarian match on epicentral distance, candidates = our families within 30 km of the subarray's
templates, null = our candidate cloud randomly rotated + shifted, 1000 draws).  Truth = Fig 4 Hyp2000
rows (NW 122, SW 54).  `frac_Xkm` = fraction of TRUTH templates matched one-to-one within X km.
Rows: every located row of each table ("all located"), plus the dedupe-kept subset where the table
has one.  The chain's own site-shift null tables get the identical scoring.
Caveat: a dense cloud of hundreds of our points inside the 30 km region matches many truth templates
by chance, which the rotation null measures -- read `frac` against `null_mean`/`p`, and real chain
against null chain.

{md(D[cols])}
"""
open(f"{OUT}/A6_prior_runs_vs_truth.md", "w").write(txt)
print(D[cols].round(3).to_string(index=False))
