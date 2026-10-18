#!/usr/bin/env python
"""Write results/WA_2010_ARRAY_RESULTS.md from whatever the queue has produced.  Prefix `o_`."""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import o_common as O                                                        # noqa: E402

R = O.OUTROOT


def jl(p, d=None):
    try:
        return json.load(open(p))
    except Exception:
        return d


def cl(p):
    try:
        return pd.read_csv(p)
    except Exception:
        return None


def md_table(df, cols=None, n=None, fl=3):
    if df is None or not len(df):
        return "_(none)_\n"
    d = df[cols] if cols else df
    if n:
        d = d.head(n)
    d = d.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda x: "" if pd.isna(x) else f"{x:.{fl}g}")
    head = "| " + " | ".join(str(c) for c in d.columns) + " |"
    sep = "|" + "|".join(["---"] * len(d.columns)) + "|"
    body = "\n".join("| " + " | ".join(str(v) for v in r) + " |" for r in d.values)
    return f"{head}\n{sep}\n{body}\n"


SEC0 = """\
## 0. What the papers actually did with the arrays (fetched and quoted)

**Chestler & Creager (2017), *JGR* 122, 3099 (10.1002/2016JB013717)** -- the 43-family Olympic
catalogue.  The methods chapter is reproduced verbatim in Chestler's UW dissertation *Using
Low-Frequency Earthquakes to Investigate Slow Slip* (sec. 2.1, open access, fetched 2026-09-21):

> "We use high signal-to-noise (SNR) level waveforms, calculated by Ghosh et al. (2012) using
> multibeam-backprojection (MBBP), to create a catalog of low-frequency earthquakes.  The data
> were collected during the Array of Arrays (AofA) experiment, which consisted of eight 1-km
> aperture arrays, each with 10-20 3-component stations ... For each array and every 1-minute
> time window Ghosh et al. (2012) filtered vertical-component seismograms between 5-9 Hz.  Then
> they apply a **delay-and-sum beamforming** approach to stack data and find the **vector
> slowness** corresponding to the direction of a plane wave that produces the stack with the most
> power.  Stacking seismograms from each window at the optimal slowness produces a continuous
> time series with excellent SNR.  This is repeated for each array and every 50-percent
> overlapping one-minute time window."

> "We use these high SNR waveforms to search for new LFE families.  For each array, j, we select
> all the peaks that are larger than **18 times the median absolute deviation** ... we create a
> **boxcar-like time series** for each array [+-1 s] ... Finally, we stack the time series for
> all arrays and use the peaks in the resulting function as potential LFEs."

> "To verify that potential LFEs repeat in time, using a set of **3 stations from each array**,
> we autocorrelate **15-second windows** of 3-component data, band-passed between **2 and 8 Hz**,
> centered on each peak through time.  We stack the autocorrelation functions for all channels
> and stations using an **nth-root stack where n=2** and look for times when the stacked
> autocorrelation function exceeds **45 times the median absolute deviation**.  For each
> potential LFE, up to **150** of the time windows that correlate best with the initial detection
> are stacked to create a template waveform ... After **three iterations**, we group the original
> LFE and all the correlating time windows into an LFE family, imposing a stricter correlation
> threshold that varies from 0.05 to 0.20 depending on the number of stations that are active."

> "To locate LFE families we first produce P- and S-wave differential times by cross correlating
> stacked templates across different stations and S minus P times by cross correlating S- (on
> horizontal channels) versus P- (on vertical channels) waves at the same stations.  We then
> perform a grid search using the **3D velocity model from Preston et al. (2003)** ... Finally we
> improve relative locations between nearby families using the **double difference** method of
> Waldhauser and Ellsworth (2000)."

> "We assess the similarity of two LFE families by the ratio R_com of the number of common
> detections ... When two LFE families have **R_com > 0.5**, the family with the lower median
> correlation coefficient is removed ... This reduced the number of LFE families from **57 to
> 43**."

So the published Olympic chain is, in order: (i) a **slowness-grid delay-and-sum BEAM per 1-km
array** (Ghosh et al. 2012, vertical component, 5-9 Hz, 1-minute windows, 50 % overlap);
(ii) 18 x MAD peaks on each beam, combined across arrays by a +-1 s boxcar stack, as the
CANDIDATE LFE times; (iii) a 15-s, 2-8 Hz, 3-component matched filter seeded by those windows on
**3 stations per array**, network-combined by a **2nd-root stack** at **45 x MAD**; (iv) templates
from the best 150 windows, 3 iterations; (v) R_com > 0.5 de-duplication; (vi) location by
differential P/S times in a 3-D model, then double difference.

**What this study can and cannot copy.**  Steps (i)-(iv) are copied here (sec. 1 and sec. 6).
Step (vi) is NOT: this project has no 3-D Olympic model and no P picks on LFE stacks, so the
location stage stays the study's own -- stack-and-pick S on the Crosson 1-D model with intraslab
station terms, a 50-realisation shuffled-station control and the wide-grid `edge` flag
(`n_locate.py`).  Every family here is therefore located less precisely than Chestler's, which is
why the Chestler check (sec. 5) is run at 5 km and 10 km, not at 1 km.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{O.RES}/WA_2010_ARRAY_RESULTS.md")
    a = ap.parse_args()
    bed = jl(f"{R}/bed_2010.json", {})
    bed11 = jl(f"{R}/bed_2011.json", {})
    stat = jl(f"{R}/statics_eq.json", {})
    smoke = jl(f"{R}/smoke_beam.json", {})
    band = cl(f"{R}/array_coherence_bands.csv")
    modes = jl(f"{R}/mode_smoke.json", {})
    L = []
    w = L.append
    w("# Washington N-Olympic 2010 -- the arrays run AS ARRAYS, scored against Chestler & Creager\n")
    w("Agent O, 2026-09-21.  PLAN.md sec.17: the goal is to **reproduce the published Washington "
      "LFE families**, so the 2010 N-Olympic null of `WA_2010_RESULTS.md` / "
      "`WA_2011_REFIRE_RESULTS.md` is treated as OUR reproduction failure, to be debugged the way "
      "the north was.  The prime suspect was the site treatment: 84 of the 96 admitted stations "
      "sit in ten clusters <= 1.7 km across, and the 2010 run gave every one of them an "
      "independent vote in the network correlation sum.  Chestler & Creager (2017), Ghosh et al. "
      "(2012) and Sweet et al. (2014) used the arrays AS ARRAYS.\n")
    w(SEC0)
    w("\n### 0b. The sites, defined by geometry\n")
    w("Single-link clustering of the 96 admitted stations at 1.5 km, partitioned by network "
      "(a borehole and a surface sensor 100 m apart are not one array).  Every (cluster, network) "
      "sub-group with >= 2 members is ONE site; singletons stay single stations.  "
      f"**{bed.get('n_member_stations','?')} stations -> {bed.get('n_sites','?')} sites "
      f"({bed.get('n_array_sites','?')} arrays covering {bed.get('n_members_in_arrays','?')} "
      f"stations + {bed.get('n_single_sites','?')} single stations).**\n")
    if bed.get("arrays"):
        w(md_table(pd.DataFrame(bed["arrays"]),
                   ["site", "net", "n", "aperture_km", "lat", "lon", "members"]))
        w(f"\nSingle stations: {', '.join(bed.get('singles', []))}\n")
    w("\n### 0c. Two controls that changed the beam, both measured on this bed\n")
    if smoke:
        w("**(i) The slowness search is not identifiable at this aperture.**  Ghosh's beam picks, "
          "per 1-minute window, the slowness that maximises 5-8 Hz stack power.  Run here, it "
          "raises the median semblance by 1.3-2.8x over the zero-slowness stack -- but it raises "
          "it by the SAME factor on the hour-matched QUIET hours, where there is no tremor to "
          "beam at.  The gain is the bias of taking a maximum over the slowness grid, not signal: "
          "a 1-km array at 6.5 Hz has a main lobe ~1/(fD) = 0.15 s/km wide, i.e. about seven "
          "resolution cells in the whole search.  Measured (median over 4 tremor and 4 quiet "
          "hours per array):\n")
        w(md_table(pd.DataFrame(smoke.get("slowness_control", [])),
                   ["site", "M", "floor", "tremor_zero", "tremor_best", "tremor_s",
                    "quiet_zero", "quiet_best", "quiet_s"]))
        w("\nThe zero-slowness semblance sits at the incoherent floor 1/M in BOTH populations.  "
          "A per-window slowness would therefore inject a random +-0.1-0.3 s jitter into the "
          "array trace -- fatal for a matched filter, which needs the beam waveform to be the "
          "same function of time every day.  **Adopted: the zero-slowness delay-and-sum beam** "
          "(the equal-weight mean of the rms-normalised members), which is the same estimator at "
          "the only slowness the data support.\n")
    if band is not None:
        w("\n**(ii) How coherent IS a 1-km array in the LFE band?**  Catalogued intraslab "
          "earthquakes (25-46 km deep, inside the box) are real deep sources at the LFE "
          "incidence.  Coherence N x <stack^2>/<trace^2> of their S window across each array "
          "(1 = incoherent, M = perfect), with the pre-event noise window as the control:\n")
        w(md_table(band, fl=3))
        w("\nSo the arrays are only ~25 % coherent at 1-2 Hz and essentially incoherent above "
          "4 Hz.  **This is the measurement that matters for the site treatment**: beam-stacking "
          "an array destroys most of the per-member information about the signal, while the 2010 "
          "run's 13-votes-per-cluster keeps it but also lets the cluster's SHARED noise vote 13 "
          "times.  Per-member statics from the same earthquakes do not cross-validate "
          "(leave-one-event-out coherence is not improved), so no statics are applied.\n")
    if stat.get("summary"):
        w("\n" + md_table(pd.DataFrame(stat["summary"]).T.reset_index().rename(
            columns={"index": "site"}),
            ["site", "n_events", "n_members_with_static", "static_sd_ms", "coh_plain",
             "coh_static_LOO", "coh_noise_plain", "coh_noise_static_LOO"]))
    if modes:
        w("\n### 0d. The site treatment, chosen by measurement (not by assumption)\n")
        w("The identical literal autocorrelation (`k_auto`, 6 s, 1-8 Hz, 8 x MAD of the block, "
          "|lag| >= 6 s, 4 s decluster) run three ways on the same tremor and quiet hours:\n")
        w(md_table(pd.DataFrame(modes.get("table", [])), fl=4))
        w(f"\n**Chosen: `{modes.get('chosen','?')}`** -- {modes.get('why','')}\n")

    # ---- pair stage --------------------------------------------------------------------------
    ps = jl(f"{R}/2010/pair_pairstats.json")
    if ps:
        w("\n## 1. Pair stage\n")
        w("```json\n" + json.dumps(ps, indent=1)[:4000] + "\n```\n")
    # ---- growth / screen ---------------------------------------------------------------------
    for tag, nm in (("2010", "real"), ("2010_null", "site-shift null")):
        sc = jl(f"{R}/{tag}/screen_summary.json")
        sj = jl(f"{R}/{tag}/screen_final_summary.json")
        if sc:
            w(f"\n## 2. Screen -- {nm} (`{tag}`)\n")
            w(md_table(pd.DataFrame(sc.get("clause_table", [])),
                       ["clause", "name", "alone", "cumulative"]))
            w(f"\nquiet-day network FA {sc.get('network_fa_per_tpl_day')}/template-day; "
              f"tremor-day rate {sc.get('network_rate_tremor_per_tpl_day')}/template-day; "
              f"pass (a)-(d) **{sc.get('n_pass_abcd')}**")
            if sj:
                w(f"; pass (e) stability **{sj.get('n_pass_e')}**; CANDIDATES "
                  f"**{sj.get('n_candidate')}**")
            w("\n")
    V = cl(f"{R}/2010/screen_final.csv")
    if V is None:
        V = cl(f"{R}/2010/screen.csv")
    if V is not None and "candidate" in V.columns and V.candidate.astype(bool).any():
        w("\n### Candidate families (real chain)\n")
        w(md_table(V[V.candidate.astype(bool)].sort_values("n_det", ascending=False),
                   ["grp", "lat", "lon", "dep", "slab2", "dz", "n_sta", "n_det", "n_days",
                    "fa_quiet_per_tpl_day", "co_timing", "co_base", "h2_km"]))
    # ---- Chestler ----------------------------------------------------------------------------
    ch = jl(f"{R}/chestler_check.json")
    if ch:
        w("\n## 5. The Chestler check\n")
        w("**A. screened families -> nearest published Olympic coordinate**\n\n```json\n"
          + json.dumps(ch["A"], indent=1) + "\n```\n")
        w("\n**B. where are THEIR families in OUR output** (any interior-located grown template "
          "within 10 km of each of the 69 derived coordinates; the site-shift null answers the "
          "same question and is the false-hit rate)\n\n```json\n"
          + json.dumps(ch["B"], indent=1) + "\n```\n")
        cov = cl(f"{R}/chestler_check_coverage_real.csv")
        if cov is not None:
            w("\n" + md_table(cov.sort_values("nearest_km"),
                              ["row", "lat", "lon", "dep", "nearest_grp", "nearest_km",
                               "hit_10km", "n_within_10km"]))
    # ---- 2011 --------------------------------------------------------------------------------
    cf = jl(f"{R}/2011/confirm_real_summary.json")
    cn = jl(f"{R}/2011/confirm_null_summary.json")
    if cf or cn:
        w("\n## 6. Confirmation on the sealed 2011 window\n")
        w("Rule, fixed before the 2011 numbers existed and identical to `WA_2011_REFIRE_RESULTS` "
          "sec.5: >= 3 distinct 2011 tremor days, quiet-day FA <= 0.1/template-day on the 7 "
          "strictly-quiet days, and an INTERIOR 2011 location within 10 km (3-D) of the 2010 "
          "location.  The null templates go through byte-identical code: their confirmation "
          "count is the false-confirmation rate.\n")
        for nm, s in (("real", cf), ("null", cn)):
            if s:
                w(f"\n**{nm}**: ```{json.dumps(s)}```\n")
        for nm in ("real", "null"):
            T = cl(f"{R}/2011/confirm_{nm}.csv")
            if T is not None and len(T):
                w(f"\n{nm}:\n" + md_table(T, None, 60))
    # ---- arm B -------------------------------------------------------------------------------
    cb = jl(f"{R}/chestler_arm/summary.json")
    if cb:
        w("\n## 7. Arm B -- Chestler's OWN detection chain, run literally\n")
        w("```json\n" + json.dumps(cb, indent=1)[:6000] + "\n```\n")
    w("\n## 8. Files\n")
    w(md_table(pd.DataFrame([
        dict(what="site definition", path=f"{R}/sites.csv"),
        dict(what="array coherence (intraslab EQ)", path=f"{R}/eq_statics_events.csv, "
             f"eq_statics_validation.csv, array_coherence_bands.csv"),
        dict(what="beam / slowness controls", path=f"{R}/smoke_beam.json"),
        dict(what="site-treatment smoke", path=f"{R}/mode_smoke.json"),
        dict(what="beam traces", path=f"{O.BEAM}/<SITE>/<year>/<doy>.npz"),
        dict(what="real chain", path=f"{R}/2010/"),
        dict(what="site-shift null chain", path=f"{R}/2010_null/"),
        dict(what="Chestler check", path=f"{R}/chestler_check*.csv|json"),
        dict(what="2011 confirmation", path=f"{R}/2011/"),
        dict(what="Arm B (Chestler chain)", path=f"{R}/chestler_arm/"),
    ])))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write("\n".join(L))
    print(f"wrote {a.out} ({sum(len(x) for x in L)} chars)", flush=True)


if __name__ == "__main__":
    main()
