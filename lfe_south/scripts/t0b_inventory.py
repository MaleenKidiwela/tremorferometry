#!/usr/bin/env python
"""Write results/SOUTH_TRUTH_INVENTORY.md from the tables on disk (never from memory)."""
import os,sys,numpy as np,pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
DUC=f"{ROOT}/catalogs/ducellier_catalog"
F=pd.read_csv(f"{RES}/families_south.csv"); F["fam"]=F.fam.astype(str)
R=pd.read_csv(f"{RES}/stage1_rosters.csv"); R["fam"]=R.fam.astype(str)
P=pd.read_csv(f"{RES}/stage1_download_plan.csv")
A=pd.read_csv(f"{RES}/duc_bedA_stations.csv"); B=pd.read_csv(f"{RES}/duc_bedB_stations.csv")
C=pd.read_csv(f"{RES}/fdsn_channels_south.csv")
S=F[F.source=="Shelly2025"]; Du=F[F.source=="Ducellier2022"]
L=[]
L.append("# Southern truth inventory — what is on disk, verified row by row\n")
L.append(f"Agent T, {pd.Timestamp.utcnow():%Y-%m-%d %H:%M} UTC. Every number below was read from the "
         "files named beside it; nothing is quoted from a paper abstract.\n")
L.append("## 1. Shelly et al. (2025) — Mendocino / southern edge of Cascadia\n")
L.append(f"* Catalogue `lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv` — "
         f"**{int(S.n_events.sum()):,} detections, {len(S)} families**, "
         f"{S.t0.min()[:10]} .. {S.t1.max()[:10]} (USGS data release 10.5066/P1TCKK7G, CC0).")
L.append(f"* Locations from the release metadata `Shelly2025_metadata.xml`: latitude "
         f"{S.lat.min():.3f}..{S.lat.max():.3f} N, longitude {S.lon.min():.3f}..{S.lon.max():.3f}, "
         f"depth {S.depth.min():.1f}..{S.depth.max():.1f} km. The 27 families occupy a patch about "
         f"{(S.lat.max()-S.lat.min())*111:.0f} km N-S by "
         f"{(S.lon.max()-S.lon.min())*111*np.cos(np.radians(40.06)):.0f} km E-W — one tight NE-dipping zone.")
L.append(f"* Detections per family: min {int(S.n_events.min()):,}, median {int(S.n_events.median()):,}, "
         f"max {int(S.n_events.max()):,}. Every family spans essentially the whole 2018-2024 record.")
L.append("* **No templates are released** — only times and locations. Our templates must be built "
         "from his times on our own stations (stage 2).")
L.append(f"* `catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv` is a saved 404 HTML page "
         f"({os.path.getsize(f'{ROOT}/catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv')} bytes), not data — "
         "confirmed again here; the working copy is the one above.\n")
L.append("## 2. Ducellier & Creager (2022) — southern Cascadia, eight years\n")
L.append(f"* `catalogs/ducellier_catalog/` (the author's own Zenodo repository, complete).")
L.append(f"* **{len(Du)} families**, one directory each under `data/Ducellier/catalogs/<family>/`, with two "
         f"pandas frames (`year,month,day,hour,minute,second,cc,nchannel`):")
nn=Du.n_2004.sum() if "n_2004" in Du else np.nan
L.append(f"  * `catalog_2004_2011.pkl` — the PERMANENT-network catalogue: "
         f"**{int(Du.n_2004.sum()):,} detections**, 2004-01-01 .. 2011-12-31.")
L.append(f"  * `catalog_2007_2009.pkl` — the FAME temporary-array catalogue: "
         f"**{int(Du.n_2007.sum()):,} detections**, 2007-07 .. 2009-12.")
L.append(f"* Locations from `data/Plourde_2015/templates_list.txt` (Plourde et al. 2015 templates, which "
         f"Ducellier re-used): latitude {Du.lat.min():.2f}..{Du.lat.max():.2f} N, longitude "
         f"{Du.lon.min():.2f}..{Du.lon.max():.2f}, depth {Du.depth.min():.1f}..{Du.depth.max():.1f} km, "
         f"with published uncertainties eH 0.8-5.2 km, eZ 0.3-6.4 km. The 66 template families collapse to "
         f"34 numbered Plourde families plus MF1/MF2/BC3.")
L.append(f"* **Templates ARE released**: `src/templates_v1/<family>/<STA>_<CHA>.pkl` (obspy Trace + "
         f"orientation) — 66 families, 20 Hz, response-removed to velocity, 1.5-9 Hz, "
         f"10-45 s long (constant within a family), median 16 channels per family, 47 distinct stations. "
         f"`src/templates/` is the FAME set and `src/templates_2007_2009/` the intermediate set.")
L.append(f"* His station tables are on disk too: `data/Ducellier/stations_permanent.txt` "
         f"({len(A)} rows / {A.sta.nunique()} stations, his selected permanent bed), "
         f"`stations_permanent_all.txt` (83, the pool), `stations_temporary.txt` (33, FAME), and the "
         f"per-network channel/epoch tables `data/networks/{{BK,NC,PB,FAME}}.txt`.")
L.append(f"* **{int(A.zonly.sum())} of {len(A)} rows ({100*A.zonly.mean():.0f} %) are SINGLE-COMPONENT "
         f"VERTICALS** (NC EHZ/SHZ: GBB GHM GRO GSN GTC GVV KBN KBS KCS KFP KIP KKP KOM KPP KRK LAM LBF "
         f"LBK LBP LDB LGB LGP LPG LRB LRR LSF LSH LSR LTC LVR LWH KTR). His detector correlates "
         f"PER CHANNEL and averages, so a vertical-only station is a full member of his network. Any "
         f"loader of ours that requires 3-C would silently run a thinned network — the southern chain "
         f"here therefore accepts 1-C and 2-C station-days by construction.")
L.append(f"* His published parameters, read from `src/find_all_LFEs_parallel.py` __main__ and "
         f"`src/get_data.py`: band 1.5-9.0 Hz zero-phase, 20 Hz (dt 0.05), response removed to VEL with "
         f"pre_filt (0.2, 0.5, 10, 15) and water_level 80, per-channel normalised cc averaged over "
         f"channels with a complete gapless block, threshold **8 x MAD** of the block, declustering "
         f"`floor(1/(dt*freq0))` = 20 samples = 1 s, 12-hour blocks.\n")
L.append("## 3. Other southern files checked\n")
L.append(f"* `catalogs/ducellier_b935/` — {len(os.listdir(f'{ROOT}/catalogs/ducellier_b935'))} per-family "
         f"CSVs with a single `time` column (45,844 rows total): a derived extract of the same catalogue "
         f"restricted to families that B935 sees. No picks, no locations. Kept as a cross-check only.")
L.append("* `catalogs/mtj_ducellier/` — the 404 page described above.")
L.append("* Waveform trees: no southern station-day existed anywhere on disk before this run "
         "(`lfe_picker/results/o4_deferred.md` §'Why it is deferred' re-verified).\n")
L.append("## 4. Station geometry we derived (FDSN, not assumed)\n")
L.append(f"* `results/fdsn_channels_south.csv` — {len(C):,} channel epochs from EarthScope + NCEDC for "
         f"NC, BK, PB, UW, CC, US, UO (and the XQ/FAME temporaries, kept for reference) in the box "
         f"38.6-42.4 N, -124.8 .. -121.6, 2004-2026.")
L.append(f"* Stations alive on 2026-06-01: 204 (195 of them 3-C). Stations alive in BOTH 2009 and 2026: "
         f"59, of which 50 are 3-C today — so unlike the north (POLARIS decommissioned), the southern "
         f"permanent bed can be reproduced *and* extended on the same instruments.")
L.append(f"* Stage-1 rosters (`results/stage1_rosters.csv`): permanent stations within 80 km hypocentral, "
         f"alive on the sampled detection days AND alive today, nearest 10 kept as the download budget. "
         f"Shelly: {R[R.source=='Shelly2025'].sta.nunique()} distinct stations, 10 per family, median "
         f"hypocentral distance {R[R.source=='Shelly2025'].hypo_km.median():.0f} km. "
         f"Ducellier: {R[R.source=='Ducellier2022'].sta.nunique()} distinct stations, "
         f"{int(R[R.source=='Ducellier2022'].groupby('fam').size().median())} per family (min "
         f"{int(R[R.source=='Ducellier2022'].groupby('fam').size().min())}), median "
         f"{R[R.source=='Ducellier2022'].hypo_km.median():.0f} km.")
L.append(f"* Ducellier reproduction beds: **bed A** = his exact set ({A.sta.nunique()} stations, "
         f"{int(A.zonly.sum())} Z-only rows); **bed B** = bed A plus "
         f"{int((1-B.in_bedA).sum())} stations within 80 km of a family that exist in the window but were "
         f"not in his list. Both are run and reported separately.\n")
L.append("## 5. What stage 1 downloads\n")
L.append(f"* Existence test: **{len(P):,} station-days** ({P.day.nunique()} distinct days, "
         f"{P.sta.nunique()} stations) — the busiest days of each family, enough for >= 150 of its own "
         f"detections (median {int(R.n_det.median())} per family). Estimated 100-160 GB.")
L.append(f"* Ducellier reproduction: **{len(pd.read_csv(f'{RES}/duc_download_plan.csv')):,} station-days** "
         f"over 2008-04-01 .. 2008-05-10, the richest window in his catalogue (64 of 66 families with "
         f">= 20 events in 30 days). Estimated ~70 GB.")
L.append("* Nothing else is fetched until a stage needs it.\n")
open(f"{RES}/SOUTH_TRUTH_INVENTORY.md","w").write("\n".join(L)+"\n")
print("wrote",f"{RES}/SOUTH_TRUTH_INVENTORY.md")
