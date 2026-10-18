# S5 setup: tidal stress on the plate interface at the 115 SVI LFE families

Built 2026-10-07, about 1:30 pm PT. This folder only provides the tidal stress. **No tidal-sensitivity analysis has been run.**
PLAN.md Revision 1 point 6 asks for ocean loading from a tide package plus an ocean model. That is what this folder delivers.

## Main product
`tidal_stress_svi_hourly.parquet` holds 12,850,560 rows: 115 families × 111,744 hours, from 2014-01-01 00:00 to 2026-09-30 23:00 **UTC**. The file is 0.86 GB and stored as float32.

| column | meaning |
|---|---|
| `fid` | family id, as in `results/s0/families.csv` |
| `time` | UTC, hourly, datetime64[ns] |
| `tau_Pa` | total shear traction in the slip direction: the 236° NA-vs-JdF convergence azimuth projected onto the Slab2 plane. **Positive means it promotes thrust slip.** |
| `sigma_n_Pa` | total normal traction, **tension (unclamping) positive** |
| `coulomb_0p1`, `coulomb_0p4` | τ + μ′σn for μ′ = 0.1 and μ′ = 0.4 |
| `body_tau_Pa`, `body_sigma_n_Pa` | solid-earth body-tide part |
| `ocean_tau_Pa`, `ocean_sigma_n_Pa` | ocean-load part |
| `tau_updip_Pa`, `tau_strike_Pa` (+ `body_`/`ocean_` parts) | shear for pure up-dip slip (rake 90) and for strike slip (rake 0, Aki-Richards strike) |

For any other rake r, compute τ(r) = cos r · τ_strike + sin r · τ_updip.
The 236° direction gives rakes of 74–104° (median 83°). The slip-model rakes in `lfe_extension/.../slip_svi_family_meta.csv` have a median of 71° and are carried as a column in `family_fault_geometry.csv`.

Other files:
- `tidal_constituents_per_family.csv`: harmonic fit over the whole record for each family and each series (total, body, ocean; τ, σn, Coulomb 0.4). It gives the amplitude in Pa and the Greenwich phase lag G, where x = A·f·cos(V+u−G). It also gives the std and the 99 % peak-to-peak range.
- `tidal_constituents_median_summary.csv`: the median of those values across families.
- `family_fault_geometry.csv`: strike and dip come from a plane fitted to Slab2 within 10 km. Dip is 7.7–13.4° and strike 310–340°. The file also holds the normal and slip unit vectors (ENU) and both depths.
- `ocean_constituent_stress.csv`: complex constituent stresses for every variant (model, fill, depth, zone).
- `body_stress_hourly.parquet`, `ocean_stress_hourly.parquet`: intermediate files. The body file also holds the horizontal stress tensor and areal strain.

## Methods
**Body tide** (`scripts/s5_tides_body.py`)
- Displacements come from pyTMD 3.0.9 `SET_displacements`, which follows IERS 2010. That includes degrees 2 and 3, frequency-dependent diurnal Love numbers, latitude terms and Montenbruck ephemerides.
- Displacements are computed at each epicentre and at ±0.02°. The surface strain tensor is taken from their first derivatives on a sphere.
- Assumptions:
  - The surface strain is assumed to be the same at 30–45 km depth (depth/radius ≈ 0.006).
  - Stress uses **plane stress** (σzz = σxz = σyz = 0).
  - Moduli are homogeneous: μ = 49.6 GPa, ν = 0.249. They come from USGS CVM v1.7 `south_VI`, averaged over 29.5–42 km (Vp 7.05, Vs 4.08 km/s, ρ 2.98 g/cc from Brocher 2005). Body stress scales linearly with μ.
- **pyTMD quirk:** the `'N'` output of `SET_displacements` in pyTMD 3.0.9 points along +colatitude, which is **southward**. Its correlation with the analytic −u_θ is −0.99995. The code accounts for this.

**Ocean load** (`scripts/s5_tides_ocean.py`)
- Each ocean cell is treated as a vertical point force ρw·g·η·A, with ρw = 1025 kg/m³. Stresses use the Boussinesq solution for a **homogeneous elastic half-space**, with ν = 0.25 (Johnson 1985, eq. 3.21). This stress does not depend on μ. They are evaluated at each family's catalogue depth.
- Cells: 0.02° cells within 330 km of 48.6N 123.6W (34,120 ocean cells), plus 0.125° cells from 330 to 700 km (3,445 cells). Land and lakes are masked with Natural Earth 10 m.
- Ten constituents are summed in the constituent domain: q1 o1 p1 k1 n2 m2 s2 k2 2n2 m4. Minor constituents are inferred by admittance, and hourly series are produced with pyTMD FES nodal corrections.
- **Elevation field:**
  - Open Pacific: **EOT20** (Hart-Davis et al. 2021).
  - Inland waters: **tide-gauge interpolation only** (IDW, k = 6, power 2). This covers Juan de Fuca Strait east of the Cape Flattery–Carmanah Pt mouth, Haro Strait, the San Juan Islands, Puget Sound, Hood Canal, the Strait of Georgia and the Discovery Islands, all inside one polygon in the script.
  - The global models fail inland. At Victoria, EOT20 gives M2 = 0.97 m against 0.37 m observed, with phase errors of 60–130°.
  - About 16 % of the inner cells come from the gauges.
- **Gauges:** 50 NOAA CO-OPS gauges, using their published constants. 167 CHS gauges, harmonically analysed from DFO IWLS hourly predictions for 2023. Gauges in restricted basins or rivers were dropped: the Victoria Gorge, Portage Inlet, Porpoise Bay and the lower Fraser.
- **Model choice:** TPXO and FES2014/2022 need registration, so they were not used. EOT20 and GOT5.5/5.6 can be downloaded freely.

## Validation (all numbers are in the CSVs here)
- **Boussinesq** (`boussinesq_validation.csv`): all 13 checks pass.
  - The force balance ∫σzz dA is −0.9999.
  - The circular-load closed forms for σz and σr match within 0.5 % for a = 10, 30 and 100 km.
  - The equilibrium residual is about 1e-7.
  - Signs are correct: horizontal compression gives positive thrust τ and negative σn.
- **Gauge analysis convention** (`gauge_convention_check.csv`): re-analysing NOAA's own Port Angeles predictions reproduces NOAA's constants. Amplitudes match within 2 mm and phases within 0.15°.
- **Body tide** (`body_check_iers_vs_analytic.csv`): compared with an independent Love-number calculation (nominal h2/l2/h3/l3, from the lunisolar potential), the correlation is 0.99 to 0.9999. The IERS areal strain is 9 % smaller. That is expected: the diurnal resonance lowers h for K1, and the nominal-number calculation leaves this out.
- **Ocean model vs gauges** (`model_vs_gauges.csv`):
  - On the open coast (La Push, Neah Bay, Tofino, Bamfield), the M2 vector misfit is 0.015–0.08 m for EOT20 and 0.02–0.08 m for GOT5.5.
  - Inland, the median M2 misfit is 1.40 m for EOT20 and 0.37 m for GOT5.5. This is why the gauge field is used inland.
- **Inland gauge field, leave-one-out** (`gauge_field_loo.csv`): the median vector misfit is 0.08 m for M2 (11 %), 0.04 m for K1 (5 %) and 0.02 m for O1 (5 %).
- **Sensitivity** (`variant_sensitivity.csv`, median over families):
  - GOT5.5 instead of EOT20: under 1 % change.
  - Slab2 depth instead of catalogue depth: 3–7 % (90th percentile 8–18 %).
  - Dropping the 330–700 km ring: 1–3 % for τ, ≤ 0.3 % for σn.
  - **EOT20 without the gauge field** is badly wrong: τ M2 falls to 0.47×, σn M2 rises to 2.4×, and the vector difference is 57 % / 234 %. The inland treatment is what this product depends on.
- **Phase sanity check:** ocean σn for M2 at family 10 has phase 149°. The local water level at Port Angeles and Victoria has phase 307–318°. So high water clamps the fault, as it should.
- **Literature:** reported tidal stresses on Cascadia tremor and LFE sources are of order kPa. Royer et al. (2015) give about 2 kPa of shear stress, as summarised in Table 4 of Zhou et al. (2026, arXiv 2602.06703). Here τ has a median std of 1.2 kPa and a 99 % range of 5.2 kPa (about ±2.6 kPa), so the size agrees.
  - I could not access the full text of Royer 2015, Houston 2015 or Hawthorne & Rubin 2010. **Phases have not been compared with published values.** That comparison is the first item for the S5 reproduction step.

## What the numbers say
Medians across the 115 families, in Pa (amplitude per constituent):

| series | M2 | S2 | N2 | K1 | O1 | P1 | std |
|---|---|---|---|---|---|---|---|
| τ total | 1440 | 429 | 287 | 514 | 331 | 180 | 1221 |
| τ body | 152 | 71 | 29 | 166 | 152 | 66 | 209 |
| τ ocean | 1344 | 366 | 278 | 424 | 215 | 137 | 1079 |
| σn total | 1637 | 424 | 370 | 2465 | 1471 | 777 | 2487 |
| Coulomb μ′ = 0.4 | 1365 | 428 | 289 | 739 | 475 | 227 | 1323 |

- **The ocean load dominates.** On the τ std, ocean/body is 1.8–8.2 with a median of 5.2. On σn, the body tide is negligible (median std 36 Pa).
  - This does **not** match the expectation in PLAN Revision 1 that the two are the same size.
  - The reason is the shallow dip (about 10°). Body-tide stress is almost all horizontal, so it projects weakly onto the plane. The water load sits directly overhead and has strong coastline edges.
- **τ is mainly semidiurnal and σn is mainly diurnal.** For τ, M2 is about 2.8× K1. For σn, K1 is about 1.5× M2. Any S5 test that separates M2 from K1/O1 will therefore be testing different stress components.
- Family-to-family range: τ M2 is 0.77–2.1 kPa and σn K1 is 1.1–4.1 kPa. Families under the strait are the highest.
- Figures:
  - `figures/s5_tidal_stress_month.png`: family 10, February 2016.
  - `figures/s5_tidal_stress_spectrum.png`: spectrum of the whole record.
  - `figures/s5_tidal_amplitude_map.png`: amplitude map.

## Known limits
- **Ocean-load model**
  - The half-space is homogeneous and elastic. There is no layering, no Earth curvature and no self-gravitation.
  - Loads beyond 700 km are left out. Loads from 330 to 700 km are included in the half-space approximation and contribute 1–3 %.
- **Body-tide model**
  - Plane stress is assumed at depth, and the 30–45 km depth correction is ignored. Body σn is therefore uncertain at the level of tens of Pa, which is small next to the ocean σn.
- **Total stress only**
  - There is no poroelastic or undrained pore-pressure response. Effective σn could be smaller by a Skempton-B factor.
- **Inland water field**
  - It is an IDW interpolation using straight-line distances that can cross land. Leave-one-out errors are about 11 % for M2 and 5 % for K1.
  - Families sit under these waters, so this is the largest error in the ocean part.
- **Missing loading sources**
  - Non-tidal ocean loading is not included: surge, seasonal cycle, S1, Sa/Ssa.
  - Hydrological loading is not included.
  - Long-period ocean tides (Mf/Mm) are not included. Body Mf/Mm are included.
- **Fault geometry**
  - The geometry is one plane per family, from Slab2.
  - Depth is the catalogue depth. The Slab2 depth changes amplitudes by about 5 %.

## Sources and downloads (under `/home/jovyan/tremorferometry/data/tides/`)
- **EOT20** (DGFI-TUM; SEANOE doi:10.17882/79489; https://www.seanoe.org/data/00683/79489/data/85762.zip).
  - The download is 2.33 GB (`eot20_dl/85762.zip`, which also contains `load_tides.zip`, not extracted).
  - The extracted `ocean_tides.zip` is 1.21 GB and sits in `eot20_dl/`.
  - The extracted ocean tides take 1.8 GB in `EOT20/ocean_tides/`.
- **GOT5.5 and GOT5.6** (R. Ray, NASA GSFC; https://earth.gsfc.nasa.gov/sites/default/files/2024-07/GOT5.5.tar%201.gz and the 5.6 equivalent).
  - The downloads are 724 MB and 108 MB (`got_dl/`).
  - The extracted GOT5.5 ocean tides and the GOT5.6 extra constituents take 371 MB in `GOT/`.
- **Gauges:** NOAA CO-OPS metadata API (harcon), and DFO/CHS IWLS API (`wlp` hourly for 2023). No login was needed for either.
- **Software:** pyTMD 3.0.9 and timescale 0.1.3, installed with `pip install --user`.

Scripts, in run order (all in `svi_science/scripts/`): `s5_tides_geometry.py`, `s5_tides_body.py`, `s5_tides_gauges.py`, `s5_tides_ocean.py`, `s5_tides_combine.py`. Shared helpers are in `s5_tides_common.py`.
