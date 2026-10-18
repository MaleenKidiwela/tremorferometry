# How well GNSS resolves slow slip under the SVI LFE families (2026-10-07)

Model: Gualandi (2025) vbICA solution `cascadia-sse-strain/data/2026-08-18`. It is rebuilt with the repo's own code, read-only, on its Slab2 fallback mesh (`data/cascadia_mesh.txt`): 742 triangles from a regular 0.25° grid, median edge 27.8 km, median area 274 km². There are 13 retained ICs, an exponential prior with λ = 21.09 km, and σ0 per IC from the CV curves.
Scripts: `lfe_extension/scripts/slip_build.py`, `slip_resolution.py`, `slip_figures.py`. Numbers come from `resolution_metrics.json`, `resolution_patch_table.csv`, `resolution_pairs.csv`, `svi_gnss_stations.csv`, `episode_slip_summary.json` and `sanity_repo_reproduction.json`.

## Verdict in one paragraph
GNSS sees SVI slow slip as one blob roughly 35 km in radius, not slip on individual families. A unit thrust spike on any patch that holds a family comes back with only **6 % of its value at that patch**, with a half-potency radius of **31–46 km (median 34 km)**. The total potency is recovered (median 0.95). The posterior variance is still **72 % of the prior** at the patch scale, so the data barely constrain individual patches. Two family patches 15–25 km apart have averaging kernels that correlate at **0.84** (0.96 below 15 km). The correlation only drops below 0.5 at about 35–50 km. In the model's own output, the first mode explains **90 %** of the variance of the 115-family × 7169-day slip matrix, and 3 modes explain 99 %. Daily slip series at families 20–30 km apart correlate at a median of 0.96.

What GNSS can do:
- Episode timing and total slip for the SVI region as a whole. Major-episode increments are 9–18 mm, against a 1σ of about 0.75 mm per episode.
- A broad gradient across the cloud, at the factor-of-~1.8 level. Within an episode, the 90th-percentile family slip is about 1.8× the 10th percentile.

What GNSS cannot do:
- Separate slip on neighbouring families (< 25–35 km apart).
- Resolve an up-dip versus down-dip contrast across the 47 km (5–95 %) down-dip width of the cloud. The kernel's down-dip RMS width is 28–57 km, so the down-dip structure in the model comes from the prior.
- Assign slip of a few mm to one episode.

## 1. GNSS stations (213 in the solution; from X.mat `llh`)
- Family-cloud centroid: 48.484 N, 123.749 W.
- Stations within 25 / 50 / 75 / 100 / 150 km of the centroid: **1 / 6 / 12 / 18 / 38**.
- Nearest station to each family: median 19 km (range 0.6–36 km).
- Closest stations: ALBH 22 km, PGC5 29 km, PTAA 45 km, CLRS 47 km, SC04 49 km, PTRF 50 km, P435 51 km, SC02 55 km, P439 67 km, P436 67 km, NEAH 68 km, SC03 74 km.
- Station availability dropped between eras:
  - PTAA: 30 % of days in 2007–2013, 90 % in 2014–2026.
  - CLRS and PTRF: 70 % in 2007–2013.
- So slip before 2014 under the western part of the cloud rests on fewer stations (see `svi_gnss_stations.csv`).
- There is no station seaward of the families on Vancouver Island except NEAH. Coverage is one-sided: stations sit mostly above and land-ward of the cloud.

## 2. Smoothing prior and patch size
- The spatial prior is exponential with correlation length **λ = 21.09 km** (`options.inversion.lambda_*`) and m0 = 0.
- σ0 per IC is 0.003–0.03, chosen by cross-validation.
- Patch edges (~28 km) are comparable to λ. The 115 families fall in **29 distinct patches**. The most crowded patch holds 16 families.
- The cloud spans 113 km along strike and 62 km down-dip; the 5–95 % ranges are 74 km and 47 km.
- Three ICs carry SVI slip:

  | IC | share of SVI slip variance |
  |---|---|
  | IC10 | 44 % |
  | IC11 | 28 % |
  | IC22 | 17 % |

- IC17 has a −88 excursion in late 2024 (minimum 2024-12-05). It carries **0.2 %** of SVI slip variance, so it is irrelevant here.

## 3. Resolution test near SVI
The repo's `fig01_resolution_test.png` is not a GNSS resolution test. It measures strain roughness against evaluation offset. So I computed the model resolution matrix of the same linear inversion:

R_k = Cm0 Gᵀ(G Cm0 Gᵀ + Cd_k)⁻¹G for each IC, weighted by the IC shares above.

Spike tests are the columns of R; averaging kernels are its rows. Results over the 29 family patches:

| metric | min | median | max |
|---|---|---|---|
| R diagonal (recovered at spike) | 0.03 | 0.06 | 0.08 |
| spike potency recovered (sum over mesh) | 0.71 | 0.95 | 1.25 |
| posterior / prior variance | 0.68 | 0.72 | 0.82 |
| spike half-potency radius (km) | 31 | 34 | 46 |
| kernel RMS width along strike (km) | 34 | 45 | 59 |
| kernel RMS width down-dip (km) | 28 | 41 | 57 |
| spike centroid shift down-dip (km) | −29 | −11 | −3 |

- The centroid shift is negative, meaning the recovered slip moves **up-dip**, towards the stations.
- Averaging-kernel correlation against patch separation:

  | separation | correlation (median) |
  |---|---|
  | < 15 km | 0.96 |
  | 15–25 km | 0.84 |
  | 25–35 km | 0.67 |
  | 35–50 km | 0.44 |
  | 50–75 km | 0.17 |
  | > 75 km | ≈ 0 |

- The representative families are 31–47 km apart. Their kernels correlate at 0.37–0.65, but their model slip series correlate at 0.95–0.97.
- **Independent "pixels" across the family cloud**: roughly 74/40 × 47/40 ≈ **2 (along strike) × 1 (down-dip)**. That agrees with the modal count: 1 mode = 90 % of variance, 3 modes = 99 %.
- **Depth**: Bostock depths run from 28.8 to 45.9 km, but the down-dip kernel is ≥ 28 km wide (RMS). Slip is therefore **not resolved in depth or down-dip** within the cloud.
  - The model's own loading on the cloud-mean slip falls with depth: Spearman −0.43, range 0.43–1.16. That gradient follows the station geometry and the prior (kernels move up-dip), so it should not be read as physics.
- The mesh itself sits 2.8 km deeper than Bostock on average (median +3.3 km, range −3.5 to +7.1 km). Mesh depth at the families is 34–45 km. This comes from Slab2 versus the Bostock hypocentres, and is small next to the resolution.

## 4. Timing resolution
- **No temporal prior.** The vbICA sources are i.i.d. mixtures of 4 Gaussians per IC (`init_parameters.HMM = 0`). Daily V values are therefore unsmoothed and each day is an independent estimate.
- The repo's slip **rate** (my `slip_rate_mm_per_day*` columns) low-passes the amplitudes with a 365-tap Hamming FIR at 1/7 per day, then applies a 7-day linear fit. Treat it as resolving about **1 week**. The `slip_rate_raw_*` columns are plain daily gradients.
- Daily noise at a middle family is **~1.0–1.15 mm** (MAD of residual against a 31-day median). The residual autocorrelation is 0.29 at a 1-day lag and about 0 by 5 days.
  - The vbICA posterior `var_V` alone gives 0.64 mm (`slip_sd_mm`). That is a lower bound, because it holds the patterns fixed.
  - Only 3.5 % of the slip-series power is at periods shorter than 7 days.
- **Episode increments**: per-episode 1σ is about 0.75 mm (5-day means either side).
  - In 9 major episodes with SNR ≥ 10, the 25→75 % rise at a family takes a median of **8 days** (6–9 days, but about 40 days for the multi-pulse 2014-11, 2021-01 and 2021-10 episodes).
  - The 10→90 % rise takes a median of 37 days, because of noise and leakage from slip along strike.
- **Spurious increments**: 14 of 39 "major" LFE episodes show **negative** model increments of up to −3 mm at the middle family. Treat |Δslip| < ~3 mm (about 3σ) as unresolved.

## 5. Circularity check (tremor in the inversion?)
- **No tremor input was found.** `options.mat` contains GNSS selection, trajectory-model, PCA/vbICA and inversion settings.
- The inversion prior is purely spatial: exponential covariance, m0 = 0.
- Components are selected by Hurst-exponent and spectral criteria (`inversion.select_comps`).
- The only tremor-related fields are plotting-only: `plot.video.tremors_size`, and `plot.hypocenters`, which is empty.
- `ICA.mat`, `X.mat`, `ind_comps.mat` and `misfit_comps.mat` hold no tremor fields.
- The paper text is not on disk, so the paper itself was not checked. The plot options suggest tremor is overlaid only in display products, and the downloaded inputs show no tremor in any prior. **Comparing LFEs with this slip is not circular** in the sense of tremor being an input.
- The repo's own `DATA_SOURCES.md` lists tremor only as an independent cross-check.

## 6. Sanity check (repo reproduction)
- **Exact reproduction** of `run_real.log`:
  - Misfit per IC: 0.9393, 0.5195, 0.672, 0.5631, 0.9064, ….
  - Peak |slip| 3.1 cm, peak slip rate 108.7 cm/yr, peak moment rate 2.338e21 N m/yr.
  - The per-family patch series equals `Solution.slip()` row for row (difference 0.0).
- **Known episodes**, north of 47 N (thrust-slip increment, μ = 40 GPa as in fig09):

  | episode | Mw | peak slip |
  |---|---|---|
  | Aug 2010 | 6.70 | 22 mm |
  | Aug 2011 | 6.61 | 18 mm |
  | Sep 2012 | 6.65 | 13 mm |
  | Sep 2013 | 6.74 | 18 mm |
  | Nov 2014 | 6.73 | 16 mm |
  | Feb 2017 | 6.68 | 19 mm |

  This matches the README's "Mw 6.6–6.8, 1–2 cm".
- **Sign check against data**: modelled thrust during the Sep 2012 episode predicts ALBH moving −2.98 mm east; ALBH was observed moving −3.78 mm. For Aug 2024 the figures are −4.01 mm predicted and −3.86 mm observed.
- **Two repo issues found.** I did not modify the repo.
  1. **Sign.** On this mesh every triangle normal points down, so cutde's +dip-slip is *normal* sense. The repo's `slip_magnitude(apply_rake_flip=True)`, which is the quantity plotted in fig07, is therefore **−thrust**: correlation −0.999 in the 47.5–49.5 N strip (`fig_sanity_fig07_strip.png`). In fig07, slow-slip episodes appear as *drops*. The colour-bar note "negative = loading faster than secular" is inverted for this mesh. All my outputs use **thrust = positive**.
  2. **Dates.** Gualandi's timeline is year + (doy − 0.5)/ndays. The repo's `decimal_year_to_datetime64` rounds instead of flooring, so half the epochs are dated one day late (3575 of 7169) and 1898 dates are duplicated. My outputs use corrected UTC days, 2007-01-01 → 2026-08-18; one day, 2009-02-24, is missing in the solution.

## Column conventions (slip_svi_family_daily.parquet, slip_wa_nw_family_daily.parquet)
- `slip_mm`: slip along the family's rake, thrust positive. The rake is the principal axis of the (strike-slip, thrust) series in the 25 km variant: SVI rakes 61–95°, and the axis explains ≥ 91 % of the variance.
  - Values are relative to the secular trend and to an ICA-centred mean. **Only differences in time are meaningful.** Inter-ETS periods show negative slope by construction.
- `ss_mm` and `ds_thrust_mm`: strike-slip in cutde's sign, and thrust-positive dip slip.
- Spatial variants:
  - no suffix: the patch containing the epicentre;
  - `_r10` and `_r25`: area × Gaussian weights (σ = R/2) over patch centroids within 10 or 25 km (3-D, from the epicentre's point on the interface). These use 1–3 and 6–10 patches respectively.
- Slip rate: `slip_rate_mm_per_day*` is the repo method (~1-week resolution); `slip_rate_raw_mm_per_day*` is the plain daily gradient. `slip_sd_mm*` is the var_V-only 1σ.
- Metadata (`*_meta.csv`): patch, mesh depth at the epicentre, patch centroid depth, Bostock (or catalogue) depth, rake, and the number of patches in each radius.
- **WA-NW**: 132 families on 45 patches. Their catalogue depths are heterogeneous: the median mesh-minus-catalogue difference is +6.7 km, and some catalogue depths are ~0. The resolution matrix was not computed for WA, so these numbers do not carry over.
