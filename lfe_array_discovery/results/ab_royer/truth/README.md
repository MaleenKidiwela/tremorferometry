# Royer & Bostock (2014) template truth table, digitized from the vector figures

Paper: Royer & Bostock 2014, EPSL 402:247-256 (`notes/1-s2.0-S0012821X13004780-main.pdf`),
supplement `notes/royer2014_supp.pdf`, and Royer's PhD thesis (UBC 2014,
`notes/royer_sources/royer2014_phd_thesis.pdf`). No catalogue was ever released. Every
figure used here is vector graphics, so the symbol coordinates come straight from the
PDF drawing commands (pymupdf `page.get_drawings()`). Nothing was clicked by hand.

Regenerate everything with `python digitize_royer2014.py` (about 2 min).

## Main output: `royer2014_templates_digitized.csv`
The file has one row per symbol, in three sets (column `locator`):

| locator | source | NVI | NW | SW | SVI* | depth | template_id |
|---|---|---|---|---|---|---|---|
| Hyp2000 | Fig 4 panels B (NVI), C (SVI), D (NW, SW) | **47** | **122** | **54** | 138 | from Fig S2 sections | thesis family ID for 88 NW and 44 SW templates |
| hypoDD | Fig S1, same panels | 43 | 95 | 44 | 127 | from Fig 5 sections | nan |
| thesis_family | thesis Fig 4.1 | - | 100 | 49 | 93 | nan | **family ID 1-242** (all rows) |

*SVI = Bostock et al. 2012 templates. They are shown on the same maps and kept as extra rows.

Columns: `subarray, source_fig, locator, lon, lat, depth_km, template_id, symbol_style,
georef_residual_km`. The next group adds detail: `tick_rms_km` (graticule fit), `model_spread_km`
(largest difference between projection models at this point), `km_per_pt`, `depth_n` (number of
sections that gave a depth), `depth_spread_km` (max minus min across those sections),
`depth_profiles`, and `min_gap2_km` (along-profile gap to the next-nearest candidate, which is the
depth-assignment ambiguity). Last: `page_x, page_y`, `thesis_match_km`.
`georef_residual_km = sqrt(tick_rms_km² + model_spread_km²)`.

**Counts compared with the paper.** The Hyp2000 set gives exactly 122 / 54 / 47, matching the text.
For hypoDD the text says 77 / 72 / 91 % retained, which predicts 94 / 39 / 43. The figure has
95 / 44 / 43. NVI matches. NW is off by 1. **SW has 44 symbols, where 72 % predicts 39.** Either
the percentage in the text is wrong, or Fig S1 plots events that hypoDD did not keep; the figure
cannot tell us which. There are no duplicate symbols.

Thesis figures 3.4 and A.1 are copies of EPSL Fig 4 and Fig S1. Their colour counts are identical
(122/54/47/138 and 95/44/43/127), so they add no new rows.

## Method
1. **Symbols.** Template symbols are filled 4-vertex diamonds. They are told apart by fill colour
   (orange NVI, red NW, green SW, yellow SVI; the RGB values differ between Fig 4 and S1) and by
   shape and size (about 3.1×4.0 pt, taller than wide). Station triangles are 3-vertex or
   12-segment paths that are wider than tall, so they are excluded. Coastlines, contours and the
   black earthquake dots are other shapes or colours. Symbol centre = centre of the bounding box.
2. **Georeferencing.** MATLAB exported every line as a thin filled polygon, so the tick marks are
   detected as 0.3×1.8 pt rectangles touching the axes box. Each axes box is found from its white
   background rectangle. Tick values come from the real text labels of Fig S1; in Fig 4 the labels
   are glyph outlines, but the tick layout is identical. The graticule converges, with top and
   bottom tick spacing differing by about 4 %, so this is a conic-type m_map projection and not an
   affine map. Each tick is fitted as a graticule-line/frame intersection, using the projection
   plus a similarity transform. Models tried: LCC with standard parallels fixed at the map limits
   (the m_map convention), tmerc, eqdc, aea and stere. The adopted model is **aea with free
   standard parallels**. The spread between aea, tmerc, lcc and eqdc is carried as
   `model_spread_km`. **Tick residuals, RMS / max:**
   - Panel A: 0.003-0.005 pt = 0.014-0.023 km RMS, 0.048 km max
   - Panel B: 0.004 km / 0.007 km
   - Panels C and D: 0.011-0.012 km / 0.023 km
   - Per tick: `royer2014_georef_ticks.csv`
   - The other models fit to at most 0.1 pt; their disagreement at the symbol positions is at most
     about 0.2 km (`model_spread_km`). This is the main georeferencing error. The full
     `georef_residual_km` is 0.01-0.23 km (median 0.16 km for NW, where D is 1.7 km/pt).
3. **Independent checks of the georeferencing.**
   - (a) The overview panel A (5 km/pt) and the zoom panels were georeferenced separately. Matched
     one-to-one, they agree to a median of 0.07-0.13 km and a maximum of 0.20 km for every
     subarray (`royer2014_overview_crosscheck.csv`).
   - (b) Thesis Fig 4.1 is an entirely different plot: a GMT Mercator map (Mercator tick RMS
     0.04 pt, against 1.05 pt for linear latitude). Against it, 88/100 NW families and 45/49 SW
     families fall within 0.3 km of an EPSL Hyp2000 template, with median nearest-neighbour
     distance 0.05 km (`royer2014_thesis_vs_epsl.csv`).
4. **Depths.** The 10 section panels are straight linear axes: −100..200 km by 0..80 km at the
   frame, with ticks checked, 1.49 km/pt, equal aspect.
   - Profile lines are the white bars in map panel A. They are identified by the X / X′ text labels
     of Fig S1: H, I, J over NVI; A, B over SVI; C, D, E over NW; F, G over SW. Endpoints are in
     `royer2014_profiles.csv`.
   - Map epicentres are projected onto each bar's great circle to give along- and cross-track km.
   - Per profile and colour, an offset c (d_section = d_map + c) is grid-searched, with a 1-D
     Hungarian match to the section diamonds.
   - Result: **c = −80.0 ± 0.5 km for H, I, J, C, D, E, F and G, in both locators**, with median
     |Δd| of 0.03-0.12 km and no match worse than 1 km. So the sections start 80 km from the SW
     end of each bar (bars are about 330 km long).
   - Depth = median over the sections where |Δd| ≤ 1 km.
   - Coverage: Hyp2000 NVI 47/47, NW 120/122, SW 54/54. hypoDD 43/43, 95/95, 43/44.
   - **Depth ranges:** Hyp2000 NW 31.7-46.6, SW 32.5-49.9, NVI 33.7-39.8 km. The text (for the
     Hyp2000 section figure S2) gives NW 31-46, SW 32-49, NVI 33-39. hypoDD: NW 36.2-43.3, SW
     34.4-40.8, NVI 33.5-38.5.
5. **Thesis Fig 4.1 family IDs.** There is one diamond per family, all 242 with unique fill
   colours. Each colour is located on the embedded colour-bar image, converted to a continuous ID
   using the bar's tick labels, then ranked to give IDs 1-242. SVI is 1-93, NW 94-193, SW 194-242.
   Validation:
   - Within each catalogue the IDs have **zero latitude inversions** (the caption says IDs run
     north to south).
   - The three labelled families in the figure, 49, 123 and 242, decode correctly.
   - Nearest-ID colour separation is ≥0.54 ID.
   - IDs were joined onto Hyp2000 NW/SW templates by one-to-one matching within 0.3 km, giving
     88 NW and 44 SW templates with an ID.
   - These are the chapter-4 "regrouped geographically independent" families: 100 NW (from 122
     templates) and 49 SW (from 54). The 12 NW and 5 SW families that did not match within 0.3 km
     are presumably merged or relocated groups. Thesis SVI families sit a median 2.2 km from the
     EPSL SVI symbols, so they were not joined.

## Ambiguities and failures
- **EPSL figures carry no template numbers.** The mechanism maps (Figs 7-10) have unlabelled
  beach balls, NVI has leader-line offsets, and the maps show only a selection of templates. They
  were not digitized. The only IDs come from thesis Fig 4.1.
- **hypoDD SW count (44 against 39 expected):** see above.
- **Depth assignment is 1-D.** Two templates at similar along-profile distance can swap depths.
  `min_gap2_km` gives the gap to the next candidate: median 0.1-0.2 km, occasionally 5-10 km.
  Templates seen in several overlapping swaths are a consistency check. Their `depth_spread_km`
  median is about 0.01 km, but some reach 4-6.6 km (NW Hyp2000), which indicates a probable swap.
  Treat depth_spread > 1 km as unreliable.
- **Swath counts do not match exactly.** The number of map symbols within ±25 km of a bar
  differs from the section count by up to 5 (for example C NW: 57 in the section, 54 in the
  swath). The swath definition or distance metric is not exactly ours. Matching used ±30 km.
- **Unmatched.** 2 NW Hyp2000 templates and 1 SW hypoDD template got no depth (nan).
- **SVI profiles A and B (not needed for NW/SW/NVI).** A fits c = −83.5 km with 15 mismatches
  over 1 km and 100 section points against 70 in the swath. B fits c = 0. These bars do not follow
  the −80 km convention, so SVI depths on those profiles are less reliable.
- **Symbol centre = data point** (the MATLAB marker convention). A bounding-box centre error of
  about 0.05 pt is negligible.

## QC figures
- `qc_map_Fig4.png`, `qc_map_FigS1.png`: rasterised panels with the fitted graticule (cyan),
  which lies on the printed ticks and labels, and the digitized lon/lat re-projected (+/o), which
  lies on every diamond.
- `qc_sections_FigS2.png`, `qc_sections_Fig5.png`: digitized diamond centres drawn over the
  sections.
- `qc_truth_table.png`: final Hyp2000 and hypoDD sets coloured by depth with the profile lines,
  plus the thesis families labelled with their IDs.

## Scoring helper
`lfe_array_discovery/scripts/ab_truth_score.py OUR.csv [--truth-set Hyp2000|hypoDD|thesis_family]`
takes our families (lon, lat, depth).
- Region: our points within 30 km of a subarray's templates. Our families carry no subarray label,
  so NW and SW candidates can overlap.
- Matching: one-to-one Hungarian on epicentral distance.
- Real scores: matched fraction within 3/5/10 km and depth residuals.
- Null: our points rotated and shifted randomly, with the same matching. Null mean, null 95th
  percentile and p-value are printed next to the real values.

Self-test (hypoDD as "ours" against Hyp2000 truth):

| | NW | SW | NVI |
|---|---|---|---|
| Real, within 5 km | 0.58 | 0.57 | 0.92 |
| Null mean | 0.11 | 0.09 | 0.22 |

Depth offsets are small: median −0.35 to +0.19 km.

Other files: `royer2014_all_panel_symbols.csv` (includes overview-panel duplicates),
`royer2014_section_points.csv`, `royer2014_section_matches.csv`, `royer2014_profile_fits.csv`,
`royer2014_summary.csv`. Code: `pdfvec.py` (path parsing), `georef.py` (projection fit),
`panels_cfg.py` (tick values), `digitize_thesis_fig41.py`.
