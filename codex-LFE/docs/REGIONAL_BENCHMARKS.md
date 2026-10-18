# Independent regional benchmark registry

Evidence checked directly against publications on 2026-09-11. A published catalog's existence is distinct from having obtained its machine-readable release. This registry must not be replaced by the old repository's ground-truth inventory.

| Region | Primary evidence | What is available in this workspace / next action |
|---|---|---|
| Southern Vancouver Island | [Bostock et al. 2015](https://doi.org/10.1002/2015JB012195); [Lin et al. 2024](https://seismica.library.mcgill.ca/article/view/1134) | Existing copies of published catalog data and local raw miniSEED. The independent CPU pilot uses the uncorrected Bostock source table and derives only training-based arrival anchors. Audit source identity and obtain original family identifiers before release. |
| Olympic Peninsula, Washington | [Chestler and Creager 2017](https://doi.org/10.1002/2016JB013717) | The paper's final sorted catalog contains 32,387 distinct LFEs in 43 families using AofA/CAFE data. Machine-readable event records have **not** been acquired by codex-LFE. Waveform availability is documented in the publication; temporary-array results do not prove permanent-network sensitivity. |
| Central Oregon | [Thomas and Bostock 2015](https://doi.org/10.1016/j.tecto.2015.07.013), [author PDF](https://pages.uoregon.edu/amthomas/pubs/thomas_tectonophysics_2015.pdf) | Nearly 16,000 events in 18 families during the 2009 episode. A concrete regional benchmark exists. Event tables/templates have **not** been acquired. Locate the data release and check surviving stations for a modern-network pilot. |
| Southern Cascadia | [Ducellier and Creager 2022](https://doi.org/10.1029/2021JB022986), [author manuscript](https://arxiv.org/abs/2212.14581) | Published eight-year extension of known southern families. Source audit and bounded southern detections are now complete. Seven templates represent five published families; the next-day check has 23 within-family candidate rows, with mixed independent-station confirmation. See `results/ducellier_validation_summary_v1`; these are not a finalized catalog. |
| Southern edge / Mendocino | [Shelly et al. 2025](https://doi.org/10.1029/2025GL116116) | Existing catalog copy is a starting point; verify release metadata, family coordinates, and clock convention before use. |

The central Oregon paper is particularly relevant: it reports a failed autocorrelation attempt followed by successful cross-station discovery and iterative matching. That supports testing a different method when one implementation fails. It does not justify importing its thresholds as universal settings.

No blanket claim that Washington or Oregon lacks published LFE examples is supported by this review. Equally, neither these papers nor their catalogs establishes complete coverage throughout the margin.

Next acquisition task: obtain machine-readable Oregon/Olympic references and waveform/template metadata from public releases. If the only remaining route is contacting authors, prepare a concrete request for the user; do not send messages without explicit authorization.

Additional Washington reference route checked on 2026-09-14: Royer, Thomas & Bostock (2015), https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2014JB011430, reports NW and SW catalogs (100/49 grouped families; 60,527/22,868 detections). Its data statement explicitly says processed waveform data are available by request from Michael Bostock; publisher supplements listed are figures/text/tables, not a linked event-time catalog. This is an alternative specific request route, not acquired phase truth. No author contacted.
