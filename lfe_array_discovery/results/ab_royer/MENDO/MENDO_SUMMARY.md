# MENDO blind Royer-chain discovery: summary

Written 2026-09-30 06:20Z by scripts/mendo_summary.py. Queue scripts/mendo_discovery_queue.sh, log logs/mendo_discovery_queue.log.

Scrambled control = the same chain with each site's day series shifted by whole days (AB_NULL=1; shift within each year's list of scan days; stations < 3 km apart shift together).

| quantity | real run | scrambled control |
|---|---|---|
| autocorrelation hours | 1920 | pending |
| autocorrelation detections | 6527 | pending |
| culled detection catalogue | 10729 | pending |
| initial templates (UPGMA) | 1989 | pending |
| growth iterations | 8 | pending |
| families grown (final iteration) | 1178 | pending |
| detections (final iteration) | 510351 | pending |
| families scored (2-8 Hz stacks) | 1178 | pending |
| families with exactly 2 sites | 1 | pending |
| REAL families (>= 3 sites, < 3 km = 1 site) | 6 | pending |
| families located (Hyp2000 headline, all) | 26 | pending |
| real families located | 5 | pending |
| real located, ERH < 5 km | 2 | pending |
| real located, median depth (km) | 23.9 | pending |

## Real families vs Shelly et al. (2025)

- real families: 6 (located 5)
- (a) epicentre within 10 km of a Shelly family: 0
- (b) waveform duplicate of a Shelly v3b it0 template: 0
- (c) shared detection times with his catalogue: 1
- known (any of a/b/c): 1
- **NEW (fail all three): 5** (located 4)

New families (approximate area = headline Hyp2000 location):

|   grp |   n_sites_3km |     lat |      lon |    dep |    erh |   near_fam |   near_km |   c_peak_frac |   n_det |
|------:|--------------:|--------:|---------:|-------:|-------:|-----------:|----------:|--------------:|--------:|
|   243 |            18 |  40.58  | -124.245 |  14.85 |   1.14 |          8 |    70.999 |         0.004 |     774 |
|   316 |            17 |  40.287 | -124.086 |  20.41 |   5.81 |          8 |    38.836 |         0.002 |     902 |
|   445 |             6 |  39.854 | -123.408 |  31.91 |   4.65 |         27 |    22.826 |         0.003 |     738 |
|   527 |            19 |  40.472 | -124.574 |  24.15 |  12.37 |          8 |    84.856 |         0.003 |    2117 |
|  4042 |            13 | nan     |  nan     | nan    | nan    |        nan |   nan     |         0.003 |     861 |
