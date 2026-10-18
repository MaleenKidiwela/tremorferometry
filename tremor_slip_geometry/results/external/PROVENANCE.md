# External data — provenance

## results_for_paper.mat
Bartlow, N. (2020), "Time series from PBO and UNR used in Bartlow (2019)", Mendeley Data, V1,
doi:10.17632/mc49zmg7n7.1, licence CC-BY-4.0. Supporting data for
Bartlow, N. M. (2020), "A Long-Term View of Episodic Tremor and Slip in Cascadia",
Geophys. Res. Lett. 47, e2019GL085303, doi:10.1029/2019GL085303.

Downloaded 2026-09-10 from
https://data.mendeley.com/public-files/datasets/mc49zmg7n7/files/66392424-e44d-4d66-b949-44eaa5b1799a/file_downloaded
size 169711 bytes; sha256 093e48cd51bbfc9db16ab8a17b5f98b4ce09259fcde58c4b3a45a42f199294d4
(matches the Mendeley manifest exactly).

The other two files in the dataset were NOT downloaded (not requested, and large):
PBO_stations.mat (787 MB), UNR_stations_new.mat (1.51 GB).

### Contents
    Lats, Lons, ID_G                (443,)      GPS station latitude, longitude, 4-char code
    Ev_SSE, Nv_SSE, Uv_SSE          (443,)      long-term-average ETS velocity field (E, N, U)
    sig_Ev_SSE, sig_Nv_SSE, sig_Uv_SSE (443,)   1-sigma uncertainties
    el                              (3607, 3)   mesh triangle connectivity (node indices)
    slip_rate                       (3607,)     LTA-ETS slip rate per triangle (min 0.000, max 0.042, mean 0.004)
    strike, rake, ds_fact, ss_fact  (3607,)     per-triangle geometry / rake factors
    Pole                            (2,)        JdF-NA Euler pole [lat, lon] (MORVEL)

### LIMITATION — slip_rate cannot be georeferenced from public data
The mesh NODE coordinates (`nd`, `nd_ll` in the published script) are not in this deposit. They live in
the Green's-function files `homog_GFs_Nov2019.mat` / `Heter_GFs_May2019`, which the author's README
states were "not computed specifically for this study" and which are not on Mendeley and not in the
code repository (https://github.com/noelbartlow/Long-Term-Cascadia-ETS-velocities contains only
Cascadia_secular_fits.m, simple_inversion_1comp_2019.m, remove_outliers.m, README.md).
To read Bartlow's Oregon slip rate we must (a) digitise her published figure, (b) request the
mesh/GF file from the author, or (c) rebuild the mesh from her source.

### What IS usable: her station distribution
443 stations, 38.83-53.91N, -129.96 to -117.00 deg. 330 west of -120 deg (the `ind_GPS` cut in her script).
Per-degree counts (lon < -120): 39-40:37, 40-41:40, 41-42:21, 42-43:21, 43-44:21, 44-45:26,
45-46:31, 46-47:38, 47-48:35, 48-49:29, 49-50:15, 50-51:4.
Central Oregon 44-46N = 57 stations; north 46-48N = 73; ratio 0.78.
