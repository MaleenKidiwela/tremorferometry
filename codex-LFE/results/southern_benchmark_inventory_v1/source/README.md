# catalog

[![DOI](https://zenodo.org/badge/244437521.svg)](https://zenodo.org/badge/latestdoi/244437521)

This repository contains the Python scripts associated to the paper:

Ducellier A., Creager K.C. An eight-year-long low-frequency earthquake catalog for Southern Cascadia, *to be submitted*.

The templates for the FAME catalog (2007-2009) are computed with src/compute_templates.py. They are stored in src/templates.

The FAME catalog is available in data/Ducellier/catalogs/XXXXXX.XX.XXX/catalog_2007_2009.pkl

The templates for the networks catalog are computed with src/compute_new_templates.py. They are stored in src/templates_2007_2009. We only kept the templates with a good signal-to-noise ratios that are stoted in templates_v1.

The networks catalog is available in data/Ducellier/catalogs/XXXXXX.XX.XXX/catalog_2004_2011.pkl

The instrument responses were downloaded using src/get_responses.py, and then the LFE detection was done using src/find_all_LFEs_parallel.py.

There is one directory data/Ducellier/catalogs/XXXXXX.XX.XXX for each LFE family. Each directory contains two files corresponding to the two catalogs. Each file contains a pandas data frame with the following columns: year, month, day, hour, minute, second (of the event), cc (average cross-correlation between waveform and template), nchannel (number of seismic channels recording at that time).

Figure 1 is done using paper/figures/map_LFEs.

Figure 5 is done using data/Ducellier/cumulative_LFEs.py.

Figures 3 and 6 are done using data/Ducellier/plot_tremor_nbLFEs.py.

Figure 7 is done using data/Ducellier/plot_zoom_ETS.py.

Figure 8 is done using data/Ducellier/plot_several_families.py.

Figure 9 is done using data/Ducellier/plot_nb_events.py and data/Ducellier/LFE_distribution_perm/nb_LFE_events.txt.

Figure 10 is done using src/plot_several_tidal_stress.py.
