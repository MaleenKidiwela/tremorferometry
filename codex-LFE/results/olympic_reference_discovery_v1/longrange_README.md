# longrange

This repository contains scripts to analyze catalogs of low-frequency earthquakes (LFEs) and explore whether there is long-range dependence in the time series.

I use different graphical methods to compute the fractional differencing parameter using synthetic time series and then I apply the methods to several LFE catalogs in different regions.

The data directory contains scripts to transform catalogs into time series by counting the number of earthquake events per time window and to visualize the catalogs. The actual data are not on Github but are available upon request.

The src directory contains the scripts to compute the fractional differencing parameters, launch the scripts for the analysis of the time series, and make plots to see the results.

test_long_range.py contains the functions to compute the Hurst parameter $H$ or the fractional differencing parameter $d$ using the following methods:
- Absolute value method
- Variance method
- Variance of residuals method
- R/S method
- Periodogram method

Please see:

Taqqu, M. and V. Teverovsky (1998). “On estimating the intensity of long-range dependence in finite and infinite variance time series”. In: A Practical Guide to Heavy Tails: Statistical Techniques and Application. Ed. by R.J. Adler, R.E. Feldman, and M.S. Taqqu. Boston, MA, USA: Birkhäuser

for more details on the methods.

test_long_range_parallel.py contains a parallelized version of these functions to accelerate the computation.

The files compute_XXX.py are scripts to use the functions from test_long_range_parallel.py on specific time series, write the Hurst parameter or fractional differencing parameter in a pkl output file and write the data to make the plots in output directories.

The files draw_XXX.py are scripts to read the results in the output directories and make the corresponding plots to scheck how well the graphical method works.

The R scripts generate synthetic ARFIMA time series to test the methods.
